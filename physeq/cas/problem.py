# -*- coding: utf-8 -*-
#
# Copyright (c) 2026, Geoffrey M. Poore
# All rights reserved.
#
# Licensed under the BSD 3-Clause License:
# http://opensource.org/licenses/BSD-3-Clause
#


from __future__ import annotations

import math
import sympy
from astropy.units import Quantity
from collections import defaultdict
from . import wrapped
from .equation import Eq, WrappedEq, checked_solveset, solveset_for_ans
from .printing import latex
from .symbol import (ConstSymbol, Symbol, WrappedConstSymbol, WrappedSymbol, SymbolWithConstraints,
                     translate_xreplace_rule)




class UnevaluatedQuantity(sympy.UnevaluatedExpr):
    # This is only intended for creating solutions equations of the form
    # `<symbol> = <quantity>`
    def __new__(cls, arg):
        if not isinstance(arg, Quantity):
            raise TypeError
        obj = sympy.Expr.__new__(cls, arg)
        obj._assumptions['commutative'] = True  # type: ignore
        return obj

    def doit(self, **hints):
        return self.args[0]




class Problem(object):
    def __init__(
        self,
        setup: dict[ConstSymbol | Symbol | WrappedConstSymbol | WrappedSymbol,
                    int | float | Quantity | ConstSymbol | WrappedConstSymbol],
        knowns: list[ConstSymbol | Symbol | WrappedConstSymbol | WrappedSymbol | SymbolWithConstraints],
        unknowns: list[Symbol | WrappedSymbol | SymbolWithConstraints],
        equations: list[WrappedEq],
        definitions: list[WrappedEq] | None = None,
        simplify: bool = False,
    ):
        self.setup: dict[ConstSymbol | Symbol, ConstSymbol | Quantity] = {}
        if not isinstance(setup, dict):
            raise TypeError
        for k, v in setup.items():
            if isinstance(k, (WrappedConstSymbol, WrappedSymbol)):
                k = k.expr
            elif not isinstance(k, (ConstSymbol, Symbol)):
                raise TypeError
            if isinstance(v, ConstSymbol):
                pass
            elif isinstance(v, WrappedConstSymbol):
                v = v.expr
            elif not isinstance(v, Quantity):
                if isinstance(k, Symbol):
                    v = k.quantity(v)
                else:
                    raise TypeError
            self.setup[k] = v

        self.constraints: dict[Symbol, SymbolWithConstraints] = {}

        self._known_symbols_set: set[ConstSymbol | Symbol] = set()
        if not isinstance(knowns, list):
            raise TypeError
        for s in knowns:
            if isinstance(s, SymbolWithConstraints):
                self.constraints[s.symbol] = s
                s = s.symbol
            elif isinstance(s, (WrappedConstSymbol, WrappedSymbol)):
                s = s.expr
            elif not isinstance(s, (ConstSymbol, Symbol)):
                raise TypeError
            self._known_symbols_set.add(s)

        self._unknown_symbols_set: set[ConstSymbol | Symbol] = set()
        if not isinstance(unknowns, list):
            raise TypeError
        for s in unknowns:
            if isinstance(s, SymbolWithConstraints):
                self.constraints[s.symbol] = s
                s = s.symbol
            elif isinstance(s, (WrappedConstSymbol, WrappedSymbol)):
                s = s.expr
            elif not isinstance(s, (ConstSymbol, Symbol)):
                raise TypeError
            self._unknown_symbols_set.add(s)

        self._equations_list: list[WrappedEq] = []
        if not isinstance(equations, list):
            raise TypeError
        for eq in equations:
            if not isinstance(eq, WrappedEq):
                raise TypeError
            self._equations_list.append(eq)
            for s in eq.free_symbols:
                if isinstance(s, ConstSymbol) and s not in self._unknown_symbols_set:
                    self.setup[s] = s
                    self._known_symbols_set.add(s)

        self._definitions_list: list[WrappedEq] = []
        if definitions is None:
            definitions = []
        elif not isinstance(definitions, list):
            raise TypeError
        for eq in definitions:
            if not isinstance(eq, WrappedEq):
                raise TypeError
            if not isinstance(eq.wrapped.lhs, Symbol):
                raise TypeError(
                    'Only definition equations with a single symbol on the left-hand side are supported'
                )
            self._definitions_list.append(eq)
            for s in eq.free_symbols:
                if isinstance(s, ConstSymbol) and s not in self._unknown_symbols_set:
                    self.setup[s] = s
                    self._known_symbols_set.add(s)

        for k, v in self.setup.items():
            if isinstance(k, Symbol) and k in self.constraints:
                if not self.constraints[k].constraints(v.value):
                    raise RuntimeError(
                        f'In setup, symbol "{k}" has value "{v}" that fails constraints "{self.constraints[k]}"'
                    )
        self._all_symbols_set = self._known_symbols_set | self._unknown_symbols_set

        if not isinstance(simplify, bool):
            raise TypeError
        self.simplify = simplify


        self._solved_from_setup: dict[ConstSymbol | Symbol, Quantity] = {}
        self._all_quantities_or_constants: dict[ConstSymbol | Symbol, ConstSymbol | Quantity] = self.setup.copy()

        self._solve_numerical()
        self.knowns: dict[ConstSymbol | Symbol, Quantity] = {}
        self.unknowns: dict[ConstSymbol | Symbol, Quantity] = {}
        self.equations: dict[WrappedEq, list[WrappedEq]] = {}
        for k in self._known_symbols_set:
            v = self._all_quantities_or_constants[k]
            if isinstance(v, ConstSymbol):
                v = v.to_quantity()
            self.knowns[k] = v
        self._generate_solutions()


    def _solve_numerical(self):
        '''
        Get numerical values for all symbols, based on `setup`.
        '''

        # This is currently designed for problems that can be solved via a
        # sequence of equations with one unknown.  In the future, it should be
        # enhanced to support systems of equations, assuming that there is a
        # straightforward way to implement step-by-step solutions for those
        # cases (or those cases don't need detailed solutions).

        # Don't need equation wrappers to get a numerical solution; don't care
        # about expression order
        eqs: list[Eq] = [eq.wrapped.num_xreplace(self.setup) for eq in self._equations_list + self._definitions_list]
        solved = {}
        solved_symbols_set = set(self.setup)
        while True:
            remaining_eqs = []
            for eq in eqs:
                if len(eq.free_symbols) == 1:
                    symbol: Symbol = next(iter(eq.free_symbols))  # type: ignore
                    solns = checked_solveset(eq, self.constraints.get(symbol, symbol))
                    if len(solns) != 1:
                        raise RuntimeError(''.join([
                            f'Failed to find a single solution for "{symbol}" using setup values.  Solutions:\n',
                            f'  {solns}\n',
                            'May need more equations or more constraints, ',
                            'or this problem may not be supported by the current solving algorithm.',
                        ]))
                    soln = next(iter(solns)).n()  # type: ignore
                    if not soln.is_Number:
                        raise RuntimeError(
                            f'Failed to find a numerical solution for "{symbol}" using setup values.  Found "{soln}".'
                        )
                    solved[symbol] = symbol.quantity(soln)  # type: ignore
                    solved_symbols_set.add(symbol)
                else:
                    remaining_eqs.append(eq)
            if len(remaining_eqs) == 0:
                break
            if not (self._all_symbols_set - solved_symbols_set):
                break
            if len(remaining_eqs) == len(eqs):
                raise RuntimeError(''.join([
                    'Cannot solve; may need more equations, ',
                    'or this problem may not be supported by the current solving algorithm\n',
                    f'Solved symbols: {solved_symbols_set}\n',
                    f'Remaining symbols: {self._all_symbols_set - solved_symbols_set}\n',

                    'Remaining equations:',
                    *['\n  * ' + latex(eq) for eq in remaining_eqs],
                ]))
            eqs = [eq.num_xreplace(solved) for eq in remaining_eqs]
        self._solved_from_setup.update(solved)
        self._all_quantities_or_constants.update(solved)


    def _generate_solutions(self):
        '''
        Generate step-by-step solutions for `equations`.

        Once an equation has only one unknown remaining, there are several
        steps to finalize solutions.  First, the equation is solved
        numerically, to make sure that `knowns` contain enough information to
        yield the same numerical solution that was previously found from
        `setup`.  (This is particularly important in cases like quadratics.)
        Then, symbolic solutions are generated, and the solution that yields
        the correct numerical value is selected.
        '''

        solutions: dict[WrappedEq, list[WrappedEq]] = defaultdict(list)
        for eq in self._equations_list:
            solutions[eq].append(eq)

        for eq in self._equations_list:
            current_eq: WrappedEq = solutions[eq][-1]
            # Keep free symbols constant through each iteration, so that there
            # can be multiple independent substitutions but no nested
            # substitutions
            current_free_symbols = current_eq.free_symbols
            for definition_eq in self._definitions_list:
                # Do replacements for zero in separate pass first
                if (definition_eq.wrapped.rhs == 0 and definition_eq.wrapped.lhs in current_free_symbols and
                        definition_eq.wrapped.lhs in current_eq.free_symbols):
                    replacement = definition_eq.wrapped_rhs or definition_eq.wrapped.rhs
                    current_eq = current_eq.xreplace({definition_eq.wrapped.lhs: replacement})  # type: ignore
            if current_eq is not solutions[eq][-1]:
                solutions[eq].append(current_eq)
                current_eq_simplified = wrapped.simplify(current_eq)
                if ((self.simplify and current_eq_simplified.wrapped != current_eq.wrapped) or
                        len(current_eq_simplified.free_symbols) < len(current_eq.free_symbols)):  # type: ignore
                    solutions[eq].append(current_eq_simplified)  # type: ignore
                    current_eq = current_eq_simplified  # type: ignore
            while True:
                current_free_symbols = current_eq.free_symbols
                for definition_eq in self._definitions_list:
                    if (definition_eq.wrapped.lhs in current_free_symbols and
                            definition_eq.wrapped.lhs in current_eq.free_symbols):
                        replacement = definition_eq.wrapped_rhs or definition_eq.wrapped.rhs
                        current_eq = current_eq.xreplace({definition_eq.wrapped.lhs: replacement})  # type: ignore
                if current_eq is solutions[eq][-1]:
                    break
                solutions[eq].append(current_eq)
                current_eq_simplified = wrapped.simplify(current_eq)
                if ((self.simplify and current_eq_simplified.wrapped != current_eq.wrapped) or
                        len(current_eq_simplified.free_symbols) < len(current_eq.free_symbols)):  # type: ignore
                    solutions[eq].append(current_eq_simplified)  # type: ignore
                    current_eq = current_eq_simplified  # type: ignore

        remaining_eqs = self._equations_list.copy()
        solved = self.knowns.copy()
        solved_symbols_set = self._known_symbols_set.copy()
        remaining_unknown_symbols = self._unknown_symbols_set.copy()
        eqs_with_multiple_solns = set()
        while True:
            if not remaining_unknown_symbols or not remaining_eqs:
                break
            for eq in remaining_eqs:
                current_eq = solutions[eq][-1]
                current_eq_remaining_symbols_set = current_eq.free_symbols - solved_symbols_set
                if len(current_eq_remaining_symbols_set) == 0:
                    # May have more equations than necessary.  Only show those
                    # that are needed.
                    del solutions[eq]
                    remaining_eqs = [req for req in remaining_eqs if req is not eq]
                    break
                if len(current_eq_remaining_symbols_set) == 1:
                    symbol: Symbol = next(iter(current_eq_remaining_symbols_set))  # type: ignore
                    solns = checked_solveset(
                        current_eq.wrapped.xreplace(translate_xreplace_rule(solved)),  # type: ignore
                        self.constraints.get(symbol, symbol),  # type: ignore
                    )
                    if len(solns) != 1:
                        if eq not in eqs_with_multiple_solns:
                            eqs_with_multiple_solns.add(eq)
                            continue
                        raise RuntimeError(''.join([
                            f'Failed to find a single solution for "{symbol}" using known values.\n',
                            f'Equation:  {latex(eq)}\n',
                            f'Known values:  {self.knowns}\n',
                            f'Solutions:  {solns}\n',
                            'May need more equations or more constraints, ',
                            'or this problem may not be supported by the current solving algorithm',
                        ]))
                    soln = next(iter(solns)).n()  # type: ignore
                    if not soln.is_Number:
                        raise RuntimeError(
                            f'Failed to find a numerical solution for "{symbol}" using known values.  Found "{soln}".'
                        )
                    if not math.isclose(soln, self._all_quantities_or_constants[symbol].value):
                        raise RuntimeError(
                            f'From setup, expected "{symbol}={self._all_quantities_or_constants[symbol].value}, '
                            f'but found "{soln}" in generating solutions'
                        )
                    symbolic_solns = solveset_for_ans(current_eq, symbol,
                                                      ans=self._all_quantities_or_constants[symbol].value,
                                                      xreplace=solved)  # type: ignore
                    if len(symbolic_solns) != 1:
                        raise RuntimeError(
                            f'Expected a single symbolic solution using known values, but found {solns}')
                    soln_eq = wrapped.Eq(symbol, next(iter(symbolic_solns)), parents=current_eq)
                    if soln_eq.wrapped != current_eq.wrapped:
                        solutions[eq].append(soln_eq)
                    simplified_soln_eq = wrapped.simplify(soln_eq)
                    if ((self.simplify and simplified_soln_eq.wrapped != soln_eq.wrapped) or
                            len(simplified_soln_eq.free_symbols) < len(soln_eq.free_symbols)):  # type: ignore
                        solutions[eq].append(simplified_soln_eq)  # type: ignore
                    solutions[eq].append(
                        wrapped.Eq(symbol, UnevaluatedQuantity(self._all_quantities_or_constants[symbol]))
                    )
                    solved_symbols_set.add(symbol)
                    remaining_unknown_symbols.remove(symbol)
                    remaining_eqs = [req for req in remaining_eqs if req is not eq]
                    # Insertion order into `self.equations` preserves order
                    # in which equations were solved
                    self.equations[eq] = solutions[eq]
                    value = self._all_quantities_or_constants[symbol]
                    if isinstance(value, ConstSymbol):
                        value = value.to_quantity()
                    solved[symbol] = value
                    self.unknowns[symbol] = value
                    eqs_with_multiple_solns.clear()
                    break
            else:
                raise RuntimeError(''.join([
                    'Cannot solve; may need more equations, ',
                    'or this problem may not be supported by the current solving algorithm.\n',
                    f'Known values:  {self.knowns}\n',
                    'Remaining equations:',
                    *['\n  * ' + latex(solutions[eq][-1]) for eq in remaining_eqs],
                ]))


    def simple_solutions(self) -> str:
        solutions = []

        solutions.append('**Known**\n')

        knowns = r'\quad '.join(f'{latex(k)} = {latex(v)}' for k, v in self.knowns.items())
        solutions.append(f'*   ${knowns}$\n')

        solutions.append(f'**Unknown**\n')
        unknowns = r'\quad '.join(latex(k) for k in self.unknowns.keys())
        solutions.append(f'*   ${unknowns}$\n')

        if self.constraints:
            solutions.append(f'**Constraints**\n')
            constraints = r'\quad '.join(f'{latex(k)}: \\text{{{v}}}' for k, v in self.constraints.items())
            solutions.append(f'*   ${constraints}$\n')

        solutions.append('**Solutions**\n')
        if len(self.equations) == 1:
            eq_list = next(iter(self.equations.values()))
            solutions.append(f'*  $\\displaystyle {latex(eq_list[0])}$\n')
            for eq in eq_list[1:-1]:
                solutions.append(f'    $\\displaystyle {latex(eq)}$\n')
            if len(eq_list) >= 2:
                num_sub_rhs = latex(eq_list[-2].rhs, symbol_replace=self._all_quantities_or_constants)
                solutions.append(f'    ${latex(eq_list[-2].lhs)} = {num_sub_rhs}$\n')
                solutions.append(f'    $\\displaystyle {latex(eq_list[-1])}$\n')
        else:
            for n, eq_list in enumerate(self.equations.values(), 1):
                solutions.append(f'{n:>2}. $\\displaystyle {latex(eq_list[0])}$\n')
                for eq in eq_list[1:-1]:
                    solutions.append(f'    $\\displaystyle {latex(eq)}$\n')
                if len(eq_list) >= 2:
                    num_sub_rhs = latex(eq_list[-2].rhs, symbol_replace=self._all_quantities_or_constants)
                    solutions.append(f'    ${latex(eq_list[-2].lhs)} = {num_sub_rhs}$\n')
                    solutions.append(f'    $\\displaystyle {latex(eq_list[-1])}$\n')

        solutions.append('')
        return '\n'.join(solutions)
