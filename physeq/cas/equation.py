# -*- coding: utf-8 -*-
#
# Copyright (c) 2025-2026, Geoffrey M. Poore
# All rights reserved.
#
# Licensed under the BSD 3-Clause License:
# http://opensource.org/licenses/BSD-3-Clause
#


from __future__ import annotations

import math
import sympy
from astropy.units import Quantity
from typing import Literal, Self
from sympy.core.assumptions import check_assumptions as sympy_check_assumptions
from . import exprorder
from .symbol import ConstSymbol, Symbol, translate_numerical_xreplace_rule, translate_xreplace_rule




def solveset_with_checked_assumptions(
    eq: sympy.Eq | exprorder.WrappedEq,
    symbol: sympy.Symbol | exprorder.WrappedExpr,
    domain=sympy.Complexes,
    **assumptions: dict[str, bool | None],
) -> sympy.FiniteSet:
    '''
    `solveset()` wrapper that only returns solutions consistent with `symbol`
    assumptions (for example, `nonnegative` or `negative`).  Additional
    optional assumptions can be provided to further narrow the solutions.
    '''
    if isinstance(eq, exprorder.WrappedEq):
        eq = eq.eq
    elif not isinstance(eq, sympy.Eq):
        raise TypeError
    if isinstance(symbol, exprorder.WrappedExpr):
        if not isinstance(symbol.expr, sympy.Symbol):
            raise TypeError
        symbol = symbol.expr
    elif not isinstance(symbol, sympy.Symbol):
        raise TypeError
    if isinstance(symbol, ConstSymbol):
        raise NotImplementedError

    solution_set = sympy.solveset(eq, symbol, domain=domain)
    if (isinstance(solution_set, sympy.Intersection) and len(solution_set.args) == 2 and
            solution_set.args[0] is domain and isinstance(solution_set.args[1], sympy.FiniteSet)):
        solution_set = solution_set.args[1]
    if not isinstance(solution_set, sympy.FiniteSet):
        raise NotImplementedError(f'Not "sympy.FiniteSet" but {type(solution_set)}:  {solution_set}')
    solution_set = sympy.FiniteSet(
        *(x for x in solution_set if sympy_check_assumptions(x, symbol) is not False),
        evaluate=False,
    )
    if assumptions:
        solution_set = sympy.FiniteSet(
            *(x for x in solution_set if sympy_check_assumptions(x, assumptions) is not False),
            evaluate=False,
        )
    return solution_set


def solveset_for_ans(
    eq: sympy.Eq | exprorder.WrappedEq,
    symbol: sympy.Symbol | exprorder.WrappedExpr,
    domain=sympy.Reals,  # expanded support to include `Complexes`?
    *,
    ans: sympy.Number | sympy.NumberSymbol | int | float,
    xreplace: dict[sympy.Symbol | exprorder.WrappedExpr,
                   sympy.Number | sympy.NumberSymbol | Quantity | ConstSymbol | int | float],
    rel_tol: float | None = None,
    abs_tol: float | None = None,
) -> sympy.FiniteSet:
    '''
    `solveset()` wrapper that returns only the solution(s) that, when modified
    with the given `.xreplace()` rule, yield a specified numerical answer, or
    solutions that cannot be reduced to a number with the `.xreplace()` rule.
    This is intended for determining which symbolic solution(s) yield a known
    numerical answer.
    '''
    if isinstance(eq, exprorder.WrappedEq):
        eq = eq.eq
    elif not isinstance(eq, sympy.Eq):
        raise TypeError
    if isinstance(symbol, exprorder.WrappedExpr):
        if not isinstance(symbol.expr, sympy.Symbol):
            raise TypeError
        symbol = symbol.expr
    elif not isinstance(symbol, sympy.Symbol):
        raise TypeError
    if isinstance(ans, Quantity):
        if isinstance(symbol, Symbol):
            ans = symbol.quantity_value_in_si_coherent_unit(ans)
        else:
            raise TypeError(
                f'"ans" can only be a Quantity when "{symbol}" is a physeq.Symbol with associated units'
            )
    elif not isinstance(ans, (sympy.Number, sympy.NumberSymbol, int, float)):
        raise TypeError
    if not isinstance(xreplace, dict):
        raise TypeError
    translated_xreplace = translate_numerical_xreplace_rule(xreplace)
    # https://docs.python.org/3/library/math.html#math.isclose
    if rel_tol is None:
        rel_tol = 1e-09
    if abs_tol is None:
        abs_tol = 0.0

    solution_set = solveset_with_checked_assumptions(eq, symbol, domain=domain)
    filtered = []
    for soln in solution_set:
        soln_ans = soln.xreplace(translated_xreplace).n()  # type: ignore
        if soln_ans.is_Number:
            if math.isclose(soln_ans, ans, rel_tol=rel_tol, abs_tol=abs_tol):
                filtered.append(soln)
        else:
            # Keep solutions that can't be reduced to a number
            filtered.append(soln)
    return sympy.FiniteSet(*filtered, evaluate=False)




class Eq(sympy.Eq):  # `sympy.Eq` is alias for `sympy.Equality`
    '''
    Subclass of `sympy.Eq` that is compatible with `exprorder.WrappedExpr` and
    `astropy.units.Quantity`.

    By default, `evaluate=False`.
    '''

    def __new__(cls, lhs, rhs, **options):
        options.setdefault('evaluate', False)
        if isinstance(lhs, exprorder.WrappedExpr):
            lhs = lhs.expr
        if isinstance(rhs, exprorder.WrappedExpr):
            rhs = rhs.expr
        obj = super().__new__(cls, lhs, rhs, **options)
        if not isinstance(obj, cls):
            raise NotImplementedError
        return obj


    def xreplace(self, raw_rule: dict, **kwargs) -> Self:
        '''
        `.xreplace()` compatible with `exprorder.WrappedExpr` and
        `astropy.units.Quantity`.
        '''
        return super().xreplace(translate_xreplace_rule(raw_rule, **kwargs))

    def num_xreplace(self, raw_rule: dict, **kwargs) -> Self:
        '''
        `.xreplace()` with purely numerical values compatible with PhysEq
        `ConstSymbol`, `exprorder.WrappedExpr`, or `astropy.units.Quantity`.
        '''
        return super().xreplace(translate_numerical_xreplace_rule(raw_rule, **kwargs))

    def subscript_with_xreplace_rule(
        self, subscript: str | int, style: Literal['normal', 'italic', 'bold'] | None = None
    ) -> tuple[Self, dict[Symbol, Symbol]]:
        '''
        Create a new `Eq` with subscripting for all PhysEq `Symbol` that
        support it, and also return the `xreplace()` rule for performing the
        subscripting replacements.
        '''
        rule = {}
        for s in self.free_symbols:
            if isinstance(s, Symbol) and s.is_subscriptable:
                rule[s] = s.subscript(subscript, style)
        return (self.xreplace(rule), rule)

    def subscript(self, subscript: str | int, style: Literal['normal', 'italic', 'bold'] | None = None) -> Self:
        '''
        Create a new `Eq` with subscripting for all PhysEq `Symbol` that
        support it.
        '''
        return self.subscript_with_xreplace_rule(subscript, style)[0]

    @property
    def nonconst_free_symbols(self) -> set[sympy.Symbol]:
        '''
        `.free_symbols` filtered to exclude PhysEq `ConstSymbol`.
        '''
        return set(x for x in self.free_symbols if not isinstance(x, ConstSymbol))  # type: ignore




class WrappedEq(exprorder.WrappedEq):
    wrapped: Eq
    eq: Eq

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        if not isinstance(self.wrapped, Eq):
            raise TypeError

    def num_xreplace(self, *args, **kwargs) -> Self:
        # This doesn't pass the xreplace rule to the new wrapped equation to
        # generate new term orders.  Order can typically be maintained
        # sufficiently by setting `unevaluated=True`, or by printing with an
        # xreplace rule.
        new_eq = self.wrapped.num_xreplace(*args, **kwargs)
        new_wrapped_eq = type(self)(new_eq, parents=self)
        if self.wrapped_lhs and self.wrapped_rhs:
            new_wrapped_eq.wrapped_lhs = self.wrapped_lhs.wrapper_class(new_eq.lhs, parents=self.wrapped_lhs)
            new_wrapped_eq.wrapped_rhs = self.wrapped_rhs.wrapper_class(new_eq.rhs, parents=self.wrapped_rhs)
        return new_wrapped_eq

    def subscript(self, *args, **kwargs) -> Self:
        new_eq, rule = self.wrapped.subscript_with_xreplace_rule(*args, **kwargs)
        new_wrapped_eq = type(self)(new_eq, parents=self, xreplace_rule=rule)
        if self.wrapped_lhs and self.wrapped_rhs:
            new_wrapped_eq.wrapped_lhs = self.wrapped_lhs.wrapper_class(new_eq.lhs, parents=self.wrapped_lhs,
                                                                        xreplace_rule=rule)
            new_wrapped_eq.wrapped_rhs = self.wrapped_rhs.wrapper_class(new_eq.rhs, parents=self.wrapped_rhs,
                                                                        xreplace_rule=rule)
        return new_wrapped_eq

    @property
    def nonconst_free_symbols(self) -> set[sympy.Symbol]:
        return self.wrapped.nonconst_free_symbols
