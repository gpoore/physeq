# -*- coding: utf-8 -*-
#
# Copyright (c) 2025, Geoffrey M. Poore
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
from .symbol import BaseSymbol
from .symbol import Symbol
from .symbol import ConstSymbol




def nsolveset(eq: Eq, symbol: Symbol, domain=sympy.Reals,
              nonnegative: bool | None = None, negative: bool | None = None) -> sympy.FiniteSet:
    '''
    `solveset()` wrapper that enforces PhysEqSymbol `nonnegative` and
    `negative` assumptions for numerical results.
    '''
    if nonnegative is None and negative is None:
        return symbol.values_from_solveset(sympy.solveset(eq, symbol, domain=domain))

    solution_set = symbol.values_from_solveset(sympy.solveset(eq, symbol, domain=domain))
    filtered = []
    if isinstance(nonnegative, bool) and negative is None:
        pass
    elif nonnegative is None and isinstance(negative, bool):
        nonnegative = not negative
    else:
        raise TypeError
    if nonnegative:
        for x in solution_set:
            if getattr(x, 'is_number', False):
                if x >= 0:
                    filtered.append(x)
            else:
                filtered.append(x)
    else:
        for x in solution_set:
            if getattr(x, 'is_number', False):
                if x < 0:
                    filtered.append(x)
            else:
                filtered.append(x)
    return sympy.FiniteSet(*filtered)




def symsolveset(eq: Eq, symbol: Symbol, domain=sympy.Complexes, *,
                ans: int | float | None = None, xreplace: dict[BaseSymbol, int | float] | None = None):
    if ans is None and xreplace is None:
        return sympy.solveset(eq, symbol, domain=domain)
    if not getattr(ans, 'is_number', False) and not isinstance(ans, (int, float)):
        raise TypeError
    if not isinstance(xreplace, dict):
        raise TypeError
    solution_set = sympy.solveset(eq, symbol, domain=domain)
    if len(solution_set) == 1:
        return solution_set
    filtered = []
    for soln in solution_set:
        soln_ans = soln.xreplace(Eq.translate_xreplace_rule(xreplace)).n()
        if getattr(soln_ans, 'is_number', False):
            if math.isclose(soln_ans, ans) or not math.isclose(soln_ans, -ans):
                filtered.append(soln)
        else:
            filtered.append(soln)
    return sympy.FiniteSet(*filtered)




class Eq(sympy.Eq):  # `sympy.Eq` is alias for `sympy.Equality`
    def __new__(cls, lhs, rhs, **options):
        options.setdefault('evaluate', False)
        try:
            lhs.unwrap_expr()
        except AttributeError:
            pass
        try:
            rhs.unwrap_expr()
        except AttributeError:
            pass
        return super().__new__(cls, lhs, rhs, **options)

    @staticmethod
    def translate_xreplace_rule(raw_rule) -> dict[BaseSymbol, int | float]:
        rule = {}
        for k, v in raw_rule.items():
            if isinstance(v, Quantity):
                rule[k] = float(v.value)
            elif isinstance(v, ConstSymbol):
                rule[k] = v.value
            else:
                rule[k] = v
        return rule

    def xreplace(self, raw_rule) -> Self:
        return super().xreplace(self.translate_xreplace_rule(raw_rule))

    def subscript(self, subscript: str | int, style: Literal['normal', 'italic', 'bold'] | None = None) -> Self:
        mapping = {s: s.subscript(subscript, style) if isinstance(s, PhysEqSymbol) else s for s in self.free_symbols}
        return self.xreplace(mapping)

    def __getitem__(self, item: int | str) -> Self:
        return self.subscript(item)

    @property
    def nonconst_free_symbols(self) -> set[PhysEqSymbol]:
        return set(x for x in self.free_symbols if isinstance(x, PhysEqSymbol))

