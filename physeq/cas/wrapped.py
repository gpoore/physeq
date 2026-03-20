# -*- coding: utf-8 -*-
#
# Copyright (c) 2025-2026, Geoffrey M. Poore
# All rights reserved.
#
# Licensed under the BSD 3-Clause License:
# http://opensource.org/licenses/BSD-3-Clause
#


from __future__ import annotations

import sympy
from astropy.units import Quantity
from . import equation, exprorder, symbol
from .symbol import WrappedExpr


# Numerical constants
pi = WrappedExpr(sympy.pi)
E = WrappedExpr(sympy.E)
I = WrappedExpr(sympy.I)


# Numerical and symbolic types
Integer = WrappedExpr.wrap_atomic_expr_class(sympy.Integer)
Float = WrappedExpr.wrap_atomic_expr_class(sympy.Float)
Rational = WrappedExpr.wrap_atomic_expr_class(sympy.Rational)
Symbol = WrappedExpr.wrap_atomic_expr_class(sympy.Symbol)


# Functions
Abs = WrappedExpr.wrap_compound_expr_class(sympy.Abs)
sqrt = WrappedExpr.wrap_expr_callable(sympy.sqrt)

sin = WrappedExpr.wrap_compound_expr_class(sympy.sin)
cos = WrappedExpr.wrap_compound_expr_class(sympy.cos)
tan = WrappedExpr.wrap_compound_expr_class(sympy.tan)
atan = WrappedExpr.wrap_compound_expr_class(sympy.atan)
atan2 = WrappedExpr.wrap_compound_expr_class(sympy.atan2)


# PhysEq
Symbol = symbol.WrappedSymbol.wrap_symbol_class(symbol.Symbol)
ConstSymbol = symbol.WrappedConstSymbol.wrap_const_symbol_class(symbol.ConstSymbol)
Eq = equation.WrappedEq.wrap_eq_class(equation.Eq)




def solveset_with_checked_assumptions(eq: sympy.Eq | exprorder.WrappedEq,
    symbol: sympy.Symbol | exprorder.WrappedExpr,
    domain=sympy.Reals,
    **assumptions: dict[str, bool | None],
) -> sympy.FiniteSet:
    if isinstance(eq, exprorder.WrappedEq):
        parents = eq
    else:
        parents = None
    solution_set = equation.solveset_with_checked_assumptions(eq, symbol, domain, **assumptions)
    return sympy.FiniteSet(
        *(WrappedExpr(x, parents=parents) for x in solution_set),  # type: ignore
        evaluate=False,
    )


def solveset_for_ans(
    eq: sympy.Eq | exprorder.WrappedEq,
    symbol: sympy.Symbol | exprorder.WrappedExpr,
    domain=sympy.Reals,  # expanded support to include `Complexes`?
    *,
    ans: sympy.Number | sympy.NumberSymbol | int | float,
    xreplace: dict[sympy.Symbol | exprorder.WrappedExpr,
                   sympy.Number | sympy.NumberSymbol | Quantity | symbol.ConstSymbol | int | float],
    # https://docs.python.org/3/library/math.html#math.isclose
    rel_tol: float | None = None,
    abs_tol: float | None = None,
) -> sympy.FiniteSet:
    if isinstance(eq, exprorder.WrappedEq):
        parents = eq
    else:
        parents = None
    solution_set = equation.solveset_for_ans(eq, symbol, domain, ans=ans, xreplace=xreplace,
                                             rel_tol=rel_tol, abs_tol=abs_tol)
    return sympy.FiniteSet(
        *(WrappedExpr(x, parents=parents) for x in solution_set),  # type: ignore
        evaluate=False,
    )
