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


def simplify(wrapped_expr: exprorder.BaseWrapped) -> exprorder.BaseWrapped:
    '''
    `sympy.simplify()` wrapper that simplifies wrapped expressions or
    equations and returns wrapped objects that preserve expression order.
    '''
    if isinstance(wrapped_expr, exprorder.WrappedExpr):
        return wrapped_expr.wrapper_class(sympy.simplify(wrapped_expr.wrapped), parents=wrapped_expr)
    if isinstance(wrapped_expr, exprorder.WrappedEq):
        simplified_expr = sympy.simplify(wrapped_expr.wrapped)
        wrapped_simplified_expr = type(wrapped_expr)(simplified_expr, parents=wrapped_expr)
        if wrapped_expr.wrapped_lhs and wrapped_expr.wrapped_rhs:
            wrapped_simplified_expr.wrapped_lhs = wrapped_expr.wrapped_lhs.wrapper_class(
                simplified_expr.lhs, parents=wrapped_expr.wrapped_lhs
            )
            wrapped_simplified_expr.wrapped_rhs = wrapped_expr.wrapped_rhs.wrapper_class(
                simplified_expr.rhs, parents=wrapped_expr.wrapped_rhs
            )
        return wrapped_simplified_expr
    if isinstance(wrapped_expr, exprorder.BaseWrapped):
        raise NotImplementedError
    raise TypeError
