# -*- coding: utf-8 -*-
#
# Copyright (c) 2025, Geoffrey M. Poore
# All rights reserved.
#
# Licensed under the BSD 3-Clause License:
# http://opensource.org/licenses/BSD-3-Clause
#


from __future__ import annotations

import sympy
from .ordering import TermOrderOperatorMixin


class WrappedExpr(TermOrderOperatorMixin, object):
    '''
    Wrapper for SymPy Expr, typically for numerical objects like `Rational`
    and `Pi`.  The wrapper allows term order to be recorded when the wrapped
    object is involved in a mathematical operation for the first time.  During
    that operation, the wrapper is removed, leaving the wrapped object.  See
    `TermOrderOperatorMixin` for implementation.

    For example, `Wrapper(Pi) + x -> Add(Py, x)`.  During the addition, the
    operation order `(Pi, x)` is recorded.
    '''
    __slots__ = ('_wrapped')

    def __init__(self, expr: sympy.Expr | WrappedExpr):
        if isinstance(expr, WrappedExpr):
            self._wrapped = expr.unwrap_expr()
        elif isinstance(expr, sympy.Expr):
            self._wrapped = expr
        else:
            raise NotImplementedError

    def unwrap_expr(self) -> sympy.Expr:
        return self._wrapped

    @property
    def _op_priority(self):
        return self._wrapped._op_priority
