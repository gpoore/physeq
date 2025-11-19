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
from .wrapper import WrappedExpr


pi = WrappedExpr(sympy.pi)
E = WrappedExpr(sympy.E)
I = WrappedExpr(sympy.I)


class Integer(WrappedExpr):
    def __init__(self, num, *args, **kwargs):
        super().__init__(sympy.Integer(num, *args, **kwargs))


class Float(WrappedExpr):
    def __init__(self, num, *args, **kwargs):
        super().__init__(sympy.Float(num, *args, **kwargs))


class Rational(WrappedExpr):
    def __init__(self, p: int, q: int):
        super().__init__(sympy.Rational(p, q))
