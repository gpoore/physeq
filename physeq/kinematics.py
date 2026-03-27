# -*- coding: utf-8 -*-
#
# Copyright (c) 2026, Geoffrey M. Poore
# All rights reserved.
#
# Licensed under the BSD 3-Clause License:
# http://opensource.org/licenses/BSD-3-Clause
#


from __future__ import annotations

from .symbols.constants import g
from .symbols.time import t_elapsed
from .symbols.space import x, x_0, Δx, y, y_0, Δy
from .symbols.kinematics import v, v_x, v_0x, v_y, v_0y, a, a_x, a_y
from .cas.wrapped import Eq, Rational, sqrt

eq_v_2d = Eq(v, sqrt(v_x**2 + v_y**2))
eq_v_0_2d = eq_v_2d.subscript(0)

eq_a_2d = Eq(a, sqrt(a_x**2 + a_y**2))
eq_a_0_2d = eq_a_2d.subscript(0)


eq_x_const_a = Eq(x, x_0 + v_0x*t_elapsed + Rational(1, 2)*a_x*t_elapsed**2)
eq_Δx_const_a = eq_delta_x_const_a = Eq(Δx, v_0x*t_elapsed + Rational(1, 2)*a_x*t_elapsed**2)
eq_v_x_const_a = Eq(v_x, v_0x + a_x*t_elapsed)
eq_sq_v_x_const_a = Eq(v_x**2, v_0x**2 + 2*a_x*Δx)

eq_y_const_a = Eq(y, y_0 + v_0y*t_elapsed + Rational(1, 2)*a_y*t_elapsed**2)
eq_Δy_const_a = eq_delta_y_const_a = Eq(Δy, v_0y*t_elapsed + Rational(1, 2)*a_y*t_elapsed**2)
eq_v_y_const_a = Eq(v_y, v_0y + a_y*t_elapsed)
eq_sq_v_y_const_a = Eq(v_y**2, v_0y**2 + 2*a_y*Δy)
