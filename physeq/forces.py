# -*- coding: utf-8 -*-
#
# Copyright (c) 2025-2026, Geoffrey M. Poore
# All rights reserved.
#
# Licensed under the BSD 3-Clause License:
# http://opensource.org/licenses/BSD-3-Clause
#


from __future__ import annotations


from .cas.wrapped import Abs, Eq
from .symbols.constants import G, g
from .symbols.force import ΣF, ΣF_x, ΣF_y, ΣF_sys, F_G, F_s, F_k, μ_s, μ_k, F_N, F_spr, F_spr_x, x_spr, Δx_spr
from .symbols.force import k as k_f
from .symbols.kinematics import a, a_x, a_y, a_sys
from .symbols.mass import m, m_1, m_2, m_sys
from .symbols.space import r


# Newton's 2nd law
eq_newtons_2nd_magnitude = Eq(ΣF, m*a)
eq_newtons_2nd_x = Eq(ΣF_x, m*a_x)
eq_newtons_2nd_y = Eq(ΣF_y, m*a_y)
eq_newtons_2nd_sys_magnitude = Eq(ΣF_sys, m_sys*a_sys)


# Gravity
eq_gravitational_force_magnitude_earth_surface = Eq(F_G, m*g)
eq_gravitational_force_magnitude = Eq(F_G, G*m_1*m_2/r**2)


# Friction
eq_static_friction = Eq(F_s, μ_s*F_N)
eq_kinetic_friction = Eq(F_k, μ_k*F_N)


# Spring
eq_spring_force_x = Eq(F_spr_x, -k_f*x_spr)
eq_spring_force_x_magnitude = Eq(F_spr, k_f*Abs(Δx_spr))
