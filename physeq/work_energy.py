# -*- coding: utf-8 -*-
#
# Copyright (c) 2025-2026, Geoffrey M. Poore
# All rights reserved.
#
# Licensed under the BSD 3-Clause License:
# http://opensource.org/licenses/BSD-3-Clause
#


from __future__ import annotations

from .cas.wrapped import Eq, Rational
from .symbols.constants import G, g
from .symbols.mass import m, M
from .symbols.force import k_spr, Δx_spr, Δx_spr_max, μ_k
from .symbols.kinematics import v, v_0, v_max
from .symbols.space import r, y, d_max
from .symbols.work_energy import KE, KE_0, ΔKE, PE, PE_0, ΔPE, PE_G, PE_spr, ΣW, W_nc


# Kinetic energy
eq_KE = eq_E_k = Eq(KE, Rational(1, 2)*m*v**2)
eq_KE_0 = Eq(KE_0, Rational(1, 2)*m*v_0**2)
eq_ΔKE = eq_delta_KE = Eq(ΔKE, KE - KE_0)


# Potential energy
eq_ΔPE = eq_delta_PE = Eq(ΔPE, PE - PE_0)
eq_PE_G_earth_surface = Eq(PE_G, m*g*y)
eq_PE_G = Eq(PE_G, -G*M*m/r)
eq_PE_spr_Δx = eq_PE_spr_delta_x = Eq(PE_spr, Rational(1, 2)*k_spr*Δx_spr**2)


# Work
eq_work_kinetic_energy_theorem = Eq(ΣW, ΔKE)
eq_nonconservative_work = Eq(W_nc, ΔKE + ΔPE)
