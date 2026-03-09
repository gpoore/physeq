# -*- coding: utf-8 -*-
#
# Copyright (c) 2025-2026, Geoffrey M. Poore
# All rights reserved.
#
# Licensed under the BSD 3-Clause License:
# http://opensource.org/licenses/BSD-3-Clause
#


from __future__ import annotations

from astropy.units.si import joule, watt  # type: ignore
from .core import m, M
from .constant import G, g
from .wrapped import Eq, Rational
from .kinematics import d, d_max, y, r, v_0, v, v_max
from .forces import k, Δx_spr, Δx_spr_0, Δx_spr_max, Δy_spr, Δy_spr_0, Δy_spr_max, μ_k
from .wrapped import Symbol


W = Symbol(r'W', 'work', joule, r'W_{<sub>}')
ΣW = sum_W = Symbol(r'\Sigma W', 'sum of work', joule, r'\Sigma W_{<sub>}')
W_nc = W.subscript(r'\text{nc}', description_prefix='sum of non-conservative work', subscriptable=True)
W_k = W.subscript(r'k', description_suffix='by kinetic friction', subscriptable=True)

KE = Symbol(r'\text{KE}', 'kinetic energy', joule, r'\text{KE}_{<sub>}')
KE_0 = KE.subscript(r'0', description_prefix='initial', subscriptable=True)
eq_KE = Eq(KE, Rational(1, 2)*m*v**2)
eq_KE_0 = Eq(KE_0, Rational(1, 2)*m*v_0**2)
ΔKE = delta_KE = Symbol(r'\Delta\text{KE}', 'change in kinetic energy', joule, r'\Delta\text{KE}_{<sub>}')
eq_ΔKE = eq_delta_KE = Eq(ΔKE, KE - KE_0)

PE = Symbol(r'\text{PE}', 'potential energy', joule, r'\text{PE}_{<sub>}')
PE_0 = PE.subscript(r'0', description_prefix='initial', subscriptable=True)
ΔPE = delta_PE = Symbol(r'\Delta\text{PE}', 'change in potential energy', joule, r'\Delta\text{PE}_{<sub>}')
eq_ΔPE = eq_delta_PE = Eq(ΔPE, PE - PE_0)

eq_work_kinetic_energy_theorem = Eq(ΣW, ΔKE)
eq_nonconservative_work = Eq(W_nc, delta_KE + delta_PE)

P = Symbol(r'P', 'power', watt, r'P_{<sub>}')

PE_G = PE.subscript('G', description_prefix='gravitational', subscriptable=True)
eq_PE_G_earth_surface = Eq(PE_G, m*g*y)
eq_PE_G = Eq(PE_G, -G*M*m/r)

PE_spr = PE.subscript(r'\text{spr}', description_prefix='spring', subscriptable=True)
eq_PE_spr = Eq(PE_spr, Rational(1, 2)*k*Δx_spr**2)
eq_E_cons_spring_mass_maxes = Eq(Rational(1, 2)*m*v_max**2, Rational(1, 2)*k*Δx_spr_max**2)
eq_spring_launch_mass_flat_friction = Eq(Rational(1, 2)*k*Δx_spr_max**2, μ_k*m*g*d_max)
