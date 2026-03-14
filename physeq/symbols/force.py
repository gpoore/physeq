# -*- coding: utf-8 -*-
#
# Copyright (c) 2025-2026, Geoffrey M. Poore
# All rights reserved.
#
# Licensed under the BSD 3-Clause License:
# http://opensource.org/licenses/BSD-3-Clause
#


from __future__ import annotations

from astropy.units import dimensionless_unscaled, meter, newton
from ..cas.wrapped import Symbol
from .space import x, y, z, Δx, Δy, Δz


# Generic forces
F = Symbol(r'F_{<sub><i>}', 'force magnitude', newton)
F_x, F_y, F_z = F.cartesian_components()
F_r, F_θ, F_φ = _, F_theta, F_phi =  F.spherical_polar_components()
F_r_polar, F_θ_polar = _, F_theta_polar = F.polar_components()

ΣF = sum_F = Symbol(r'\sum F_{<sub><i>}', 'sum of forces magnitude', newton)
ΣF_x, ΣF_y, ΣF_z = sum_F_x, sum_F_y, sum_F_z = ΣF.cartesian_components()
ΣF_r, ΣF_θ, ΣF_φ = sum_F_r, sum_F_theta, sum_F_phi = ΣF.spherical_polar_components()
ΣF_r_polar, ΣF_θ_polar = sum_F_r_polar, sum_F_theta_polar = ΣF.polar_components()

ΣF_sys = sum_F_sys = ΣF.subscript('sys', style='normal')
ΣF_sys_x, ΣF_sys_y, ΣF_sys_z = sum_F_sys_x, sum_F_sys_y, sum_F_sys_z = ΣF_sys.cartesian_components()
ΣF_sys_r, ΣF_sys_θ, ΣF_sys_φ = sum_F_sys_r, sum_F_sys_theta, sum_F_sys_phi = ΣF_sys.spherical_polar_components()
ΣF_sys_r_polar, ΣF_sys_θ_polar = sum_F_sys_r_polar, sum_F_sys_theta_polar = ΣF_sys.polar_components()


# Gravitational force
F_G = F.subscript('G', description_prefix='gravitational')
F_G_x, F_G_y, F_G_z = F_G.cartesian_components()
F_G_r, F_G_θ, F_G_φ = _, F_G_theta, F_G_phi = F_G.spherical_polar_components()
F_G_r_polar, F_G_θ_polar = _, F_G_theta_polar = F_G.polar_components()


# Normal force
F_N = F.subscript('N', description_prefix='normal')
F_N_x, F_N_y, F_N_z = F_N.cartesian_components()


# Tension force
F_T = F.subscript('T', description_prefix='tension')
F_T_x, F_T_y, F_T_z = F_T.cartesian_components()


# Static and kinetic friction
F_s = F.subscript('s', description_prefix='static friction')
F_s_x, F_s_y, F_s_z = F_s.cartesian_components()
F_k = F.subscript('k', description_prefix='kinetic friction')
F_k_x, F_k_y, F_k_z = F_k.cartesian_components()
μ_s = mu_s = Symbol(r'\mu_{s<sub>}', 'coefficient of static friction', dimensionless_unscaled)
μ_k = mu_k = Symbol(r'\mu_{k<sub>}', 'coefficient of kinetic friction', dimensionless_unscaled)


# Spring
F_spr = F.subscript(r'spr', description_prefix='spring', style='normal')
F_spr_x, F_spr_y, F_spr_z = F_spr.cartesian_components()
F_spr_r, F_spr_θ, F_spr_φ = _, F_spr_theta, F_spr_phi = F_spr.spherical_polar_components()
F_spr_r_polar, F_spr_θ_polar = _, F_spr_theta_polar = F_spr.polar_components()
k = Symbol(r'k', 'spring constant', newton/meter, r'k_{<sub>}')
x_spr, y_spr, z_spr = [
    s.subscript('spr', description_prefix='spring', description_suffix='(origin is equilibrium)', style='normal')
    for s in (x, y, z)
]
Δx_spr, Δy_spr, Δz_spr = delta_x_spr, delta_y_spr, delta_z_spr = [
    s.subscript('spr', description_prefix='spring', description_suffix='from equilibrium', style='normal')
    for s in (x, y, z)
]
Δx_spr_0, Δy_spr_0, Δz_spr_0 = delta_x_spr_0, delta_y_spr_0, delta_z_spr_0 = [
    s.subscript(0, description_prefix='initial', style='normal') for s in (Δx_spr, Δy_spr, Δz_spr)
]
Δx_spr_max, Δy_spr_max, Δz_spr_max = delta_x_spr_max, delta_y_spr_max, delta_z_spr_max = [
    s.subscript('max', description_prefix='max', style='normal') for s in (Δx_spr, Δy_spr, Δz_spr)
]


# Pulling or pushing force
F_P = F.subscript('P', description_prefix='pulling/pushing')
F_P_x, F_P_y, F_P_z = F_P.cartesian_components()
