# -*- coding: utf-8 -*-
#
# Copyright (c) 2025, Geoffrey M. Poore
# All rights reserved.
#
# Licensed under the BSD 3-Clause License:
# http://opensource.org/licenses/BSD-3-Clause
#


from __future__ import annotations

from astropy.units import dimensionless_unscaled
from astropy.units.si import meter, newton  # type: ignore
from sympy import Abs
from .equation import Eq
from .constant import G, g
from .core import m, m_1, m_2, m_sys
from .kinematics import a_x, a_y, a, a_sys, r, θ
from .symbol import Symbol


# Symbols ####################################################################

# generic forces and sums of forces
F_x = Symbol(r'F_x', 'force in x', newton, r'F_{<n>\,x}')
F_y = Symbol(r'F_y', 'force in y', newton, r'F_{<n>\,y}')
F_r = Symbol(r'F_r', 'force in r', newton, r'F_{<n>\,r}')
F_θ = F_theta = Symbol(r'F_\theta', 'force in theta', newton, r'F_{<n>\,\theta}')
F = Symbol(r'F', 'force magnitude', newton, r'F_{<n>}')
ΣF_x = sum_F_x = Symbol(r'\sum F_x', 'sum of forces in x', newton, r'\sum F_{<n>\,x}')
ΣF_y = sum_F_y = Symbol(r'\sum F_y', 'sum of forces in y', newton, r'\sum F_{<n>\,y}')
ΣF_r = sum_F_r = Symbol(r'\sum F_r', 'sum of forces in r', newton, r'\sum F_{<n>\,r}')
ΣF_θ = sum_F_θ= ΣF_theta = Symbol(r'\sum F_\theta', 'sum of forces in theta', newton, r'\sum F_{<n>\,\theta}')
ΣF = sum_F = Symbol(r'\sum F', 'sum of forces magnitude', newton, r'\sum F_{<n>}')
ΣF_sys_x = ΣF_x.subscript('sys', style='normal')
ΣF_sys_y = ΣF_y.subscript('sys', style='normal')
ΣF_sys = ΣF.subscript('sys', style='normal')
# gravity
F_Gx = F_x.subscript('G', description_prefix='gravitational', subscriptable=True)
F_Gy = F_y.subscript('G', description_prefix='gravitational', subscriptable=True)
F_G = F.subscript('G', description_prefix='gravitational', subscriptable=True)
# normal
F_Nx = F_x.subscript('N', description_prefix='normal', subscriptable=True)
F_Ny = F_y.subscript('N', description_prefix='normal', subscriptable=True)
F_N = F.subscript('N', description_prefix='normal', subscriptable=True)
# tension
F_Tx = F_x.subscript('T', description_prefix='tension', subscriptable=True)
F_Ty = F_y.subscript('T', description_prefix='tension', subscriptable=True)
F_T = F.subscript('T', description_prefix='tension', subscriptable=True)
# kinetic and static friction
F_sx = F_x.subscript('s', description_prefix='tension', subscriptable=True)
F_sy = F_y.subscript('s', description_prefix='tension', subscriptable=True)
F_s = F.subscript('s', description_prefix='tension', subscriptable=True)
F_kx = F_x.subscript('k', description_prefix='tension', subscriptable=True)
F_ky = F_y.subscript('k', description_prefix='tension', subscriptable=True)
F_k = F.subscript('k', description_prefix='tension', subscriptable=True)
# spring
F_spr = F.subscript(r'\text{spr}', description_prefix='spring', subscriptable=True)
F_spr_x = F_spr.subscript(r'\text{spr}\,x', description_suffix='in x', subscriptable=True)
F_spr_y = F_spr.subscript(r'\text{spr}\,y', description_suffix='in y', subscriptable=True)
k = Symbol(r'k', 'spring constant', newton/meter, r'k_{<n>}')
x_spr = Symbol(r'x_\text{spr}', 'spring position in x (x=0 is equilibrium)', meter, r'x_{\text{spr}\,<n>}')
Δx_spr = delta_x_spr = Symbol(
    r'\Delta x_\text{spr}', 'displacement of spring from equilibrium length in x',
    meter, r'\Delta x_{\text{spr}\,<n>}'
)
Δx_spr_0 = delta_x_spr_0 = Δx_spr.subscript(r'0', description_prefix='initial', subscriptable=True)
Δx_spr_max = delta_x_spr_max = Δx_spr.subscript(r'\text{max}', description_prefix='maximum', subscriptable=True)
y_spr = Symbol(r'y_\text{spr}', 'spring position in y (y=0 is equilibrium)', meter, r'y_{\text{spr}\,<n>}')
Δy_spr = delta_y_spr = Symbol(
    r'\Delta y_\text{spr}', 'displacement of spring from equilibrium length in y',
    meter, r'\Delta y_{\text{spr}\,<n>}'
)
Δy_spr_0 = delta_y_spr_0 = Δy_spr.subscript(r'0', description_prefix='initial', subscriptable=True)
Δy_spr_max = delta_y_spr_max = Δy_spr.subscript(r'\text{max}', description_prefix='maximum', subscriptable=True)
# pulling or pushing
F_Px = F_x.subscript('P', description_prefix='tension', subscriptable=True)
F_Py = F_y.subscript('P', description_prefix='tension', subscriptable=True)
F_P = F.subscript('P', description_prefix='tension', subscriptable=True)




# Equations ##################################################################

# Newton's 2nd law
eq_newtons_2nd_x = Eq(ΣF_x, m*a_x)
eq_newtons_2nd_y = Eq(ΣF_y, m*a_y)
eq_newtons_2nd_magnitude = Eq(ΣF, m*a)
eq_newtons_2nd_sys = Eq(ΣF_sys, m_sys*a_sys)
# Gravity
eq_gravitational_force_magnitude_earth_surface = Eq(F_G, m*g)
eq_gravitational_force_magnitude = Eq(F_G, G*m_1*m_2/r**2)
# Friction
μ_s = mu_s = Symbol(r'\mu_s', 'coefficient of static friction', dimensionless_unscaled, r'\mu_{s\,<n>}')
μ_k = mu_k = Symbol(r'\mu_k', 'coefficient of kinetic friction', dimensionless_unscaled, r'\mu_{k\,<n>}')
eq_static_friction = Eq(F_s, μ_s*F_N)
eq_kinetic_friction = Eq(F_k, μ_k*F_N)
# Spring
eq_spring_force_magnitude = Eq(F_spr, k*Abs(Δx_spr))
eq_spring_force_x = Eq(F_spr_x, -k*x_spr)
