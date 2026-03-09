# -*- coding: utf-8 -*-
#
# Copyright (c) 2026, Geoffrey M. Poore
# All rights reserved.
#
# Licensed under the BSD 3-Clause License:
# http://opensource.org/licenses/BSD-3-Clause
#


from __future__ import annotations

from astropy.units.si import newton, volt, coulomb  # type: ignore
from .constant import k
from .wrapped import Abs, atan, atan2, Eq, Symbol, sin, cos
from .kinematics import r, r_1, r_2, x, x_1, x_2, y, θ
from .work_energy import PE




Q = Symbol(r'Q', 'Electric source charge', coulomb, r'Q_{<sub>}')
Q_1 = Q.subscript(1, description_suffix=1)
Q_2 = Q.subscript(2, description_suffix=2)

q = Symbol(r'q', 'Electric test charge', coulomb, r'q_{<sub>}')

F_E = Symbol(r'F_E', 'Coulomb force magnitude', newton, r'F_{E\,<sub><i>}', is_vector_magnitude=True)
F_E_x, F_E_y, F_E_z = F_E.cartesian_components()
F_E_1 = F_E.subscript(1, description='Coulomb force 1 magnitude')
F_E_1_x, F_E_1_y, F_E_1_z = F_E_1.cartesian_components()
F_E_2 = F_E.subscript(2, description='Coulomb force 2 magnitude')
F_E_2_x, F_E_2_y, F_E_2_z = F_E_2.cartesian_components()

eq_coulombs_law = eq_F_E = Eq(F_E, k*Abs(q*Q)/r**2)
eq_pythagorean_theorem = Eq(r**2, x**2 + y**2)
#eq_θ_E = Eq(θ, atan2(y, x))
eq_θ_E = Eq(θ, atan(y/x))
eq_F_E_x = Eq(F_E_x, F_E*cos(θ))
eq_F_E_y = Eq(F_E_y, F_E*sin(θ))

E = Symbol(r'E', 'Electric field magnitude', newton/coulomb, r'E_{<sub><i>}', is_vector_magnitude=True)
E_x, E_y, E_z = E.cartesian_components()
eq_E = Eq(E, k*Abs(Q)/r**2)

E_1 = E.subscript(1, description='Electric field 1 magnitude')
eq_E_1 = Eq(E_1, k*Abs(Q_1)/r_1**2)
E_1_x, E_1_y, E_1_z = E_1.cartesian_components()

E_2 = E.subscript(2, description='Electric field 2 magnitude')
eq_E_2 = Eq(E_2, k*Abs(Q_2)/r_2**2)
E_2_x, E_2_y, E_2_z = E_2.cartesian_components()

V = Symbol(r'V', 'Electric potential', volt, r'V_{<sub>}')
PE_E = PE.subscript('E', description_prefix='electric', subscriptable=True)
eq_electric_potential_energy = Eq(PE_E, q*V)
