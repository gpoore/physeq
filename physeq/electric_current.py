# -*- coding: utf-8 -*-
#
# Copyright (c) 2026, Geoffrey M. Poore
# All rights reserved.
#
# Licensed under the BSD 3-Clause License:
# http://opensource.org/licenses/BSD-3-Clause
#


from __future__ import annotations

from astropy.units.si import amp, ohm  # type: ignore
from .constant import k
from .wrapped import Abs, atan, atan2, Eq, Symbol, sin, cos
from .kinematics import r, r_1, r_2, x, x_1, x_2, y, θ
from .work_energy import P
from .coulomb_force import V

I = Symbol(r'I', 'Electric current', amp, r'I_{<sub>}')
R = Symbol(r'R', 'Electric resistance', ohm, r'R_{<sub>}')

eq_ohms_law = Eq(V, I*R)
eq_electric_power = Eq(P, I*V)
