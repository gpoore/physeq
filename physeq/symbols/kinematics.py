# -*- coding: utf-8 -*-
#
# Copyright (c) 2025-2026, Geoffrey M. Poore
# All rights reserved.
#
# Licensed under the BSD 3-Clause License:
# http://opensource.org/licenses/BSD-3-Clause
#


from __future__ import annotations

from astropy.units.si import meter, second  # type: ignore
from ..cas.wrapped import Symbol


# speed and velocity
v = Symbol(r'v_{<sub><i>}', 'speed', meter/second)
v_x, v_y, v_z = v.cartesian_components()
v_min = v.subscript('min', style='normal')
v_max = v.subscript('max', style='normal')
v_sys = v.subscript('sys', style='normal')

v_0 = v.subscript(0, description_prefix='initial')
v_0x, v_0y, v_0z = v_0.cartesian_components()
v_0_sys = v_0.subscript('sys', style='normal')


# acceleration
a = Symbol(r'a_{<sub><i>}', 'acceleration', meter/second**2)
a_x, a_y, a_z = a.cartesian_components()
a_sys = a.subscript('sys', style='normal')
