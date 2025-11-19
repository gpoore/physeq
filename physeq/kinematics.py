# -*- coding: utf-8 -*-
#
# Copyright (c) 2025, Geoffrey M. Poore
# All rights reserved.
#
# Licensed under the BSD 3-Clause License:
# http://opensource.org/licenses/BSD-3-Clause
#


from __future__ import annotations

from astropy.units.si import meter, radian, second  # type: ignore
from .symbol import Symbol


# distance
d = Symbol(r'd', 'distance', meter, r'd_{<n>}', nonnegative=True)
d_max = d.subscript(r'\text{max}', description_prefix='maximum', subscriptable=True)


# x
x_0 = Symbol(r'x_0', 'initial position in x', meter, r'x_{0\,<n>}')
x = Symbol(r'x', 'position in x', meter, r'x_{<n>}')
Δx = delta_x = Symbol(r'\Delta x', 'displacement in x', meter, r'\Delta x_{<n>}')
v_0x = Symbol(r'v_{0\,x}', 'initial velocity in x', meter/second, r'v_{0\,<n>\,x}')
v_x = Symbol(r'v_x', 'velocity in x', meter/second, r'v_{<n>\,x}')
a_x = Symbol(r'a_x', 'acceleration in x', meter/second**2, r'a_{<n>\,x}')


# y
y_0 = Symbol(r'y_0', 'initial position in y', meter, r'y_{0\,<n>}')
y = Symbol(r'y', 'position in y', meter, r'y_{<n>}')
Δy = delta_y = Symbol(r'\Delta y', 'displacement in y', meter, r'\Delta y_{<n>}')
v_0y = Symbol(r'v_{0\,y}', 'initial velocity in y', meter/second, r'v_{0\,<n>\,y}')
v_y = Symbol(r'v_y', 'velocity in y', meter/second, r'v_{<n>\,y}')
a_y = Symbol(r'a_y', 'acceleration in y', meter/second**2, r'a_{<n>\,y}')


# polar
# r_0 = Symbol('r_0', meter, nonnegative=True)
r = Symbol(r'r', 'distance from origin', meter, r'r_{<n>}', nonnegative=True)
# θ_0 = theta_0 = Symbol('\\theta_0', radian)
θ = theta = Symbol(r'\theta', 'angle', radian, r'\theta_{<n>}')
# v_0r = Symbol('v_{0r}', meter/second)
# v_0θ = v_0theta = Symbol('v_{0\\theta}', meter/second)
# ω_0 = omega_0 = Symbol('\\omega_0', radian/second)
# v_r = Symbol('v_r', meter/second)
# v_theta = Symbol('v_\\theta', meter/second)
# ω = omega = Symbol('\\omega', radian/second)
# a_r = Symbol('a_r', meter/second**2)
# a_θ = a_theta = Symbol('a_\\theta', radian/second**2)
# α = alpha = Symbol('\\alpha', radian/second**2)


# magnitudes
v_0 = Symbol(r'v_0', 'initial speed', meter/second, r'v_{0\,<n>}', nonnegative=True)
v = Symbol(r'v', 'speed', meter/second, r'v_{<n>}', nonnegative=True)
v_min = Symbol(r'v_\text{min}', 'minimum speed', meter/second, r'v_{\text{min}\,<n>}', nonnegative=True)
v_max = Symbol(r'v_\text{max}', 'maximum speed', meter/second, r'v_{\text{max}\,<n>}', nonnegative=True)
a = Symbol(r'a', 'acceleration magnitude', meter/second**2, r'a_{<n>}', nonnegative=True)
v_0_sys = v_0.subscript('sys', style='normal')
v_sys = v.subscript('sys', style='normal')
a_sys = a.subscript('sys', style='normal')


# multiple objects
#### optimize so these can be derived via a `.subscript()` style method?
v_01x = Symbol(r'v_{0\,1\,x}', '1: initial velocity in x', meter/second)
v_1x = Symbol(r'v_{1\,x}', '1: velocity in x', meter/second)
v_02x = Symbol(r'v_{0\,2\,x}', '2: initial velocity in x', meter/second)
v_2x = Symbol(r'v_{2\,x}', '2: velocity in x', meter/second)

# const accel equations
# x

# y
