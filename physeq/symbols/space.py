# -*- coding: utf-8 -*-
#
# Copyright (c) 2025-2026, Geoffrey M. Poore
# All rights reserved.
#
# Licensed under the BSD 3-Clause License:
# http://opensource.org/licenses/BSD-3-Clause
#


from __future__ import annotations

from astropy.units import meter
from ..cas.wrapped import Symbol


# distance
d = Symbol(r'd_{<sub>}', 'distance', meter)
d_min = d.subscript(r'min', description_prefix='min', style='normal')
d_max = d.subscript(r'max', description_prefix='max', style='normal')


# Cartesian
x, y, z = Symbol.cartesian_symbols_from_template('<i>_{<sub>}', 'position in <i>', meter)
x_0, y_0, z_0 = (s.subscript(0, description_prefix='initial') for s in (x, y, z))
Δx, Δy, Δz = delta_x, delta_y, delta_z = Symbol.cartesian_symbols_from_template(r'\Delta <i>_{<sub>}',
                                                                                'displacement in <i>', meter)

# spherical polar
r, θ, φ = _, theta, phi = Symbol.spherical_polar_symbols_from_template('<i>_{<sub>}', 'position in <i>', meter)
r_0, θ_0, φ_0 = _, theta_0, phi_0 = tuple(s.subscript(0, description_prefix='initial') for s in (r, θ, φ))


# polar
r_polar, θ_polar = _, theta_polar = Symbol.polar_symbols_from_template('<i>_{<sub>}', 'position in <i>', meter)
r_0_polar, θ_0_polar = _, theta_0_polar = tuple(s.subscript(0, description_prefix='initial') for s in (r_polar, θ_polar))
