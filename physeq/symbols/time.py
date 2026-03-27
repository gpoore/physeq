# -*- coding: utf-8 -*-
#
# Copyright (c) 2026, Geoffrey M. Poore
# All rights reserved.
#
# Licensed under the BSD 3-Clause License:
# http://opensource.org/licenses/BSD-3-Clause
#


from __future__ import annotations

from astropy.units import second
from ..cas.wrapped import Symbol


t = Symbol(r't_{<sub>}', 'time', second)
t_0 = t.subscript(0, description_prefix='initial')
t_elapsed = Symbol(r't_{<sub>}', 'elapsed time', second)
