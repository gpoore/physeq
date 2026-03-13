# -*- coding: utf-8 -*-
#
# Copyright (c) 2026, Geoffrey M. Poore
# All rights reserved.
#
# Licensed under the BSD 3-Clause License:
# http://opensource.org/licenses/BSD-3-Clause
#


from __future__ import annotations

from astropy.units import kilogram
from .wrapped import Symbol


m = Symbol(r'm_{<sub>}', 'mass', kilogram)
m_1 = m.subscript(1)
m_2 = m.subscript(2)

M = Symbol(r'M_{<sub>}', 'Mass', kilogram)

m_sys = m.subscript('sys', description_prefix='system', style='normal')

