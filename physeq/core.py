# -*- coding: utf-8 -*-
#
# Copyright (c) 2025, Geoffrey M. Poore
# All rights reserved.
#
# Licensed under the BSD 3-Clause License:
# http://opensource.org/licenses/BSD-3-Clause
#


from __future__ import annotations

from astropy.units.si import kilogram, second  # type: ignore
from .symbol import Symbol


# Time
t = Symbol(r't', 'time', second, r't_{<n>}')


# Mass
m = Symbol(r'm', 'mass', kilogram, r'm_{<n>}')
m_1 = m.subscript('1')
m_2 = m.subscript('2')
M = Symbol(r'M', 'mass', kilogram, r'M_{<n>}')
m_sys = m.subscript('sys', style='normal')
