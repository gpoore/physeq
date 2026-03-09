# -*- coding: utf-8 -*-
#
# Copyright (c) 2025-2026, Geoffrey M. Poore
# All rights reserved.
#
# Licensed under the BSD 3-Clause License:
# http://opensource.org/licenses/BSD-3-Clause
#


from __future__ import annotations

import astropy.constants as ac
from astropy.units import Quantity
from astropy.units.si import coulomb, meter, newton, second  # type: ignore
from math import pi
from .wrapped import ConstSymbol


c = ConstSymbol(r'c', ac.c.name.lower(), ac.c)
k = ConstSymbol(r'k', 'Coulomb constant', Quantity(1/(4*pi*ac.eps0), newton*meter**2/coulomb**2))
g = ConstSymbol(r'g', 'gravitational acceleration at surface of Earth', Quantity(9.80, meter/second**2))
G = ConstSymbol(r'G', ac.G.name.lower(), ac.G)
