# -*- coding: utf-8 -*-
#
# Copyright (c) 2025, Geoffrey M. Poore
# All rights reserved.
#
# Licensed under the BSD 3-Clause License:
# http://opensource.org/licenses/BSD-3-Clause
#


from __future__ import annotations

import astropy.constants as ac
from astropy.units import Quantity
from astropy.units.si import meter, second  # type: ignore
from .symbol import ConstSymbol


c = ConstSymbol('c', ac.c.name.lower(), ac.c)
g = ConstSymbol('g', 'gravitational acceleration at surface of Earth', Quantity(9.80, meter/second**2))
G = ConstSymbol('G', ac.G.name.lower(), ac.G)
