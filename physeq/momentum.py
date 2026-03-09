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
from .core import m_1, m_2
from .wrapped import Eq
from .kinematics import v_x, v_01x, v_1x, v_02x, v_2x
from .wrapped import Symbol

p_01x = Symbol(r'p_{0\,1\,x}', '1: initial momentum in x', meter/second)
p_1x = Symbol(r'p_{1\,x}', '1: momentum in x', meter/second)
p_02x = Symbol(r'p_{0\,2\,x}', '2: initial momentum in x', meter/second)
p_2x = Symbol(r'p_{2\,x}', '2: momentum in x', meter/second)

eq_p_01x = Eq(p_01x, m_1*v_01x)
eq_p_1x = Eq(p_1x, m_1*v_1x)
eq_p_02x = Eq(p_02x, m_2*v_02x)
eq_p_2x = Eq(p_2x, m_2*v_2x)

eq_p_cons_mv_2_objects_collide_bounce_apart = Eq(m_1*v_1x + m_2*v_2x, m_1*v_01x + m_2*v_02x)
eq_p_cons_mv_2_objects_collide_stick_together = Eq((m_1 + m_2)*v_x, m_1*v_01x + m_2*v_02x)
