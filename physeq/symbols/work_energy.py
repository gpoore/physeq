# -*- coding: utf-8 -*-
#
# Copyright (c) 2025-2026, Geoffrey M. Poore
# All rights reserved.
#
# Licensed under the BSD 3-Clause License:
# http://opensource.org/licenses/BSD-3-Clause
#


from __future__ import annotations

from astropy.units import joule, watt
from ..cas.wrapped import Symbol


# Work
W = Symbol(r'W_{<sub>}', 'work', joule)
ΣW = sum_W = Symbol(r'\sum W_{<sub>}', 'sum of work', joule)
W_nc = W.subscript(r'nc', description_prefix='non-conservative work', style='normal')
W_k = W.subscript(r'k', description_suffix='by kinetic friction')


# Kinetic energy
KE = E_k = Symbol(r'<config.kinetic_energy>_{<sub>}', 'kinetic energy', joule, nonnegative=True)
KE_0 = E_k_0 = KE.subscript(0, description_prefix='initial')
ΔKE = ΔE_k = delta_KE = delta_E_k = Symbol(r'\Delta <config.kinetic_energy>_{<sub>}',
                                           'change in kinetic energy', joule)


# Potential energy
PE = E_p = Symbol(r'<config.potential_energy>_{<sub>}', 'potential energy', joule)
PE_0 = E_p_0 = PE.subscript(0, description_prefix='initial')
ΔPE = ΔE_p = delta_PE = delta_E_p = Symbol(r'\Delta <config.potential_energy>_{<sub>}',
                                           'change in potential energy', joule)

PE_G = PE.subscript('G', description_prefix='gravitational')
PE_spr = PE.subscript('spr', description_prefix='spring', style='normal')


# Power
P = Symbol(r'P_{<sub>}', 'power', watt)
