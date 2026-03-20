# -*- coding: utf-8 -*-
#
# Copyright (c) 2026, Geoffrey M. Poore
# All rights reserved.
#
# Licensed under the BSD 3-Clause License:
# http://opensource.org/licenses/BSD-3-Clause
#


from __future__ import annotations

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from .cas.symbol import BaseSymbol


class Config(object):
    def __init__(self):
        self._kinetic_energy = r'\text{KE}'
        self._kinetic_energy_symbols: set[BaseSymbol] = set()
        self._potential_energy = r'\text{PE}'
        self._potential_energy_symbols: set[BaseSymbol] = set()
        self._spring_constant = 'k'
        self._spring_constant_symbols: set[BaseSymbol] = set()
        self._coulombs_constant = 'k'
        self._coulombs_constant_symbols: set[BaseSymbol] = set()


    def register_symbol(self, symbol: BaseSymbol, field_name: str):
        try:
            registered_symbols: set[BaseSymbol] = getattr(self, f'_{field_name}_symbols')
        except AttributeError:
            raise AttributeError(f'PhysEq config does not have attribute "{field_name}"')
        registered_symbols.add(symbol)

    def _update_registered_symbols(self, field_name: str, symbols: set[BaseSymbol]):
        outdated = None
        for s in symbols:
            if field_name in s.name_template:
                s.rename(s.name_template)
            else:
                if outdated is None:
                    outdated = []
                outdated.append(s)
        if outdated is not None:
            for s in outdated:
                symbols.remove(s)


    @property
    def kinetic_energy(self) -> str:
        return self._kinetic_energy

    @kinetic_energy.setter
    def kinetic_energy(self, value: str):
        if not isinstance(value, str):
            raise TypeError
        self._kinetic_energy = value
        self._update_registered_symbols('kinetic_energy', self._kinetic_energy_symbols)


    @property
    def potential_energy(self) -> str:
        return self._potential_energy

    @potential_energy.setter
    def potential_energy(self, value: str):
        if not isinstance(value, str):
            raise TypeError
        self._potential_energy = value
        self._update_registered_symbols('potential_energy', self._potential_energy_symbols)


    @property
    def spring_constant(self) -> str:
        return self._spring_constant

    @spring_constant.setter
    def spring_constant(self, value: str):
        if not isinstance(value, str):
            raise TypeError
        self._spring_constant = value
        self._update_registered_symbols('spring_constant', self._spring_constant_symbols)


    @property
    def coulombs_constant(self) -> str:
        return self._coulombs_constant

    @coulombs_constant.setter
    def coulombs_constant(self, value: str):
        if not isinstance(value, str):
            raise TypeError
        self._coulombs_constant = value
        self._update_registered_symbols('coulombs_constant', self._coulombs_constant_symbols)
