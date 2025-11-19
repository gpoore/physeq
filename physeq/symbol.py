# -*- coding: utf-8 -*-
#
# Copyright (c) 2025, Geoffrey M. Poore
# All rights reserved.
#
# Licensed under the BSD 3-Clause License:
# http://opensource.org/licenses/BSD-3-Clause
#


from __future__ import annotations

import sympy
from astropy.units import Quantity
from astropy.units.core import UnitBase
from typing import Any, Literal, Self
from .ordering import TermOrderOperatorMixin




class BaseSymbol(TermOrderOperatorMixin, sympy.Symbol):
    __slots__ = ('description', 'si_coherent_unit')
    _symbol_cache: dict[tuple, Self] = {}
    _symbol_collision_cache: dict[tuple, Self] = {}


    def __new__(cls, name: str, description: str, si_coherent_unit: UnitBase, **assumptions):
        if cls is BaseSymbol:
            raise NotImplementedError
        return cls.__new_inner__(name, description, si_coherent_unit, **assumptions)

    @classmethod
    def __new_inner__(cls, name: str, description: str, si_coherent_unit: UnitBase,
                      subclass_attr: dict[str, Any] | None, **assumptions):
        # https://docs.sympy.org/latest/guides/assumptions.html
        assumptions.setdefault('extended_real', True)
        assumptions_tuple = tuple(sorted(assumptions.items()))
        cache_key = (name, description, si_coherent_unit, assumptions_tuple)
        cls._sanitize(assumptions, cls)
        try:
            return cls._symbol_cache[cache_key]
        except KeyError:
            if not all(isinstance(x, str) for x in (name, description)):
                raise TypeError
            if not isinstance(si_coherent_unit, UnitBase):
                raise TypeError
            collision_key = (name, description)
            if collision_key in cls._symbol_collision_cache:
                raise ValueError(
                    'Cannot define new Symbol that is identical to existing Symbol except for units or assumptions'
                )

            obj = super().__xnew__(cls, name, **assumptions)
            obj.description = description
            obj.si_coherent_unit = si_coherent_unit
            if subclass_attr:
                for k, v in subclass_attr.items():
                    setattr(obj, k, v)

            cls._symbol_cache[cache_key] = obj
            cls._symbol_collision_cache[collision_key] = obj
            return cls._symbol_cache[cache_key]


    is_phys_const: bool


    def _hashable_content(self):
        return (self.name, self.description, self.si_coherent_unit) + self._assumptions0




class Symbol(BaseSymbol):
    __slots__ = ('subscript_template',)

    def __new__(cls, name: str, description: str, si_coherent_unit: UnitBase,
                subscript_template: str | None = None, **assumptions):
        if subscript_template is not None and not isinstance(subscript_template, str):
            raise TypeError
        return cls.__new_inner__(name, description, si_coherent_unit, dict(subscript_template=subscript_template),
                                 **assumptions)

    is_phys_const = False


    def subscript(self, subscript: str | int, description: str | None = None,
                    description_prefix: str | None = None, description_suffix: str | None = None,
                    subscriptable: bool = False, style: Literal['normal', 'italic', 'bold'] | None = None) -> Self:
        if self.subscript_template is None:
            raise TypeError('Subscripting is only supported when a subscript template is defined')

        if isinstance(subscript, int):
            if subscript < 1:
                raise ValueError
            subscript_formatted = str(subscript)
        elif isinstance(subscript, str):
            if not subscript:
                raise ValueError
            if style is None:
                subscript_formatted = subscript
            elif style == 'normal':
                subscript_formatted = rf'\text{{{subscript}}}'
            elif style == 'italic':
                subscript_formatted = rf'\text{{\textit{{{subscript}}}}}'
            elif style == 'bold':
                subscript_formatted = rf'\text{{\textbf{{{subscript}}}}}'
            else:
                raise TypeError
        else:
            raise TypeError

        name = self.subscript_template.replace('<n>', subscript_formatted)
        if not all(x is None or isinstance(x, str) for x in (description, description_prefix, description_suffix)):
            raise TypeError
        if all(x is None for x in (description, description_prefix, description_suffix)):
            description = f'{subscript}: {self.description}'
        elif isinstance(description, str) and all(x is None for x in (description_prefix, description_suffix)):
            pass
        elif description is None and any(isinstance(x, str) for x in (description_prefix, description_suffix)):
            if isinstance(description_prefix, str):
                if isinstance(description_suffix, str):
                    description = f'{description_prefix} {self.description} {description_suffix}'
                else:
                    description = f'{description_prefix} {self.description}'
            elif isinstance(description_suffix, str):
                description = f'{self.description} {description_suffix}'
            else:
                raise TypeError
        else:
            raise TypeError

        if not isinstance(subscriptable, bool):
            raise TypeError
        if subscriptable:
            subscript_template = self.subscript_template.replace('<n>', rf'{subscript}\,<n>')
        else:
            subscript_template = None

        return type(self)(name, description, self.si_coherent_unit, subscript_template, **self._assumptions_orig)


    # https://en.wikipedia.org/wiki/Coherence_(units_of_measurement)
    def quantity_value_in_si_coherent_unit(self, quantity: str | Quantity) -> float:
        if isinstance(quantity, str):
            quantity = Quantity(quantity)
        elif not isinstance(quantity, Quantity):
            raise TypeError
        return float(quantity.to(self.si_coherent_unit).value)

    def quantity_from_value_in_si_coherent_unit(self, value: float) -> Quantity:
        return Quantity(value, self.si_coherent_unit)

    def values_from_solveset(self, solveset: sympy.FiniteSet) -> sympy.FiniteSet:
        # `solveset()` doesn't take `negative` and `nonnegative` into account.
        if self.assumptions0.get('nonnegative'):
            filtered = []
            for x in solveset:
                if getattr(x, 'is_number', False):
                    if x >= 0:
                        filtered.append(x)
                else:
                    filtered.append(x)
            return sympy.FiniteSet(*filtered)
        if self.assumptions0.get('negative'):
            filtered = []
            for x in solveset:
                if getattr(x, 'is_number', False):
                    if x < 0:
                        filtered.append(x)
                else:
                    filtered.append(x)
            return sympy.FiniteSet(*filtered)
        return solveset




class ConstSymbol(BaseSymbol):
    __slots__ = ('quantity',)

    def __new__(cls, name: str, description: str, quantity: str | Quantity, **assumptions):
        if isinstance(quantity, str):
            quantity = Quantity(quantity)
        elif not isinstance(quantity, Quantity):
            raise TypeError
        return cls.__new_inner__(name, description, quantity.unit, dict(quantity=quantity),  # type: ignore
                                 **assumptions)

    is_phys_const = True

    @property
    def value(self):
        return float(self.quantity.value)

    def to_expr(self) -> sympy.UnevaluatedExpr:
        return sympy.UnevaluatedExpr(self.value)

    def as_quantity(self) -> Quantity:
        return self.quantity
