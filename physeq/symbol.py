# -*- coding: utf-8 -*-
#
# Copyright (c) 2025-2026, Geoffrey M. Poore
# All rights reserved.
#
# Licensed under the BSD 3-Clause License:
# http://opensource.org/licenses/BSD-3-Clause
#


from __future__ import annotations

import sympy
from astropy.units import Quantity
from astropy.units.core import UnitBase
from typing import Any, Callable, Literal, Self
from sympy.core.assumptions import check_assumptions as sympy_check_assumptions
from . import exprorder




def translate_xreplace_rule(
    raw_rule: dict[sympy.Symbol | exprorder.WrappedExpr, sympy.Expr | Quantity | int | float],
    check_assumptions: bool = True,
    strict_quantities: bool = True,
) -> dict[sympy.Symbol, sympy.Expr | int | float]:
    '''
    Translate an `.xreplace()` rule that contains `exprorder.WrappedExpr`
    or `astropy.units.Quantity` into a form compatible with SymPy.

    If `check_assumptions`:  check that a symbol's assumptions are
    compatible with its replacement value.

    If `strict_quantities`:  when a replacement value is an
    `astropy.units.Quantity`, require the symbol to be a PhysEq `Symbol`
    with compatible units.
    '''
    rule: dict[sympy.Symbol, sympy.Expr | int | float] = {}
    for k, v in raw_rule.items():
        if isinstance(k, exprorder.WrappedExpr):
            k = k.expr
        if not isinstance(k, sympy.Symbol):
            raise TypeError
        if isinstance(v, exprorder.WrappedExpr):
            v = v.expr
        elif isinstance(v, (sympy.Expr, int, float)):
            pass
        elif isinstance(v, Quantity):
            if isinstance(k, Symbol):
                v = k.quantity_value_in_si_coherent_unit(v)
            elif not strict_quantities:
                v = float(v.value)
            else:
                raise TypeError(
                    f'Symbol "{k}" does not have associated units, '
                    f'so cannot extract value from Quantity "{v}" when "strict_quantities = True"'
                )
        else:
            raise TypeError
        if check_assumptions and sympy_check_assumptions(v, k) is False:
            raise TypeError(f'SymPy assumptions for "{k}" are incompatible with replacement value "{v}"')
        rule[k] = v
    return rule


def translate_numerical_xreplace_rule(
    raw_rule: dict[sympy.Symbol | exprorder.WrappedExpr,
                    sympy.Number | sympy.NumberSymbol | Quantity | ConstSymbol | int | float],
    check_assumptions: bool = True,
    strict_quantities: bool = True,
    evalf_number_symbol: bool = False,
    unevaluated: bool = False,
) -> dict[sympy.Symbol, sympy.Number | int | float]:
    '''
    Translate an `.xreplace()` rule with purely numerical values including
    PhysEq `ConstSymbol`, `exprorder.WrappedExpr`, or
    `astropy.units.Quantity` into a form compatible with SymPy.

    If `check_assumptions`:  check that a symbol's assumptions are
    compatible with its replacement value.

    If `strict_quantities`:  when a replacement value is an
    `astropy.units.Quantity`, require the symbol to be a PhysEq `Symbol`
    with compatible units.
    '''
    rule: dict[sympy.Symbol, sympy.Number | int | float] = {}
    for k, v in raw_rule.items():
        if isinstance(k, exprorder.WrappedExpr):
            k = k.expr
        if not isinstance(k, sympy.Symbol):
            raise TypeError
        if isinstance(v, exprorder.WrappedExpr):
            v = v.expr
        if isinstance(v, (sympy.Number, int, float)):
            pass
        elif isinstance(v, Quantity):
            if isinstance(k, Symbol):
                v = k.quantity_value_in_si_coherent_unit(v)
            elif not strict_quantities:
                v = float(v.value)
            else:
                raise TypeError(
                    f'Symbol "{k}" does not have associated units, '
                    f'so cannot extract value from Quantity "{v}" when "strict_quantities = True"'
                )
        elif isinstance(v, ConstSymbol):
            if k is v:
                v = v.value
            elif isinstance(k, Symbol):
                v = k.quantity_value_in_si_coherent_unit(v.quantity)
            elif not strict_quantities:
                v = v.value
            else:
                raise TypeError(
                    f'Symbol "{k}" does not have associated units, '
                    f'so cannot extract value from ConstSymbol "{v}" when "strict_quantities = True"'
                )
        elif isinstance(v, sympy.NumberSymbol):
            if evalf_number_symbol:
                v = v.n()
        else:
            raise TypeError
        if check_assumptions and sympy_check_assumptions(v, k) is False:
            raise TypeError(f'SymPy assumptions for "{k}" are incompatible with replacement value "{v}"')
        if unevaluated:
            rule[k] = exprorder.OrderUnevaluatedExpr(sympy.sympify(v))
        else:
            rule[k] = v
    return rule




class BaseSymbol(sympy.Symbol):
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
        cls._sanitize(assumptions, cls)  # type: ignore
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
        return (self.name, self.description, self.si_coherent_unit) + self._assumptions0  # type: ignore




class Symbol(BaseSymbol):
    __slots__ = ('is_vector_magnitude', 'is_vector_component', 'subscript_template',)

    def __new__(cls, name: str, description: str, si_coherent_unit: UnitBase, subscript_template: str | None = None, *,
                is_vector_magnitude: bool = False, is_vector_component: bool = False, **assumptions):
        if not isinstance(is_vector_magnitude, bool):
            raise TypeError
        if not isinstance(is_vector_component, bool):
            raise TypeError
        if is_vector_magnitude and is_vector_component:
            raise ValueError
        # https://docs.sympy.org/latest/guides/assumptions.html
        if is_vector_magnitude:
            assumptions['extended_real'] = True
            assumptions['nonnegative'] = True
        elif is_vector_component:
            assumptions['extended_real'] = True
        if subscript_template is not None:
            if not isinstance(subscript_template, str):
                raise TypeError
            if '<sub>' not in subscript_template:
                raise ValueError
        subclass_attr = dict(
            is_vector_magnitude = is_vector_magnitude,
            is_vector_component = is_vector_component,
            subscript_template = subscript_template,
        )
        return cls.__new_inner__(name, description, si_coherent_unit, subclass_attr,
                                 **assumptions)

    is_phys_const = False


    def subscript(self, subscript: str | int, description: str | None = None,
                  description_prefix: str | int | None = None, description_suffix: str | int | None = None,
                  subscriptable: bool = True, style: Literal['normal', 'italic', 'bold'] | None = None) -> Self:
        if self.subscript_template is None:
            raise TypeError('Subscripting is only supported when a subscript template is defined')

        if isinstance(subscript, int):
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

        name = self.subscript_template.replace('<sub>', subscript_formatted).replace('<i>', '')
        if all(x is None for x in (description, description_prefix, description_suffix)):
            description = f'{subscript}: {self.description}'
        elif isinstance(description, str) and all(x is None for x in (description_prefix, description_suffix)):
            pass
        elif description is None and any(isinstance(x, (str, int)) for x in (description_prefix, description_suffix)):
            if isinstance(description_prefix, (str, int)):
                if isinstance(description_suffix, (str, int)):
                    description = f'{description_prefix} {self.description} {description_suffix}'
                elif description_suffix is None:
                    description = f'{description_prefix} {self.description}'
                else:
                    raise TypeError
            elif description_prefix is None and isinstance(description_suffix, (str, int)):
                description = f'{self.description} {description_suffix}'
            else:
                raise TypeError
        else:
            raise TypeError

        if not isinstance(subscriptable, bool):
            raise TypeError
        if subscriptable:
            subscript_template = self.subscript_template.replace('<sub>', rf'{subscript}\,<sub>')
        else:
            subscript_template = None

        return type(self)(name, description, self.si_coherent_unit, subscript_template,
                          is_vector_magnitude=self.is_vector_magnitude, is_vector_component=self.is_vector_component,
                          **self._assumptions_orig)  # type: ignore


    # https://en.wikipedia.org/wiki/Coherence_(units_of_measurement)
    def quantity_value_in_si_coherent_unit(self, quantity: str | Quantity) -> float:
        if isinstance(quantity, str):
            quantity = Quantity(quantity)
        elif not isinstance(quantity, Quantity):
            raise TypeError
        return float(quantity.to(self.si_coherent_unit).value)

    def quantity_from_value_in_si_coherent_unit(self, value: float) -> Quantity:
        value = float(value)
        return Quantity(value, self.si_coherent_unit)

    def compatible_values_from_set(self, solnset: sympy.Set) -> sympy.FiniteSet:
        '''
        Filter a set to include only values that are compatible with SymPy
        assumptions.
        '''
        if not isinstance(solnset, sympy.Set):
            raise TypeError
        if not isinstance(solnset, sympy.FiniteSet):
            raise NotImplementedError
        # Use `.is_number` instead of `.is_Number` to handle `NumberSymbol`
        # objects like `pi`.
        return sympy.FiniteSet(*(x for x in solnset if sympy_check_assumptions(x, self) is not False), evaluate=False)


    _component_description_symbol_map: dict[str, str] = {
        r'\theta': 'θ',
        r'\phi': 'φ',
    }

    def _components(self, coords: tuple[str, str, str]) -> tuple[Self, Self, Self]:
        if not self.is_vector_magnitude:
            raise TypeError('Can only derive vector component symbols from a vector magnitude')
        if not self.subscript_template or '<i>' not in self.subscript_template:
            raise ValueError(
                'Can only derive vector component symbols from a vector magnitude '
                'with a `subscript_template` that includes an <i> field'
            )
        components = []
        for i in coords:
            i_description = self._component_description_symbol_map.get(i, i)
            if self.description.endswith(' magnitude'):
                description = f'{self.description.rsplit(' magnitude', 1)[0]} in {i_description}'
            else:
                description = f'{self.description} in {i_description}'
            components.append(type(self)(
                self.subscript_template.replace('<sub>','').replace('<i>', i),
                description,
                self.si_coherent_unit,
                is_vector_component=True,
                subscript_template=self.subscript_template.replace('<i>', ''),
            ))
        return tuple(components)

    def cartesian_components(self) -> tuple[Self, Self, Self]:
        return self._components(('x', 'y', 'z'))

    def spherical_polar_components(self) -> tuple[Self, Self, Self]:
        return self._components((r'r', r'\theta', r'\phi'))




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
    def value(self) -> float:
        return float(self.quantity.value)

    @property
    def unit(self) -> UnitBase:
        return self.quantity.unit

    def to_expr(self) -> sympy.UnevaluatedExpr:
        return sympy.UnevaluatedExpr(self.value)

    def as_quantity(self) -> Quantity:
        return self.quantity




class WrappedExpr(exprorder.WrappedExpr):
    @staticmethod
    def wrapper_class(*args, **kwargs) -> WrappedExpr:
        return WrappedExpr(*args, **kwargs)

    def xreplace(self, raw_rule: dict, **kwargs) -> Self:
        '''
        `.xreplace()` compatible with `exprorder.WrappedExpr` and
        `astropy.units.Quantity`.
        '''
        return type(self)(self.wrapped.xreplace(translate_xreplace_rule(raw_rule, **kwargs)), parents=self)

    def num_xreplace(self, raw_rule: dict, **kwargs) -> Self:
        '''
        `.xreplace()` with purely numerical values compatible with PhysEq
        `ConstSymbol`, `exprorder.WrappedExpr`, or `astropy.units.Quantity`.
        '''
        return type(self)(self.wrapped.xreplace(translate_numerical_xreplace_rule(raw_rule, **kwargs)), parents=self)




class WrappedSymbol(WrappedExpr):
    wrapped: Symbol
    expr: Symbol

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        if not isinstance(self.wrapped, Symbol):
            raise TypeError

    @classmethod
    def wrap_symbol_class(cls, symbol_class: type[Symbol]) -> Callable[..., Self]:
        if not issubclass(symbol_class, Symbol):
            raise TypeError
        def wrapper(*args, **kwargs):
            return cls(symbol_class(*args, **kwargs))
        return wrapper

    @classmethod
    def wrap_atomic_expr_class(cls, *args, **kwargs):
        raise NotImplementedError

    @classmethod
    def wrap_compound_expr_class(cls, *args, **kwargs):
        raise NotImplementedError

    @staticmethod
    def _unwrap_args(args) -> tuple:
        return tuple(arg.wrapped if isinstance(arg, WrappedExpr) else arg for arg in args)

    def subscript(self, *args, **kwargs) -> Self:
        return type(self)(self.expr.subscript(*args, **kwargs))

    def quantity_value_in_si_coherent_unit(self, *args, **kwargs) -> float:
        args = self._unwrap_args(args)
        return self.expr.quantity_value_in_si_coherent_unit(*args, **kwargs)

    def quantity_from_value_in_si_coherent_unit(self, *args, **kwargs) -> Quantity:
        args = self._unwrap_args(args)
        return self.expr.quantity_from_value_in_si_coherent_unit(*args, **kwargs)

    def compatible_values_from_set(self, *args, **kwargs) -> sympy.FiniteSet:
        args = self._unwrap_args(args)
        return self.expr.compatible_values_from_set(*args, **kwargs)

    def cartesian_components(self) -> tuple[Self, Self, Self]:
        return tuple(type(self)(x) for x in self.expr.cartesian_components())  # type: ignore

    def spherical_polar_components(self) -> tuple[Self, Self, Self]:
        return tuple(type(self)(x) for x in self.expr.spherical_polar_components())  # type: ignore




class WrappedConstSymbol(WrappedExpr):
    wrapped: ConstSymbol
    expr: ConstSymbol

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        if not isinstance(self.wrapped, ConstSymbol):
            raise TypeError

    @classmethod
    def wrap_const_symbol_class(cls, const_symbol_class: type[ConstSymbol]) -> Callable[..., Self]:
        if not issubclass(const_symbol_class, ConstSymbol):
            raise TypeError
        def wrapper(*args, **kwargs):
            return cls(const_symbol_class(*args, **kwargs))
        return wrapper

    @classmethod
    def wrap_atomic_expr_class(cls, *args, **kwargs):
        raise NotImplementedError

    @classmethod
    def wrap_compound_expr_class(cls, *args, **kwargs):
        raise NotImplementedError
