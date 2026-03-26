# -*- coding: utf-8 -*-
#
# Copyright (c) 2025-2026, Geoffrey M. Poore
# All rights reserved.
#
# Licensed under the BSD 3-Clause License:
# http://opensource.org/licenses/BSD-3-Clause
#


from __future__ import annotations

import random
import re
import sympy
from astropy.units import meter, radian, second, Quantity
from astropy.units.core import UnitBase
from sympy.core.assumptions import check_assumptions as sympy_check_assumptions
from typing import Any, Callable, Literal, Self
from .. import config
from . import exprorder




def translate_xreplace_rule(
    raw_rule: dict[sympy.Symbol | exprorder.WrappedExpr, sympy.Expr | Quantity | int | float],
    check_assumptions: bool = True,
    check_value_constraints: bool = True,
    strict_symbols: bool = True,
    strict_quantities: bool = True,
) -> dict[sympy.Symbol, sympy.Expr | int | float]:
    '''
    Translate an `.xreplace()` rule that contains `exprorder.WrappedExpr`
    or `astropy.units.Quantity` into a form compatible with SymPy.

    If `check_assumptions`:  check that a symbol's assumptions are
    compatible with its replacement value.

    If `check_value_constraints`:  check that a symbol's value constraints
    are satisfied (this goes beyond assumptions).

    If `strict_quantities`:  when a replacement value is an
    `astropy.units.Quantity`, require the symbol to be a PhysEq `Symbol`
    with compatible units.
    '''
    rule: dict[sympy.Symbol, sympy.Expr | int | float] = {}
    for k, v in raw_rule.items():
        if isinstance(k, exprorder.WrappedExpr):
            k = k.expr
        if strict_symbols:
            if not isinstance(k, Symbol) and not isinstance(k, ConstSymbol):
                raise TypeError(
                    'Keys must be instances of physeq.Symbol or physeq.ConstSymbol; '
                    'sympy.Symbol is not accepted when "strict_symbols = True"  (default)'
                )
        elif not isinstance(k, sympy.Symbol):
            raise TypeError
        if isinstance(v, exprorder.WrappedExpr):
            v = v.expr
        if isinstance(v, (sympy.Expr, int, float)):
            pass
        elif isinstance(v, Quantity):
            if isinstance(k, Symbol):
                v = k.quantity_value_in_si_coherent_unit(v)
            elif isinstance(k, ConstSymbol):
                if v != k.to_quantity():
                    raise ValueError(
                        f'Constant symbol "{k}" has value "{k.to_quantity()}", but "{v}" was given'
                    )
                v = float(v.value)
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
        if check_value_constraints and isinstance(k, Symbol) and k.value_constraints is not None:
            if ((isinstance(v, (int, float)) and not k.value_constraints(v)) or
                    (isinstance(v, sympy.Expr) and v.is_number and not k.value_constraints(v.n()))):
                raise TypeError(f'"value_constraints" for "{k}" are incompatible with replacement value "{v}"')
        if isinstance(v, int):
            v = sympy.Integer(v)
        elif isinstance(v, float):
            v = sympy.Float(v)
        rule[k] = v
    return rule


def translate_numerical_xreplace_rule(
    raw_rule: dict[sympy.Symbol | exprorder.WrappedExpr,
                    sympy.Number | sympy.NumberSymbol | Quantity | ConstSymbol | int | float],
    check_assumptions: bool = True,
    check_value_constraints: bool = True,
    strict_symbols: bool = True,
    strict_quantities: bool = True,
    evalf_number_symbol: bool = False,
    unevaluated: bool = False,
) -> dict[sympy.Symbol, sympy.Number | int | float | exprorder.OrderUnevaluatedExpr]:
    '''
    Translate an `.xreplace()` rule with purely numerical values including
    PhysEq `ConstSymbol`, `exprorder.WrappedExpr`, or
    `astropy.units.Quantity` into a form compatible with SymPy.

    If `check_assumptions`:  check that a symbol's assumptions are
    compatible with its replacement value.

    If `check_value_constraints`:  check that a symbol's value constraints
    are satisfied (this goes beyond assumptions).

    If `strict_quantities`:  when a replacement value is an
    `astropy.units.Quantity`, require the symbol to be a PhysEq `Symbol`
    with compatible units.
    '''
    rule: dict[sympy.Symbol, sympy.Number | int | float | exprorder.OrderUnevaluatedExpr] = {}
    for k, v in raw_rule.items():
        if isinstance(k, exprorder.WrappedExpr):
            k = k.expr
        if strict_symbols:
            if not isinstance(k, Symbol) and not isinstance(k, ConstSymbol):
                raise TypeError(
                    'Keys must be instances of physeq.Symbol or physeq.ConstSymbol; '
                    'sympy.Symbol is not accepted when "strict_symbols = True" (default)'
                )
        elif not isinstance(k, sympy.Symbol):
            raise TypeError
        if isinstance(v, exprorder.WrappedExpr):
            v = v.expr
        if isinstance(v, (sympy.Number, int, float)):
            pass
        elif isinstance(v, Quantity):
            if isinstance(k, Symbol):
                v = k.quantity_value_in_si_coherent_unit(v)
            elif isinstance(k, ConstSymbol):
                if v != k.to_quantity():
                    raise ValueError(
                        f'Constant symbol "{k}" has value "{k.to_quantity()}", but "{v}" was given'
                    )
                v = float(v.value)
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
                v = k.quantity_value_in_si_coherent_unit(v.to_quantity())
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
        if check_value_constraints and isinstance(k, Symbol) and k.value_constraints is not None:
            if ((isinstance(v, (int, float)) and not k.value_constraints(v)) or
                    (isinstance(v, sympy.Expr) and v.is_number and not k.value_constraints(v.n()))):
                raise TypeError(f'"value_constraints" for "{k}" are incompatible with replacement value "{v}"')
        if isinstance(v, int):
            v = sympy.Integer(v)
        elif isinstance(v, float):
            v = sympy.Float(v)
        if unevaluated:
            rule[k] = exprorder.OrderUnevaluatedExpr(sympy.sympify(v))
        else:
            rule[k] = v  # type: ignore
    return rule




class BaseSymbol(sympy.Symbol):
    # SymPy hashes/distinguishes symbols using
    # `(self.name,) + self._assumptions0`.  This doesn't work for PhysEq,
    # because it must be possible to modify/customize names.  `description`
    # is used instead of `name` to uniquely identify symbols.
    #
    # Reference:
    #   https://github.com/sympy/sympy/blob/master/sympy/core/symbol.py

    __slots__ = ('description', 'name_template', 'si_coherent_unit')
    description: str
    name_template: str
    si_coherent_unit: UnitBase
    _symbol_cache: dict[tuple, Self] = {}
    _symbol_description_cache: dict[str, Self] = {}


    is_phys_const: bool


    def __new__(cls, name: str, description: str, si_coherent_unit: UnitBase, **assumptions) -> Self:
        raise NotImplementedError
        # <code>
        # return cls.__new_inner__(...)

    @classmethod
    def __new_inner__(cls, name: str, name_template: str, description: str, si_coherent_unit: UnitBase,
                      subclass_attr: dict[str, Any] | None, **assumptions) -> Self:
        # https://docs.sympy.org/latest/guides/assumptions.html
        assumptions.setdefault('extended_real', True)
        assumptions_tuple = tuple(sorted(assumptions.items()))
        cache_key = (description, si_coherent_unit, assumptions_tuple)
        cls._sanitize(assumptions, cls)  # type: ignore
        try:
            obj = cls._symbol_cache[cache_key]
        except KeyError:
            if not all(isinstance(x, str) for x in (name, name_template, description)):
                raise TypeError
            if not isinstance(si_coherent_unit, UnitBase):
                raise TypeError
            if description in cls._symbol_description_cache:
                existing = cls._symbol_description_cache[description]
                raise ValueError(
                    f'Cannot define new Symbol that has same description as existing Symbol "{existing}"'
                )

            obj = super().__xnew__(cls, name, **assumptions)
            obj.description = description
            obj.name_template = name_template
            obj.si_coherent_unit = si_coherent_unit
            if subclass_attr:
                for k, v in subclass_attr.items():
                    setattr(obj, k, v)
            obj._register_config_template_fields(name_template)

            cls._symbol_cache[cache_key] = obj
            cls._symbol_description_cache[description] = obj
            return cls._symbol_cache[cache_key]
        else:
            if obj.name != name or obj.name_template != name_template:
                raise ValueError(
                    'To rename Symbol, use ".rename()"; cannot rename existing, cached object using "__new__()"'
                )
            if subclass_attr:
                for k, v in subclass_attr.items():
                    if getattr(obj, k) != v:
                        raise ValueError(
                            f'Cannot modify "{k}" attribute of existing, cached object using "__new__()"; '
                            f'for existing, cached object "{k}" = "{v}"'
                        )
            return obj


    def _hashable_content(self):
        return (self.description,) + self._assumptions0  # type: ignore


    def rename(self, name: str):
        raise NotImplementedError


    _latex_empty_formatting_patterns = (
        r'\text{}',
        r'\textrm{}', r'\textit{}', r'\textsf{}',
        r'\mathrm{}', r'\mathit{}', r'\mathsf{}',
    )

    _config_template_field_re = re.compile(r'<config\.(.+?)>(.*)')
    _template_field_re = re.compile(r'<(.*?)>')
    _config_subscript_space = r'\,'

    @classmethod
    def _config_template_field_replacer(cls, match: re.Match) -> str:
        try:
            replacement = getattr(config, match.group(1))
            if not isinstance(replacement, str):
                raise TypeError
        except (AttributeError, TypeError):
            raise AttributeError(f'PhysEq config does not have attribute "{match.group(1)}"')
        after = match.group(2)
        if '_' in replacement and after.startswith('_{') and after.endswith('}'):
            replacement_main, replacement_sub = replacement.split('_', 1)
            if replacement_main.count('{') == replacement_main.count('}'):
                return rf'{replacement_main}_{{{replacement_sub}{cls._config_subscript_space}{after[2:-1]}}}'
        return f'{replacement}{after}'

    @classmethod
    def _process_template_fields_in_name(cls, name: str) -> str:
        while True:
            if name.endswith(r'_{}'):
                name = name[:-3]
                continue
            for pattern in cls._latex_empty_formatting_patterns:
                if pattern in name:
                    name = name.replace(pattern, '')
                    break
            else:
                break
        if '<' in name and '>' in name:
            name = cls._config_template_field_re.sub(lambda x: cls._config_template_field_replacer(x), name)
            unknown_match = cls._template_field_re.search(name)
            if unknown_match is not None:
                raise ValueError(f'Symbol name "{name}" contains unknown template field "{unknown_match.group(1)}"')
        return name

    def _register_config_template_fields(self, name_with_template_fields: str):
        for match in self._config_template_field_re.finditer(name_with_template_fields):
            config.register_symbol(self, match.group(1))




class Symbol(BaseSymbol):
    __slots__ = ('is_vector_magnitude', 'is_vector_component', 'parent', '_children_naming_data', 'value_constraints')
    is_vector_magnitude: bool
    is_vector_component: bool
    parent: Self | None
    _children_naming_data: list[tuple[Self, Callable[..., str], tuple]]
    value_constraints: Callable[[int | float], bool] | None


    is_phys_const = False


    def __new__(cls, name: str, description: str, si_coherent_unit: UnitBase, *,
                is_vector_magnitude: bool | None = None, is_vector_component: bool | None = None,
                parent: Self | None = None, value_constraints: Callable[[int | float], bool] | None = None,
                **assumptions) -> Self:
        if not isinstance(name, str):
            raise TypeError
        name_without_template_fields = cls._process_template_fields_in_name(name)

        if is_vector_magnitude is not None and not isinstance(is_vector_magnitude, bool):
            raise TypeError
        if is_vector_component is not None and not isinstance(is_vector_component, bool):
            raise TypeError
        if is_vector_magnitude and is_vector_component:
            raise ValueError

        if cls._subscript_template_field in name:
            if name.count(cls._subscript_template_field) > 1:
                raise ValueError(f'Subscript template field "{cls._subscript_template_field}" can only be used once')
            if name.startswith(cls._subscript_template_field) or name.endswith(cls._subscript_template_field):
                raise ValueError(
                    f'Symbol name cannot start or end with subscript template field "{cls._subscript_template_field}"'
                )
        if cls._component_template_field in name:
            if name.count(cls._component_template_field) > 1:
                raise ValueError(
                    f'Vector component template field "{cls._component_template_field}" can only be used once'
                )
            if name.startswith(cls._component_template_field) or name.endswith(cls._component_template_field):
                raise ValueError(
                    'Symbol name cannot start or end with vector component template field '
                    f'"{cls._component_template_field}"'
                )
            if name.split(cls._component_template_field, 1)[1].rstrip('}'):
                raise ValueError(
                    f'Vector component template field "{cls._component_template_field}" can only be followed '
                    'by right curly braces "}"'
                )
            if is_vector_magnitude is None:
                is_vector_magnitude = True
            elif not is_vector_magnitude:
                raise TypeError(
                    'Only vector magnitudes can have a vector component '
                    f'template field "{cls._component_template_field}"'
                )
            if is_vector_component is None:
                is_vector_component = False
            elif is_vector_component:
                raise TypeError(
                    'Only vector magnitudes can have a vector component '
                    f'template field "{cls._component_template_field}"'
                )
        else:
            if is_vector_magnitude is None:
                is_vector_magnitude = False
            if is_vector_component is None:
                is_vector_magnitude = False

        # https://docs.sympy.org/latest/guides/assumptions.html
        if is_vector_magnitude:
            if assumptions.get('extended_real', True) is not True:
                raise ValueError
            assumptions['extended_real'] = True
            if assumptions.get('nonnegative', True) is not True:
                raise ValueError
            assumptions['nonnegative'] = True
        elif is_vector_component:
            if assumptions.get('extended_real', True) is not True:
                raise ValueError
            assumptions['extended_real'] = True

        if value_constraints is not None and not callable(value_constraints):
            raise TypeError

        subclass_attr = dict(
            is_vector_magnitude = is_vector_magnitude,
            is_vector_component = is_vector_component,
            parent = parent,
            value_constraints = value_constraints,
        )

        return cls.__new_inner__(name_without_template_fields, name, description, si_coherent_unit, subclass_attr,
                                 **assumptions)


    _component_template_field = '<i>'
    _subscript_template_field = '<sub>'
    _latex_math_spaces = (r'\;', r'\:', r'\,', r'\!')
    _subscript_space = r'\,'
    _pre_vector_component_space = r'\,'

    @classmethod
    def _process_template_fields_in_name(cls, name: str) -> str:
        name = name.replace(cls._subscript_template_field, '').replace(cls._component_template_field, '')
        return super()._process_template_fields_in_name(name)


    def _register_config_template_fields(self, name_with_template_fields: str):
        if self.parent is None:
            super()._register_config_template_fields(name_with_template_fields)


    @property
    def is_subscriptable(self) -> bool:
        return self._subscript_template_field in self.name_template

    def subscript(self, subscript: str | int, description: str | None = None,
                  description_prefix: str | int | None = None, description_suffix: str | int | None = None,
                  subscriptable: bool = True, style: Literal['normal', 'italic', 'bold'] | None = None) -> Self:
        if not self.is_subscriptable:
            raise TypeError(
                'Subscripting is only supported for Symbols that are defined with names including '
                f'a template field "{self._subscript_template_field}"'
            )

        if isinstance(subscript, int):
            subscript = str(subscript)
        elif not isinstance(subscript, str):
            raise TypeError
        elif not subscript:
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
            if isinstance(style, str):
                raise ValueError
            raise TypeError

        if not isinstance(subscriptable, bool):
            raise TypeError
        name = self._derive_subscripted_name(self.name_template, subscript_formatted, subscriptable)
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

        obj = type(self)(name, description, self.si_coherent_unit,
                        is_vector_magnitude=self.is_vector_magnitude, is_vector_component=self.is_vector_component,
                        parent=self, value_constraints=self.value_constraints,
                        **self._assumptions_orig)  # type: ignore
        try:
            _children_naming_data = self._children_naming_data
        except AttributeError:
            self._children_naming_data = []
            _children_naming_data = self._children_naming_data
        _children_naming_data.append((obj, self._derive_subscripted_name, (subscript_formatted, subscriptable)))
        return obj

    @classmethod
    def _derive_subscripted_name(cls, name_template: str, subscript_formatted: str, subscriptable: bool) -> str:
        before, after = name_template.split(cls._subscript_template_field)
        if not before.endswith('{') and not before.endswith(cls._latex_math_spaces):
            before += cls._subscript_space
        if (not after.startswith('}') and not after.startswith(cls._latex_math_spaces) and
                not after.startswith(cls._component_template_field)):
            after = cls._subscript_space + after
        if subscriptable:
            name = before + subscript_formatted + cls._subscript_template_field + after
        else:
            name = before + subscript_formatted + after
        return name


    _cartesian_coords = (r'x', r'y', r'z')
    _spherical_polar_coords = (r'r', r'\theta', r'\phi')
    _polar_coords = (r'r', r'\theta')
    _coord_description_symbol_map: dict[str, str] = {
        r'\theta': 'θ',
        r'\phi': 'φ',
    }

    def _components(self, coords: tuple[str, ...],
                    description_template: str | None = None) -> tuple[Self, ...]:
        if description_template is None:
            pass
        elif isinstance(description_template, str):
            if description_template.count(self._component_template_field) != 1:
                raise ValueError
        else:
            raise TypeError
        if not self.is_vector_magnitude:
            raise TypeError('Can only derive vector components from a vector magnitude')
        if not self.name_template or self._component_template_field not in self.name_template:
            raise ValueError(
                'Can only derive vector components from a vector magnitude '
                f'whose name contains a vector component template field "{self._component_template_field}"'
            )
        if self.value_constraints is not None:
            raise NotImplementedError
        components = []
        for coord in coords:
            coord_description = self._coord_description_symbol_map.get(coord, coord)
            if coords == self._spherical_polar_coords:
                coord_description = f'{coord_description} (spherical polar)'
            elif coords == self._polar_coords:
                coord_description = f'{coord_description} (2D polar)'
            if description_template is not None:
                description = description_template.replace(self._component_template_field, coord_description)
            elif self.description.endswith(' magnitude'):
                description = f'{self.description[:-len(" magnitude")]} in {coord_description}'
            elif self.description == 'speed' or self.description.endswith(' speed'):
                description = f'{self.description.rsplit("speed", 1)[0]}velocity in {coord_description}'
            else:
                description = f'{self.description} in {coord_description}'
            name = self._derive_component_name(self.name_template, coord)
            obj = type(self)(name, description, self.si_coherent_unit, is_vector_component=True, parent=self)
            components.append(obj)
            try:
                _children_naming_data = self._children_naming_data
            except AttributeError:
                self._children_naming_data = []
                _children_naming_data = self._children_naming_data
            _children_naming_data.append((obj, self._derive_component_name, (coord,)))
        return tuple(components)

    @classmethod
    def _derive_component_name(cls, name_template: str, coord: str) -> str:
        before, after = name_template.split(cls._component_template_field)
        if (not before.endswith('{') and not any(before.endswith(space) for space in cls._latex_math_spaces)
                and not before.endswith(cls._subscript_template_field)):
            before += cls._subscript_space
        # Don't need to check `after` for space insertion, since it can only
        # be braces `}`; vector component must be final element of name
        return before + coord + after

    def cartesian_components(self, description_template: str | None = None) -> tuple[Self, Self, Self]:
        return self._components(self._cartesian_coords, description_template)  # type: ignore

    def spherical_polar_components(self, description_template: str | None = None) -> tuple[Self, Self, Self]:
        return self._components(self._spherical_polar_coords, description_template)  # type: ignore

    def polar_components(self, description_template: str | None = None) -> tuple[Self, Self]:
        return self._components(self._polar_coords, description_template)  # type: ignore

    @classmethod
    def _coord_symbols_from_template(cls, coords: tuple[str, ...], name_template: str, description_template: str,
                                     unit: UnitBase) -> tuple[Self, ...]:
        if not isinstance(name_template, str) or not isinstance(description_template, str):
            raise TypeError
        if name_template.count(cls._component_template_field) != 1:
            raise ValueError
        if description_template.count(cls._component_template_field) != 1:
            raise ValueError
        if not isinstance(unit, UnitBase):
            raise TypeError
        if unit not in (meter, meter/second, meter/second**2):
            raise NotImplementedError
        coord_symbols = []
        for coord in coords:
            coord_description = cls._coord_description_symbol_map.get(coord, coord)
            if coords == cls._spherical_polar_coords:
                coord_description = f'{coord_description} (spherical polar)'
            elif coords == cls._polar_coords:
                coord_description = f'{coord_description} (2D polar)'
            name = name_template.replace(cls._component_template_field, coord)
            description = description_template.replace(cls._component_template_field, coord_description)
            if coords == cls._cartesian_coords:
                obj = cls(name, description, unit)
            elif coords in (cls._spherical_polar_coords, cls._polar_coords):
                if len(coord_symbols) == 0:
                    obj = cls(name, description, unit, nonnegative=True)
                else:
                    obj = cls(name, description, unit/meter*radian)
            else:
                raise NotImplementedError
            coord_symbols.append(obj)
        return tuple(coord_symbols)

    @classmethod
    def cartesian_symbols_from_template(cls, name_template: str, description_template: str,
                                        unit: UnitBase) -> tuple[Self, Self, Self]:
        return cls._coord_symbols_from_template(cls._cartesian_coords, name_template, description_template,
                                                unit)  # type: ignore

    @classmethod
    def spherical_polar_symbols_from_template(cls, name_template: str, description_template: str,
                                              unit: UnitBase) -> tuple[Self, Self, Self]:
        return cls._coord_symbols_from_template(cls._spherical_polar_coords, name_template, description_template,
                                                unit)  # type: ignore

    @classmethod
    def polar_symbols_from_template(cls, name_template: str, description_template: str,
                                    unit: UnitBase) -> tuple[Self, Self]:
        return cls._coord_symbols_from_template(cls._polar_coords, name_template, description_template,
                                                unit)  # type: ignore


    def rename(self, name: str):
        if self.parent is not None:
            raise TypeError(
                f'Cannot rename a symbol whose name is based on another symbol (name based on "{self.parent}")'
            )
        if self.name_template.count(self._subscript_template_field) != name.count(self._subscript_template_field):
            raise ValueError(
                'Cannot rename a symbol unless the new name and old name either both contain or both do not contain '
                f'the subscript template field "{self._subscript_template_field}"'
            )
        if self.name_template.count(self._component_template_field) != name.count(self._component_template_field):
            raise ValueError(
                'Cannot rename a symbol unless the new name and old name either both contain or both do not contain '
                f'the vector component template field "{self._component_template_field}"'
            )
        self._rename(name)

    def _rename(self, name: str):
        self.name = self._process_template_fields_in_name(name)
        self.name_template = name
        self._register_config_template_fields(name)
        try:
            _children_naming_data = self._children_naming_data
        except AttributeError:
            return
        for symbol, func, other_args in _children_naming_data:
            symbol._rename(func(name, *other_args))


    def _is_value_compatible(self, value: int | float | Quantity | sympy.Expr) -> bool:
        if isinstance(value, sympy.Expr):
            if value.is_number:
                value = float(value.n())
            else:
                return sympy_check_assumptions(value, self) is not False
        elif isinstance(value, Quantity):
            try:
                value.to(self.si_coherent_unit)
            except Exception:
                return False
            value = float(value.value)
        if sympy_check_assumptions(value, self) is False:
            return False
        if self.value_constraints is not None and not self.value_constraints(value):
            return False
        return True

    def _check_value(self, value: int | float | Quantity | sympy.Expr):
        if isinstance(value, sympy.Expr):
            if value.is_number:
                value = float(value.n())
            else:
                if sympy_check_assumptions(value, self) is False:
                    raise ValueError(f'Symbol "{self}" cannot have value "{value}"; this violates symbol assumptions')
                return
        elif isinstance(value, Quantity):
            value.to(self.si_coherent_unit)
            value = float(value.value)
        if sympy_check_assumptions(value, self) is False:
            raise ValueError(f'Symbol "{self}" cannot have value "{value}"; this violates symbol assumptions')
        if self.value_constraints is not None and not self.value_constraints(value):
            raise ValueError(f'Symbol "{self}" cannot have value "{value}"; this violates symbol value constraints')

        # https://en.wikipedia.org/wiki/Coherence_(units_of_measurement)
    def quantity_value_in_si_coherent_unit(self, quantity: str | Quantity) -> float:
        if isinstance(quantity, str):
            quantity = Quantity(quantity)
        elif not isinstance(quantity, Quantity):
            raise TypeError
        value = float(quantity.to(self.si_coherent_unit).value)
        self._check_value(value)
        return value

    def quantity_from_value_in_si_coherent_unit(self, value: float) -> Quantity:
        value = float(value)
        self._check_value(value)
        return Quantity(value, self.si_coherent_unit)

    quantity = quantity_from_value_in_si_coherent_unit

    def compatible_values_from_set(self, solnset: sympy.Set) -> sympy.FiniteSet:
        '''
        Filter a set to include only values that are compatible with SymPy
        assumptions plus value constraints.
        '''
        if not isinstance(solnset, sympy.Set):
            raise TypeError
        if not isinstance(solnset, sympy.FiniteSet):
            raise NotImplementedError
        return sympy.FiniteSet(*(x for x in solnset if self._is_value_compatible(x)), evaluate=False)  # type: ignore

    def randrange_quantity(self, *args, **kwargs) -> Quantity:
        value = random.randrange(*args, **kwargs)
        self._check_value(value)
        return Quantity(value, self.si_coherent_unit)

    def randint_quantity(self, *args) -> Quantity:
        value = random.randint(*args)
        self._check_value(value)
        return Quantity(value, self.si_coherent_unit)

    def choice_quantity(self, seq) -> Quantity:
        value = random.choice(seq)
        self._check_value(value)
        return Quantity(value, self.si_coherent_unit)

    def uniform_quantity(self, *args) -> Quantity:
        value = random.uniform(*args)
        self._check_value(value)
        return Quantity(value, self.si_coherent_unit)




class ConstSymbol(BaseSymbol):
    __slots__ = ('_quantity',)

    is_phys_const = True

    def __new__(cls, name: str, description: str, quantity: str | Quantity, **assumptions) -> Self:
        if isinstance(quantity, str):
            quantity = Quantity(quantity)
        elif not isinstance(quantity, Quantity):
            raise TypeError
        name_without_template_fields = cls._process_template_fields_in_name(name)
        return cls.__new_inner__(name_without_template_fields, name, description, quantity.unit,  # type: ignore
                                 dict(_quantity=quantity), **assumptions)

    @property
    def value(self) -> float:
        return float(self._quantity.value)

    @property
    def unit(self) -> UnitBase:
        return self._quantity.unit

    def to_expr(self) -> sympy.UnevaluatedExpr:
        return sympy.UnevaluatedExpr(self.value)

    def to_quantity(self) -> Quantity:
        return self._quantity

    def rename(self, name: str):
        self.name = self._process_template_fields_in_name(name)
        self.name_template = name
        self._register_config_template_fields(name)





class WrappedExpr(exprorder.WrappedExpr):
    @staticmethod
    def wrapper_class(*args, **kwargs) -> WrappedExpr:
        return WrappedExpr(*args, **kwargs)

    def xreplace(self, raw_rule: dict, **kwargs) -> WrappedExpr:
        '''
        `.xreplace()` compatible with `exprorder.WrappedExpr` and
        `astropy.units.Quantity`.
        '''
        return self.wrapper_class(self.wrapped.xreplace(translate_xreplace_rule(raw_rule, **kwargs)), parents=self)

    def num_xreplace(self, raw_rule: dict, **kwargs) -> WrappedExpr:
        '''
        `.xreplace()` with purely numerical values compatible with PhysEq
        `ConstSymbol`, `exprorder.WrappedExpr`, or `astropy.units.Quantity`.
        '''
        return self.wrapper_class(self.wrapped.xreplace(translate_numerical_xreplace_rule(raw_rule, **kwargs)), parents=self)




class WrappedSymbol(WrappedExpr):
    wrapped: Symbol
    expr: Symbol

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        if not isinstance(self.wrapped, Symbol):
            raise TypeError

    @classmethod
    def wrap_symbol_class(cls, symbol_class: type[Symbol]):
        if not issubclass(symbol_class, Symbol):
            raise TypeError
        class Wrapper(cls):
            def __new__(sub_cls, *args, **kwargs):
                if len(args) == 1 and not kwargs and isinstance(args[0], symbol_class):
                    return cls(args[0])
                return cls(symbol_class(*args, **kwargs))
        return Wrapper

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

    def cartesian_components(self, *args, **kwargs) -> tuple[Self, Self, Self]:
        return tuple(type(self)(x) for x in self.expr.cartesian_components(*args, **kwargs))  # type: ignore

    def spherical_polar_components(self, *args, **kwargs) -> tuple[Self, Self, Self]:
        return tuple(type(self)(x) for x in self.expr.spherical_polar_components(*args, **kwargs))  # type: ignore

    def polar_components(self, *args, **kwargs) -> tuple[Self, Self]:
        return tuple(type(self)(x) for x in self.expr.polar_components(*args, **kwargs))  # type: ignore

    @classmethod
    def cartesian_symbols_from_template(cls, *args, **kwargs) -> tuple[Self, Self, Self]:
        return tuple(cls(x) for x in Symbol.cartesian_symbols_from_template(*args, **kwargs))  # type: ignore

    @classmethod
    def spherical_polar_symbols_from_template(cls, *args, **kwargs) -> tuple[Self, Self, Self]:
        return tuple(cls(x) for x in Symbol.spherical_polar_symbols_from_template(*args, **kwargs))  # type: ignore

    @classmethod
    def polar_symbols_from_template(cls, *args, **kwargs) -> tuple[Self, Self]:
        return tuple(cls(x) for x in Symbol.polar_symbols_from_template(*args, **kwargs))  # type: ignore

    def rename(self, *args, **kwargs):
        self.expr.rename(*args, **kwargs)

    @property
    def name_template(self) -> str:
        return self.expr.name_template

    @property
    def si_coherent_unit(self) -> UnitBase:
        return self.expr.si_coherent_unit

    def quantity_value_in_si_coherent_unit(self, *args, **kwargs) -> float:
        args = self._unwrap_args(args)
        return self.expr.quantity_value_in_si_coherent_unit(*args, **kwargs)

    def quantity_from_value_in_si_coherent_unit(self, *args, **kwargs) -> Quantity:
        args = self._unwrap_args(args)
        return self.expr.quantity_from_value_in_si_coherent_unit(*args, **kwargs)

    quantity = quantity_from_value_in_si_coherent_unit

    def compatible_values_from_set(self, *args, **kwargs) -> sympy.FiniteSet:
        args = self._unwrap_args(args)
        return self.expr.compatible_values_from_set(*args, **kwargs)

    def randrange_quantity(self, *args, **kwargs) -> Quantity:
        return self.expr.randrange_quantity(*args, **kwargs)

    def randint_quantity(self, *args) -> Quantity:
        return self.expr.randint_quantity(*args)

    def choice_quantity(self, seq) -> Quantity:
        return self.expr.choice_quantity(seq)

    def uniform_quantity(self, *args) -> Quantity:
        return self.expr.uniform_quantity(*args)




class WrappedConstSymbol(WrappedExpr):
    wrapped: ConstSymbol
    expr: ConstSymbol

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        if not isinstance(self.wrapped, ConstSymbol):
            raise TypeError

    @classmethod
    def wrap_const_symbol_class(cls, const_symbol_class: type[ConstSymbol]):
        if not issubclass(const_symbol_class, ConstSymbol):
            raise TypeError
        class Wrapper(cls):
            def __new__(sub_cls, *args, **kwargs):
                if len(args) == 1 and not kwargs and isinstance(args[0], const_symbol_class):
                    return cls(args[0])
                return cls(const_symbol_class(*args, **kwargs))
        return Wrapper

    @classmethod
    def wrap_atomic_expr_class(cls, *args, **kwargs):
        raise NotImplementedError

    @classmethod
    def wrap_compound_expr_class(cls, *args, **kwargs):
        raise NotImplementedError

    def rename(self, *args, **kwargs):
        self.expr.rename(*args, **kwargs)

    @property
    def name_template(self) -> str:
        return self.expr.name_template

    @property
    def si_coherent_unit(self) -> UnitBase:
        return self.expr.si_coherent_unit
