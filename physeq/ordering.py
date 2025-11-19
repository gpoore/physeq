# -*- coding: utf-8 -*-
#
# Copyright (c) 2025, Geoffrey M. Poore
# All rights reserved.
#
# Licensed under the BSD 3-Clause License:
# http://opensource.org/licenses/BSD-3-Clause
#
# ----------------------------------------------------------------------------
#
# Some methods of the `TermOrderOperatorMixin` class are adapted from SymPy
# (https://github.com/sympy/sympy).  SymPy license:
#
#
# Copyright (c) 2006-2023 SymPy Development Team
#
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
#   a. Redistributions of source code must retain the above copyright notice,
#      this list of conditions and the following disclaimer.
#   b. Redistributions in binary form must reproduce the above copyright
#      notice, this list of conditions and the following disclaimer in the
#      documentation and/or other materials provided with the distribution.
#   c. Neither the name of SymPy nor the names of its contributors
#      may be used to endorse or promote products derived from this software
#      without specific prior written permission.
#
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS BE LIABLE FOR
# ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
# OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
# DAMAGE.
#


from __future__ import annotations

import sympy
from sympy.core.decorators import call_highest_priority
from sympy.core.sympify import _sympify, SympifyError
from sympy.core.singleton import S
from typing import Iterable




term_order_slots = ('_first_term', '_last_term')


class TermOrderCompoundMixin(object):
    __slots__ = ()

    def __new__(cls, *args, evaluate: bool = True) -> sympy.Expr:
        obj = super().__new__(cls, *args, evaluate=evaluate)
        cls._set_new_obj_order_attr(obj, args)
        return obj

    @classmethod
    def _set_new_obj_order_attr(cls, obj, args):
        if isinstance(obj, (Add, Mul)):
            left = args[0]
            right = args[-1]
            if isinstance(left, cls):
                obj._first_term = left._first_term
            else:
                obj._first_term = left
            if isinstance(right, cls):
                obj._last_term = right._last_term
            else:
                obj._last_term = right
            if len(args) == 2:
                try:
                    cls._update_term_order(left, right)
                except Exception:
                    pass

    @classmethod
    def _update_term_order(cls, left, right):
        if left.is_Number or right.is_Number:
            return
        order = (left, right)
        key = tuple(sorted(order, key=id))
        try:
            term_order = cls.term_order
            registered_terms = cls.registered_terms
        except AttributeError:
            raise NotImplementedError
        if key in term_order:
            return
        term_order[key] = order
        registered_terms.update(order)

        try:
            left_last_term = left._last_term
        except AttributeError:
            left_last_term = None
        try:
            right_first_term = right._first_term
        except AttributeError:
            right_first_term = None
        if left_last_term is not None:
            cls._update_term_order(left_last_term, right)
        if right_first_term is not None:
            cls._update_term_order(left, right_first_term)

    @classmethod
    def sort_terms(cls, args: Iterable[sympy.Expr]) -> list[sympy.Expr]:
        try:
            term_order = cls.term_order
            registered_terms = cls.registered_terms
        except AttributeError:
            raise NotImplementedError
        not_sortable_terms = []
        sortable_terms = []
        for arg in args:
            if arg in registered_terms:
                sortable_terms.append(arg)
            else:
                print(f'{arg} ({type(arg)})')
                not_sortable_terms.append(arg)
        if not sortable_terms:
            return list(args)
        sorted_term_seqs = []
        while sortable_terms:
            seq = []
            for n, term in enumerate(sortable_terms):
                for other_n, other_term in enumerate(sortable_terms):
                    if n == other_n:
                        continue
                    key = tuple(sorted((term, other_term), key=id))
                    try:
                        order = term_order[key]
                    except KeyError:
                        continue
                    else:
                        seq.extend(order)
                        break
                if seq:
                    break
            if not seq:
                break
            for x in seq:
                sortable_terms.remove(x)
            last_len_sortable_terms = len(sortable_terms)
            while True:
                for n, term in enumerate(sortable_terms):
                    start_key = tuple(sorted((term, seq[0]), key=id))
                    try:
                        start_order = term_order[start_key]
                    except KeyError:
                        pass
                    else:
                        if term is start_order[0]:
                            seq.insert(0, term)
                            sortable_terms.remove(term)
                            break
                    end_key = tuple(sorted((term, seq[-1]), key=id))
                    try:
                        end_order = term_order[end_key]
                    except KeyError:
                        pass
                    else:
                        if term is end_order[1]:
                            sorted.append(term)
                            sortable_terms.remove(term)
                            break
                if len(sortable_terms) == last_len_sortable_terms:
                    break
                last_len_sortable_terms = len(sortable_terms)
            sorted_term_seqs.append(seq)
        return not_sortable_terms + sortable_terms + [x for seq in sorted_term_seqs for x in seq]


class TermOrderOperatorMixin(object):
    __slots__ = ()

    @call_highest_priority('__radd__')
    def __add__(self, other) -> sympy.Expr:
        try:
            self = self.unwrap_expr()
        except AttributeError:
            pass
        try:
            other = other.unwrap_expr()
        except AttributeError:
            try:
                other = _sympify(other)
            except SympifyError:
                return NotImplemented
        return Add(self, other)

    @call_highest_priority('__add__')
    def __radd__(self, other) -> sympy.Expr:
        try:
            self = self.unwrap_expr()
        except AttributeError:
            pass
        try:
            other = other.unwrap_expr()
        except AttributeError:
            try:
                other = _sympify(other)
            except SympifyError:
                return NotImplemented
        return Add(other, self)

    def __neg__(self) -> sympy.Expr:
        try:
            self = self.unwrap_expr()
        except AttributeError:
            pass
        c = self.is_commutative
        return Mul._from_args((S.NegativeOne, self), c)

    @call_highest_priority('__rsub__')
    def __sub__(self, other) -> sympy.Expr:
        try:
            self = self.unwrap_expr()
        except AttributeError:
            pass
        try:
            other = other.unwrap_expr()
        except AttributeError:
            try:
                other = _sympify(other)
            except SympifyError:
                return NotImplemented
        return Add(self, -other)

    @call_highest_priority('__sub__')
    def __rsub__(self, other) -> sympy.Expr:
        try:
            self = self.unwrap_expr()
        except AttributeError:
            pass
        try:
            other = other.unwrap_expr()
        except AttributeError:
            try:
                other = _sympify(other)
            except SympifyError:
                return NotImplemented
        return Add(other, -self)


    @call_highest_priority('__rmul__')
    def __mul__(self, other) -> sympy.Expr:
        try:
            self = self.unwrap_expr()
        except AttributeError:
            pass
        try:
            other = other.unwrap_expr()
        except AttributeError:
            try:
                other = _sympify(other)
            except SympifyError:
                return NotImplemented
        return Mul(self, other)

    @call_highest_priority('__mul__')
    def __rmul__(self, other) -> sympy.Expr:
        try:
            self = self.unwrap_expr()
        except AttributeError:
            pass
        try:
            other = other.unwrap_expr()
        except AttributeError:
            try:
                other = _sympify(other)
            except SympifyError:
                return NotImplemented
        return Mul(other, self)

    @call_highest_priority('__rtruediv__')
    def __truediv__(self, other) -> sympy.Expr:
        try:
            self = self.unwrap_expr()
        except AttributeError:
            pass
        try:
            other = other.unwrap_expr()
        except AttributeError:
            try:
                other = _sympify(other)
            except SympifyError:
                return NotImplemented
        denom = Pow(other, S.NegativeOne)
        if self is S.One:
            return denom
        return Mul(self, denom)

    @call_highest_priority('__truediv__')
    def __rtruediv__(self, other) -> sympy.Expr:
        try:
            self = self.unwrap_expr()
        except AttributeError:
            pass
        try:
            other = other.unwrap_expr()
        except AttributeError:
            try:
                other = _sympify(other)
            except SympifyError:
                return NotImplemented
        denom = Pow(self, S.NegativeOne)
        if other is S.One:
            return denom
        return Mul(other, denom)


    @call_highest_priority('__rpow__')
    def _pow(self, other) -> sympy.Expr:
        try:
            self = self.unwrap_expr()
        except AttributeError:
            pass
        try:
            other = other.unwrap_expr()
        except AttributeError:
            try:
                other = _sympify(other)
            except SympifyError:
                return NotImplemented
        return Pow(self, other)

    @call_highest_priority('__pow__')
    def __rpow__(self, other) -> sympy.Expr:
        try:
            self = self.unwrap_expr()
        except AttributeError:
            pass
        try:
            other = other.unwrap_expr()
        except AttributeError:
            try:
                other = _sympify(other)
            except SympifyError:
                return NotImplemented
        return Pow(other, self)




class Add(TermOrderCompoundMixin, TermOrderOperatorMixin, sympy.Add):
    __slots__ = term_order_slots
    term_order: dict[tuple[int, int], tuple] = {}
    registered_terms: set[sympy.Expr] = set()

    @classmethod
    def _from_args(cls, args, is_commutative=None) -> sympy.Expr:
        obj = super()._from_args(args, is_commutative)
        cls._set_new_obj_order_attr(obj, args)
        return obj


class Mul(TermOrderCompoundMixin, TermOrderOperatorMixin, sympy.Mul):
    __slots__ = term_order_slots
    term_order: dict[tuple[int, int], tuple] = {}
    registered_terms: set[sympy.Expr] = set()

    @classmethod
    def _from_args(cls, args, is_commutative=None) -> sympy.Expr:
        obj = super()._from_args(args, is_commutative)
        cls._set_new_obj_order_attr(obj, args)
        return obj


class Pow(TermOrderCompoundMixin, TermOrderOperatorMixin, sympy.Pow):
    __slots__ = ()
