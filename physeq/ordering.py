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
from sympy.core.sympify import _sympify, SympifyError
from sympy.core.singleton import S
from typing import Any, Callable, Self




class OrderUnevaluatedExpr(sympy.UnevaluatedExpr):
    '''
    Subclass of SymPy's `UnevaluatedExpr` used to wrap numbers and other
    expressions to prevent automatic simplification.  Automatic simplification
    would prevent user-defined expression order from being preserved.
    Typically, operations such as addition and multiplication are performed
    with a raw expression, the result is checked for automatic simplification,
    and then if simplification took place, the operation is repeated using
    `OrderUnevaluatedExpr`.  This ensures that `OrderUnevaluatedExpr` is only
    used when necessary.

    To prevent automatic simplification, it is necessary to wrap a given
    expression in an additional layer of `OrderUnevaluatedExpr` each time it
    recurs.  (An alternative would be to create a separate class for each
    expression, but that would add overhead through class creation.  That
    would also interfere with using `type(expr) is OrderUnevaluatedExpr` to
    efficiently identify expressions that need to be unwrapped during
    expression ordering.)
    '''

    _order_unevaluated_expr_cache: dict[sympy.Expr, Self] = {}

    def __new__(cls, expr: sympy.Expr) -> Self:
        # Every time that a given expression is wrapped, add an additional
        # `OrderUnevaluatedExpr` layer to prevent automatic simplification
        try:
            existing_unevaluated_expr = cls._order_unevaluated_expr_cache[expr]
        except KeyError:
            obj = sympy.Expr.__new__(cls, expr)
        else:
            obj = sympy.Expr.__new__(cls, existing_unevaluated_expr)
        cls._order_unevaluated_expr_cache[expr] = obj
        return cls._order_unevaluated_expr_cache[expr]

    @classmethod
    def from_any(cls, expr: Any) -> Self:
        return cls(_sympify(expr))


    @property
    def raw_expr(self) -> sympy.Expr:
        expr = self.args[0]
        if type(expr) is OrderUnevaluatedExpr:
            while type(expr) is OrderUnevaluatedExpr:
                expr = self.args[0]
        return expr  # type: ignore

    @classmethod
    def unwrap_expr(cls, expr: sympy.Expr) -> sympy.Expr:
        if type(expr) is cls:
            return expr.raw_expr
        return expr

    @classmethod
    def unwrap_exprs(cls, exprs: tuple[sympy.Expr, ...] | list[sympy.Expr]) -> tuple[sympy.Expr, ...]:
        return tuple(expr.raw_expr if type(expr) is cls else expr for expr in exprs)

    @classmethod
    def is_Numberlike(cls, expr: sympy.Expr) -> bool:
        return expr.is_Number or (type(expr) is cls and expr.args[0].is_Number)

    @classmethod
    def is_Numberlike_or_Mul_of_Numberlikes(cls, expr: sympy.Expr) -> bool:
        if type(expr) is cls:
            expr = expr.raw_expr
        if expr.is_Number:
            return True
        if expr.is_Mul and all(cls.is_Numberlike(x) for x in expr.args):  # type: ignore
            return True
        return False

    @classmethod
    def count_factors_Numberlike(cls, factors: tuple[sympy.Expr, ...] | list[sympy.Expr]) -> int:
        return sum(1 for f in factors if cls.is_Numberlike(f))

    @classmethod
    def count_terms_Numberlike(cls, terms: tuple[sympy.Expr, ...] | list[sympy.Expr]) -> int:
        return sum(1 for t in terms if cls.is_Numberlike_or_Mul_of_Numberlikes(t))




_normalize_expr_cache: dict[sympy.Expr, tuple[sympy.Expr, ...]] = {}
def normalize_expr(expr: sympy.Expr) -> tuple[sympy.Expr, ...]:
    '''
    Recursively replace SymPy `Pow` and `Function` (but not `Piecewise`) with
    their arguments.  Also replace `OrderUnevaluatedExpr`.  This provides
    a way of identifying exprs when they are manipulated so that the arguments
    remain the same (or similar), while the operation acting on the arguments
    changes.
    '''
    if expr.is_Atom:
        return (expr,)
    try:
        return _normalize_expr_cache[expr]
    except KeyError:
        if (expr.is_Function and not expr.is_Piecewise):
            normalized = normalize_expr(expr.args[0])  # type: ignore
        elif expr.is_Pow:
            normalized = normalize_expr(expr.args[0]) + normalize_expr(expr.args[1])  # type: ignore
        elif type(expr) is OrderUnevaluatedExpr:
            normalized = normalize_expr(expr.raw_expr)
        else:
            normalized = (expr,)
        _normalize_expr_cache[expr] = normalized
        return _normalize_expr_cache[expr]

_normalize_expr_to_frozenset_cache: dict[sympy.Expr, frozenset[sympy.Expr]] = {}
def normalize_expr_to_frozenset(expr: sympy.Expr) -> frozenset[sympy.Expr]:
    if expr.is_Atom:
        return frozenset((expr,))
    try:
        return _normalize_expr_to_frozenset_cache[expr]
    except KeyError:
        _normalize_expr_to_frozenset_cache[expr] = frozenset(normalize_expr(expr))
        return _normalize_expr_to_frozenset_cache[expr]

_unsorted_expand_normalize_Mul_cache: dict[sympy.Expr, tuple[sympy.Expr, ...]] = {}
def unsorted_expand_normalize_Mul(expr: sympy.Expr) -> tuple[sympy.Expr, ...]:
    '''
    Same as `normalize_expr()` with an additional first step of replacing
    `sympy.Mul` instances with their arguments `expr.args` before performing
    normalization.  This does not use a `OrderCollection` instance to order
    the `Mul` factors; it is intended for cases when order is not important,
    such as in `unsorted_expand_normalize_Mul_to_frozenset()`.
    '''
    try:
        return _unsorted_expand_normalize_Mul_cache[expr]
    except KeyError:
        if not expr.is_Mul:
            raise TypeError
        expanded_normalized = tuple(n_t for arg in expr.args for n_t in normalize_expr(arg))  # type: ignore
        _unsorted_expand_normalize_Mul_cache[expr] = expanded_normalized
        return _unsorted_expand_normalize_Mul_cache[expr]

_unsorted_expand_normalize_Mul_to_frozenset_cache: dict[sympy.Expr, frozenset[sympy.Expr]] = {}
def unsorted_expand_normalize_Mul_to_frozenset(expr: sympy.Expr) -> frozenset[sympy.Expr]:
    try:
        return _unsorted_expand_normalize_Mul_to_frozenset_cache[expr]
    except KeyError:
        expanded_normalized_frozenset = frozenset(unsorted_expand_normalize_Mul(expr))
        _unsorted_expand_normalize_Mul_to_frozenset_cache[expr] = expanded_normalized_frozenset
        return _unsorted_expand_normalize_Mul_to_frozenset_cache[expr]

_normalize_frozenset_cache: dict[frozenset[sympy.Expr], frozenset[sympy.Expr]] = {}
def normalize_frozenset(exprs_frozenset: frozenset[sympy.Expr]) -> frozenset[sympy.Expr]:
    '''
    Create a new `frozenset` containing the normalized expressions from an
    existing `frozenset`.  The new `frozenset` may contain more elements than
    the original, since each `Pow` instance will be replaced with two
    expressions.
    '''
    try:
        return _normalize_frozenset_cache[exprs_frozenset]
    except KeyError:
        normalized_frozenset = frozenset(n_t for t in exprs_frozenset for n_t in normalize_expr(t))
        _normalize_frozenset_cache[exprs_frozenset] = normalized_frozenset
        return _normalize_frozenset_cache[exprs_frozenset]




class BaseExprOrder(object):
    '''
    Expression order for addition or multiplication.  Order is represented in
    two ways, as well as summarized in a set.

      * `sorted_sympy_args: tuple[sympy.Expr, ...]` is the complete, literal,
        correctly ordered sequence of SymPy expressions that makes up an
        addition or multiplication operation.  It is equivalent to `expr.args`
        with correct ordering, except that all instances of
        `OrderUnevaluatedExpr` are replaced with the expressions that they
        wrap.

      * `<term|factor>_tuple` is a more abstract representation of expression
        order designed for detecting similarities with expressions that lack a
        defined order.  For `sympy.Add`, the order is represented as
        `tuple[frozenset[sympy.Expr], ...]`.  Each term in the addition is
        represented by a `frozenset` of the expressions that make up the term
        (or of simplified expressions derived from them).  For `sympy.Mul`,
        the order is represented as `tuple[sympy.Expr, ...]` (where the
        expressions may be simplified expressions derived from those that are
        literally present).

      * `<term|factor>_frozenset: frozenset[sympy.Expr]` summarizes all
        expressions that are present (or simplified expressions derived from
        them).  It is useful in quickly determining whether the expression
        order is relevant for a given unordered expression.

    When an expression is manipulated, its form may change in ways that make
    it difficult to map the manipulated expression to the original expression
    order.  This is addressed with two methods for simplifying expressions
    into a form that is better for detecting similarities:

     1. `normalize`:  Replace certain operations or functions with their
        arguments or expressions derived from their arguments.
          * Replace `sympy.Pow(arg, exp)` with the sequence `(arg, exp)`.
          * Replace `sympy.Function` instances `f(arg)` with `arg` (except
            ignore piecewise functions).
          * Replace `OrderUnevaluatedExpr` with the expression it wraps.
          * Perform normalization recursively, but without breaking the
            current `arg` (or `exp`) into subexpressions.

     2. `flatten`:  Convert a sequence of normalized expressions into a
        sequence of `sympy.Atom` or other expressions that are irreducible by
        `normalize` and `flatten` operations.  This involves breaking `arg`
        (or `exp`) from `normalize` into separate expressions and potentially
        performing additional, recursive `normalize` and `flatten` operations.
        For example, for multiplication, `normalize` converts `f(x*y)` into
        `x*y`, resulting in an expression order `(..., x*y, ...)`, which is
        flattened into `(..., x, y, ...)`.  As part of `flatten`,
        `sympy.Derivative(arg, ...)`, `sympy.Integral(arg, ...)`,
        `sympy.Product(arg, ...)`, and `sympy.Sum(arg, ...)` are replaced by
        `arg`, and then `arg` is flattened, with bound/dummy variables being
        omitted from the resulting sequence.

    Normalized/flattened expressions are combined with the expressions from
    which they are derived.  Continuing the previous `flatten` example for
    function multiplication, normalizing and flattening `f(x*y)` results in
    the final expression order `(..., x, y, x*y, f(x*y), ...)`.  The original
    and intermediate expressions are retained at the end of the sequence of
    `sympy.Atom` so that order is maintained more precisely in cases when
    flattening or normalization are not necessary for a given expression.
    This also helps in cases when a given `sympy.Atom` occurs multiple times
    within an expression order before deduplication.

    `inferred`:  Whether the expression order is literally present or instead
    rather inferred from other relationships.  For example, `sympy.Eq` results
    in the inferred order `eq.lhs - eq.rhs = 0`.  Inferred orders in a
    `ExprOrderCollection` can be overwritten by explicit orders defined later.
    '''

    __slots__ = ('_parent_collection', 'expr', 'inferred', 'sorted_sympy_args')

    # Class attributes set in subclasses
    is_Add_order: bool = False
    is_Mul_order: bool = False
    normalize: bool = False
    flatten: bool = False

    def __init__(self, *args, **kwargs):
        raise NotImplementedError


class BaseAddExprOrder(BaseExprOrder):
    __slots__ = ('terms_frozenset', 'terms_tuple')
    is_Add_order = True

    def __init__(self, expr: sympy.Add, sorted_sympy_args: tuple[sympy.Expr, ...],
                 parent_collection: ExprOrderCollection, inferred: bool = False):
        self.expr = expr
        self.sorted_sympy_args = OrderUnevaluatedExpr.unwrap_exprs(sorted_sympy_args)
        self._parent_collection = parent_collection
        self.inferred = inferred
        self.terms_tuple: tuple[frozenset[sympy.Expr], ...]
        self.terms_frozenset: frozenset[sympy.Expr]
        terms_set = set()
        terms_list = []
        if not self.normalize:
            for arg in self.sorted_sympy_args:
                if arg.is_Mul:
                    # Order doesn't matter, only which terms are present.
                    try:
                        arg_args = parent_collection.Mul_orders[arg].sorted_sympy_args
                    except KeyError:
                        arg_args = OrderUnevaluatedExpr.unwrap_exprs(arg.args)  # type: ignore
                    terms_set.update(arg_args)
                    terms_list.append(frozenset(arg_args))
                else:
                    unwrapped_arg = OrderUnevaluatedExpr.unwrap_expr(arg)
                    terms_set.add(unwrapped_arg)
                    terms_list.append(frozenset((unwrapped_arg,)))
            self.terms_frozenset = frozenset(terms_set)
            self.terms_tuple = tuple(terms_list)
            return
        for arg in self.sorted_sympy_args:
            if arg.is_Mul:
                try:
                    arg_args = parent_collection.Mul_orders[arg].sorted_sympy_args
                except KeyError:
                    arg_args = OrderUnevaluatedExpr.unwrap_exprs(arg.args)  # type: ignore
                normalized_arg_tuple = tuple(n_a for arg_args_n in arg_args for n_a in normalize_expr(arg_args_n))
            else:
                normalized_arg_tuple = normalize_expr(arg)
            normalized_arg_set = set()
            if self.flatten:
                for normalized_arg in normalized_arg_tuple:
                    normalized_arg_set.update(parent_collection.normalize_flatten_expr(normalized_arg))
            normalized_arg_set.update(normalized_arg_tuple)
            normalized_arg_set.add(arg)
            terms_set.update(normalized_arg_set)
            terms_list.append(frozenset(normalized_arg_set))
        self.terms_frozenset = frozenset(terms_set)
        self.terms_tuple = tuple(terms_list)

    def term_to_index(self, term: sympy.Expr, reference_collection: ExprOrderCollection,
                      ignore_terms_shared_subexprs: tuple[sympy.Expr, ...] | None = None,
                      ignore_Number_coeffs: bool = True,
                      skip_literal: bool = False, skip_normalize_without_flatten: bool = False) -> int | float:
        '''
        Return the index of `term` within `.terms_tuple`, or -1 if `term` is
        not found.  When the expression order has `normalize` and `flatten`,
        process `term` with these in searching for an index.

        If `term` is a `Mul` instance, it must be replaced by its arguments in
        searching.

        `ignore_terms_shared_subexprs`:  For a sequence of terms (typically
        the sequence that `term` is a member of), determine which
        subexpressions of these terms are shared by two or more terms.  Ignore
        these shared subexpressions in determining the index for `term`.  This
        can help in finding an order for an expression that has few
        similarities to defined orders, but it can also cause `term` to match
        multiple orders and thereby introduce ambiguity.

        `ignore_Number_coeffs`:  Ignore numbers unless `term` (or the relevant
        normalized/flattened sequence derived from it) consists only of
        numbers.
        '''
        if type(term) is OrderUnevaluatedExpr:
            term = term.raw_expr
        if skip_literal:
            if not self.normalize:
                return -1
        else:
            term_frozenset = reference_collection.term_to_frozenset(
                term, ignore_terms_shared_subexprs, ignore_Number_coeffs,  normalize=False, flatten=False
            )
            if not term_frozenset:
                # Can be empty set due to `ignore_terms_shared_subexprs`
                return -1
            for n, terms_n_frozenset in enumerate(self.terms_tuple):
                if not (term_frozenset - terms_n_frozenset):
                    return n
            else:
                if not self.normalize:
                    return -1
        # Try search with `normalize` but without `flatten`, since that can
        # give more accurate results in cases when flattened term components
        # duplicate unflattened terms.
        if skip_normalize_without_flatten:
            if not self.flatten:
                return -1
        else:
            term_frozenset = reference_collection.term_to_frozenset(
                term, ignore_terms_shared_subexprs, ignore_Number_coeffs, normalize=True, flatten=False
            )
            if not term_frozenset:
                return -1
            for n, terms_n_frozenset in enumerate(self.terms_tuple):
                if not (term_frozenset - terms_n_frozenset):
                    return n
            else:
                if not self.flatten:
                    return -1
        term_frozenset = reference_collection.term_to_frozenset(
                term, ignore_terms_shared_subexprs, ignore_Number_coeffs, normalize=True, flatten=True
            )
        if not term_frozenset:
                return -1
        for n, terms_n_frozenset in enumerate(self.terms_tuple):
            if not (term_frozenset - terms_n_frozenset):
                return n
        else:
            return -1

    def terms_to_indices(self, terms: tuple[sympy.Expr, ...] | list[sympy.Expr],
                         reference_collection: ExprOrderCollection,
                         ignore_shared_subexprs: bool = False, ignore_Number_coeffs: bool = True,
                         allow_repeated_indices: bool = False) -> tuple[tuple[int | float, ...], int, int, int]:
        '''
        For the `Add` operation represented by `terms`, determine the relative
        ordering of `terms` based on `.terms_tuple`.
        '''
        if ignore_shared_subexprs:
            indices = tuple(
                self.term_to_index(term, reference_collection, None, ignore_Number_coeffs) for term in terms
            )
        else:
            terms = tuple(terms)
            indices = tuple(
                self.term_to_index(term, reference_collection, terms, ignore_Number_coeffs) for term in terms
            )
        index_set = set()
        count_repeated = 0
        repeated_index_set = set()
        for index in indices:
            if index == -1:
                continue
            if index in index_set:
                count_repeated += 1
                repeated_index_set.add(index)
            else:
                index_set.add(index)
        if not allow_repeated_indices and repeated_index_set:
            indices = tuple(index if index not in repeated_index_set else -1 for index in indices)
        count_found = 0
        count_found_not_Numberlike = 0
        for term, index in zip(terms, indices):
            if index != -1:
                count_found += 1
                if not OrderUnevaluatedExpr.is_Numberlike_or_Mul_of_Numberlikes(term):
                    count_found_not_Numberlike += 1
        return (indices, count_found, count_found_not_Numberlike, count_repeated)

class AddExprOrder(BaseAddExprOrder):
    __slots__ = ('_normalized', '_normalized_flattened')

    _normalized: NormalizedAddExprOrder
    _normalized_flattened: NormalizedFlattenedAddExprOrder

    def normalize(self):
        try:
            return self._normalized
        except AttributeError:
            self._normalized = NormalizedAddExprOrder(
                self.expr, self.sorted_sympy_args, self._parent_collection, self.inferred
            )
            return self._normalized

    def normalize_flatten(self):
        try:
            return self._normalized_flattened
        except AttributeError:
            self._normalized_flattened = NormalizedFlattenedAddExprOrder(
                self.expr, self.sorted_sympy_args, self._parent_collection, self.inferred
            )
            return self._normalized_flattened


class NormalizedAddExprOrder(BaseAddExprOrder):
    __slots__ = ()
    normalize = True

class NormalizedFlattenedAddExprOrder(BaseAddExprOrder):
    __slots__ = ()
    normalize = True
    flatten = True


class BaseMulExprOrder(BaseExprOrder):
    __slots__ = ('factors_frozenset', 'factors_tuple')
    is_Mul_order = True

    def __init__(self, expr: sympy.Mul, sorted_sympy_args: tuple[sympy.Expr, ...],
                 parent_collection: ExprOrderCollection, inferred: bool = False):
        self.expr = expr
        self.sorted_sympy_args = OrderUnevaluatedExpr.unwrap_exprs(sorted_sympy_args)
        self._parent_collection = parent_collection
        self.inferred = inferred
        self.factors_tuple: tuple[sympy.Expr, ...]
        self.factors_frozenset: frozenset[sympy.Expr]
        if not self.normalize:
            self.factors_frozenset = frozenset(self.sorted_sympy_args)
            self.factors_tuple = self.sorted_sympy_args
            return
        sympy_args_frozenset = frozenset(self.sorted_sympy_args)
        factors_set = set()
        factors_list = []
        for arg in self.sorted_sympy_args:
            for normalized_arg in normalize_expr(arg):
                if self.flatten:
                    for flattened_arg in parent_collection.normalize_flatten_expr(normalized_arg):
                        if flattened_arg not in sympy_args_frozenset and flattened_arg not in factors_set:
                            factors_set.add(flattened_arg)
                            factors_list.append(flattened_arg)
                if normalized_arg not in sympy_args_frozenset and normalized_arg not in factors_set:
                    factors_set.add(normalized_arg)
                    factors_list.append(normalized_arg)
            if arg not in factors_set:
                factors_set.add(arg)
                factors_list.append(arg)
        self.factors_frozenset = frozenset(factors_set)
        self.factors_tuple = tuple(factors_list)

    def factor_to_index(self, factor: sympy.Expr, reference_collection: ExprOrderCollection,
                        ignore_Number_coeffs: bool = True,
                        skip_literal: bool = False, skip_normalize_without_flatten: bool = False) -> int | float:
        '''
        Return the index of `factor` within `.factors_tuple`, or -1 if
        `factor` is not found.  When the expression order has `normalize` and
        `flatten`, process `factor` with these in searching for an index.

        `ignore_Number_coeffs`:  For a normalized/flattened sequence derived
        from a factor, ignore numbers unless the normalized/flattened sequence
        consists only of numbers.
        '''
        if type(factor) is OrderUnevaluatedExpr:
            factor = factor.raw_expr
        if skip_literal:
            if not self.normalize:
                return -1
        else:
            if factor in self.factors_frozenset:
                return self.factors_tuple.index(factor)
            if not self.normalize:
                return -1
        if skip_normalize_without_flatten:
            if not self.flatten:
                return -1
        else:
            normalized_factor_tuple = normalize_expr(factor)
            if ignore_Number_coeffs and any(not n_f.is_Number for n_f in normalized_factor_tuple):
                if all(n_f in self.factors_frozenset or n_f.is_Number for n_f in normalized_factor_tuple):
                    indices = tuple(self.factors_tuple.index(n_f) for n_f in normalized_factor_tuple
                                    if n_f in self.factors_frozenset)
                    return sum(indices) / len(indices)
            elif all(n_f in self.factors_frozenset for n_f in normalized_factor_tuple):
                indices = tuple(self.factors_tuple.index(n_f) for n_f in normalized_factor_tuple)
                return sum(indices) / len(indices)
            if not self.flatten:
                return -1
        normalized_flattened_factor_tuple = reference_collection.normalize_flatten_expr(factor)
        if ignore_Number_coeffs and any(not n_f_f.is_Number for n_f_f in normalized_flattened_factor_tuple):
            if all(n_f_f in self.factors_frozenset or n_f_f.is_Number for n_f_f in normalized_flattened_factor_tuple):
                indices = tuple(self.factors_tuple.index(n_f_f) for n_f_f in normalized_flattened_factor_tuple
                                if n_f_f in self.factors_frozenset)
                return sum(indices) / len(indices)
        elif all(n_f_f in self.factors_frozenset for n_f_f in normalized_flattened_factor_tuple):
            indices = tuple(self.factors_tuple.index(n_f_f) for n_f_f in normalized_flattened_factor_tuple)
            return sum(indices) / len(indices)
        return -1

    def factors_to_indices(self, factors: tuple[sympy.Expr, ...] | list[sympy.Expr],
                           reference_collection: ExprOrderCollection,
                           ignore_Number_coeffs: bool = True,
                           allow_repeated_indices: bool = False) -> tuple[tuple[int | float, ...], int, int, int]:
        '''
        For the `Mul` operation represented by `factors`, determine the
        relative ordering of `factors` based on `.factors_tuple`.
        '''
        indices = tuple(self.factor_to_index(factor, reference_collection, ignore_Number_coeffs) for factor in factors)
        index_set = set()
        count_repeated = 0
        repeated_index_set = set()
        for index in indices:
            if index == -1:
                continue
            if index in index_set:
                count_repeated += 1
                repeated_index_set.add(index)
            else:
                index_set.add(index)
        if not allow_repeated_indices and repeated_index_set:
            indices = tuple(index if index not in repeated_index_set else -1 for index in indices)
        count_found = 0
        count_found_not_Numberlike = 0
        for factor, index in zip(factors, indices):
            if index != -1:
                count_found += 1
                if not OrderUnevaluatedExpr.is_Numberlike(factor):
                    count_found_not_Numberlike += 1
        return (indices, count_found, count_found_not_Numberlike, count_repeated)

class MulExprOrder(BaseMulExprOrder):
    __slots__ = ('_normalized', '_normalized_flattened')

    _normalized: NormalizedMulExprOrder
    _normalized_flattened: NormalizedFlattenedMulExprOrder

    def normalize(self):
        try:
            return self._normalized
        except AttributeError:
            self._normalized = NormalizedMulExprOrder(
                self.expr, self.sorted_sympy_args, self._parent_collection, self.inferred
            )
            return self._normalized

    def normalize_flatten(self):
        try:
            return self._normalized_flattened
        except AttributeError:
            self._normalized_flattened = NormalizedFlattenedMulExprOrder(
                self.expr, self.sorted_sympy_args, self._parent_collection, self.inferred
            )
            return self._normalized_flattened

class NormalizedMulExprOrder(BaseMulExprOrder):
    __slots__ = ()
    normalize = True

class NormalizedFlattenedMulExprOrder(BaseMulExprOrder):
    __slots__ = ()
    normalize = True
    flatten = True




class ExprOrderCollection(object):
    '''
    Collection of expression orders that summarizes the ordered operations
    involved in constructing one or more expressions.  This can be used to
    attempt to order the subexpressions within an expression that lacks a
    defined order.

    Ordering an expression that lacks a defined order can involve a number of
    `flatten` operations and the creation of many sets of expressions.  There
    is extensive caching to optimize this.
    '''
    def __init__(self, parents: Self | tuple[Self, ...] | list[Self] | None = None):
        self.Add_orders: dict[sympy.Expr, AddExprOrder] = {}
        self.Mul_orders: dict[sympy.Expr, MulExprOrder] = {}
        if isinstance(parents, ExprOrderCollection):
            parents = (parents,)
        if isinstance(parents, (tuple, list)) and all(isinstance(x, ExprOrderCollection) for x in parents):
            # `parents` precedence is left-to-right, so update in reverse
            for parent in reversed(parents):
                for expr, order in parent.Add_orders.items():
                    if not order.inferred or expr not in self.Add_orders or self.Add_orders[expr].inferred:
                        self.Add_orders[expr] = order
                for expr, order in parent.Mul_orders.items():
                    if not order.inferred or expr not in self.Mul_orders or self.Mul_orders[expr].inferred:
                        self.Mul_orders[expr] = order
        elif parents is not None:
            raise TypeError
        # Normalize/flatten caches aren't copied from parents.  To re-use
        # caches, would need to track whether `normalize_flatten_expr()` could
        # be modified by additional term orders, and would also need to merge
        # caches from multiple parents in a way that maintains parent order
        # precedence except for cases when one parent can provide a more
        # complete ordering.
        self._normalize_flatten_expr_cache: dict[sympy.Expr, tuple[sympy.Expr, ...]] = {}
        self._normalize_flatten_expr_to_frozenset_cache: dict[sympy.Expr, frozenset[sympy.Expr]] = {}
        self._normalize_flatten_frozenset_cache: dict[frozenset[sympy.Expr], frozenset[sympy.Expr]] = {}
        self._term_to_frozenset_cache: dict[tuple, frozenset[sympy.Expr]] = {}
        self._terms_to_shared_subexprs_frozenset_cache: dict[tuple[sympy.Expr, ...], frozenset[sympy.Expr]] = {}


    def normalize_flatten_expr(self, expr: sympy.Expr) -> tuple[sympy.Expr, ...]:
        '''
        Normalize a term while flattening `Add` and `Mul` into sequences of
        normalized terms.
        '''
        # No need to check for `OrderUnevaluatedExpr`, since that is
        # handled within `normalize_expr()`.
        if expr.is_Atom:
            return (expr,)
        try:
            return self._normalize_flatten_expr_cache[expr]
        except KeyError:
            normalized_flattened: list[sympy.Expr] = []
            if expr.is_Add:
                try:
                    args = self.Add_orders[expr].sorted_sympy_args
                except KeyError:
                    args = expr.args
                for arg in args:
                    normalized_flattened.extend(self.normalize_flatten_expr(arg))  # type: ignore
            elif expr.is_Mul:
                try:
                    args = self.Mul_orders[expr].sorted_sympy_args
                except KeyError:
                    args = expr.args
                for arg in args:
                    normalized_flattened.extend(self.normalize_flatten_expr(arg))  # type: ignore
            elif expr.is_Derivative:
                normalized_flattened.extend(self.normalize_flatten_expr(expr.expr))  # type: ignore
            elif hasattr(expr, 'has_finite_limits'):
                # `ExprWithLimits` is the only SymPy class with
                # `.has_finite_limits`.  It is the superclass of `Integral`,
                # `Sum`, and `Product`.
                function = expr.function  # type: ignore
                bound_symbols = frozenset(expr.bound_symbols)  # type: ignore
                for n_f in self.normalize_flatten_expr(function):  # type: ignore
                    if not (n_f.free_symbols & bound_symbols):
                        normalized_flattened.append(n_f)
                for _, lower_lim, upper_lim in expr.limits:  # type: ignore
                    # This can skip some valid expressions in cases when one
                    # limit involves a dummy variable for a subsequent limit.
                    # Situations when this prevents finding a good term order
                    # should be rare.  Handling the exceptional cases would
                    # probably add significant complexity.
                    if not (lower_lim.free_symbols & bound_symbols):
                        normalized_flattened.extend(self.normalize_flatten_expr(lower_lim))
                    if not (upper_lim.free_symbols & bound_symbols):
                        normalized_flattened.extend(self.normalize_flatten_expr(upper_lim))
            else:
                for n_t in normalize_expr(expr):
                    normalized_flattened.extend(self.normalize_flatten_expr(n_t))
            self._normalize_flatten_expr_cache[expr] = tuple(normalized_flattened)
            return self._normalize_flatten_expr_cache[expr]

    def normalize_flatten_expr_to_frozenset(self, expr: sympy.Expr) -> frozenset[sympy.Expr]:
        if expr.is_Atom:
            return frozenset((expr,))
        try:
            return self._normalize_flatten_expr_to_frozenset_cache[expr]
        except KeyError:
            self._normalize_flatten_expr_to_frozenset_cache[expr] = frozenset(self.normalize_flatten_expr(expr))
            return self._normalize_flatten_expr_to_frozenset_cache[expr]

    def normalize_flatten_frozenset(self, exprs_frozenset: frozenset[sympy.Expr]) -> frozenset[sympy.Expr]:
        '''
        Create a `frozenset` containing normalized, flattened expressions from
        an existing `frozenset`.
        '''
        try:
            return self._normalize_flatten_frozenset_cache[exprs_frozenset]
        except KeyError:
            normalized_flattened_frozenset = frozenset(
                n_f_e for e in exprs_frozenset for n_f_e in self.normalize_flatten_expr(e)
            )
            self._normalize_flatten_frozenset_cache[exprs_frozenset] = normalized_flattened_frozenset
            return self._normalize_flatten_frozenset_cache[exprs_frozenset]


    def term_to_frozenset(self, term: sympy.Expr,
                          ignore_terms_shared_subexprs: tuple[sympy.Expr, ...] | None, ignore_Number_coeffs: bool,
                          normalize: bool, flatten: bool) -> frozenset[sympy.Expr]:
        '''
        Convert a term of an `Add` operation into a `frozenset` of its
        subexpressions.  There are options to narrow which subexpressions are
        included in the `frozenset`, and to determine whether subexpressions
        are normalized and flattened.
        '''
        if ignore_terms_shared_subexprs is None:
            skipped_frozenset = None
        else:
            skipped_frozenset = self._terms_to_shared_subexprs_frozenset(
                ignore_terms_shared_subexprs, normalize, flatten
            )
        cache_key = (term, skipped_frozenset, ignore_Number_coeffs, normalize, flatten)
        try:
            return self._term_to_frozenset_cache[cache_key]
        except KeyError:
            if normalize:
                if flatten:
                    base_term_frozenset = self.normalize_flatten_expr_to_frozenset(term)
                elif term.is_Mul:
                    base_term_frozenset = unsorted_expand_normalize_Mul_to_frozenset(term)
                else:
                    base_term_frozenset = normalize_expr_to_frozenset(term)
            elif term.is_Mul:
                base_term_frozenset = frozenset(OrderUnevaluatedExpr.unwrap_exprs(term.args))  # type: ignore
            else:
                base_term_frozenset = frozenset((OrderUnevaluatedExpr.unwrap_expr(term),))

            if skipped_frozenset is None:
                if not ignore_Number_coeffs or all(t.is_Number for t in base_term_frozenset):
                    term_frozenset = base_term_frozenset
                else:
                    term_frozenset = frozenset(t for t in base_term_frozenset if not t.is_Number)
            else:
                if normalize:
                    if flatten:
                        skipped_frozenset = self.normalize_flatten_frozenset(skipped_frozenset)
                    else:
                        skipped_frozenset = normalize_frozenset(skipped_frozenset)
                if not ignore_Number_coeffs or all(t.is_Number for t in base_term_frozenset):
                    term_frozenset = base_term_frozenset - skipped_frozenset
                else:
                    term_frozenset = frozenset(
                        t for t in base_term_frozenset if not t.is_Number and t not in skipped_frozenset
                    )

            self._term_to_frozenset_cache[cache_key] = term_frozenset
            return self._term_to_frozenset_cache[cache_key]


    def _terms_to_shared_subexprs_frozenset(self, terms: tuple[sympy.Expr, ...],
                                            normalize: bool, flatten: bool) -> frozenset[sympy.Expr]:
        '''
        Create a frozenset of subexpressions that are shared by the terms in a
        sequence of terms.
        '''
        try:
            return self._terms_to_shared_subexprs_frozenset_cache[terms]
        except KeyError:
            union_set = set()
            repeated_set = set()
            if normalize:
                if flatten:
                    for term in terms:
                        n_f_t_frozenset = self.normalize_flatten_expr_to_frozenset(term)
                        repeated_set.update(n_f_t for n_f_t in n_f_t_frozenset if n_f_t in union_set)
                        union_set |= n_f_t_frozenset
                else:
                    for term in terms:
                        n_t_frozenset = normalize_expr_to_frozenset(term)
                        repeated_set.update(n_t for n_t in n_t_frozenset if n_t in union_set)
                        union_set |= n_t_frozenset
            else:
                for term in terms:
                    if term.is_Mul:
                        t_frozenset = frozenset(OrderUnevaluatedExpr.unwrap_exprs(term.args))  # type: ignore
                    else:
                        t_frozenset = frozenset((OrderUnevaluatedExpr.unwrap_expr(term),))
                    repeated_set.update(t for t in t_frozenset if t in union_set)
                    union_set |= t_frozenset
            self._terms_to_shared_subexprs_frozenset_cache[terms] = frozenset(union_set - repeated_set)
            return self._terms_to_shared_subexprs_frozenset_cache[terms]


    def append_Add_order(self, expr: sympy.Add, order: tuple[sympy.Expr, ...], inferred: bool = False):
        if expr not in self.Add_orders or (not inferred and self.Add_orders[expr].inferred):
            self.Add_orders[expr] = AddExprOrder(expr, order, self, inferred)

    def append_Mul_order(self, expr: sympy.Mul, order: tuple[sympy.Expr, ...], inferred: bool = False):
        if expr not in self.Mul_orders or (not inferred and self.Mul_orders[expr].inferred):
            self.Mul_orders[expr] = MulExprOrder(expr, order, self, inferred)


    # Sorting methods take an expression as an optional argument.  Printing
    # functions sometimes only operate on part of an expression, so more
    # context is available by having both the full expression and the separate
    # sequence of subexpressions to be sorted.

    def sort_terms(self, terms: list[sympy.Expr] | tuple[sympy.Expr, ...],
                   expr: sympy.Expr | None = None) -> tuple[sympy.Expr, ...]:
        if len(terms) == 1:
            return tuple(terms)
        if expr is not None:
            try:
                expr_order = self.Add_orders[expr]
            except KeyError:
                pass
            else:
                indices, _, _, _ = expr_order.terms_to_indices(terms, self)
                terms_to_keys = {t: k for t, k in zip(terms, indices)}
                return tuple(sorted(terms, key=lambda t: terms_to_keys[t]))

        remaining_terms: tuple[sympy.Expr, ...] = tuple(terms)
        ordered_terms: list[sympy.Expr] = []
        ignore_shared_subexprs: bool = False
        flatten: bool = False
        while True:
            count_remaining_Numberlike: int = OrderUnevaluatedExpr.count_terms_Numberlike(remaining_terms)
            count_remaining_not_Numberlike: int = len(remaining_terms) - count_remaining_Numberlike
            best_indices: tuple[int | float, ...] | None = None
            best_count_found: int = 0
            best_count_found_not_Numberlike: int = 0
            best_count_repeated: int = 0

            if not flatten:
                for term_orders, normalize in ((self.Add_orders.values(), False),
                                               ((t.normalize() for t in self.Add_orders.values()), True)):
                    if ignore_shared_subexprs:
                        if not self._terms_to_shared_subexprs_frozenset(remaining_terms, normalize, flatten):
                            continue
                    for term_order in term_orders:
                        indices_data = term_order.terms_to_indices(remaining_terms, self, ignore_shared_subexprs)
                        indices, count_found, count_found_not_Numberlike, count_repeated = indices_data
                        if count_found < 2:
                            continue
                        if (count_found_not_Numberlike > best_count_found_not_Numberlike or
                                (count_found_not_Numberlike == best_count_found_not_Numberlike and
                                 (count_found > best_count_found or count_repeated < best_count_repeated))):
                            best_indices = indices
                            best_count_found = count_found
                            best_count_found_not_Numberlike = count_found_not_Numberlike
                            best_count_repeated = count_repeated
                            if best_count_found == len(remaining_terms):
                                break
                    if best_indices and best_count_found_not_Numberlike == count_remaining_not_Numberlike:
                        # Don't normalize if that could only improve numbers
                        break
                if count_remaining_not_Numberlike >= 2 and best_count_found_not_Numberlike < 2:
                    if ignore_shared_subexprs:
                        # Only switch to flatten when there are multiple
                        # non-Numberlike that still need sorting.
                        flatten = True
                        ignore_shared_subexprs = False
                    else:
                        ignore_shared_subexprs = True
                        continue
            # Only normalize and flatten when absolutely necessary, since
            # flattening can be less precise and introduce ambiguity.
            if flatten:
                if ignore_shared_subexprs:
                    if not self._terms_to_shared_subexprs_frozenset(remaining_terms, normalize, flatten):
                        break
                for term_order in (t.normalize_flatten() for t in self.Add_orders.values()):
                    indices_data = term_order.terms_to_indices(remaining_terms, self, ignore_shared_subexprs)
                    indices, count_found, count_found_not_Numberlike, count_repeated = indices_data
                    if count_found < 2:
                        continue
                    if (count_found_not_Numberlike > best_count_found_not_Numberlike or
                            (count_found_not_Numberlike == best_count_found_not_Numberlike and
                                (count_found > best_count_found or count_repeated < best_count_repeated))):
                        best_indices = indices
                        best_count_found = count_found
                        best_count_found_not_Numberlike = count_found_not_Numberlike
                        best_count_repeated = count_repeated
                        if best_count_found == len(remaining_terms):
                            break
                if count_remaining_not_Numberlike >= 2 and best_count_found_not_Numberlike < 2:
                    if not ignore_shared_subexprs:
                        ignore_shared_subexprs = True
                        continue

            if best_indices is None or (count_remaining_not_Numberlike >= 2 and best_count_found_not_Numberlike < 2):
                break

            terms_to_keys = {t: k for t, k in zip(remaining_terms, best_indices)}
            partially_ordered_terms = sorted(remaining_terms, key=lambda t: terms_to_keys[t])
            remaining_terms = tuple(t for t in partially_ordered_terms if terms_to_keys[t] == -1)
            ordered_terms.extend(t for t in partially_ordered_terms if terms_to_keys[t] != -1)
            if len(remaining_terms) < 2:
                break

        remaining_Numberlike_terms: list[sympy.Expr] = []
        for t in remaining_terms:
            if OrderUnevaluatedExpr.is_Numberlike(t):
                remaining_Numberlike_terms.append(t)
            else:
                ordered_terms.append(t)

        return tuple(remaining_Numberlike_terms + ordered_terms)


    def sort_factors(self, factors: list[sympy.Expr] | tuple[sympy.Expr, ...],
                     expr: sympy.Expr | None = None) -> tuple[sympy.Expr, ...]:
        if len(factors) == 1:
            return tuple(factors)
        if expr is not None:
            try:
                expr_order = self.Mul_orders[expr]
            except KeyError:
                pass
            else:
                indices, _, _, _ = expr_order.factors_to_indices(factors, self)
                factors_to_keys = {f: k for f, k in zip(factors, indices)}
                return tuple(sorted(factors, key=lambda f: factors_to_keys[f]))

        remaining_factors: tuple[sympy.Expr, ...] = tuple(factors)
        ordered_factors: list[sympy.Expr] = []
        flatten: bool = False
        while True:
            count_remaining_Numberlike: int = OrderUnevaluatedExpr.count_factors_Numberlike(remaining_factors)
            count_remaining_not_Numberlike: int = len(remaining_factors) - count_remaining_Numberlike
            best_indices: tuple[int | float, ...] | None = None
            best_count_found: int = 0
            best_count_found_not_Numberlike: int = 0
            best_count_repeated: int = 0

            if not flatten:
                for factor_orders, normalize in ((self.Mul_orders.values(), False),
                                               ((f.normalize() for f in self.Mul_orders.values()), True)):
                    for factor_order in factor_orders:
                        indices_data = factor_order.factors_to_indices(remaining_factors, self)
                        indices, count_found, count_found_not_Numberlike, count_repeated = indices_data
                        if count_found < 2:
                            continue
                        if (count_found_not_Numberlike > best_count_found_not_Numberlike or
                                (count_found_not_Numberlike == best_count_found_not_Numberlike and
                                 (count_found > best_count_found or count_repeated < best_count_repeated))):
                            best_indices = indices
                            best_count_found = count_found
                            best_count_found_not_Numberlike = count_found_not_Numberlike
                            best_count_repeated = count_repeated
                            if best_count_found == len(remaining_factors):
                                break
                    if best_indices and best_count_found_not_Numberlike == count_remaining_not_Numberlike:
                        # Don't normalize if that could only improve numbers
                        break
                if count_remaining_not_Numberlike >= 2 and best_count_found_not_Numberlike < 2:
                    # Only switch to flatten when there are multiple
                    # non-Numberlike that still need sorting.
                    flatten = True
            # Only normalize and flatten when absolutely necessary, since
            # flattening can be less precise and introduce ambiguity.
            if flatten:
                for factor_order in (f.normalize_flatten() for f in self.Mul_orders.values()):
                    indices_data = factor_order.factors_to_indices(remaining_factors, self)
                    indices, count_found, count_found_not_Numberlike, count_repeated = indices_data
                    if count_found < 2:
                        continue
                    if (count_found_not_Numberlike > best_count_found_not_Numberlike or
                            (count_found_not_Numberlike == best_count_found_not_Numberlike and
                             (count_found > best_count_found or count_repeated < best_count_repeated))):
                        best_indices = indices
                        best_count_found = count_found
                        best_count_found_not_Numberlike = count_found_not_Numberlike
                        best_count_repeated = count_repeated
                        if best_count_found == len(remaining_factors):
                            break

            if best_indices is None or (count_remaining_not_Numberlike >= 2 and best_count_found_not_Numberlike < 2):
                break

            factors_to_keys = {f: k for f, k in zip(remaining_factors, best_indices)}
            partially_ordered_factors = sorted(remaining_factors, key=lambda f: factors_to_keys[f])
            remaining_factors = tuple(f for f in partially_ordered_factors if factors_to_keys[f] == -1)
            ordered_factors.extend(f for f in partially_ordered_factors if factors_to_keys[f] != -1)
            if len(remaining_factors) < 2:
                break

        remaining_Numberlike_factors = []
        for f in remaining_factors:
            if OrderUnevaluatedExpr.is_Numberlike(f):
                remaining_Numberlike_factors.append(f)
            else:
                ordered_factors.append(f)

        return tuple(remaining_Numberlike_factors + ordered_factors)




class BaseWrapped(object):
    '''
    Wrapper for SymPy `Basic` that records expression order during operations.
    This allows order to be reconstructed, at least partially, during
    printing.
    '''

    __slots__ = ('wrapped', 'expr_order_collection')

    __sympy__ = True  # So that `_sympify()` accepts instances

    wrapped: sympy.Basic
    expr_order_collection: ExprOrderCollection

    def __new__(cls, *args, **kwargs) -> Self:
        if cls is BaseWrapped:
            raise NotImplementedError
        return super().__new__(cls)

    def _init_expr_order_collection(self, parents: BaseWrapped | list[BaseWrapped] | tuple[BaseWrapped, ...] | None):
        if isinstance(parents, BaseWrapped):
            self.expr_order_collection = ExprOrderCollection(parents.expr_order_collection)
        elif isinstance(parents, (list, tuple)) and all(isinstance(x, BaseWrapped) for x in parents):
            self.expr_order_collection = ExprOrderCollection(tuple(p.expr_order_collection for p in parents))
        elif parents is None:
            self.expr_order_collection = ExprOrderCollection()
        else:
            raise TypeError



class WrappedExpr(BaseWrapped):
    __slots__ = ()
    wrapped: sympy.Expr

    def __init__(self, expr: sympy.Expr, *,
                 parents: BaseWrapped | list[BaseWrapped] | tuple[BaseWrapped, ...] | None = None):
        if not isinstance(expr, sympy.Expr):
            raise TypeError
        self.wrapped = expr
        self._init_expr_order_collection(parents)

    @staticmethod
    def wrapper_class(*args, **kwargs) -> WrappedExpr:
        return WrappedExpr(*args, **kwargs)

    @property
    def expr(self) -> sympy.Expr:
        return self.wrapped

    # Max `_op_priority` in SymPy is 20.  Use a much larger value here to
    # force wrapped expr methods.  Then the wrapped expr methods will invoke
    # the expr methods, which will use the normal expr `_op_priority`.
    _op_priority = 1_000_000

    @property
    def has_order(self) -> bool:
        if self.expr.is_Add and self.expr in self.expr_order_collection.Add_orders:
            return True
        if self.expr.is_Mul and self.expr in self.expr_order_collection.Mul_orders:
            return True
        return False

    @property
    def needs_order(self) -> bool:
        if self.expr.is_Add:
            try:
                return self.expr_order_collection.Add_orders[self.expr].inferred
            except KeyError:
                return True
        if self.expr.is_Mul:
            try:
                return self.expr_order_collection.Mul_orders[self.expr].inferred
            except KeyError:
                return True
        return False


    @classmethod
    def wrap_atomic_expr_class(cls, expr_class: type[sympy.Expr]) -> Callable[..., Self]:
        if not issubclass(expr_class, sympy.Expr) or not expr_class.is_Atom:
            raise TypeError
        def wrapper(*args, **kwargs):
            filtered_args = []
            parents = []
            for arg in args:
                if isinstance(arg, WrappedExpr):
                    filtered_args.append(arg.expr)
                    parents.append(arg)
                else:
                    filtered_args.append(arg)
            return cls(expr_class(*filtered_args, **kwargs), parents=parents)
        return wrapper

    @classmethod
    def wrap_compound_expr_class(cls, expr_class: type[sympy.Expr]) -> Callable[..., Self]:
        if not issubclass(expr_class, sympy.Expr) or expr_class.is_Atom:
            raise TypeError
        def wrapper(*args, **kwargs):
            filtered_args = []
            parents = []
            for arg in args:
                if isinstance(arg, WrappedExpr):
                    filtered_args.append(arg.expr)
                    parents.append(arg)
                else:
                    filtered_args.append(arg)
            expr = expr_class(*filtered_args, **kwargs)
            if not isinstance(expr, expr_class):
                if len(filtered_args) == 1:
                    arg_0 = filtered_args[0]
                    if (isinstance(arg_0, (sympy.Expr, float, complex)) or
                            (isinstance(arg_0, int) and not isinstance(arg_0, bool))):
                        filtered_args[0] = OrderUnevaluatedExpr.from_any(arg_0)
                        expr = expr_class(*filtered_args, **kwargs)
                elif filtered_args:
                    for n, arg in enumerate(filtered_args):
                        if (n < 1 or isinstance(arg, (sympy.Expr, float, complex)) or
                                (isinstance(arg, int) and not isinstance(arg, bool))):
                            continue
                        filtered_args[n-1] = OrderUnevaluatedExpr.from_any(filtered_args[n-1])
                        expr = expr_class(*filtered_args, **kwargs)
                        break
            return cls(expr, parents=parents)
        return wrapper

    @classmethod
    def wrap_expr_callable(cls, callable: Callable[..., sympy.Expr]) -> Callable[..., Self]:
        def wrapper(*args, **kwargs):
            kwargs.setdefault('evaluate', False)
            filtered_args = []
            parents = []
            for arg in args:
                if isinstance(arg, WrappedExpr):
                    filtered_args.append(arg.expr)
                    parents.append(arg)
                else:
                    filtered_args.append(arg)
            expr = callable(*filtered_args, **kwargs)
            return cls(expr, parents=parents)
        return wrapper


    def n(self) -> sympy.Expr:
        return self.wrapped.n()

    @property
    def is_Symbol(self) -> bool:
        return False

    def as_dummy(self):
        raise TypeError


    def __pos__(self) -> WrappedExpr:
        self_expr = self.expr
        new_expr = self_expr.__pos__()
        if new_expr is self_expr:
            return self
        raise NotImplementedError

    def __neg__(self) -> WrappedExpr:
        self_expr = self.expr
        if not self_expr.is_Add and not self_expr.is_Mul:
            new_expr = self_expr.__neg__()
            if new_expr.is_Mul and new_expr.args == (S.NegativeOne, self_expr):
                new_wrapped = self.wrapper_class(new_expr, parents=self)
                new_wrapped.expr_order_collection.append_Mul_order(new_expr, new_expr.args)
                return new_wrapped
        # For `Add()` and `Mul()`:  Negated expression must be in parentheses,
        # be created explicitly via `Add()` or `Mul()`, or be the output of an
        # evaluated function or operation.  Preserve expression order by
        # preventing an `Add` from being converted into a new `Add` with all
        # terms negated, and by preventing a `Mul` from being converted into a
        # new `Mul` prepended with `-1`.
        unevaluated_self_expr = OrderUnevaluatedExpr(self_expr)
        new_expr = unevaluated_self_expr.__neg__()
        new_wrapped = self.wrapper_class(new_expr, parents=self)
        if new_expr.is_Mul and new_expr.args == (S.NegativeOne, unevaluated_self_expr):
            new_wrapped.expr_order_collection.append_Mul_order(new_expr, new_expr.args)
        return new_wrapped


    def __abs__(self) -> WrappedExpr:
        self_expr = self.expr
        new_expr = self_expr.__abs__()
        if not (new_expr.is_Function and new_expr.func is sympy.Abs and new_expr.args[0] is self_expr):
            new_expr = OrderUnevaluatedExpr(self_expr).__abs__()
        return self.wrapper_class(new_expr, parents=self)


    def __add__(self, other) -> WrappedExpr:
        self_expr = self.expr
        if isinstance(other, WrappedExpr):
            other_expr = other.expr
            parents = (self, other)
        else:
            try:
                other_expr = _sympify(other)
            except SympifyError:
                return NotImplemented
            parents = self
        return self._add_exprs(self_expr, other_expr, parents)

    def __radd__(self, other) -> WrappedExpr:
        self_expr = self.expr
        if isinstance(other, WrappedExpr):
            # This shouldn't happen unless `__radd__()` is invoked explicitly
            other_expr = other.expr
            parents = (other, self)
        else:
            try:
                other_expr = _sympify(other)
            except SympifyError:
                return NotImplemented
            parents = self
        return self._add_exprs(other_expr, self_expr, parents)

    def _add_exprs(self, left_expr: sympy.Expr, right_expr: sympy.Expr,
                   parents: WrappedExpr | tuple[WrappedExpr, ...]) -> WrappedExpr:
        if right_expr.is_Add:
            # Must be in parentheses, be created explicitly, or be the output
            # of an evaluated function or operation
            right_expr = OrderUnevaluatedExpr(right_expr)
        if left_expr is S.Zero:
            return self.wrapper_class(right_expr, parents=parents)
        # Try regular addition regardless of whether this was invoked by
        # `__add__()` or `__radd()__`, so that SymPy precedence is invoked.
        new_expr = left_expr + right_expr  # type: ignore
        need_unevaluated: bool = False
        if not new_expr.is_Add:
            need_unevaluated = True
        elif left_expr.is_Add:
            if len(new_expr.args) != len(left_expr.args) + 1:
                need_unevaluated = True
        elif len(new_expr.args) != 2:
            need_unevaluated = True
        if need_unevaluated:
            right_expr = OrderUnevaluatedExpr(right_expr)
            new_expr = left_expr + right_expr  # type: ignore
        if not new_expr.is_Add:
            raise NotImplementedError(
                f'Addition that tracks expression order is not supported for "{left_expr}" + "{right_expr}"'
            )
        new_wrapped = self.wrapper_class(new_expr, parents=parents)
        if not new_wrapped.needs_order:
            return new_wrapped
        # Length checks below should be redundant with those above, but
        # safeguard against custom classes that evaluate `UnevaluatedExpr`
        # instances.
        if left_expr.is_Add:
            if len(new_expr.args) == len(left_expr.args) + 1:
                left_key = OrderUnevaluatedExpr.unwrap_expr(left_expr)
                try:
                    left_sorted_args = new_wrapped.expr_order_collection.Add_orders[left_key].sorted_sympy_args
                except KeyError:
                    return new_wrapped
                new_sorted_args = left_sorted_args + (right_expr,)
                new_wrapped.expr_order_collection.append_Add_order(new_expr, new_sorted_args)  # type: ignore
        elif len(new_expr.args) == 2:
            new_wrapped.expr_order_collection.append_Add_order(new_expr, (left_expr, right_expr))
        return new_wrapped


    def __sub__(self, other) -> WrappedExpr:
        self_expr = self.expr
        if isinstance(other, WrappedExpr):
            other_expr = other.expr
            parents = (self, other)
        else:
            try:
                other_expr = _sympify(other)
            except SympifyError:
                return NotImplemented
            parents = self
        return self._sub_exprs(self_expr, other_expr, parents)

    def __rsub__(self, other) -> WrappedExpr:
        self_expr = self.expr
        if isinstance(other, WrappedExpr):
            # This shouldn't happen unless `__rsub__()` is invoked explicitly
            other_expr = other.expr
            parents = (other, self)
        else:
            try:
                other_expr = _sympify(other)
            except SympifyError:
                return NotImplemented
            parents = self
        return self._sub_exprs(other_expr, self_expr, parents)

    def _sub_exprs(self, left_expr: sympy.Expr, right_expr: sympy.Expr,
                   parents: WrappedExpr | tuple[WrappedExpr, ...]) -> WrappedExpr:
        if right_expr.is_Add:
            # Must be in parentheses, be created explicitly, or be the output
            # of an evaluated function or operation
            right_expr = OrderUnevaluatedExpr(right_expr)
        if left_expr is S.Zero:
            return -self.wrapper_class(right_expr, parents=parents)
        new_expr = left_expr - right_expr  # type: ignore
        need_unevaluated: bool = False
        if not new_expr.is_Add:
            need_unevaluated = True
        elif left_expr.is_Add:
            if len(new_expr.args) != len(left_expr.args) + 1:
                need_unevaluated = True
        elif len(new_expr.args) != 2:
            need_unevaluated = True
        if need_unevaluated:
            right_expr = OrderUnevaluatedExpr(right_expr)
            new_expr = left_expr - right_expr  # type: ignore
        if not new_expr.is_Add:
            raise NotImplementedError(
                f'Subtraction that tracks expression order is not supported for "{left_expr}" - "{right_expr}"'
            )
        new_wrapped = self.wrapper_class(new_expr, parents=parents)
        if not new_wrapped.needs_order:
            return new_wrapped
        # Length checks below should be redundant with those above, but
        # safeguard against custom classes that evaluate `UnevaluatedExpr`
        # instances.
        if left_expr.is_Add:
            if len(new_expr.args) == len(left_expr.args) + 1:
                left_key = OrderUnevaluatedExpr.unwrap_expr(left_expr)
                try:
                    left_sorted_args = new_wrapped.expr_order_collection.Add_orders[left_key].sorted_sympy_args
                except KeyError:
                    return new_wrapped
                new_sorted_args = left_sorted_args + (-right_expr,)
                new_wrapped.expr_order_collection.append_Add_order(new_expr, new_sorted_args)  # type: ignore
        elif len(new_expr.args) == 2:
            new_wrapped.expr_order_collection.append_Add_order(new_expr, (left_expr, -right_expr))
        return new_wrapped


    def __mul__(self, other) -> WrappedExpr:
        self_expr = self.expr
        if isinstance(other, WrappedExpr):
            other_expr = other.expr
            parents = (self, other)
        else:
            try:
                other_expr = _sympify(other)
            except SympifyError:
                return NotImplemented
            parents = self
        return self._mul_exprs(self_expr, other_expr, parents)

    def __rmul__(self, other) -> WrappedExpr:
        self_expr = self.expr
        if isinstance(other, WrappedExpr):
            # This shouldn't happen unless `__rmul__()` is invoked explicitly
            other_expr = other.expr
            parents = (other, self)
        else:
            try:
                other_expr = _sympify(other)
            except SympifyError:
                return NotImplemented
            parents = self
        return self._mul_exprs(other_expr, self_expr, parents)

    def _mul_exprs(self, left_expr: sympy.Expr, right_expr: sympy.Expr,
                   parents: WrappedExpr | tuple[WrappedExpr, ...]) -> WrappedExpr:
        if right_expr.is_Mul:
            # Must be in parentheses, be created explicitly, or be the output
            # of an evaluated function or operation
            right_expr = OrderUnevaluatedExpr(right_expr)
        if left_expr is S.One:
            return self.wrapper_class(right_expr, parents=parents)
        new_expr = left_expr * right_expr  # type: ignore
        need_unevaluated: bool = False
        if not new_expr.is_Mul:
            need_unevaluated = True
        elif left_expr.is_Mul:
            if len(new_expr.args) != len(left_expr.args) + 1:
                need_unevaluated = True
        elif len(new_expr.args) != 2:
            need_unevaluated = True
        if need_unevaluated:
            right_expr = OrderUnevaluatedExpr(right_expr)
            new_expr = left_expr * right_expr  # type: ignore
        if not new_expr.is_Mul:
            raise NotImplementedError(
                f'Multiplication that tracks expression order is not supported for "{left_expr}" * "{right_expr}"'
            )
        new_wrapped = self.wrapper_class(new_expr, parents=parents)
        if not new_wrapped.needs_order:
            return new_wrapped
        # Length checks below should be redundant with those above, but
        # safeguard against custom classes that evaluate `UnevaluatedExpr`
        # instances.
        if left_expr.is_Mul:
            if len(new_expr.args) == len(left_expr.args) + 1:
                left_key = OrderUnevaluatedExpr.unwrap_expr(left_expr)
                try:
                    left_sorted_args = new_wrapped.expr_order_collection.Mul_orders[left_key].sorted_sympy_args
                except KeyError:
                    return new_wrapped
                new_sorted_args = left_sorted_args + (right_expr,)
                new_wrapped.expr_order_collection.append_Mul_order(new_expr, new_sorted_args)  # type: ignore
        elif len(new_expr.args) == 2:
            new_wrapped.expr_order_collection.append_Mul_order(new_expr, (left_expr, right_expr))
        return new_wrapped


    def __truediv__(self, other) -> WrappedExpr:
        self_expr = self.expr
        if isinstance(other, WrappedExpr):
            other_expr = other.expr
            parents = (self, other)
        else:
            try:
                other_expr = _sympify(other)
            except SympifyError:
                return NotImplemented
            parents = self
        return self._truediv_exprs(self_expr, other_expr, parents)

    def __rtruediv__(self, other) -> WrappedExpr:
        self_expr = self.expr
        if isinstance(other, WrappedExpr):
            # This shouldn't happen unless `__rtruediv__()` is invoked
            # explicitly
            other_expr = other.expr
            parents = (other, self)
        else:
            try:
                other_expr = _sympify(other)
            except SympifyError:
                return NotImplemented
            parents = self
        return self._truediv_exprs(other_expr, self_expr, parents)

    def _truediv_exprs(self, left_expr: sympy.Expr, right_expr: sympy.Expr,
                       parents: WrappedExpr | tuple[WrappedExpr, ...]) -> WrappedExpr:
        if right_expr.is_Mul:
            # Must be in parentheses, be created explicitly, or be the output
            # of an evaluated function or operation
            right_expr = OrderUnevaluatedExpr(right_expr)
        new_expr = left_expr / right_expr  # type: ignore
        if left_expr is S.One:
            if not new_expr.is_Pow or new_expr.args != (right_expr, S.NegativeOne):
                right_expr = OrderUnevaluatedExpr(right_expr)
                new_expr = left_expr / right_expr  # type: ignore
            return self.wrapper_class(new_expr, parents=parents)
        need_unevaluated: bool = False
        if not new_expr.is_Mul:
            need_unevaluated = True
        elif left_expr.is_Mul:
            if len(new_expr.args) != len(left_expr.args) + 1:
                need_unevaluated = True
        elif len(new_expr.args) != 2:
            need_unevaluated = True
        if need_unevaluated:
            right_expr = OrderUnevaluatedExpr(right_expr)
            new_expr = left_expr / right_expr  # type: ignore
        if not new_expr.is_Mul:
            raise NotImplementedError(
                f'Division that tracks expression order is not supported for "{left_expr}" / "{right_expr}"'
            )
        new_wrapped = self.wrapper_class(new_expr, parents=parents)
        if not new_wrapped.needs_order:
            return new_wrapped
        # Length checks below should be redundant with those above, but
        # safeguard against custom classes that evaluate `UnevaluatedExpr`
        # instances.
        if left_expr.is_Mul:
            if len(new_expr.args) == len(left_expr.args) + 1:
                left_key = OrderUnevaluatedExpr.unwrap_expr(left_expr)
                try:
                    left_sorted_args = new_wrapped.expr_order_collection.Mul_orders[left_key].sorted_sympy_args
                except KeyError:
                    return new_wrapped
                new_sorted_args = left_sorted_args + (sympy.Pow(right_expr, -1),)
                new_wrapped.expr_order_collection.append_Mul_order(new_expr, new_sorted_args)  # type: ignore
        elif len(new_expr.args) == 2:
            new_wrapped.expr_order_collection.append_Mul_order(new_expr, (left_expr, sympy.Pow(right_expr, -1)))
        return new_wrapped


    def __pow__(self, other) -> WrappedExpr:
        self_expr = self.expr
        if isinstance(other, WrappedExpr):
            other_expr = other.expr
            parents = (self, other)
        else:
            try:
                other_expr = _sympify(other)
            except SympifyError:
                return NotImplemented
            parents = self
        return self._pow_exprs(self_expr, other_expr, parents)

    def __rpow__(self, other) -> WrappedExpr:
        self_expr = self.expr
        if isinstance(other, WrappedExpr):
            # This shouldn't happen unless `__rpow__()` is invoked explicitly
            other_expr = other.expr
            parents = (other, self)
        else:
            try:
                other_expr = _sympify(other)
            except SympifyError:
                return NotImplemented
            parents = self
        return self._pow_exprs(other_expr, self_expr, parents)

    def _pow_exprs(self, left_expr: sympy.Expr, right_expr: sympy.Expr,
                   parents: WrappedExpr | tuple[WrappedExpr, ...]) -> WrappedExpr:
        new_expr = left_expr ** right_expr
        if not new_expr.is_Pow or new_expr.args != (left_expr, right_expr):
            unevaluated_right_expr = OrderUnevaluatedExpr(right_expr)
            new_expr = left_expr ** unevaluated_right_expr
            if not new_expr.is_Pow or new_expr.args != (left_expr, unevaluated_right_expr):
                unevaluated_left_expr = OrderUnevaluatedExpr(left_expr)
                new_expr = unevaluated_left_expr ** right_expr
                if not new_expr.is_Pow or new_expr.args != (unevaluated_left_expr, right_expr):
                    new_expr = unevaluated_left_expr ** unevaluated_right_expr
        return self.wrapper_class(new_expr, parents=parents)




class WrappedEq(BaseWrapped):
    wrapped: sympy.Eq

    __slots__ = ('wrapped_lhs', 'wrapped_rhs')

    def __init__(self, eq: sympy.Eq, *,
                 parents: BaseWrapped | tuple[BaseWrapped, ...] | list[BaseWrapped] | None = None):
        if not isinstance(eq, sympy.Eq):
            raise TypeError
        self.wrapped = eq
        self._init_expr_order_collection(parents)
        self.wrapped_lhs = None
        self.wrapped_rhs = None

        left_key = OrderUnevaluatedExpr.unwrap_expr(self.wrapped.lhs)  # type: ignore
        if left_key.is_Add:
            try:
                left_sorted_args = self.expr_order_collection.Add_orders[left_key].sorted_sympy_args
            except KeyError:
                return
        else:
            left_sorted_args = (left_key,)
        right_key = OrderUnevaluatedExpr.unwrap_expr(self.wrapped.rhs)  # type: ignore
        if right_key.is_Add:
            try:
                right_sorted_args = self.expr_order_collection.Add_orders[right_key].sorted_sympy_args
            except KeyError:
                return
        else:
            right_sorted_args = (right_key,)
        # These are only inferred term orders, not actual operations ordered
        # by user input, so don't prevent simplification with
        # `OrderUnevaluatedExpr`.  Instead, just skip the inferred order if
        # simplification occurs.
        len_lhs_and_rhs_args = 0
        if self.wrapped.lhs.is_Add:
            len_lhs_and_rhs_args += len(self.wrapped.lhs.args)
        else:
            len_lhs_and_rhs_args += 1
        if self.wrapped.rhs.is_Add:
            len_lhs_and_rhs_args += len(self.wrapped.rhs.args)
        else:
            len_lhs_and_rhs_args += 1
        lhs_minus_rhs_inferred_expr = self.wrapped.lhs - self.wrapped.rhs  # type: ignore
        if lhs_minus_rhs_inferred_expr.is_Add and len(lhs_minus_rhs_inferred_expr.args) == len_lhs_and_rhs_args:
            self.expr_order_collection.append_Add_order(lhs_minus_rhs_inferred_expr,
                                                        left_sorted_args + tuple(-x for x in right_sorted_args),
                                                        inferred=True)
        rhs_minus_lhs_inferred_expr = self.wrapped.rhs - self.wrapped.lhs  # type: ignore
        if rhs_minus_lhs_inferred_expr.is_Add and len(rhs_minus_lhs_inferred_expr.args) == len_lhs_and_rhs_args:
            self.expr_order_collection.append_Add_order(rhs_minus_lhs_inferred_expr,
                                                        right_sorted_args + tuple(-x for x in left_sorted_args),
                                                        inferred=True)


    @property
    def eq(self) -> sympy.Eq:
        return self.wrapped

    @property
    def lhs(self) -> WrappedExpr:
        if self.wrapped_lhs is not None:
            return self.wrapped_lhs
        return WrappedExpr(self.wrapped.lhs, parents=self)

    @property
    def rhs(self) -> WrappedExpr:
        if self.wrapped_rhs is not None:
            return self.wrapped_rhs
        return WrappedExpr(self.wrapped.rhs, parents=self)


    _wrap_eq_class_cache: dict[type[sympy.Eq], Callable[..., Self]] = {}

    @classmethod
    def wrap_eq_class(cls, eq_class: type[sympy.Eq]) -> Callable[..., Self]:
        try:
            return cls._wrap_eq_class_cache[eq_class]
        except KeyError:
            pass
        if not issubclass(eq_class, sympy.Eq):
            raise TypeError
        def wrapper(lhs: sympy.Expr | WrappedExpr, rhs: sympy.Expr | WrappedExpr,
                    parents: BaseWrapped | tuple[BaseWrapped, ...] | list[BaseWrapped] | None = None,
                    **options) -> Self:
            if parents is None:
                parents = []
            elif isinstance(parents, BaseWrapped):
                parents = [parents]
            elif isinstance(parents, (list, tuple)):
                parents = list(parents)
            else:
                raise TypeError
            # Parents precedence (left to right): lhs, rhs, other parents
            if isinstance(rhs, WrappedExpr):
                wrapped_rhs = rhs
                parents.insert(0, wrapped_rhs)
                rhs = rhs.expr
            else:
                wrapped_rhs = None
            if isinstance(lhs, WrappedExpr):
                wrapped_lhs = lhs
                parents.insert(0, wrapped_lhs)
                lhs = lhs.expr
            else:
                wrapped_lhs = None
            options['evaluate'] = False
            eq: sympy.Eq = eq_class(lhs, rhs, **options)  # type: ignore
            wrapped_eq = cls(eq, parents=parents)
            wrapped_eq.wrapped_lhs = wrapped_lhs
            wrapped_eq.wrapped_rhs = wrapped_rhs
            return wrapped_eq
        cls._wrap_eq_class_cache[eq_class] = wrapper
        return cls._wrap_eq_class_cache[eq_class]



    @property
    def free_symbols(self) -> set[sympy.Expr]:
        return self.wrapped.free_symbols  # type: ignore


    def xreplace(self, raw_rule: dict) -> Self:
        # Might consider also processing key and values with
        # `OrderUnevaluatedExpr.unwrap_expr()`
        rule = {}
        for k, v in raw_rule.items():
            if isinstance(k, WrappedExpr):
                k = k.expr
            if isinstance(v, WrappedExpr):
                v = v.expr
            rule[k] = v
        if self.wrapped_lhs and self.wrapped_rhs:
            new_eq = self.wrapped.xreplace(rule)
            new_wrapped_eq = type(self)(new_eq, parents=self)
            new_wrapped_eq.wrapped_lhs = self.wrapped_lhs.wrapper_class(new_eq.lhs, parents=self.wrapped_lhs)
            new_wrapped_eq.wrapped_rhs = self.wrapped_rhs.wrapper_class(new_eq.rhs, parents=self.wrapped_rhs)
            return new_wrapped_eq
        return type(self)(self.wrapped.xreplace(rule), parents=self)
