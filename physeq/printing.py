# -*- coding: utf-8 -*-
#
# Copyright (c) 2025-2026, Geoffrey M. Poore
# All rights reserved.
#
# Licensed under the BSD 3-Clause License:
# http://opensource.org/licenses/BSD-3-Clause
#
# ----------------------------------------------------------------------------
#
# Some methods of the `LatexPrinter` class are adapted from SymPy
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

from astropy.units import Quantity, Unit
from sympy import Basic, Expr, Float, Number, Mul, Pow, Symbol
from sympy.core.singleton import S
from sympy.printing.latex import _between_two_numbers_p
from sympy.printing.latex import LatexPrinter as SymPyLatexPrinter
from sympy.utilities.iterables import sift
from . import exprorder, symbol




class LatexPrinter(SymPyLatexPrinter):
    def __init__(self, settings):
        separate_numerical_frac = settings.pop('separate_numerical_frac', False)
        if not isinstance(separate_numerical_frac, bool):
            raise TypeError
        float_fmtspec = settings.pop('float_fmtspec', None)
        if float_fmtspec is not None and not isinstance(float_fmtspec, str):
            raise TypeError
        symbol_replace = settings.pop('symbol_replace', None)
        if symbol_replace is not None and not isinstance(symbol_replace, dict):
            raise TypeError
        if symbol_replace is not None:
            symbol_replace = symbol.translate_numerical_xreplace_rule(symbol_replace)
        symbol_replace_all = settings.pop('symbol_replace_all', None)
        if symbol_replace_all is not None and not isinstance(symbol_replace_all, dict):
            raise TypeError
        if symbol_replace_all is not None:
            symbol_replace_all = symbol.translate_numerical_xreplace_rule(symbol_replace_all)
        if symbol_replace and symbol_replace_all:
            raise TypeError('cannot use "symbol_replace" and "symbol_replace_all" at the same time')
        mul_symbol_latex_numbers = settings.pop('mul_symbol_latex_numbers', r' \times ')
        if mul_symbol_latex_numbers is not None and not isinstance(mul_symbol_latex_numbers, str):
            raise TypeError
        # Default `\,` follows `siunitx`; Astropy uses `\;`
        unit_sep = settings.pop('unit_sep', r'\,')
        if not isinstance(unit_sep, str):
            raise TypeError

        super().__init__(settings)

        self._settings['separate_numerical_frac']  = separate_numerical_frac
        self._settings['float_fmtspec'] = float_fmtspec
        self._settings['symbol_replace'] = symbol_replace
        self._settings['symbol_replace_all'] = symbol_replace_all
        if mul_symbol_latex_numbers:
            self._settings['mul_symbol_latex_numbers'] = mul_symbol_latex_numbers
        self._settings['unit_sep'] = unit_sep
        self._wrapped: exprorder.BaseWrapped | None = None
        self._wrapped_lhs: exprorder.WrappedExpr | None = None
        self._wrapped_rhs: exprorder.WrappedExpr | None = None


    def doprint(self, expr: Basic | exprorder.BaseWrapped | Quantity | Unit) -> str:
        if isinstance(expr, exprorder.BaseWrapped):
            self._wrapped = expr
            if isinstance(expr, exprorder.WrappedEq):
                self._wrapped_lhs = expr.wrapped_lhs
                self._wrapped_rhs = expr.wrapped_rhs
            printed = super().doprint(expr.wrapped)
            self._wrapped = None
            self._wrapped_lhs = None
            self._wrapped_rhs = None
        else:
            printed = super().doprint(expr)
        return printed

    def _print(self, expr, **kwargs) -> str:
        if self._print_level == 1:
            if self._wrapped_lhs and self._wrapped_lhs.expr is expr:
                _wrapped = self._wrapped
                self._wrapped = self._wrapped_lhs
                try:
                    string = super()._print(expr, **kwargs)
                finally:
                    self._wrapped = _wrapped
                return string
            if self._wrapped_rhs and self._wrapped_rhs.expr is expr:
                _wrapped = self._wrapped
                self._wrapped = self._wrapped_rhs
                try:
                    string = super()._print(expr, **kwargs)
                finally:
                    self._wrapped = _wrapped
                return string
        return super()._print(expr, **kwargs)


    def _print_Float(self, expr: Float) -> str:
        if not self._settings['float_fmtspec']:
            return super()._print_Float(expr)
        str_real = f'{expr:{self._settings['float_fmtspec']}}'
        #### Begin code from `sympy/printing/latex.py`, from `_print_Float()`
        #### https://github.com/sympy/sympy/blob/37615e9bcac2e0938c6957ec3484a182d223e086/sympy/printing/latex.py
        # Must always have a mul symbol (as 2.5 10^{20} just looks odd)
        # thus we use the number separator
        separator = self._settings['mul_symbol_latex_numbers']

        if 'e' in str_real:
            (mant, exp) = str_real.split('e')

            if exp[0] == '+':
                exp = exp[1:]
            if self._settings['decimal_separator'] == 'comma':
                mant = mant.replace('.','{,}')

            return r"%s%s10^{%s}" % (mant, separator, exp)
        elif str_real == "+inf":
            return r"\infty"
        elif str_real == "-inf":
            return r"- \infty"
        else:
            if self._settings['decimal_separator'] == 'comma':
                str_real = str_real.replace('.','{,}')
            return str_real
        #### End code from `sympy/printing/latex.py`, from `_print_Float()`


    _dimensionless_unit = Unit()

    def _print_Unit(self, expr: Unit) -> str:
        if self._settings['mode'] == 'inline':
            return expr.to_string(format='latex_inline').strip('$')
        return expr.to_string(format='latex').replace(r'\frac', r'\tfrac').strip('$')

    _print_CompositeUnit = _print_Unit

    def _print_Quantity(self, expr: Quantity) -> str:
        # The current implementation uses SymPy's `Float` to handle numbers,
        # so that number handling is uniform.  Possible alternative
        # implementation that leaves number handling to Astropy's `Quantity`:
        # --------------------------------------------------------------------
        # if self._settings['mode'] == 'inline':
        #     return expr.to_string(format='latex', subfmt='inline', precision=precision).strip("$")
        # string = expr.to_string(format='latex', subfmt='inline', precision=precision).strip("$")
        # if string.startswith(r'\displaystyle'):
        #     string = string.replace(r'\displaystyle', '').lstrip()
        # return string
        # --------------------------------------------------------------------
        value = self._print(float(expr.value))
        if expr.unit == self._dimensionless_unit:
            return value
        unit = self._print_Unit(expr.unit)
        return rf'{value} {self._settings['unit_sep']} {unit}'


    def _print_Symbol(self, expr: Symbol) -> str:
        if self._settings['symbol_replace']:
            try:
                replacement = self._settings['symbol_replace'][expr]
            except KeyError:
                pass
            else:
                return rf'\left({self._print(replacement)}\right)'
        if self._settings['symbol_replace_all']:
            try:
                replacement = self._settings['symbol_replace_all'][expr]
            except KeyError:
                raise KeyError(f'"symbol_replace_all" is missing replacement for "{expr}"')
            else:
                return rf'\left({self._print(replacement)}\right)'
        return super()._print_Symbol(expr)


    #### Function from `sympy/printing/latex.py`
    #### https://github.com/sympy/sympy/blob/37615e9bcac2e0938c6957ec3484a182d223e086/sympy/printing/latex.py
    def _print_Add(self, expr, order=None):
        terms = self._as_ordered_terms(expr, order=order)
        #### Begin modification
        if self._wrapped:
            terms = self._wrapped.expr_order_collection.sort_terms(terms, expr)
        #### End modification

        tex = ""
        for i, term in enumerate(terms):
            if i == 0:
                pass
            elif term.could_extract_minus_sign():
                tex += " - "
                term = -term
            else:
                tex += " + "
            term_tex = self._print(term)
            if self._needs_add_brackets(term):
                term_tex = r"\left(%s\right)" % term_tex
            tex += term_tex

        return tex
    #### End function from `sympy/printing/latex.py`


    #### Function from `sympy/printing/latex.py`
    #### https://github.com/sympy/sympy/blob/37615e9bcac2e0938c6957ec3484a182d223e086/sympy/printing/latex.py
    def _print_Mul(self, expr: Expr):
        from sympy.simplify import fraction
        separator: str = self._settings['mul_symbol_latex']
        numbersep: str = self._settings['mul_symbol_latex_numbers']

        def convert(expr) -> str:
            if not expr.is_Mul:
                return str(self._print(expr))
            else:
                if self.order not in ('old', 'none'):
                    args = expr.as_ordered_factors()
                else:
                    args = list(expr.args)

                # If there are quantities or prefixes, append them at the back.
                units, nonunits = sift(args, lambda x: (hasattr(x, "_scale_factor") or hasattr(x, "is_physical_constant")) or
                              (isinstance(x, Pow) and
                               hasattr(x.base, "is_physical_constant")), binary=True)
                prefixes, units = sift(units, lambda x: hasattr(x, "_scale_factor"), binary=True)
                #### Begin modification
                #### # return convert_args(nonunits + prefixes + units)
                if self._wrapped is None:
                    return convert_args(nonunits + prefixes + units)
                return convert_args(
                    list(self._wrapped.expr_order_collection.sort_factors(nonunits, expr)) + prefixes + units
                )
                #### End modification

        def convert_args(args) -> str:
            _tex = last_term_tex = ""

            for i, term in enumerate(args):
                term_tex = self._print(term)
                if not (hasattr(term, "_scale_factor") or hasattr(term, "is_physical_constant")):
                    if self._needs_mul_brackets(term, first=(i == 0),
                                                last=(i == len(args) - 1)):
                        term_tex = r"\left(%s\right)" % term_tex

                    if  _between_two_numbers_p[0].search(last_term_tex) and \
                        _between_two_numbers_p[1].match(term_tex):
                        # between two numbers
                        _tex += numbersep
                    elif _tex:
                        _tex += separator
                elif _tex:
                    _tex += separator

                _tex += term_tex
                last_term_tex = term_tex
            return _tex

        # Check for unevaluated Mul. In this case we need to make sure the
        # identities are visible, multiple Rational factors are not combined
        # etc so we display in a straight-forward form that fully preserves all
        # args and their order.
        # XXX: _print_Pow calls this routine with instances of Pow...
        if isinstance(expr, Mul):
            args = expr.args
            if args[0] is S.One or any(isinstance(arg, Number) for arg in args[1:]):
                #### Begin modification
                #### # return convert_args(args)
                if self._wrapped is None:
                    return convert_args(args)
                return convert_args(self._wrapped.expr_order_collection.sort_factors(args, expr))
                #### End modification

        include_parens = False
        if expr.could_extract_minus_sign():
            expr = -expr
            tex = "- "
            if expr.is_Add:
                tex += "("
                include_parens = True
        else:
            tex = ""

        numer, denom = fraction(expr, exact=True)

        if denom is S.One and Pow(1, -1, evaluate=False) not in expr.args:
            # use the original expression here, since fraction() may have
            # altered it when producing numer and denom
            tex += convert(expr)

        else:
            snumer = convert(numer)
            sdenom = convert(denom)
            ldenom = len(sdenom.split())
            ratio = self._settings['long_frac_ratio']
            if self._settings['fold_short_frac'] and ldenom <= 2 and \
                    "^" not in sdenom:
                # handle short fractions
                if self._needs_mul_brackets(numer, last=False):
                    tex += r"\left(%s\right) / %s" % (snumer, sdenom)
                else:
                    tex += r"%s / %s" % (snumer, sdenom)
            elif ratio is not None and \
                    len(snumer.split()) > ratio*ldenom:
                # handle long fractions
                if self._needs_mul_brackets(numer, last=True):
                    tex += r"\frac{1}{%s}%s\left(%s\right)" \
                        % (sdenom, separator, snumer)
                elif numer.is_Mul:
                    # split a long numerator
                    a = S.One
                    b = S.One
                    for x in numer.args:
                        if self._needs_mul_brackets(x, last=False) or \
                                len(convert(a*x).split()) > ratio*ldenom or \
                                (b.is_commutative is x.is_commutative is False):
                            b *= x
                        else:
                            a *= x
                    if self._needs_mul_brackets(b, last=True):
                        tex += r"\frac{%s}{%s}%s\left(%s\right)" \
                            % (convert(a), sdenom, separator, convert(b))
                    else:
                        tex += r"\frac{%s}{%s}%s%s" \
                            % (convert(a), sdenom, separator, convert(b))
                else:
                    tex += r"\frac{1}{%s}%s%s" % (sdenom, separator, snumer)
            else:
                #### Begin modification
                #### # tex += r"\frac{%s}{%s}" % (snumer, sdenom)
                separate_numerical_frac = self._settings['separate_numerical_frac']
                if separate_numerical_frac and denom.is_Integer:
                    if self._settings['mode'] == 'inline':
                        frac = r"\tfrac"
                    else:
                        frac = r"\frac"
                    if numer.is_Integer:
                        tex += r"%s{%s}{%s}" % (frac, snumer, sdenom)
                    elif numer.is_Mul and numer.args[0].is_Integer:
                        snumer_num, snumer_rest = snumer.split(maxsplit=1)
                        if self._needs_mul_brackets(numer, last=False):
                            tex += r"%s{%s}{%s} \left(%s\right)" % (frac, snumer_num, sdenom, snumer_rest)
                        else:
                            tex += r"%s{%s}{%s} %s" % (frac, snumer_num, sdenom, snumer_rest)
                    else:
                        if self._needs_mul_brackets(numer, last=False):
                            tex += r"%s{1}{%s} \left(%s\right)" % (frac, sdenom, snumer)
                        else:
                            tex += r"%s{1}{%s} %s" % (frac, sdenom, snumer)
                else:
                    tex += r"\frac{%s}{%s}" % (snumer, sdenom)
                #### End modification

        if include_parens:
            tex += ")"
        return tex
    #### End function from `sympy/printing/latex.py`


def latex(expr: Basic | exprorder.BaseWrapped | Quantity | Unit, **settings) -> str:
    #settings.setdefault('mul_symbol_latex_numbers', r'\times')
    return LatexPrinter(settings).doprint(expr)
