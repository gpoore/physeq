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

import sympy
from sympy.core.singleton import S
from sympy.printing.latex import LatexPrinter as SymPyLatexPrinter
from sympy.printing.latex import _between_two_numbers_p
from sympy.utilities.iterables import sift
from .ordering import Add, Mul




class LatexPrinter(SymPyLatexPrinter):
    def _print_Mul(self, expr: sympy.Expr):
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
                              (isinstance(x, sympy.Pow) and
                               hasattr(x.base, "is_physical_constant")), binary=True)
                prefixes, units = sift(units, lambda x: hasattr(x, "_scale_factor"), binary=True)
                return convert_args(Mul.sort_terms(nonunits) + prefixes + units)

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
        if isinstance(expr, sympy.Mul):
            args = expr.args
            if args[0] is S.One or any(isinstance(arg, sympy.Number) for arg in args[1:]):
                return convert_args(Mul.sort_terms(args))

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

        if denom is S.One and sympy.Pow(1, -1, evaluate=False) not in expr.args:
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
                tex += r"\frac{%s}{%s}" % (snumer, sdenom)

        if include_parens:
            tex += ")"
        return tex


def latex(expr: sympy.Expr, **settings):
    return LatexPrinter(**settings).doprint(expr)
