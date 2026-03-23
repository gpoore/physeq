# -*- coding: utf-8 -*-
#
# Copyright (c) 2026, Geoffrey M. Poore
# All rights reserved.
#
# Licensed under the BSD 3-Clause License:
# http://opensource.org/licenses/BSD-3-Clause
#


from physeq.symbols.constants import g
from physeq.symbols.mass import m
from physeq.symbols.space import d
from physeq.symbols.kinematics import v, v_0
from physeq.symbols.force import μ_k, F_k, F_N
from physeq.force import eq_kinetic_friction
from physeq.symbols.work_energy import W_nc, ΔPE

from physeq.work_energy import eq_nonconservative_work, eq_ΔKE, eq_KE, eq_KE_0
from physeq.cas.wrapped import Eq

from physeq.cas.problem import Problem
from physeq.cas.printing import latex


problem = Problem(
    setup={
        m: m.randint_quantity(15, 25),
        v_0: v_0.uniform_quantity(5, 10),
        v: v.quantity(0),
        μ_k: μ_k.uniform_quantity(0.1, 0.4),
        g: g,
    },
    knowns=[m, v_0, v, μ_k, g],
    unknowns=[d],
    equations=[
        eq_nonconservative_work,
    ],
    definitions=[
        Eq(ΔPE, 0),
        eq_ΔKE,
        eq_KE,
        eq_KE_0,
        Eq(W_nc, -F_k*d),
        eq_kinetic_friction,
        Eq(F_N, m*g),
    ],
)

print(f'''\
A mass ${latex(m)}$ is initially traveling at speed ${latex(v_0)}$ across a
flat surface with coefficient of kinetic friction ${latex(μ_k)}$.  It comes to
rest after traveling a distance ${latex(d)}$.
''')
print(problem.simple_solutions())
