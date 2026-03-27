# -*- coding: utf-8 -*-
#
# Copyright (c) 2026, Geoffrey M. Poore
# All rights reserved.
#
# Licensed under the BSD 3-Clause License:
# http://opensource.org/licenses/BSD-3-Clause
#

import random
from physeq.kinematics import Δy, v_y, v_0y, a_y, t_elapsed, g, eq_Δy_const_a, Eq, eq_v_y_const_a, eq_sq_v_y_const_a
from physeq.cas.problem import Problem
from physeq.cas.printing import latex


symbols = [Δy, v_0y, v_y, t_elapsed]
knowns = []
while len(knowns) < 2:
    choice = random.choice(symbols)
    if choice not in knowns:
        knowns.append(choice)
unknowns = [s for s in symbols if s not in knowns]
if v_0y in unknowns:
    unknowns = [u if u is not v_0y else v_0y.constrain_positive() for u in unknowns]
if v_y in unknowns:
    unknowns = [u if u is not v_y else v_y.constrain_nonnegative() for u in unknowns]

problem = Problem(
    setup={
        v_0y: v_0y.uniform_quantity(25, 45),
        t_elapsed: t_elapsed.uniform_quantity(0.5, 2.5),
    },
    knowns=knowns,
    unknowns=unknowns,
    equations=[
        eq_Δy_const_a,
        eq_v_y_const_a,
        eq_sq_v_y_const_a,
    ],
    definitions=[
        Eq(a_y, -g),
    ],
)

print(f'''\
A projectile is launched directly upward at initial velocity ${latex(v_0y)}$.
After time ${latex(t_elapsed)}$, it has a displacement ${latex(Δy)}$ relative
to the launch point.
''')
print(problem.simple_solutions())
