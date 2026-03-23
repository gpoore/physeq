# PhysEq – Generate physics problems with realistic values and step-by-step solutions


PhysEq simplifies the process of creating realistic physics problems with
step-by-step solutions at the high school and introductory college levels.
You specify the setup and the relevant equations, and PhysEq solves for
unknowns and creates solutions.  There is built-in support for randomizing
values, with safeguards against unphysical values.

Currently, there is not built-in support for calculus.

## Features

  * Symbols, equations, and algebra for solutions are based on
    [SymPy](https://www.sympy.org/).

  * Equations and algebra track expression order.  When you define equations,
    you are also defining your preferred order for symbols, terms, and
    factors.  For example, if you define an equation as $F=ma$, then it will
    appear as $F=ma$ in solutions.  The appearance of equations is not based
    on the ordering used internally within the computer algebra system
    (similar to alphabetical order), so $F=ma$ is **not** transformed into
    $F=am$ in solutions.  This makes solutions much easier to read.

  * Physical values are represented with `Quantity` from
    [Astropy](https://docs.astropy.org/en/stable/api/astropy.units.Quantity.html).
    `Quantity` combines a numerical value with a unit and has built-in
    support for unit conversions.  PhysEq checks inputs to ensure that they
    have appropriate units.  Inputs are automatically converted into standard
    SI units before numerical calculations.

  * When numerical values are substituted into symbolic equations, values are
    checked against mathematical assumptions and any additional value
    constraints that have been defined.  For example, masses and vector
    magnitudes cannot be negative, so attempting to replace them with negative
    numbers results in an error.  Similarly, functions for generating
    randomized quantities raise an error for unphysical values.


## Examples

Examples are under `examples/` in the project repository.  Here is example
output for a friction problem:

*A mass $m$ is initially traveling at speed $v_{0}$ across a
flat surface with coefficient of kinetic friction $\mu_{k}$.  It comes to
rest after traveling a distance $d$.*

**Known**

*   $g = 9.80 \, \mathrm{\tfrac{m}{s^{2}}}\quad v = 0 \, \mathrm{\tfrac{m}{s}}\quad m = 22.0 \, \mathrm{kg}\quad \mu_{k} = 0.154\quad v_{0} = 8.10 \, \mathrm{\tfrac{m}{s}}$

**Unknown**

*   $d$

**Solutions**

*  $\displaystyle W_{\text{nc}} = \Delta \text{KE} + \Delta \text{PE}$

    $\displaystyle W_{\text{nc}} = \Delta \text{KE}$

    $\displaystyle - F_{k} d = \text{KE} - \text{KE}_{0}$

    $\displaystyle - \mu_{k} F_{N} d = \frac{1}{2} m v^{2} - \frac{1}{2} m v_{0}^{2}$

    $\displaystyle - \mu_{k} m g d = \frac{1}{2} m v^{2} - \frac{1}{2} m v_{0}^{2}$

    $\displaystyle d = \frac{- \frac{1}{2} m v^{2} + \frac{1}{2} m v_{0}^{2}}{\mu_{k} m g}$

    $\displaystyle d = \frac{- v^{2} + v_{0}^{2}}{2 \mu_{k} g}$

    $d = \frac{- \left(0\right)^{2} + \left(8.10\right)^{2}}{2 \left(0.154\right) \left(9.80\right)}$

    $\displaystyle d = 21.7 \, \mathrm{m}$



## Creating problems

The examples under `examples/` in the project repository may be a good place
to start in creating problems.

Problems are created with the `Problem` class:
```
problem = Problem(
    setup: dict[<symbols>, <numbers or quantities>],
    knowns: list[<symbols>],
    unknowns: list[<symbols>],
    equations: list[<equation>],
    definitions: list[<equation>] | None = None,
    simplify: bool = False,
)
```
`setup` provides a mapping of symbols to numerical values that allows a
complete problem to be generated.  It may include symbols that will ultimately
be used as unknowns.  The goal at this stage is to generate a complete set of
values, regardless of which ones will ultimately treated as unknowns.  In many
cases, it is simpler to create a physically realistic problem by starting with
quantities such as velocity and mass, even if one or more of these are
ultimately treated as unknowns.

`setup` values may be numbers, in which case they are interpreted as values in
the standard SI units for the corresponding symbols.  Alternatively, values may
be provided as Astropy `Quantity`, in which case they are automatically
converted into standard SI units before use.

`knowns` are the symbols whose values will be given in the problem.
`unknowns` are the symbols that must be solved for.

`equations` are the primary equations to be used in solving.  This list would
typically include fundamental equations such as $\sum F_x=ma_x$ or
$W_{\text{nc}} = \Delta \text{KE} + \Delta \text{PE}$, plus important
problem-specific relationships.  These equations are shown in detail in
step-by-step solutions.

Currently, there is only support for problems that can be solved as a sequence
of solving equations with one unknown.  For example, an equation might have
two unknowns, but it is possible to solve for one of these using a separate
equation with only one unknown, so ultimately the problem only involves
solving equations with one unknown.  Limited support for problems involving
systems of equations with multiple unknowns may be added in the future.

`definitions` are substituted into `equations` as needed, but typically would
not appear independently in step-by-step solutions.  They may be simple
equations that define terms in `equations` (such as $F_G=mg$) or equations
that provide information about initial or final conditions (such as $v_0=0$).
These equations are limited to the form `<symbol>=<expression>`, where
`<symbol>` appears in `equations`.

By default, step-by-step solutions use minimal simplifications, so that they
will be easier to follow and maintain a form closer to the starting equations.
This can be changed by setting `simplify=True`.

When a problem is created with `problem = Problem(...)`, there are three
relevant attributes for converting `problem` into a form that can be used in
teaching.

  * `problem.knowns: dict[<symbol>, <quantity>]`:  This maps known symbols to
    the quantities that are given.  Quantities are instances of
    [Astropy `Quantity`](https://docs.astropy.org/en/stable/api/astropy.units.Quantity.html),
    so they have both a numerical value and a unit.

  * `problem.unknowns: dict[<symbol>, <quantity>]`:  This maps unknown symbols
    to the quantities that should be found.

  * `problem.equations: dict[<equation>, list[<equation>]]`:  This maps
    equations to step-by-step solutions.  Solutions are represented as lists
    of equations.

There is also a method `problem.simple_solutions()` that creates simple
solutions in Markdown format.

`physeq.cas.printing` provides a `latex()` function for converting symbols,
equations, and quantities into strings in LaTeX format, compatible with LaTeX
documents and HTML libraries like MathJax and KaTeX.


### Customizing symbols and equations

PhysEq includes an expanding library of pre-defined physics symbols
(`physeq.symbols.<physics_category>`) and equations
(`physeq.<physics_category>`).

Some physical quantities do not have completely standardized representations
in introductory physics materials.  For example, kinetic energy commonly
appears as either $\text{KE}$ or $E_k$.  By default, PhysEq uses $\text{KE}$,
but this can be changed for all built-in kinetic energy symbols by setting
`physeq.config.kinetic_energy = 'E_k'` (or using another definition).  It is
also possible to customize `config.potential_energy` (default $\text{PE}$),
`config.coulombs_constant` (default $k$), and `config.spring_constant`
(default $k$).

Most symbols support subscripting.  For example, for a problem involving
multiple masses, additional mass symbols can be defined from the predefined
mass symbol `m` via `m_1 = m.subscript(1)`, `m_2 = m.subscript(2)`, and so
forth.  This is typically simpler than creating new symbols from scratch.

Additional custom symbols should be created using `physeq.cas.wrapped.Symbol`.
Additional equations should be created using `physeq.cas.wrapped.Eq`.  When
possible, SymPy representations of numbers (for example,
`Rational`) and functions (for example, `sqrt` and `cos`) should not be used
directly; the wrappers in `physeq.cas.wrapped` should be used instead.  These
`Symbol` and `Eq` classes along with the wrappers track expression order, so
that solutions will be easier to follow.  Without expression order tracking,
equations would appear with symbols sorted in something like alphabetical
order (for example, $F=ma$ would appear as $F=am$).


### Creating randomized values

Symbols include built-in support for randomization.  (This does not apply to
symbols representing physical constants.)  For example,
`m.randint_quantity(1, 10)` might return the Astropy `Quantity` 2 kg.

There are four methods for generating randomized quantities.  All of these
build on Python's
[`random` module](https://docs.python.org/3/library/random.html).  These
methods select a random number, and then create a `Quantity` by combining that
number with the symbol's standard SI unit.

  * `randrange_quantity(start, stop[, step])`:  Create a `Quantity` using a
    random integer from `range(start, stop, step)`.

  * `randint_quantity(first, last)`:  Create a `Quantity` using a random
    integer `N` selected so that `first <= N <= last`.

  * `choice_quantity(seq)`:  Randomly select a number from the sequence
    `seq`, and create a `Quantity` from it.

  * `uniform_quantity(a, b)`:  Randomly select a floating-point number
    between `a` and `b` (inclusive of both), and create a `Quantity` from it.

When quantities are randomly generated using these methods, the resulting
numerical values are checked against mathematical assumptions and any
additional value constraints that have been defined.  An error is raised for
unphysical values.
