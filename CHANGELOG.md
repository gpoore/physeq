# Changelog


## v0.3.0 (dev)

 * `Symbol` now has constraint methods for greater/less than (or equal to) a
   specified value:  `constrain_gt()`, `constrain_ge()`, `constrain_lt()`,
   `constrain_le(<value>)`.



## v0.2.0 (2026-03-26)

*  Added kinematics equations.

*  In `Problem`, more robust generation of values from `setup` and generation
   of solutions from `knowns`.

   Solutions are now generated from known values only and do not rely on setup
   values once known values are generated.  There are now errors when known
   values alone (without setup) cannot produce a single solution (for example,
   quadratics with insufficient constraints).

   When solutions are being generated, an error is no longer raised
   immediately if a given equation yields multiple solutions.  Instead, the
   remaining equations are tried first, in case they have sufficient
   constraints.

   Definitions equal to zero are now substituted all at once in solutions,
   instead of one per line.

   Fixed bug that caused a duplicate line in solutions when an equation was
   already solved for the desired symbol.

   Fixed bugs that could prevent `Problem` generation for certain combinations
   of `setup`, `knowns`, and `unknowns`.  Improved error messages for problem
   generation.

*  In `Problem`, added support for specifying constraints for `knowns` and
   `unknowns`.  These are now displayed in `Problem.simple_solutions()`.

*  In `Problem`, constants no longer have to be specified explicitly in
   `setup` or `knowns`.  They are automatically detected in `equations` and
   `definitions`, and then automatically added to the collection of known
   symbols.

*  Fixed bugs that could prevent correct term order for addition/subtraction.

*  `translate_xreplace_rule()` and `translate_numerical_xreplace_rule()`
   are now compatible with mappings from `ConstSymbol` -> `Quantity`.

*  Renamed `solveset_with_checked_assumptions()` to `checked_solveset()`,
   to reflect constraint checking capabilities.


## v0.1.0 (2026-03-23)

*  Initial release.
