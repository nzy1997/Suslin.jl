# Issue 7 Laurent Ring Design

## Context

Suslin currently exposes `suslin_polynomial_ring(F, names)` as a thin wrapper
around Oscar's polynomial-ring constructor. The Laurent roadmap needs a matching
entry point for exact Laurent polynomial rings before later normalization,
elementary-matrix, factorization-driver, and fixture-validation issues can share
strict parent checks.

Issue #7 depends on the #21 test contract, which documents the default package
test entry point and the full-suite command. Sibling issues #6, #8, #12, #13,
#14, and #17 repeatedly require Laurent inputs to be checked for a single common
parent and rejected when mixed with ordinary polynomial elements.

## Approach Options

1. Add a minimal constructor and internal validator helpers in
   `src/core/rings.jl`. This keeps the public surface narrow while giving later
   code reusable parent checks. This is the selected approach.
2. Export all validator helpers immediately. This would make downstream tests
   convenient but commits the package to helper names before any public Laurent
   workflow has exercised them.
3. Create a new Laurent-focused core file. This is unnecessary for the first
   issue because the existing ring module is small and already owns ring
   constructors.

## Selected Design

Add `suslin_laurent_polynomial_ring(F, names::Vector{String})` next to the
ordinary polynomial constructor. It will call `Oscar.laurent_polynomial_ring(F,
names)` and return `(R, collect(vars))`, matching the existing package pattern
of returning the parent ring plus a concrete generator vector. Oscar's
names-vector Laurent constructor returns a multivariate Laurent parent, including
for one variable, which is acceptable for this roadmap because later ToricBuilder
inputs are multivariate Laurent-polynomial matrices.

Export only `suslin_laurent_polynomial_ring`. Keep validators internal with
underscore-prefixed names so later Suslin internals can reuse them without
prematurely defining public validator API names.

## Validator Contract

Add these helpers in `src/core/rings.jl`:

- `_is_laurent_polynomial_ring(R)`: returns `true` for Oscar/AbstractAlgebra
  univariate or multivariate Laurent polynomial parents and `false` otherwise.
- `_require_laurent_polynomial_ring(R; label="ring")`: returns `R` when it is a
  Laurent polynomial parent and throws `ArgumentError` otherwise.
- `_require_laurent_element(value; label="value")`: returns `value` when its
  parent is a Laurent polynomial ring and throws `ArgumentError` otherwise.
- `_require_laurent_element(value, R; label="value")`: also checks that
  `parent(value) === R`, rejecting mixed Laurent parents even when the variable
  names and coefficient field match.
- `_require_same_laurent_parent(values; label="values")`: validates a nonempty
  iterable of Laurent elements and returns the shared Laurent parent.

The helpers intentionally use identity (`===`) for parent matching. Later exact
algebra should not silently coerce between different cached or non-cached
parents.

## Tests

Add `test/internal/laurent_rings.jl` and register it in the `internal` group in
`test/runtests.jl`.

The focused internal tests must verify:

- the constructor builds a Laurent parent with the requested generator names,
- returned generators belong to the reported parent,
- `x^-1` and `x^-1*y` can be constructed and stay in the same parent,
- the Laurent ring validator accepts the constructed parent and rejects an
  ordinary polynomial parent,
- element validators accept Laurent values from the expected parent,
- passing an element from a different Laurent parent throws `ArgumentError`,
- passing an ordinary polynomial element where a Laurent element is required
  throws `ArgumentError`,
- collection validation rejects mixed-parent Laurent values.

Update `test/public/api_surface.jl` only to cover the new exported constructor.
Do not change the ordinary polynomial constructor or existing algorithm tests.

## Verification

Run the issue-specific command:

```bash
julia --project=. -e 'include(test/internal/laurent_rings.jl)'
```

Then run the package entry point required by the Agent Desk task:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Then run the full suite documented by #21:

```bash
julia --project=. test/runtests.jl all
```

## Automatic Approval

This Agent Desk run is non-interactive. Under the standing answer policy, the
design and written spec are approved automatically because the selected approach
is the narrowest one that satisfies issue #7, follows existing repository
patterns, and avoids committing validator helper names as public API before the
later Laurent issues need that surface.
