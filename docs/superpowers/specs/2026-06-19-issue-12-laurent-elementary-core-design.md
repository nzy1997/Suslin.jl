# Issue 12 Laurent Elementary Core Design

## Goal

Lift the existing elementary matrix, Cohn-type realization, conjugated
elementary realization, and factorization verification routines so they work
over Oscar Laurent polynomial parents with exact parent checks and unchanged
ordinary polynomial behavior.

## Context

The repository already provides `suslin_laurent_polynomial_ring` and internal
Laurent parent validators. `verify_factorization` and the elementary routines
are public API and currently operate through generic Oscar matrices and
`_coerce_into_ring`. The requested scope is limited to the elementary/Cohn/
normality layer; determinant normalization and the full Laurent factorization
driver remain out of scope.

No `AGENTS.md`, `CLAUDE.md`, or `GEMINI.md` file is present in this checkout.
The README documents the package test entry point and the full-suite command.

## Approaches Considered

1. Reuse the existing generic routines and harden their Laurent parent checks.
   This is the recommended approach because it preserves the public API,
   avoids parallel Laurent-only entry points, and matches the issue guidance.
2. Add separate Laurent-specific realization functions. This would isolate the
   new behavior but would expand API surface and duplicate algebra that should
   remain shared.
3. Broaden the full `elementary_factorization` driver. This is too large for
   this issue and overlaps with later Laurent roadmap work.

The chosen design is approach 1.

## Design

`elementary_matrix(n, i, j, a, R)` remains the single constructor for ordinary
polynomial and Laurent parents. It must coerce `a` into `R`, reject diagonal
elementary requests, and return an exact `n x n` identity matrix with only the
requested off-diagonal entry changed.

`realize_cohn_type(n, i, j, a, v, R)` remains parent-parametric. It validates
dimension and indexing as before, coerces every vector entry and scalar into
`R`, and builds the existing Cohn commutator factor sequence from those
coerced entries. For Laurent inputs, every factor must have base ring exactly
`R`.

`realize_conjugate_elementary(B, i, j, a)` derives `R = base_ring(B)` and
continues to use the Cohn-type routine for all nonzero terms in the normality
decomposition. Oscar's generic `inv(B)` path is not implemented for Laurent
matrix parents in the supported dependency set, so the Laurent path needs a
small exact inverse helper based on determinant, minors, and the adjugate while
leaving ordinary polynomial matrices on the existing `inv(B)` path. For
Laurent matrices, the returned factors must multiply back to
`B * elementary_matrix(n, i, j, a, R) * inverse(B)` exactly.

`verify_factorization(A, factors)` keeps exact multiplication semantics. It
must throw `ArgumentError` when any factor size differs from `A` or when a
factor's base ring is not equal or identical to `base_ring(A)`. This is the
negative control for distinct Laurent parents.

## Tests

Add `test/expert/laurent_elementary_core.jl`, modeled on the existing
`cohn_type` and `normality` expert tests and using
`suslin_laurent_polynomial_ring`. The test file covers:

- Laurent elementary matrix construction, including a negative exponent
  off-diagonal entry.
- Laurent Cohn-type reconstruction for `n = 3` and `n = 4`.
- Laurent conjugated elementary reconstruction of `B * E * inv(B)`.
- `verify_factorization` returning true for the Laurent Cohn and conjugated
  factor lists.
- `verify_factorization` throwing `ArgumentError` when a factor belongs to a
  different Laurent parent than the target matrix.

Register the new expert file in `test/runtests.jl` so the documented
expert/full-suite commands exercise it.

## Verification

Run these commands:

```bash
julia --project=. -e 'include("test/expert/laurent_elementary_core.jl")'
julia --project=. -e 'using Pkg; Pkg.test()'
julia --project=. test/runtests.jl all
```

The first command is the issue-specific check. The second is the package test
entry point required by Agent Desk. The third is the documented full suite from
issue #21.

## Scope Boundaries

This design does not add determinant normalization, does not enable the full
Laurent `elementary_factorization` driver, and does not change public API names.
