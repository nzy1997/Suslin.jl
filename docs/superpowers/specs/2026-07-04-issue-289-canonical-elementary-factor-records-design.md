# Issue 289 Canonical Elementary Factor Records Design

Issue #289 adds a private canonical record layer for ordinary-polynomial
elementary matrix factors used by upcoming Steinberg rewrites in #188. The
helpers let rewrite code inspect either an identity factor or a single
elementary factor `E_ij(a)`, including the zero-coefficient case
`E_ij(0) = I`, without exporting a public API or changing factorization routes.

## Context

There is no repository `AGENTS.md`; the README, existing Superpowers docs, and
the elementary-matrix core patterns apply. `gh issue view 289` failed in this
Agent Desk sandbox because the configured GitHub proxy is not reachable, so the
supplied issue body and the merged #288 fixture catalog are the source of
truth. `gh pr view 302` returned Codecov context for the #288 dependency, and
the local `main` merge commit includes the Steinberg optimization catalog.

The existing `_elementary_factor_data` helper in
`src/algorithm/laurent_column_peel.jl` accepts only nonzero elementary factors.
That is correct for Laurent column-peel inversion, but too narrow for Steinberg
identity removal and verifier-friendly logs.

## Approach Options

Recommended: add private helpers in `src/core/elementary_matrices.jl` returning
named-tuple records with `kind = :identity` or `kind = :elementary`. This keeps
the canonical layer close to `elementary_matrix`, reuses `_same_base_ring` and
shape checks, and avoids committing to public struct names before #188 consumes
the layer.

Alternative: widen `_elementary_factor_data` in `laurent_column_peel.jl`. That
would couple Steinberg rewrite needs to a Laurent-specific algorithm helper and
would require all current callers to handle identity records they do not need.

Alternative: export a public record type now. That is premature because the
issue explicitly asks to keep helpers internal and #188 has not yet finalized
the optimizer-facing API.

## Chosen Design

Add `_canonical_elementary_factor_record(factor)` as a private helper. It
requires a square matrix, records `n = nrows(factor)` and `ring = base_ring(factor)`,
requires every diagonal entry to equal `one(ring)`, and scans off-diagonal
entries.

If every off-diagonal entry is zero, it returns:

```julia
(; kind = :identity, n, ring)
```

If exactly one off-diagonal entry is nonzero, it returns:

```julia
(; kind = :elementary, n, ring, row, col, coefficient)
```

Matrices with non-one diagonal entries or more than one nonzero off-diagonal
entry throw `ArgumentError`. Nonsquare matrices throw `DimensionMismatch` via
the shared `_require_square_matrix` shape helper.

Add `_elementary_factor_record_matrix(record)` as the private round-trip
constructor. Identity records reconstruct `identity_matrix(record.ring, record.n)`.
Elementary records reconstruct `elementary_matrix(record.n, record.row,
record.col, record.coefficient, record.ring)`, so the coefficient is coerced by
the same path as ordinary elementary matrix construction.

The helper accepts ordinary-polynomial matrix factors used by #188. It does
not special-case Laurent rings, does not inspect factor sequences, and does not
replace `_elementary_factor_data`.

## Tests

Create `test/expert/steinberg_factor_count_optimization.jl` with a focused
testset for the new private helpers:

- nonzero `E_12(a)` round-trips as `kind = :elementary`;
- `elementary_matrix(n, i, j, zero(R), R)` is recognized as an identity record;
- `_elementary_factor_record_matrix(record)` reconstructs the original matrix
  for each accepted record;
- nonsquare matrices, matrices with a non-one diagonal entry, and matrices with
  two nonzero off-diagonal entries throw before any rewrite can consume them.

Register the file in the expert group of `test/runtests.jl`.

## Verification

Focused verification:

```bash
julia --project=. -e 'include("test/expert/steinberg_factor_count_optimization.jl")'
```

Full verification:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

The focused command must exit 0. The default package suite must remain green.

## Out Of Scope

Do not optimize factor sequences. Do not export a public API. Do not change
current factorization routes. Do not make the Laurent column-peel helper accept
identity factors.

## Self-Review

This design has no placeholders, keeps the feature internal, and limits code
changes to the elementary-matrix core plus focused expert tests. It explicitly
covers the issue's negative controls and preserves existing factorization and
Laurent inversion behavior.
