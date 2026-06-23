# Issue 109 Park-Woodburn Polynomial Fixture Catalog Design

## Context

Issue 64 needs ordinary-polynomial `SL_n` examples before the public
Park-Woodburn driver starts routing through local `SL_3`, block-local
`SL_n`, recursive column peel, and Quillen patch layers. The repository already
uses checked-in Julia fixture catalogs for Laurent boundaries, ECP columns,
Murthy-Gupta local data, and Quillen patch data. Issue 109 should add the same
kind of catalog for ordinary polynomial driver cases without changing
`elementary_factorization`.

The relevant current examples are:

- `test/public/factorization_driver_shell.jl` for the supported univariate
  `SL_3` fast local case and existing staged failures.
- `test/expert/sln_to_sl3_reduction.jl` for supported univariate `n > 3`
  disjoint local-block reductions.
- `test/fixtures/quillen_patch_cases.jl` for #99 Quillen fixture ids that
  later multivariate routing should reuse.

## Design Choice

Add a checked-in fixture catalog plus an internal validator.

The fixture file will define ordinary polynomial rings and return named entries
with literal matrices, expected route tags, current staged/pass status,
dependency issue ids, provenance, and negative-control metadata. The validator
will include the fixture file and check exact shape, determinant-one claims,
ordinary field-backed ring metadata, expected route/status compatibility, and
that catalog negative controls fail validation.

Alternative considered: put the examples directly in future public driver
tests. That would make issue 109 less reusable and would not give later route
issues a single shared source of truth.

Alternative considered: add package APIs for loading these fixtures. That is
out of scope because the catalog is test support, not a public Suslin surface.

## Catalog Scope

The positive catalog entries are:

- `pw-poly-univariate-sl3-fast-local-qq`: `QQ[X]`, `3 x 3`, supported now by
  the local `SL_3` path.
- `pw-poly-univariate-sln-disjoint-blocks-qq`: `QQ[X]`, `6 x 6`, supported
  now by the disjoint local-block `SL_n` reduction from the expert test.
- `pw-poly-recursive-column-peel-gf2`: `GF(2)[x, y]`, `4 x 4`, determinant
  one, staged until #113.
- A multivariate Quillen-backed case reusing the #99 id
  `quillen-patched-substitution-witness-qq`, over `QQ[X, r, g]`, staged until
  #105 and #115.

The negative controls are:

- A determinant-not-one ordinary polynomial matrix.
- A determinant-one matrix outside the implemented witness families.
- A known matrix whose metadata claims the wrong route.

## Validation Rules

The internal validator will require:

- Every entry has a non-empty id, role, route, expected status, ordinary
  polynomial ring constructor metadata, ring metadata, matrix, determinant
  expectation, dependency issue ids, source refs, and consumer issue ids.
- Every positive matrix is square of size at least 3 and lives over its recorded
  ring.
- Positive entries with determinant expectation `:one` have
  `det(matrix) == one(R)`.
- Supported statuses are only allowed for current implemented witness families:
  `:fast_local_sl3` and `:disjoint_local_blocks`.
- The Quillen entry must list #105 as a blocking dependency and must not claim
  current support.
- Negative controls are part of the catalog and each must be rejected by the
  same validator.

## Files

- Create `test/fixtures/park_woodburn_polynomial_cases.jl`: fixture catalog
  module and helpers such as `catalog()` and `cases_by_id()`.
- Create `test/internal/park_woodburn_polynomial_fixtures.jl`: validators,
  required id checks, and mutation/negative-control tests.
- Modify `test/runtests.jl`: register the validator in the `internal` group.

## Verification

Focused command:

```bash
julia --project=. -e 'include("test/internal/park_woodburn_polynomial_fixtures.jl")'
```

Required package command:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

The focused command must validate at least four positive entries and three
negative controls, report the multivariate Quillen fixture as blocked by #105,
and prove metadata is enforced by rejecting a determinant-one wrong-route
mutation and a determinant-not-one `:supported` mutation.

## Spec Self-Review

- No public API or driver behavior changes are included.
- The design reuses #99 Quillen fixture ids rather than inventing a second
  Quillen schema.
- The staged route/status model is explicit enough for #64, #113, #115, and
  final acceptance tests to consume later.
- Negative controls validate both algebraic facts and metadata claims.
