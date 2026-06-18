# Issue 6 Laurent Fixture Catalog Design

## Context

Suslin now has a Laurent ring constructor and a first ToricBuilder contract
fixture. Later Laurent algebra issues need a shared, exact, small fixture
catalog instead of copying matrices and vectors into each test file. The
catalog is test support, not a public Suslin API.

Issue 19 provides the first real ToricBuilder boundary: exact GF(2)
two-variable Laurent matrices from `factor_toric_block(3, x, y, R)` with
checked inverse-style relations. Issue 21 defines the package entry point
`julia --project=. -e 'using Pkg; Pkg.test()'` and the full-suite command
`julia --project=. test/runtests.jl all`.

## Design Choice

Add a dedicated shared catalog at `test/fixtures/laurent_cases.jl` and a
validator test at `test/internal/laurent_fixtures.jl`.

The catalog returns structured named tuples. Each fixture entry records:

- `id`
- `kind`
- `ring_constructor`
- `ring`
- `dimensions`
- `inputs`
- `expected_relation`
- `provenance`
- `determinant_profile`
- `consumer_test_ids`

The validator lives in internal tests because the catalog is not a supported
library interface. It checks each entry by exact Oscar arithmetic and includes
negative controls for wrong expected relations, determinant classifications,
and ToricBuilder provenance-derived relations.

Alternatives considered:

- Add a public Suslin fixture API. This is too broad for a test catalog issue.
- Extend the existing ToricBuilder fixture only. That would not serve synthetic
  Laurent tests that need reusable solvable, unsolvable, and normalization
  cases.
- Duplicate ToricBuilder matrices into the shared catalog. Reusing the existing
  Issue 19 fixture is less error-prone and keeps provenance in one place.

## Fixture Scope

The initial catalog stays intentionally small:

- A solvable Laurent linear-system case with a triangular `2 x 2` matrix,
  column vector right-hand side, and exact expected solution.
- An unsolvable linear-system case with a zero `1 x 1` matrix and a nonzero
  right-hand side, witnessed by a simple zero-matrix/nonzero-rhs certificate.
- A negative-exponent normalization case that records a Laurent vector, a
  monomial multiplier, and the expected nonnegative-exponent vector.
- A ToricBuilder-motivated relation case using the smaller `Pinv` entry from
  the Issue 19 `factor_toric_block(3, x, y, R)` fixture.

## Validation Rules

Every fixture must have the required metadata fields and at least one consumer
test id.

Linear-system fixtures must have matrix/vector dimensions that agree. Solvable
fixtures must satisfy `A * expected_solution == rhs`. Unsolvable fixtures must
carry the supported zero-matrix/nonzero-rhs certificate and the validator must
check that certificate directly.

Negative-exponent normalization fixtures must satisfy
`normalization_unit * input_vector == normalized_vector`, and the recorded
normalized vector must not include the original negative exponent entries.

Matrix-relation fixtures must verify the claimed exact multiplication relation.
The ToricBuilder case verifies `source_matrix * matrix == identity_matrix(R, n)`.

Determinant profiles are checked when a matrix is present and the fixture marks
the profile as relevant. The initial classes are `one`, `laurent_monomial_unit`,
and `non-unit`.

## Files

- Create `test/fixtures/laurent_cases.jl`: loadable catalog module that builds
  exact Oscar/Suslin Laurent fixtures.
- Create `test/internal/laurent_fixtures.jl`: validator, focused tests, and
  negative controls.
- Modify `test/runtests.jl`: include the new internal validator in the default
  internal group.

## Verification

Focused validator:

```bash
julia --project=. -e 'include("test/internal/laurent_fixtures.jl")'
```

Package entry point:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Documented full suite from Issue 21:

```bash
julia --project=. test/runtests.jl all
```

## Spec Self-Review

- No incomplete markers remain.
- Scope is limited to test fixtures and internal validation.
- The ToricBuilder case is shaped by the Issue 19 contract.
- Full-suite verification uses the Issue 21 command.
