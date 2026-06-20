# Issue 39 ToricBuilder Issue 38 Fixture Design

## Context

Issue #38 records a real `6 x 6` ToricBuilder `Q` block over
`GF(2)[u^+/-1, v^+/-1]`. Its determinant is the Laurent monomial unit `u*v`.
The existing Laurent `GL_n` normalization boundary can divide by that unit, but
both determinant-one cores still fail in the staged `SL_n` to local `SL_3`
reduction with the same local obligation error.

Issue #39 asks for an offline, checked-in fixture for that case. The fixture
must not depend on a local ToricBuilder checkout, and it must preserve both row
and column normalization variants because they fail for the same deeper reason.

## Design Choice

Add a dedicated fixture module at
`test/fixtures/toricbuilder_issue38_cases.jl` and a focused internal validator
at `test/internal/toricbuilder_issue38_fixture.jl`.

This follows the repository's existing fixture pattern: fixture files define a
small module that constructs Oscar matrices and returns named-tuple entries;
internal tests validate metadata, exact algebraic relations, and negative
controls. A dedicated fixture keeps the Issue #38 matrix easy to cite without
overloading the broader Laurent fixture catalog.

Alternatives considered:

- Add the case to `test/fixtures/laurent_cases.jl`. This would reuse the
  general catalog, but the Issue #38 entry has richer normalization and failure
  metadata than the current catalog fields.
- Put the matrix directly in the internal test. This would verify the current
  behavior, but it would not give later algorithm work a reusable source of
  truth.

## Fixture Interface

The fixture module returns a catalog with one stable entry:

- `id = "toricbuilder-issue-38-q-block"`.
- `kind = :toricbuilder_issue38_q_block`.
- `ring` metadata containing the Oscar Laurent parent and generators `(u, v)`.
- `dimensions = (; matrix = (6, 6))`.
- `inputs.matrix`: the original `Q` block copied from Issue #38.
- `determinant_profile`: expected determinant `u*v`, class
  `:laurent_monomial_unit`, and exponent metadata `(1, 1)`.
- `normalizations.row`: the existing left row normalization returned by
  `normalize_laurent_gl_matrix(Q)`.
- `normalizations.column`: an explicit right diagonal column normalization
  with `Dcol = diag(inv(u*v), 1, ..., 1)` and `core = Q * Dcol`.
- `expected_current_status`: current failure metadata for the row and column
  determinant-one cores.
- `provenance`: Issue #38 URL, source description, and the commit reported in
  the Issue #38 MWE.

The module constructs the matrix directly with Oscar. It does not read cache
files and does not import ToricBuilder.

## Validation

The internal test registers in the `internal` group and exposes validator
helpers for the fixture entry. Validation checks:

- Required metadata fields exist and IDs are unique.
- The original matrix is exactly `6 x 6` over the fixture ring.
- `det(Q) == u*v`.
- `classify_laurent_determinant(Q).classification ==
  :laurent_monomial_unit`.
- The row-normalized core has determinant `1`, reconstructs through
  `verify_laurent_gl_normalization`, and records the expected left correction.
- The column-normalized core has determinant `1` and reconstructs as
  `Q * Dcol`.
- Both normalized cores still throw `ArgumentError` from
  `elementary_factorization`, with messages containing
  `staged SL_n to local SL_3 reduction failure` and
  `failed to solve local SL_3 obligation`.

The negative control corrupts one entry of the original matrix while leaving
the expected determinant and normalization metadata unchanged. The validator
must reject that corrupted entry, proving the fixture is not accepted by
metadata alone.

## Out of Scope

This issue does not make the Issue #38 matrix factorize. It preserves the
current staged failure as expected status metadata so later issues can update
the expected status when the algorithm grows.

## Verification

Focused verification:

```bash
julia --project=. -e 'include("test/internal/toricbuilder_issue38_fixture.jl")'
```

Package verification:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Spec Self-Review

- No placeholders or incomplete requirements remain.
- The design is limited to fixture construction, validation, registration, and
  tests.
- The fixture is independent of a ToricBuilder checkout.
- The negative control has an exact observable failure mode.
