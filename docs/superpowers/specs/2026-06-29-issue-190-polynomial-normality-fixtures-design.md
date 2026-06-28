# Issue 190 Polynomial Normality Fixture Catalog Design

## Context

Issue 190 needs a small ordinary-polynomial fixture catalog for the three
Park-Woodburn Section 2 normality layers: Cohn-type matrices, orthogonal
rank-one matrices, and conjugated elementary matrices. The repository already
has direct expert tests for `realize_cohn_type` and
`realize_conjugate_elementary`, plus reusable checked-in fixture catalogs under
`test/fixtures`.

This issue should not add a normality certificate API. The deliverable is test
support for later #181 work: named, source-grounded cases with enough metadata
for later certificate tests to consume exact inputs and expected targets.

## Design Choice

Add a dedicated `test/fixtures/polynomial_normality_cases.jl` catalog and a
focused expert validator at `test/expert/polynomial_normality_fixtures.jl`.
Register the expert validator in `test/runtests.jl`, matching the repository's
explicit expert test list.

The catalog will expose one positive case per Section 2 layer:

- `pw-section2-cohn-type-qq`: `I + a*v*(v_j*e_i - v_i*e_j)`.
- `pw-section2-orthogonal-rank-one-qq`: `I + v*w` with `w*v == 0` and an
  explicit Bezout row `g` with `g*v == 1`.
- `pw-section2-conjugated-elementary-qq`: `B * E_ij(a) * inv(B)`.

Each case records ring constructor metadata, the ring object and generators,
inputs, target matrix, expected convention, source metadata, and consumer issue
ids. The catalog also carries negative controls rejected by the same verifier:
a tampered Cohn target, a rank-one row with invalid orthogonality, and a
tampered conjugated-elementary target.

Alternative considered: extend `park_woodburn_polynomial_cases.jl`. That file
tracks driver-route matrices and status metadata for polynomial factorization,
while this issue needs normality-layer input conventions. A separate catalog is
clearer and avoids overloading the route schema.

Alternative considered: only add assertions to `test/expert/cohn_type.jl` and
`test/expert/normality.jl`. That would test current helpers, but it would not
give later #181 certificate tests a reusable fixture source.

## Validation Rules

The expert validator will require:

- Every positive entry has a non-empty id, Section 2 layer symbol, source,
  ordinary exact field-backed polynomial ring metadata, target matrix,
  expected convention, and consumer issue ids.
- Cohn-type targets equal the literal convention
  `I + a*v*(v_j*e_i - v_i*e_j)`.
- Rank-one targets equal `I + v*w`, with `w*v == 0` and `g*v == 1`.
- Conjugated elementary targets equal `B * E_ij(a) * inv(B)`, with extracted
  `v`, `w`, and `g` matching the stored presentation.
- Every positive target has determinant `one(R)`.
- Every catalog negative control throws `ArgumentError` under the same
  validator.

## Files

- Create `test/fixtures/polynomial_normality_cases.jl`.
- Create `test/expert/polynomial_normality_fixtures.jl`.
- Modify `test/runtests.jl` to include the expert validator.

## Verification

Focused command:

```bash
julia --project=. -e 'include("test/expert/polynomial_normality_fixtures.jl")'
```

Required package command:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

The focused command must print a fixture-catalog testset, validate the three
named positive cases by exact equality to stored target matrices, and reject at
least one invalid orthogonality/unimodularity or tampered target control.

## Spec Self-Review

- The design stays inside test fixtures and expert validation.
- It uses ordinary polynomial rings over `QQ`.
- It does not add certificate APIs, Laurent/ToricBuilder/Murthy/Quillen logic,
  or factor-count behavior.
- The negative controls exercise both algebraic invalidity and target
  tampering.
