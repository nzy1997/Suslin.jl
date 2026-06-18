# ToricBuilder Integration Contract

## Consumer Boundary

The first ToricBuilder boundary Suslin records is the small output of
`factor_toric_block(3, x, y, R)` over `GF(2)[x^+/-1, y^+/-1]`.
ToricBuilder produces exact Laurent transformation matrices; Suslin records
them as contract fixtures before expanding factorization support.

## Fixture Entries

| Entry | Role | Ring | Size | Determinant classification | Current Suslin status |
| --- | --- | --- | --- | --- | --- |
| `factor_toric_block_3_qinv` | `Qinv` | `GF(2)[x^+/-1, y^+/-1]` | `16 x 16` | `one` | staged unsupported input |
| `factor_toric_block_3_pinv` | `Pinv` | `GF(2)[x^+/-1, y^+/-1]` | `8 x 8` | `one` | staged unsupported input |

## Provenance

- ToricBuilder path: `/Users/nzy/jcode/topological-code-decoupling/julia_code/ToricBuilder`
- ToricBuilder commit: `fa7f82252d42fdc0b2726bc48af24ac4c70a8d73`
- Source function: `src/toric_form/toric_factorization.jl:factor_toric_block`
- Generation command: `factor_toric_block(3, x, y, R)` with
  `R, (x, y) = laurent_polynomial_ring(GF(2), ["x", "y"])`
- Returned block coarse-graining size: `(2, 2)`

## Exact Relations

The checked-in fixture stores the source transformation matrix for each
inverse-style entry. Tests verify these equalities by exact multiplication:

- `column_transformation * Qinv == I_16`
- `row_transformation * Pinv == I_8`

## Determinant Classes

Fixture metadata records one of:

- `one`: determinant is exactly `1`.
- `laurent_monomial_unit`: determinant is an invertible Laurent monomial other
  than `1`.
- `other_unit`: determinant is a non-monomial unit if a future coefficient
  ring admits one.
- `non-unit`: determinant is not a unit.

The current fixture entries both use `one`.

## Suslin Output Contract

ToricBuilder ultimately needs a verified transformation certificate. For this
first contract, that certificate is the exact inverse relation plus determinant
classification. Raw elementary factors and normalized factor metadata are not
required until later implementation issues expand Laurent support.

Current Suslin behavior is deliberately unsupported for these matrices:
`elementary_factorization` accepts only a narrow `3 x 3` univariate polynomial
`SL_3` slice, while these fixtures are larger two-variable Laurent matrices.
