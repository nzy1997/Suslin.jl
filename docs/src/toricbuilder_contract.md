# ToricBuilder Integration Contract

## Consumer Boundary

The first ToricBuilder boundary Suslin records is the small output of
`factor_toric_block(3, x, y, R)` over `GF(2)[x^+/-1, y^+/-1]`.
ToricBuilder produces exact Laurent transformation matrices; Suslin records
them as contract fixtures before expanding factorization support.

## Fixture Entries

| Entry | Role | Ring | Size | Determinant classification | Current Suslin status |
| --- | --- | --- | --- | --- | --- |
| `factor_toric_block_3_qinv` | `Qinv` | `GF(2)[x^+/-1, y^+/-1]` | `16 x 16` | `one` | supported by Laurent column peel |
| `factor_toric_block_3_pinv` | `Pinv` | `GF(2)[x^+/-1, y^+/-1]` | `8 x 8` | `one` | supported by Laurent column peel |

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

## GL_n Laurent Normalization Boundary

Suslin exposes `normalize_laurent_gl_matrix(A)` for exact Laurent `GL_n`
inputs before any staged `SL_n` factorization attempt. The boundary computes
and classifies the determinant, then either returns a determinant-one core
with explicit correction metadata or throws a staged `ArgumentError`.

Supported corrections are determinant `1`, permutation/sign determinant `-1`
where the coefficient ring distinguishes it from `1`, and Laurent monomial
unit determinants such as `x^-1*y` over `GF(2)`. The correction metadata stores
a left diagonal factor `D` and verifies exact reconstruction as
`D * normalized_matrix == A`.

For supported Laurent monomial-unit `GL_n` inputs, the public decomposition is
the verified `laurent_gl_factorization_certificate(A)` relation
`A == certificate.correction.factor * prod(certificate.core_factors)`. The
core factors are elementary factors of the determinant-one normalized core; the
diagonal correction is recorded separately because it is not elementary when
its determinant is nontrivial.

## Suslin Output Contract

ToricBuilder ultimately needs a verified transformation certificate. For
supported Laurent monomial-unit `GL_n` inputs, it consumes the public
Laurent-`GL_n` certificate contract above rather than a core-only elementary
factor sequence.
