# Issue 20 GL_n Laurent Normalization Design

## Context

Suslin's current elementary factorization core only handles a narrow
univariate polynomial `SL_3` slice. ToricBuilder-style transformation matrices
live one layer earlier in the pipeline: they are exact Laurent matrices in
`GL_n`, and their determinants may be `1`, a Laurent monomial unit, a
permutation sign, or a non-unit. The factorization driver must not silently
pretend every invertible Laurent input is already an `SL_n` input.

Issues #19 and #6 already provide exact GF(2) Laurent fixtures and determinant
metadata for ToricBuilder-style matrices. Issue #21 defines the test command
contract. This issue adds the named boundary that those later consumers should
call before the staged `SL_n` factorization path.

## Design Choice

Add a small public normalization layer:

- `classify_laurent_determinant(A)`: compute `det(A)` and classify it.
- `normalize_laurent_gl_matrix(A)`: return an `SL_n` core plus correction
  metadata when the determinant is supported.
- `verify_laurent_gl_normalization(A, normalization)`: check exact
  reconstruction and determinant-one normalization.

The normalization will construct a left diagonal correction. For determinant
`d`, the correction factor is `D = diag(d, 1, ..., 1)`, and the normalized
core is `D^-1 * A`. This gives `det(D^-1 * A) == 1` and reconstructs exactly as
`D * normalized_matrix == A`.

Alternative considered: reject every Laurent `GL_n` determinant except `1`.
That would be safe, but it would not implement the issue's requested
monomial-unit normalization boundary.

Alternative considered: add determinant correction directly inside
`elementary_factorization`. That obscures the boundary the issue asks to make
explicit. The driver may call the named layer before continuing to the current
staged `SL_3` checks, but the layer itself remains separately testable.

## Determinant Classification

The classifier returns structured named-tuple metadata:

- `classification = :one` when the determinant is exactly `1`.
- `classification = :permutation_sign_unit` when the determinant is exactly
  `-1` in a ring where `-1 != 1`.
- `classification = :laurent_monomial_unit` when the determinant is a single
  Laurent monomial with an invertible coefficient, including GF(2) monomial
  units such as `x^-1*y`.
- `classification = :other_unit` when Oscar reports a unit but it is not a
  single Laurent monomial under the available inspection API.
- `classification = :non_unit` when the determinant is not a unit.

The classifier records the exact determinant and, for monomial cases, the
single exponent vector and coefficient. A permutation matrix over GF(2) whose
sign collapses to `1` is classified as `:one`; over a ring where `-1 != 1`, a
simple transposition is classified as `:permutation_sign_unit`.

## Normalization Output

`normalize_laurent_gl_matrix(A)` validates that `A` is square over a Laurent
polynomial ring. It returns a named tuple with:

- `input_size`: `(n, n)`.
- `ring`: the Laurent parent.
- `determinant`: exact determinant.
- `determinant_classification`: the classifier symbol.
- `normalized_matrix`: the determinant-one core.
- `correction`: named-tuple metadata containing `kind`, `side`, `factor`,
  `inverse_factor`, and `determinant`.

For `:one`, the correction factor is the identity matrix and the normalized
matrix is `A`.

For `:permutation_sign_unit` and `:laurent_monomial_unit`, the correction is
the left diagonal determinant correction described above.

For `:other_unit` and `:non_unit`, normalization throws `ArgumentError` with a
staged message explaining that the determinant is outside the currently
supported `SL_n` path.

## Driver Boundary

`elementary_factorization(A)` will call the named normalization layer before
the current algorithm-specific checks when `A` is over a Laurent polynomial
ring. It will then continue to reject unsupported matrix sizes and multivariate
Laurent rings with the existing staged errors. This preserves the current
algorithm scope while ensuring Laurent `GL_n` inputs cross an explicit
determinant boundary first.

Polynomial inputs keep the existing `SL_3` behavior and determinant check.

## Tests

Add `test/internal/gl_laurent_normalization.jl` and register it in the
`internal` group. The focused test must verify:

- A determinant-one Laurent matrix passes unchanged and verifies exactly.
- A Laurent monomial-unit determinant normalizes to an `SL_n` core, records a
  left diagonal correction, and reconstructs exactly.
- A permutation/sign determinant over a characteristic-not-two Laurent ring is
  classified as `:permutation_sign_unit` and normalizes exactly.
- A non-unit determinant throws an `ArgumentError` with the staged unsupported
  determinant message.
- Tampering with the correction metadata makes exact reconstruction
  verification return `false`.
- `elementary_factorization` reaches the normalization boundary for Laurent
  inputs before continuing to the current staged `SL_3` limits.

Final verification uses the issue #21 package command and the full-suite
command:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
julia --project=. test/runtests.jl all
```

## Spec Self-Review

- No incomplete markers remain.
- The scope is limited to determinant classification, normalization metadata,
  exact verification, driver boundary wiring, docs, and focused tests.
- The design does not claim full Laurent elementary factorization support.
- Reconstruction verification is explicit and testable, including a negative
  control for tampered metadata.
