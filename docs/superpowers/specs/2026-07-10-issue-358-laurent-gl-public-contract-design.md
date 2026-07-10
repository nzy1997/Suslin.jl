# Issue 358 Laurent GL Public Contract Design

## Goal

Finalize the public Laurent `GL_n` decomposition contract for supported
monomial-unit determinant inputs without changing the elementary factor API.
`elementary_factorization(A)` remains an elementary-only `SL_n` factor sequence
API. `laurent_gl_factorization_certificate(A)` is the public API for supported
Laurent `GL_n` inputs whose original determinant is `1` or a supported Laurent
monomial unit.

## Selected Approach

Use the existing eager `LaurentGLFactorizationCertificate` as the final public
contract. The certificate records the original matrix, determinant profile,
normalization, left diagonal correction, inverse correction, determinant-one
core, core factorization certificate, ordered elementary core factors,
reconstructed product, and detailed verification status.

The supported reconstruction relation is:

```julia
A == certificate.correction.factor * prod(certificate.core_factors)
```

The monomial-unit correction is not an elementary factor. Every elementary
matrix has determinant one, so a product of elementary matrices also has
determinant one. A Laurent matrix with determinant `u*v` or another nontrivial
monomial unit therefore cannot be represented by an elementary-only factor
sequence for the original input.

## Behavior

- `elementary_factorization(A)` keeps returning only matrices whose product is
  exactly `A`.
- Laurent determinant-one inputs continue through the Laurent `SL_n` fallback.
- Supported Laurent monomial-unit `GL_n` inputs rejected by
  `elementary_factorization(A)` get an intentional determinant-contract error
  that points callers to `laurent_gl_factorization_certificate(A)`.
- `laurent_gl_factorization_certificate(A)` independently verifies the
  determinant classification, correction, inverse correction, normalized core,
  elementary core factor product, and exact original reconstruction.

## Verification Design

The verifier recomputes `normalize_laurent_gl_matrix(A)` from the original
matrix and compares it with stored certificate fields. It then checks that the
normalized core has determinant one, the stored core factorization replay
verifies, the stored `core_factors` match the core factorization, every core
factor is elementary, their product equals the normalized core, and the left
correction times that product reconstructs the original matrix exactly.

Negative controls must fail for a modified correction, incorrect correction
side, reordered or modified core factor, copied reconstructed product that does
not match replay, malformed certificate, and non-unit determinant input.

## Documentation And Example

Update README and Documenter scope text to state the split contract:

- `elementary_factorization(A)` is an `SL_n` elementary factor API;
- `laurent_gl_factorization_certificate(A)` is the supported Laurent `GL_n`
  decomposition API;
- the diagonal monomial-unit correction is deliberately outside the elementary
  factor sequence.

Update `example/toric_decoupling/issue38_unit_correction.jl` so the issue-38
`case_001` walkthrough demonstrates the final certificate contract and prints
verified reconstruction status instead of describing the original input as a
temporary staged boundary.

## Tests

Extend the focused issue-38 expert test and public API surface checks. Run the
issue verification commands plus package tests:

```bash
julia --project=. example/toric_decoupling/issue38_unit_correction.jl
julia --project=. -e 'include("test/expert/issue38_laurent_gl_certificate.jl")'
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Out Of Scope

Do not introduce heterogeneous factor types, do not pretend the determinant
correction is elementary, do not add support for arbitrary non-monomial units,
and do not bundle LaurentToPoly or d14 conversion work.
