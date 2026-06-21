# Issue 59 Laurent GL Certificate Design

## Goal

Add a narrow Laurent `GL_n` certificate path for the original Issue #38 `Q`
block. The certificate must reconstruct the original `Q` exactly by combining
the Laurent determinant correction from `normalize_laurent_gl_matrix(Q)` with
elementary factors for the determinant-one core from the Issue #57 column-peel
route.

## Context

Issue #38's original `6 x 6` matrix over `GF(2)[u^+/-1, v^+/-1]` has
determinant `u*v`. It is therefore not a pure elementary product, because every
elementary factor has determinant `1`.

Issue #20 already added the Laurent `GL_n` normalization boundary. For a
supported monomial-unit determinant `d`, it builds a left diagonal correction
`D = diag(d, 1, ..., 1)`, a normalized determinant-one core `D^-1 * A`, and
metadata that verifies `D * core == A`. Issue #57 then added the recursive
column-peel factorization path for the Issue #38 determinant-one row and column
cores.

Issue #59 should compose those two existing capabilities without changing the
meaning of `verify_factorization(A, factors)`.

## Approaches Considered

1. Add a focused certificate object for the Issue #38 fixture and expose a
   verifier.

   This is the selected approach. It is a direct composition of the existing
   normalization and column-peel layers. It keeps the non-elementary determinant
   correction explicit and returns elementary factors only for the normalized
   `SL_n` core.

2. Teach `elementary_factorization(Q)` to return the core factors for the
   original `Q`.

   This is rejected. Returning core-only factors for a `det(Q)=u*v` input would
   make `verify_factorization(Q, factors)` false and would violate the existing
   pure elementary-factor API.

3. Add a broad Laurent `GL_n` factorization API for arbitrary supported
   monomial-unit matrices.

   This is out of scope. Issue #59 asks for the original Issue #38 fixture path
   only, and broader support would need additional acceptance cases and API
   design.

## Design

Create `src/algorithm/laurent_gl_certificate.jl` and include it after
`laurent_column_peel.jl`, so it can call both `normalize_laurent_gl_matrix(A)`
and `_factor_laurent_sl_column_peel(core)`.

The main public entry point is:

```julia
laurent_gl_factorization_certificate(A)
```

It accepts the original Issue #38 matrix shape through the existing exact
Laurent normalization and supported column-peel core path. It returns a
`LaurentGLFactorizationCertificate` with:

- `original_matrix`: the input matrix.
- `determinant_profile`: the exact profile from
  `classify_laurent_determinant(A)` as stored by normalization.
- `normalization`: the named tuple returned by `normalize_laurent_gl_matrix(A)`.
- `correction`: the normalization correction metadata.
- `inverse_correction`: the inverse correction factor for direct access.
- `normalized_core`: the determinant-one core.
- `core_factorization`: the Issue #57 column-peel certificate for the core.
- `core_factors`: the elementary factor sequence for the core.
- `reconstructed_product`: the exact product `correction.factor *
  product(core_factors)`.
- `verification`: verifier metadata.

The selected variant is the left-correction variant already implemented by the
normalization layer:

```text
Q = correction.factor * normalized_core
normalized_core = product(core_factors)
Q = correction.factor * product(core_factors)
```

Export:

```julia
laurent_gl_factorization_certificate
verify_laurent_gl_factorization_certificate
LaurentGLFactorizationCertificate
```

`elementary_factorization(Q)` remains unchanged for the original `Q`: it should
continue to throw the existing Laurent `GL_n` normalization boundary error
instead of returning factors that do not reconstruct `Q`.

## Verification Rules

`verify_laurent_gl_factorization_certificate(certificate)::Bool` must recompute
and check the certificate instead of trusting stored fields. It should return
`false` on malformed or tampered certificates rather than throwing for ordinary
bad metadata.

The verifier checks:

- The original matrix is square over a Laurent polynomial ring.
- Recomputed `normalize_laurent_gl_matrix(original_matrix)` matches the stored
  determinant profile, correction, inverse correction, and normalized core.
- The determinant classification is supported by normalization and the
  normalized core has determinant `1`.
- The core factors verify against `normalized_core` with
  `verify_factorization(normalized_core, core_factors)`.
- The stored core column-peel certificate verifies with
  `_verify_laurent_column_peel_replay`.
- The stored core factors exactly match the core certificate's factors.
- The exact reconstructed product is recomputed from the correction factor and
  core factors and equals both the stored `reconstructed_product` and the
  original matrix.

Tampering with either the determinant correction factor or one core factor must
make the verifier return `false`.

## Tests

Add `test/expert/issue38_laurent_gl_certificate.jl` and register it in the
expert group. The focused test will:

- Load `test/fixtures/toricbuilder_issue38_cases.jl`.
- Assert the original `Q` determinant classification is
  `:laurent_monomial_unit` and the determinant is `u*v`.
- Build `laurent_gl_factorization_certificate(Q)`.
- Assert the normalized core has determinant `1`.
- Assert `verify_factorization(core, core_factors)` is true.
- Assert `verify_laurent_gl_factorization_certificate(certificate)` is true.
- Assert `certificate.reconstructed_product == Q`.
- Assert the left-correction reconstruction equals `Q`.
- Assert `verify_factorization(Q, core_factors)` is false.
- Tamper with the correction factor and one core factor, and assert the
  certificate verifier returns false for each case.
- Assert `elementary_factorization(Q)` still rejects the original `Q` with the
  Laurent `GL_n` normalization boundary message.

Focused verification command:

```bash
julia --project=. -e 'include("test/expert/issue38_laurent_gl_certificate.jl")'
```

Package verification command:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Out Of Scope

Do not claim that the original Issue #38 `Q` is a pure elementary product. Do
not broaden support to arbitrary Laurent `GL_n` inputs beyond the existing
normalization and Issue #38 core factorization path. Do not change
`verify_factorization(A, factors)`.

## Spec Self-Review

- No placeholders or incomplete sections remain.
- The selected approach directly composes the merged Issue #20 and Issue #57
  layers.
- The original `Q` remains outside the pure elementary-factor API.
- The verifier requirements explicitly reject trusted stored products and cover
  determinant-correction and core-factor tampering.
