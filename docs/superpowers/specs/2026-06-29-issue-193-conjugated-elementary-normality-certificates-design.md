# Issue 193 Conjugated Elementary Normality Certificates Design

## Context

Issue 193 builds on the #190 ordinary-polynomial normality fixture catalog, the
#191 Cohn-type certificate API, and the #192 orthogonal rank-one normality
certificate API. The current `realize_conjugate_elementary(B, i, j, a)` helper
already rewrites selected conjugates into elementary factors using the
Park-Woodburn rank-one route, but it returns only raw factors and supports a
Laurent path for older callers.

The new layer should certify the ordinary-polynomial `SL_n` case only. Given
`A`, distinct `i` and `j`, and coefficient `a`, the certificate records the
explicit convention

```julia
A * elementary_matrix(n, i, j, a, R) * inv(A)
```

and proves that the final elementary factors multiply exactly to that target.
The GitHub CLI could not fetch live issue comments in this sandbox because the
configured proxy was unavailable. The issue body supplied by Agent Desk,
merged #190/#191/#192 code, and checked-in specs were used as the binding
context.

## Design Choice

Add a concrete `ConjugatedElementaryNormalityCertificate` in
`src/algorithm/normality.jl` with:

- `n`, `A`, `i`, `j`, coerced `a`, and the ordinary polynomial `ring`.
- `determinant`, `inverse_A`, and the stored `elementary_matrix`.
- `conjugation_convention = :A_E_invA` and `conjugation_target`, with field
  names making the route explicit.
- Derived vectors `v = A[:, i]`, `w = a * inv(A)[j, :]`, and
  `g = inv(A)[i, :]`.
- A child `RankOneNormalityCertificate` for `I + v*w`.
- The final elementary `factors`, replayed `product`, and verification
  metadata.

Export `ConjugatedElementaryNormalityCertificate`,
`realize_conjugate_elementary_certificate(A, i, j, a)`, and
`verify_conjugate_elementary_certificate(cert)`. Keep
`realize_conjugate_elementary(B, i, j, a)` as the factor-only helper for
current callers; it must retain its Laurent behavior and zero-coefficient
behavior.

Alternative considered: route `realize_conjugate_elementary` through the new
certificate. Rejected because this issue is polynomial-only while the helper
currently supports Laurent matrices through its inverse helper.

Alternative considered: duplicate the rank-one coefficient and child Cohn
certificate replay inside the conjugated certificate. Rejected because #192
already owns that replay contract, and duplication would make later normality
consumers compare two independently maintained conventions.

Alternative considered: support both `A * E * inv(A)` and `inv(A) * E * A`.
Rejected for this issue because existing callers and fixtures already use
`A * E * inv(A)`, and the issue asks that later modules not infer the route.
The single convention is explicit in stored metadata and verification flags.

## Verification Rules

The constructor rejects:

- non-square matrices, `n < 3`, out-of-range indices, or `i == j`;
- Laurent or non-polynomial certificate rings;
- coefficients that cannot be coerced into `base_ring(A)`; and
- ordinary-polynomial matrices with determinant not equal to `one(R)`.

The verifier returns `false` unless a fresh replay confirms:

- the stored matrix/ring/indices/coefficient pass the constructor checks;
- the stored determinant is `one(R)`;
- `inverse_A` is a two-sided inverse of `A`;
- the stored elementary matrix equals `E_ij(a)`;
- the stored convention is exactly `:A_E_invA`;
- the stored target equals `A * E_ij(a) * inverse_A`;
- the stored `v`, `w`, and `g` match the selected convention;
- the child rank-one certificate verifies and matches `v`, `w`, `g`, `R`, and
  the conjugation target;
- the stored factors are exactly the child rank-one factors;
- replaying stored factors reproduces the stored product; and
- stored verification metadata equals fresh core verification.

## Files

- Modify `src/algorithm/normality.jl` for the certificate type, checked-data
  helper, constructor, core verifier, and public verifier.
- Modify `src/Suslin.jl` to export the new type and public functions.
- Modify `test/expert/normality.jl` to cover non-fixture `SL_3` and `SL_4`
  examples, fixture-backed metadata, raw helper compatibility, and negative
  controls.
- Modify `test/public/api_surface.jl` to cover the new public exports.

## Verification

Focused command:

```bash
julia --project=. -e 'include("test/expert/normality.jl")'
```

Required package command:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

The focused command must include non-fixture ordinary-polynomial `SL_3` and
`SL_4` cases where certificate factors multiply exactly to the stored
`A * E_ij(a) * inv(A)` target. Negative controls must reject determinant-not-one
matrices, equal indices, coefficients from the wrong ring, and a certificate
whose factor sequence has been tampered with.

## Spec Self-Review

- The design is limited to ordinary-polynomial conjugated elementary
  certificates.
- The convention is named and verified as `:A_E_invA`.
- The implementation reuses #192 rank-one certificates and does not duplicate
  Cohn-type child replay.
- Existing raw factor-only behavior remains available, including Laurent
  callers.
- It does not implement Quillen patching, Murthy local solving, ECP reducers,
  recursive `SL_n` factorization, Laurent certificates, ToricBuilder support,
  or Steinberg factor-count optimization.
