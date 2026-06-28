# Issue 157 Deferred Laurent Submatrix Normalization Design

## Context

Issue #156 added `LaurentDeterminantDeferredPeelCertificate`, an internal
certificate that records elementary Laurent column peel steps and stops at a
smaller `deferred_submatrix`. Its replay metadata proves
`left_product * original_matrix * right_product == blockdiag(deferred_submatrix, I)`.
Because the recorded peel operations are elementary, the determinant of the
deferred submatrix is the determinant still needing classification.

The current lazy GL peel route still treats every non-one deferred determinant
as an unimplemented error. Issue #157 changes that boundary. The route should
classify only the deferred submatrix determinant, normalize supported Laurent
monomial-unit deferred cores to determinant one, and return explicit structured
metadata for unsupported determinant classes. It must not compute the original
full matrix determinant and must not hoist the deferred correction back to the
original matrix-level certificate.

## Approach Options

Recommended: add a small internal enrichment function that consumes a
`LaurentDeterminantDeferredPeelCertificate` and returns named metadata for the
deferred determinant. This keeps the #156 replay certificate intact, reuses
`classify_laurent_determinant` plus existing diagonal correction helpers, and
gives tests a direct interface matching the issue statement.

Alternative: extend `LaurentDeterminantDeferredPeelCertificate` with
normalization fields. This would make every certificate construction perform
determinant classification, including tests that only want replay metadata, and
would mix replay proof state with determinant correction state.

Alternative: only change `_factor_laurent_gl_lazy_determinant_peel` to return
different values for every determinant class. That would address the immediate
failure but would hide the certificate-to-metadata interface requested by the
issue and make unsupported boundary tests less direct.

## Chosen Design

Add `_normalize_laurent_determinant_deferred_submatrix(certificate; ...)` in
`src/algorithm/laurent_column_peel.jl`. The function accepts a #156 deferred
peel certificate and classifies `certificate.deferred_submatrix` using the
provided `determinant_probe`, defaulting to `classify_laurent_determinant`.
The probe is called on the deferred submatrix only.

The function returns a named tuple with these fields:

- `peel_certificate`;
- `deferred_submatrix`;
- `determinant_source`;
- `determinant_profile`;
- `overall_determinant`;
- `determinant_classification`;
- `supported`;
- `deferred_correction`;
- `deferred_diagonal_correction`;
- `normalized_deferred_core`;
- `staged_boundary`;
- `verification`.

For `:one`, the metadata is supported. `overall_determinant` is one,
`deferred_correction` is an identity correction scoped to the deferred
submatrix, `deferred_diagonal_correction` is `nothing`, and
`normalized_deferred_core == deferred_submatrix`.

For `:laurent_monomial_unit`, the metadata is supported.
`deferred_correction` and `deferred_diagonal_correction` hold the existing
left diagonal determinant correction scoped to the deferred submatrix, and
`normalized_deferred_core = correction.inverse_factor * deferred_submatrix`.
Verification checks `det(normalized_deferred_core) == 1` and
`correction.factor * normalized_deferred_core == deferred_submatrix`.

For `:non_unit` and all other non-supported classes, the metadata is an
explicit staged boundary. `supported` is false, both correction fields are
`nothing`, `normalized_deferred_core` is `nothing`, and `staged_boundary`
contains the determinant source, determinant value, classification, deferred
submatrix size, and a reason such as `:non_unit_deferred_determinant` or
`:unsupported_deferred_unit_class`.

Update `_factor_laurent_gl_lazy_determinant_peel` at the existing deferred
determinant probe point. It should build a one-step deferred certificate,
normalize/classify that certificate, and return the metadata for non-one
deferred determinant classes instead of throwing. The determinant-one branch
continues the existing SL peel recursion and returns the current
`LaurentColumnPeelFactorization` type.

## Error Handling

The enrichment function verifies the incoming peel certificate replay before
trusting its deferred submatrix. Invalid replay metadata throws an internal
error, matching the existing column-peel verification style.

Unsupported determinant classes are not errors in the enrichment layer. They
return structured staged boundaries so callers can preserve the exact boundary
without coercing, inverting, or silently dropping the determinant.

## Testing

Add `test/expert/laurent_lazy_submatrix_normalization.jl` and register it in
the expert group. The test covers:

- a determinant-one deferred submatrix that remains unchanged and verifies as a
  supported SL core;
- a visible monomial-unit fixture with deferred determinant `u*v`, proving
  `overall_determinant == u*v`,
  `determinant_classification == :laurent_monomial_unit`, and
  `det(normalized_deferred_core) == 1`;
- a probe guard proving classification is invoked on the deferred submatrix,
  not on the original full matrix;
- a #154 non-unit determinant fixture embedded as the deferred submatrix,
  proving the result is a staged boundary with no correction and no normalized
  core;
- the lazy GL peel route no longer throws on the existing monomial-unit lazy
  fixture and returns the enriched metadata.

Verification commands:

```bash
julia --project=. -e 'include("test/expert/laurent_lazy_submatrix_normalization.jl")'
julia --project=. -e 'include("test/expert/laurent_lazy_peel_no_initial_det.jl")'
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Scope

Do not export new APIs. Do not add original matrix-level correction replay.
Do not claim elementary factors for the original matrix when the deferred
submatrix required a non-one correction.
