# Issue 158 Deferred Laurent Correction Hoist Design

## Context

Issue #156 records a determinant-deferred Laurent peel certificate with the
exact relation

```text
left_product * original_matrix * right_product == blockdiag(deferred_submatrix, I)
```

Issue #157 enriches that certificate by classifying only the deferred
submatrix determinant. For supported monomial-unit determinants it builds a
left diagonal correction `D_deferred` and a determinant-one
`normalized_deferred_core` such that

```text
D_deferred * normalized_deferred_core == deferred_submatrix
```

The missing step is to turn that deferred correction into original
matrix-level certificate data. The correction cannot be moved by only storing
`blockdiag(D_deferred, I)`: the elementary factors to its left must be rewritten
by exact diagonal conjugation/scaling.

## Approach Options

Recommended: add a narrow internal certificate constructor for the supported
row/left deferred correction path. It consumes #157 metadata, embeds the
deferred diagonal correction into the original dimension, factors the
normalized deferred core, rewrites the inverse recorded left peel factors while
moving the diagonal correction to the front, appends the embedded core factors
and inverse right peel factors, and verifies exact reconstruction.

Alternative: add the embedded diagonal correction to the current metadata
without rewriting elementary coefficients. This is rejected because it is the
bug the issue calls out: `P * D` is not generally equal to `D * P` for an
elementary product `P`.

Alternative: expose a public user-facing row/column hoist API. This is out of
scope because the issue explicitly says not to expose row/column options and to
implement one correction side first.

## Chosen Design

Add `LaurentLazyGLHoistCertificate` in
`src/algorithm/laurent_gl_certificate.jl`. It is internal for now and stores:

- `original_matrix`;
- `deferred_metadata`;
- `overall_determinant`;
- `correction`;
- `inverse_correction`;
- `normalized_deferred_core`;
- `normalized_deferred_factorization`;
- `normalized_deferred_factors`;
- `rewritten_left_factors`;
- `elementary_factors`;
- `elementary_product`;
- `reconstructed_product`;
- `verification`.

Add `_laurent_gl_lazy_deferred_correction_certificate(metadata; ...)`, plus
`_laurent_gl_lazy_deferred_correction_certificate(A; ...)` as a convenience
constructor that obtains #157 metadata through
`_factor_laurent_gl_lazy_determinant_peel`. The constructor supports metadata
where `supported == true`, `determinant_source == :deferred_submatrix`, and
the correction side is `:left`. Determinant-one metadata can produce the same
shape with an identity correction. Unsupported metadata remains a staged
boundary and is rejected by this constructor.

For a monomial-unit left correction, the constructor computes:

```text
D = blockdiag(D_deferred, I)
D^-1 = blockdiag(D_deferred^-1, I)
C = blockdiag(normalized_deferred_core, I)
core factors = factorization(normalized_deferred_core)
A = left_product^-1 * D * C * right_product^-1
```

The hoisted certificate stores the original-level correction `D` and an
elementary product `P` such that:

```text
A == D * P
P == rewritten(left_product^-1 across D) *
     blockdiag(core factors, I) *
     right_product^-1
```

The rewrite uses the issue's row/left algebraic relation. If `D` has diagonal
entries `d_i`, then for an elementary row factor `E_ij(a)`:

```text
E_ij(a) * D == D * E_ij(d_i^-1 * a * d_j)
```

The helper layer extracts diagonal entries from the correction matrix and
rewrites elementary factors one at a time. Tests cover the helper directly, so
the main constructor does not open-code coefficient scaling.

## Verification

Add `_laurent_gl_lazy_deferred_correction_certificate_verification(certificate)`
and `_verify_laurent_gl_lazy_deferred_correction_certificate(certificate)`.
The verifier recomputes and checks:

- the #157 metadata verifies;
- the correction and inverse correction are square original-dimension matrices
  over the same ring and multiply to identity in both orders;
- the embedded correction has determinant `overall_determinant`;
- the normalized deferred factorization verifies exactly;
- `rewritten_left_factors` are exactly the factors obtained by hoisting
  `left_product^-1` across the embedded diagonal correction;
- `elementary_product == product(elementary_factors)`;
- `reconstructed_product == correction.factor * elementary_product`;
- `reconstructed_product == original_matrix`.

The negative control builds a deliberately wrong certificate that embeds the
diagonal correction but keeps the inverse left factors unrewritten. The
verifier must reject it, proving the hoist checks the algebra instead of only
trusting a stored correction matrix.

## Testing

Add `test/expert/laurent_lazy_correction_hoist.jl` and register it in the
expert group. The focused test covers:

- direct coefficient rewriting for the helper relation
  `E_ij(a) * D == D * E_ij(d_i^-1 * a * d_j)`;
- a supported monomial-unit deferred fixture where
  `certificate.overall_determinant` equals the #157 deferred determinant;
- exact reconstruction at original level:
  `certificate.correction.factor * certificate.elementary_product == A`;
- verification of the embedded correction and inverse correction;
- rejection of the wrong certificate that moves `D` without rewriting the
  left elementary coefficients.

Focused verification command:

```bash
julia --project=. -e 'include("test/expert/laurent_lazy_correction_hoist.jl")'
```

Package verification command:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Scope

Do not export new APIs. Do not expose user-facing row/column options. Do not
change the pure elementary-factor contract for non-`SL_n` Laurent inputs. Do
not claim support for unsupported deferred determinant classes.

## Automatic Decisions

- Clarifying questions were resolved by the Standing Answer Policy because this
  is a non-interactive Agent Desk run.
- The visual companion was skipped because no visual decision would clarify
  the algebraic certificate path.
- The recommended narrow row/left hoist approach was selected because #157
  already produces left deferred corrections and the issue asks to implement
  one correction side first.
