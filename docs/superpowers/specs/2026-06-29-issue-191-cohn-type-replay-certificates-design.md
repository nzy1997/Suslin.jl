# Issue 191 Cohn-Type Replay Certificates Design

## Context

Issue 191 builds on the #190 ordinary-polynomial fixture catalog and the
existing `realize_cohn_type(n, i, j, a, v, R)` helper. The current helper
returns only elementary factors. Later #181 normality work needs a reusable
certificate object that records the exact Cohn-type target, the replayed factor
product, and verification metadata instead of trusting a raw vector of
matrices.

No repository instruction file is present. The relevant merged #190 PR added
`test/fixtures/polynomial_normality_cases.jl`, including
`pw-section2-cohn-type-qq`, which this issue should consume in a focused expert
test.

## Design Choice

Add a small certificate layer in `src/algorithm/cohn_type.jl`:

- `CohnTypeRealizationCertificate` stores `n`, `i`, `j`, coerced `a`, coerced
  one-based `v`, the ring, the chosen auxiliary index, the Cohn-type target,
  elementary factors, replayed product, and verification metadata.
- `realize_cohn_type_certificate(n, i, j, a, v, R)` validates the same shape
  and indexing constraints as `realize_cohn_type`, rejects Laurent and
  non-polynomial certificate rings, coerces inputs into `R`, builds the exact
  target `I + a*v*(v_j*e_i - v_i*e_j)`, reuses the existing Cohn-type factor
  construction, replays the product, and stores the resulting metadata.
- `verify_cohn_type_certificate(cert)` recomputes the target from stored
  coerced inputs, replays the stored factors, checks the stored product, and
  checks that stored verification metadata matches a fresh verification pass.
- `realize_cohn_type` remains available as the factor-only API and continues to
  use the existing factor construction path. The ordinary-polynomial rejection
  is limited to the new certificate route.

Alternative considered: make `realize_cohn_type` a wrapper over the certificate
constructor. That would be shorter, but it would narrow the existing factor API
to ordinary polynomial rings. Keeping factor generation as the primitive
preserves compatibility while still enforcing the issue's ordinary-polynomial
boundary for certificates.

Alternative considered: use an anonymous `NamedTuple` certificate. A concrete
struct matches neighboring certificate APIs such as
`PolynomialFactorizationRouteCertificate` and makes tampered certificates easier
to test and reason about.

## Verification Rules

The verifier returns `false` rather than accepting corrupted proof data when:

- a factor is tampered and no longer replays to the target,
- the stored target no longer equals the Cohn-type formula from stored inputs,
- the stored product no longer equals the replayed factor product,
- the stored `a` or `v` no longer matches the stored target and factors, or
- the stored verification metadata differs from fresh verification metadata.

The constructor raises `ArgumentError` for invalid `n`, indices, vector shape,
zero-based vectors, Laurent rings, or rings outside the ordinary-polynomial
certificate route.

## Files

- Modify `src/algorithm/cohn_type.jl` for the certificate struct, constructor,
  helper routines, and verifier.
- Modify `src/Suslin.jl` to export the new public certificate type and API.
- Modify `test/expert/cohn_type.jl` with fixture-backed certificate tests and
  negative controls.
- Modify `test/expert/polynomial_normality_fixtures.jl` only to avoid redefining
  the #190 fixture module when multiple expert tests include it in one session.

## Verification

Focused command:

```bash
julia --project=. -e 'include("test/expert/cohn_type.jl")'
```

Required package command:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

The focused command must include the #190 Cohn-type fixture, assert that
`cert.target == cert.product`, assert exact equality to the fixture target, and
reject tampered factor, target, and stored-input controls.

## Spec Self-Review

- The design is limited to Cohn-type replay certificates.
- It does not implement rank-one normality, conjugated elementary normality,
  Murthy, Quillen, ECP, Laurent, ToricBuilder, or factor-count optimization.
- The ordinary-polynomial boundary is explicit for the certificate route.
- The negative controls cover the issue's required corruption modes.
