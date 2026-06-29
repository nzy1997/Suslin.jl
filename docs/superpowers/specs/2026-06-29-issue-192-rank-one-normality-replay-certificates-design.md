# Issue 192 Rank-One Normality Replay Certificates Design

## Context

Issue 192 builds on the #190 ordinary-polynomial Park-Woodburn fixture catalog
and the #191 Cohn-type replay certificate API. The merged #191 code already
records exact Cohn-type targets, factors, products, and verification metadata
for ordinary polynomial rings. The existing `realize_conjugate_elementary`
helper in `src/algorithm/normality.jl` also contains the rank-one-to-Cohn
coefficient formula, but it returns only raw factors.

The new layer should certify a presentation `I + v*w` where `w*v == 0` and a
Bezout row `g` satisfies `g*v == 1`. It must keep the Cohn-type factor formula
owned by the Cohn certificate constructor and only assemble child certificates
using the deterministic Park-Woodburn coefficients
`a_ij = w_i*g_j - w_j*g_i` for `i < j`.

The GitHub CLI could not fetch live issue comments in this sandbox because the
configured proxy was unavailable. The issue body supplied by Agent Desk,
merged #190/#191 commits, and checked-in #190/#191 specs were used as the
binding context.

## Design Choice

Add a concrete `RankOneNormalityCertificate` in `src/algorithm/normality.jl`
with:

- `n`, coerced one-based `v`, `w`, and `g`, plus the ordinary polynomial ring.
- Stored orthogonality and Bezout replay values, namely `w*v` and `g*v`.
- A deterministic `cohn_coefficients` table containing every lexicographic
  pair `i < j`, including zero coefficients, so replay metadata is stable.
- `child_certificates` for the nonzero coefficient entries only, built by
  calling `realize_cohn_type_certificate(n, i, j, a_ij, v, R)`.
- The concatenated elementary `factors`, the rank-one `target = I + v*w`, the
  replayed `product`, and verification metadata.

Export `RankOneNormalityCertificate`,
`realize_rank_one_normality_certificate(v, w, g, R)`, and
`verify_rank_one_normality_certificate(cert)`.

Alternative considered: include zero-coefficient child Cohn certificates.
That would make the child list mirror the table exactly, but it would add
identity/canceling elementary factors and diverge from the existing
`realize_conjugate_elementary` convention, which skips zero coefficients.
The table still records zero entries; only child certificates and factors skip
them.

Alternative considered: expose a factor-only rank-one helper. The issue asks
for replayable certificates, and existing callers can already obtain factors
through conjugated elementary normality. Keeping the public surface certificate
focused avoids adding an unrequested factor route.

## Verification Rules

The constructor rejects:

- `n < 3`, vectors with different lengths, or non-one-based vectors.
- Laurent or non-polynomial rings in the certificate route.
- Inputs that fail `w*v == 0` or `g*v == 1` after coercion into `R`.

The verifier returns `false` unless a fresh replay confirms:

- the stored ring and coerced vectors are valid;
- stored orthogonality equals `zero(R)` and stored Bezout equals `one(R)`;
- the stored coefficient table matches lexicographic recomputation exactly;
- every child Cohn-type certificate verifies and matches the corresponding
  nonzero table entry;
- the stored factor sequence is exactly the concatenation of child factors;
- replaying the stored factors reproduces the stored product; and
- the stored target equals both the rank-one formula and the replayed product,
  with stored verification metadata equal to the fresh core verification.

## Files

- Modify `src/algorithm/normality.jl` for the rank-one certificate type,
  constructor, replay helpers, and verifier.
- Modify `src/Suslin.jl` to export the new type and public functions.
- Create `test/expert/normality_rank_one.jl` with fixture-backed positive
  coverage and negative controls for bad orthogonality, bad Bezout data, and
  tampered child Cohn-type factors.
- Modify `test/runtests.jl` to register the focused expert test for the
  `all` expert suite.

## Verification

Focused command:

```bash
julia --project=. -e 'include("test/expert/normality_rank_one.jl")'
```

Required package command:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

The focused command must consume
`pw-section2-orthogonal-rank-one-qq` from
`test/fixtures/polynomial_normality_cases.jl`, assert exact equality to
`I + v*w`, verify every child Cohn certificate, and reject the requested
negative controls.

## Spec Self-Review

- The design is limited to ordinary-polynomial orthogonal rank-one
  certificates.
- It reuses #191 Cohn-type certificates for all child factor formulas.
- Deterministic metadata is explicit: all coefficient table entries are
  lexicographic `i < j`, and child certificates preserve nonzero table order.
- It does not implement conjugated elementary normality certificates, Murthy,
  Quillen, ECP, Laurent, ToricBuilder, or Steinberg factor-count
  optimization.
