# Issue 214 Quillen Local Factor Sequence Design

## Context

Issue #214 is the representation step between the #213 Park-Woodburn mainline
fixture catalog and later #183 patch assembly work. The current
`QuillenLocalRealizationCertificate` stores one normalized
`QuillenElementaryCorrection` and proves one weighted elementary factor. A real
local realization can be an ordered sequence of localized elementary factors,
and later cover extraction needs to inspect each factor's denominator and
provenance without trusting fixture ids or opaque matrices.

The local code already has replayable denominator-cover certificates, local
realization certificates, contribution normalization, and deterministic patch
assembly in `src/algorithm/quillen_induction.jl`. The #213 catalog now exposes
factor evidence and raw denominator provenance in
`test/fixtures/quillen_mainline_cases.jl`.

## Approaches Considered

Recommended: add a new internal `QuillenLocalFactorSequenceCertificate` layer
with structured factor records. Each record carries `row`, `col`, `numerator`,
`denominator`, `coverage_multiplier`, `provenance`, and a derived
`LocalCertificate`. The verifier replays elementary matrices from those fields,
preserves order, extracts raw denominators and a product denominator, and stores
normalized contribution records plus weighted global elementary factors for
later assembly.

Alternative: extend `QuillenLocalRealizationCertificate` so its `factors` field
accepts sequence metadata. This would blur the existing single-correction
contract and risks breaking the normalization code that assumes one
`QuillenElementaryCorrection`.

Alternative: store only a vector of already materialized matrices and annotate
it with denominator metadata. This would not satisfy the issue: the verifier
would still be trusting opaque matrices rather than replaying from
`row`/`col`/`numerator`/`denominator`.

The chosen design is the recommended internal sequence certificate because it
keeps existing behavior intact, makes denominator-bearing entries explicit, and
creates the representation that #215/#218 can consume without solving local
data again.

## Certificate Shape

Add these internal records in `src/algorithm/quillen_induction.jl`:

- `QuillenLocalElementaryFactor`: one localized elementary factor record with
  factor-level provenance and optional metadata.
- `QuillenLocalFactorSequenceCertificate`: the ordered sequence certificate.
- `QuillenLocalFactorSequenceVerification`: replay output used by the verifier.

The sequence certificate records:

- `factors`: ordered `QuillenLocalElementaryFactor` records;
- `raw_denominators`: one raw denominator per factor;
- `product_denominator`: the exact product of raw denominators;
- `local_product`: the product replayed from the ordered structured factors;
- `local_correction`: the expected local product/correction matrix;
- `normalized_local_contributions`: `QuillenLocalContribution` records derived
  from each factor;
- `normalized_global_elementary_factors`: weighted elementary matrices from the
  existing normalization layer;
- optional patched-substitution and chain witness metadata;
- replay metadata containing ordered provenance, denominator data, and witness
  metadata.

Existing `QuillenLocalRealizationCertificate` values convert to a length-one
sequence certificate. The conversion uses the stored correction entry as the
factor numerator, the stored denominator and coverage multiplier as factor
metadata, the stored local certificate, and the existing local correction as
the expected sequence product.

## Validation

Construction normalizes every factor over the selected ordinary polynomial
ring, validates elementary indices, coerces numerator and denominator entries,
derives or checks a local certificate, and requires nonempty provenance. If a
provenance record supplies `factor_index`, `sequence_index`, or `local_index`,
the verifier requires it to match the ordered position. This makes factor-order
tampering visible even for commuting elementary factors in the small fixture
catalog.

Replay never trusts stored matrices. It rebuilds normalized
`QuillenLocalContribution` records from factor data, rebuilds elementary
matrices with the existing `_quillen_factors` normalization path, multiplies
them in order, recomputes denominator data and replay metadata, and compares
the result to the stored certificate fields.

Negative controls must return `false` from the verifier or throw
`ArgumentError` during construction when any of the factor entry, denominator,
selected variable, provenance, or factor order is corrupted.

## Files

- Modify `src/algorithm/quillen_induction.jl`.
- Create `test/expert/quillen_local_factor_sequence.jl`.
- Modify `test/expert/quillen_local_certificate.jl` for length-one conversion
  coverage.
- Modify `test/runtests.jl` to register the new expert test.
- Add `docs/superpowers/plans/2026-06-29-issue-214-quillen-local-factor-sequence.md`.

No exported toy `QuillenPatch` API, public `elementary_factorization` routing,
denominator-cover discovery, ideal-membership solving, Murthy solving, or
public patch API behavior changes are included.

## Verification

Focused sequence verifier:

```bash
julia --project=. -e 'include("test/expert/quillen_local_factor_sequence.jl")'
```

Existing local-certificate regression:

```bash
julia --project=. -e 'include("test/expert/quillen_local_certificate.jl")'
```

Package entry point:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Automatic Approval

This Agent Desk run is non-interactive. Under the Standing Answer Policy, the
design is approved automatically because it chooses the structured internal
certificate recommended by the issue, avoids public API churn, preserves
existing single-correction behavior, and keeps solver/cover/Murthy work out of
scope.

## Spec Self-Review

- No placeholders remain.
- The design distinguishes ordered factors, raw denominators, product
  denominator, replayed product/correction, and normalized contribution data.
- The verifier replays from structured factor fields rather than trusting a
  vector of matrices.
- The conversion path from existing single-correction certificates is explicit.
- Negative controls cover factor data, denominator data, variable metadata,
  provenance, and order.
