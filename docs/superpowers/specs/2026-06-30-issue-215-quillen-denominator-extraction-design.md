# Issue 215 Quillen Denominator Extraction Design

## Context

Issue #215 sits between the #214 local factor-sequence certificate and later
cover solving. The repository already has checked denominator-cover
certificates, local sequence certificates, and replay helpers in
`src/algorithm/quillen_induction.jl`. The missing stage is an internal record
that extracts the raw local support denominators from verified local sequence
certificates without proving that those denominators cover the ring.

The input is a nonempty collection of verified
`QuillenLocalFactorSequenceCertificate` values for the same original input,
ring, matrix size, and selected variable. The output is a denominator-cover
candidate containing the raw denominators `r_i`, the per-local factor entries
that justify each one, replay metadata, and the original input context.

## Approaches Considered

Recommended: add a new internal `QuillenDenominatorCoverCandidate` layer. Each
local sequence certificate is verified first, then a
`QuillenLocalDenominatorSupport` is extracted from its structured factors. The
support records the factor denominators, factor entries, factor provenance, the
support kind, and a replayed equality showing that the local support
denominator is justified by the factor data. The candidate verifier replays the
same extraction and compares it to the stored candidate.

Alternative: construct a `QuillenDenominatorCoverCertificate` directly from the
extracted denominators. This is out of scope because #216 must choose an
exponent and prove an ideal-membership cover identity later.

Alternative: infer denominators from the materialized local matrices. This
does not meet the issue because #215 must use factor-level denominator
provenance from #214, not opaque matrices.

The chosen design is the recommended internal candidate layer. It keeps raw
support denominators separate from cover powers and makes tampering with stored
candidate denominator data or factor-level denominator data visible to the
verifier.

## Candidate Shape

Add these internal records in `src/algorithm/quillen_induction.jl`:

- `QuillenLocalDenominatorSupport`: one local realization's denominator
  support, including support kind, support denominator, factor denominators,
  factor entries, factor provenance, replayed denominator, replay equality, and
  replay status.
- `QuillenDenominatorCoverCandidateVerification`: replay metadata for a full
  candidate.
- `QuillenDenominatorCoverCandidate`: the extracted candidate for later cover
  solving.

The candidate records:

- `original_input`, `ring`, `size`, and `selected_variable`;
- the verified local sequence certificates used for extraction;
- `raw_denominators`, one local support denominator per input certificate;
- `local_supports`, one `QuillenLocalDenominatorSupport` per input
  certificate;
- replay metadata with factor-level provenance and local verification
  summaries;
- a verification object whose booleans explain whether context alignment,
  local certificate verification, support replay, raw-denominator replay, and
  metadata replay all succeeded.

The extraction function is `extract_quillen_denominator_cover_candidate`.
Verification is done by `replay_quillen_denominator_cover_candidate` and
`verify_quillen_denominator_cover_candidate`.

## Support Semantics

For #215, the local support denominator is the #214 product denominator:

```julia
prod(factor.denominator for factor in certificate.factors) ==
    certificate.product_denominator
```

The support kind is recorded as `:product`. This is the conservative choice
because #214 explicitly exposes a per-local product denominator. No common
denominator simplification, supplied support denominator, cover power, or
ideal-membership proof is introduced in this issue.

If later code supplies common-denominator or supplied-support evidence, it can
extend the support-kind replay without changing the candidate-level contract:
stored support data must still replay from structured factor entries.

## Validation

Construction rejects:

- an empty local sequence collection;
- any local sequence certificate that does not verify;
- mixed original inputs, rings, sizes, or selected variables;
- local support data whose product equality does not replay.

Verification rejects tampering after construction. In particular, changing a
candidate raw denominator, changing a stored support denominator, changing a
factor-level denominator in the stored support, dropping a local certificate
from a candidate, or mixing a certificate with a different selected variable
causes replay to disagree with the stored candidate or context checks to fail.

## Files

- Modify `src/algorithm/quillen_induction.jl`.
- Create `test/expert/quillen_denominator_extraction.jl`.
- Modify `test/runtests.jl` to register the new expert test after
  `expert/quillen_local_factor_sequence.jl`.

No cover solving, exponent choice, global factor assembly, Murthy solving, or
public `QuillenPatch` API change is included.

## Verification

Focused extraction verifier:

```bash
julia --project=. -e 'include("test/expert/quillen_denominator_extraction.jl")'
```

Package entry point:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Diff hygiene:

```bash
git diff --check
```

## Automatic Approval

This Agent Desk run is non-interactive. Under the Standing Answer Policy, the
design is approved automatically because it follows the issue's recommended
structured extraction, keeps raw denominators separate from later cover powers,
and avoids out-of-scope cover solving.

## Spec Self-Review

- No placeholders remain.
- The design reads #214 sequence records and never infers from opaque matrices.
- Raw denominators are separate from later powers and cover multipliers.
- Factor-level denominator and provenance data are stored and replayed.
- Negative controls cover dropped certificates, edited candidate denominators,
  edited factor-level denominators, and mixed selected variables.
