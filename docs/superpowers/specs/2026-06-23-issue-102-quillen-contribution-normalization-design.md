# Issue 102 Quillen Contribution Normalization Design

## Context

Issue #102 is the normalization child of the constructive Quillen patching
parent issue #63. Its direct dependencies are merged in this checkout:

- #100 added replayable `QuillenLocalRealizationCertificate` records and
  `verify_quillen_local_certificate`.
- #101 added exact `QuillenDenominatorCoverCertificate` records and
  `verify_quillen_denominator_cover`.

The existing public `QuillenLocalContribution` type is already used by exact
patching tests. This issue needs a richer replayable normalization record, but
it does not need to change final global patch construction.

No `AGENTS.md`, `CLAUDE.md`, or `GEMINI.md` file is present in this checkout.

## Approaches Considered

1. Add a richer expert/internal normalization record in
   `src/algorithm/quillen_induction.jl`. The record consumes a replayable local
   realization certificate and a verified cover certificate, binds the matching
   denominator/coverage multiplier pair, recomputes patched-substitution data,
   local products, and weighted global elementary factors, and exposes a replay
   verifier.
2. Extend the existing public `QuillenLocalContribution` type with replay
   fields. This risks breaking exact patching compatibility and expands a
   public shape before final patch assembly exists.
3. Route normalization through `elementary_factorization` or a fresh local
   solver. This is out of scope and would make normalization depend on new
   solving behavior instead of replaying supplied deterministic certificates.

The chosen design is approach 1. It reuses the #100 and #101 schemas, preserves
existing patching compatibility, and keeps the new surface expert/internal.

## API Surface

The new names remain expert/internal and are intentionally not exported from
`src/Suslin.jl`:

- `QuillenLocalContributionNormalization`
- `QuillenLocalContributionNormalizationVerification`
- `normalize_quillen_local_contribution(certificate, cover; original_input = certificate.original_input, selected_variable = certificate.selected_variable, patched_substitution_witness = certificate.patched_substitution_witness)`
- `normalize_quillen_local_contributions(certificates, cover; original_input = nothing, selected_variable = nothing)`
- `replay_quillen_local_contribution_normalization(normalized)`
- `verify_quillen_local_contribution_normalization(normalized)`

Expert tests use qualified `Suslin.<name>` access. `test/public/api_surface.jl`
does not change.

## Normalization Semantics

`normalize_quillen_local_contribution` accepts one
`QuillenLocalRealizationCertificate` and one
`QuillenDenominatorCoverCertificate`. Construction requires both inputs to
verify. It normalizes the local contribution against the certificate ring and
size, finds the exact cover pair whose denominator and coverage multiplier
match the local certificate, recomputes the weighted global elementary factor
from `_quillen_factors`, and records:

- the denominator,
- the coverage multiplier,
- the denominator/coverage multiplier index in the cover,
- the selected Quillen substitution variable,
- the patched-substitution witness data and replay summary,
- the local product and local correction from the local certificate replay,
- the weighted global elementary factor data,
- replay metadata tying the record back to the local certificate and cover.

`normalize_quillen_local_contributions` applies the same normalization to a
collection of certificates and returns the normalized records in input order.
When an explicit `selected_variable` is supplied, every local certificate must
use that generator.

`replay_quillen_local_contribution_normalization` recomputes the local
certificate replay, cover replay, cover pairing, patched substitution, local
product, weighted global elementary factor, and stored metadata. The verifier
returns `true` only when the replay matches the stored normalization and all
component verifiers succeed.

Construction throws `ArgumentError` for unverified local certificates,
unverified cover certificates, ring/coverage mismatches, selected-variable
mismatches, cover-pair mismatches, or patched-substitution witness mismatches.
Verification catches the same tampering and returns `false`.

## Tests

Add `test/expert/quillen_contribution_normalization.jl`. The test includes the
#99 fixture catalog and builds replayable local certificates from #100 against
a verified cover from #101. Positive coverage normalizes both local entries
from `quillen-patched-substitution-witness-qq`; each normalized record must
replay:

- the patched substitution,
- the local factor product,
- the weighted global elementary factor,
- the exact denominator pairing from the cover certificate.

Negative controls tamper the selected substitution variable, substitution
exponent, substitution shift, coverage multiplier, and weighted global factor.
Each tampered record must fail verification or construction.

Register the new expert test in `test/runtests.jl`.

## Files

- Modify `src/algorithm/quillen_induction.jl`: add normalization structs,
  constructor, batch helper, replay helper, and verifier.
- Create `test/expert/quillen_contribution_normalization.jl`: focused positive
  and negative normalization coverage.
- Modify `test/runtests.jl`: register the expert test.
- Do not modify `src/Suslin.jl` or `test/public/api_surface.jl`.

## Verification

Focused command:

```bash
julia --project=. -e 'include("test/expert/quillen_contribution_normalization.jl")'
```

Expert group command:

```bash
julia --project=. test/runtests.jl expert
```

Package command required by Agent Desk:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Automatic Approval

This Agent Desk run is non-interactive. Under the Standing Answer Policy, the
design is approved automatically because it chooses the conservative internal
normalization API, reuses merged dependency schemas, avoids public API churn,
and avoids final global patch construction.

## Spec Self-Review

- No incomplete markers remain.
- Scope is limited to contribution normalization and replay.
- The API surface decision is explicit and avoids exported names.
- Negative controls cover variable, exponent, shift, coverage multiplier, and
  weighted-factor tampering.
