# Issue 219 Murthy Adapter Consumption Design

## Context

Issue #219 connects the landed Murthy-to-Quillen adapter from #211 to the
supplied-evidence Quillen patch assembler from #218. The codebase already has
these pieces:

- `MurthyQuillenLocalAdapter`, built from a verified
  `SL3LocalRealizationCertificate`;
- `_verify_murthy_quillen_local_adapter`, which replays the Murthy certificate,
  local factor replay, materialized factor alignment, selected variable, and
  adapter metadata;
- `QuillenLocalFactorSequenceCertificate`, the #214 local sequence certificate;
- `assemble_quillen_patch_from_local_evidence`, the #218 constructor that
  assembles a replayable global patch from supplied local sequence evidence.

The missing layer is a consumer that treats #211 adapter records as the input
evidence. It must not call the Murthy solver or accept raw factor vectors. It
must verify that every adapter belongs to the requested matrix and selected
variable before handing sequence certificates to the #218 assembler.

## Approaches Considered

Recommended: add an internal adapter-consumption certificate in
`src/algorithm/quillen_induction.jl`. It accepts `MurthyQuillenLocalAdapter`
records, verifies each adapter, extracts its #214 sequence certificate, records
the Murthy provenance for replay, and delegates assembly to
`assemble_quillen_patch_from_local_evidence`. This keeps #219 additive and
evidence-backed.

Alternative: let `assemble_quillen_patch_from_local_evidence` accept either
sequence certificates or Murthy adapters. This would make the supplied-evidence
assembler responsible for a Murthy-specific contract and would blur ownership
between #211 and #218.

Alternative: reconstruct Murthy certificates from matrices inside the Quillen
patch route. This duplicates Murthy branch logic and violates the issue's
evidence-backed boundary.

Chosen approach: the recommended additive internal consumer.

## Consumer Shape

Add two internal records:

- `QuillenMurthyAdapterConsumptionVerification`, a replay summary for the
  adapter-consumption boundary.
- `QuillenMurthyAdapterConsumption`, the consumer output. It records the
  original input matrix, selected variable, verified Murthy adapters, converted
  local sequence certificates, assembled supplied-evidence patch, replay
  metadata, and verification object.

Add internal entry points:

- `quillen_local_sequences_from_murthy_adapters(A, X, adapters)`;
- `assemble_quillen_patch_from_murthy_adapters(A, X, adapters; kwargs...)`;
- `replay_quillen_murthy_adapter_consumption(consumption)`;
- `verify_quillen_murthy_adapter_consumption(consumption)::Bool`.

These names remain internal by convention and are not exported from
`src/Suslin.jl`. Expert tests call them through `Suslin.`.

## Validation

The consumer verifies:

- the input is a supported ordinary-polynomial matrix and `X` is a generator of
  the matrix ring;
- every supplied item is a `MurthyQuillenLocalAdapter`;
- `_verify_murthy_quillen_local_adapter(adapter)` is true;
- adapter ring, size, original input, and selected variable match the requested
  `A` and `X`;
- adapter replay metadata matches `_murthy_quillen_local_replay_metadata`;
- adapter local factor replay target, selected variable, factor order,
  denominator product, cleared product, and materialized factors align with the
  Murthy certificate and adapter fields;
- the adapter contains an ordinary #214
  `QuillenLocalFactorSequenceCertificate`; localized denominator-cleared
  handoffs throw a staged `ArgumentError`;
- each converted sequence certificate verifies, keeps the same original input
  and selected variable, exposes denominator data through its replay metadata,
  and carries factor provenance whose source is the Murthy adapter path.

After validation, the consumer calls
`assemble_quillen_patch_from_local_evidence` with the converted sequence
certificates and passes through the existing #218 keyword arguments for cover
search, substitution-chain override, base-term policy, base-term factors, and
metadata.

## Replay Metadata

The consumer replay metadata records:

- `source = :quillen_murthy_adapter_consumption`;
- adapter count;
- selected variable;
- per-adapter Murthy replay metadata from #211;
- per-sequence replay metadata from #214;
- assembled patch metadata from #218;
- caller metadata.

This lets a reviewer trace each local sequence factor back to the Murthy
certificate and local factor replay that produced it.

## Error Handling

Throw `ArgumentError` before assembly when:

- no adapters are supplied;
- a supplied item is not a `MurthyQuillenLocalAdapter`;
- an adapter fails its #211 replay verifier;
- original target matrix, ring, size, or selected variable do not match;
- factor order, local denominators, materialized factors, or products no longer
  match adapter replay metadata;
- an adapter contains localized denominator-cleared replay instead of ordinary
  sequence evidence;
- a converted sequence certificate lacks denominator provenance or fails
  replay.

Rely on #218 to throw its existing staged error when the `A(0)` base-term policy
or supplied base-term evidence is missing or contradictory.

## Tests

Add `test/expert/quillen_murthy_adapter_consumption.jl`.

Positive coverage:

- build at least one verified #211 ordinary Murthy adapter record from the
  Murthy fixture catalog;
- convert it to #214 local sequence evidence through the new consumer helper;
- assemble a #218 patch through `assemble_quillen_patch_from_murthy_adapters`;
- verify exact factor multiplication, selected variable, raw denominators,
  local product, replay metadata, sequence factor provenance, and adapter
  provenance.

Negative controls:

- tamper a Murthy adapter local factor;
- tamper selected variable;
- tamper denominator/local-witness metadata;
- tamper original target;
- pass a localized denominator-cleared adapter and assert the staged conversion
  error;
- omit base-term evidence/policy and assert the #218 staged error is raised
  after adapter validation.

Update `test/expert/park_woodburn_quillen_route_adapter.jl` with a focused
regression that existing #211 adapters remain accepted by the new consumer, and
register the new expert test in `test/runtests.jl`.

## Files

- Modify `src/algorithm/quillen_induction.jl`.
- Modify `test/expert/park_woodburn_quillen_route_adapter.jl`.
- Create `test/expert/quillen_murthy_adapter_consumption.jl`.
- Modify `test/runtests.jl`.
- Add `docs/superpowers/plans/2026-06-30-issue-219-murthy-adapter-consumption.md`.

Do not broaden Murthy solving, add a public factorization route, implement the
general `SL_3` path, ECP, recursive `SL_n`, Laurent support, or Steinberg
optimization.

## Verification

Focused Murthy adapter consumption:

```bash
julia --project=. -e 'include("test/expert/quillen_murthy_adapter_consumption.jl")'
```

Existing route adapter:

```bash
julia --project=. -e 'include("test/expert/park_woodburn_quillen_route_adapter.jl")'
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
additive internal consumer design is approved automatically because it follows
the landed #211 and #218 contracts, keeps Murthy solving out of the Quillen
consumer, and rejects unsupported localized handoffs before assembly.

## Spec Self-Review

- No placeholders remain.
- The design matches the landed #211 adapter contract rather than an abandoned
  draft shape.
- The consumer validates every issue-required alignment before assembly.
- Base-term handling remains delegated to the #218 assembler with staged
  errors preserved.
- Negative controls cover local factors, selected variables, denominator
  metadata, original targets, localized adapters, and base-term policy.
