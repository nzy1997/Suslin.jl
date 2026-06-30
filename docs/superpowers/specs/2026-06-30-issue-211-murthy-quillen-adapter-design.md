# Issue 211 Murthy-Quillen Adapter Design

## Context

Issue #211 asks for an internal bridge from verified Murthy local `SL_3`
certificates to the data shape that Quillen local patching can consume later.
The existing Murthy local solver records `SL3LocalElementaryFactorReplay` for
local q(0)-unit and q(0)-nonunit branches. That replay already distinguishes
ordinary denominator-one factors (`mode == :ordinary`) from localized
denominator-cleared records (`mode == :denominator_cleared`).

The Quillen side currently has two relevant certificate shapes:

- `QuillenLocalRealizationCertificate`, whose `factors` field contains ordinary
  matrix factors and whose replay multiplies them directly.
- `QuillenLocalFactorSequenceCertificate`, whose factors are
  `QuillenLocalElementaryFactor` records with denominator and provenance data,
  but whose normalized replay still expects proven Quillen-compatible ordinary
  global factors.

This issue is a handoff layer only. It must not discover denominator covers,
assemble global patches, or change public factorization routing.

## Approaches Considered

1. Materialize every Murthy local replay into Quillen local factors.

   This is too broad. Localized Murthy records can carry nontrivial denominators
   whose replay is only denominator-cleared. Passing those records to an
   ordinary-factor-only Quillen API would make the localized proof look like an
   ordinary matrix factorization.

2. Add an internal adapter record that materializes only denominator-one Murthy
   replays and otherwise preserves localized replay metadata.

   This matches the issue boundary. Denominator-one certificates can become
   `QuillenLocalRealizationCertificate` values because `sl3_local_materialize_*`
   proves each record is an ordinary elementary matrix. Nontrivial localized
   replays become a guarded adapter record for #183, and ordinary-only conversion
   fails with a clear message.

3. Defer all conversion until #183.

   This avoids risk but misses the requested bridge and leaves #183 without a
   stable Murthy handoff contract.

Chosen approach: option 2.

## Design

Add an internal `MurthyQuillenLocalAdapter` record in the algorithm layer. It
will store:

- the original polynomial input and selected variable supplied by the caller;
- the verified `SL3LocalRealizationCertificate`;
- the Murthy local factor replay extracted from the certificate or rebuilt from
  its factors;
- the adapter mode, either `:ordinary_quillen_local` or
  `:localized_replay_handoff`;
- optional `QuillenLocalRealizationCertificate` data when ordinary conversion is
  proven;
- witness metadata and a replay summary that #183 can inspect without
  re-solving local `SL_3`.

The constructor will accept only `SL3LocalRealizationCertificate`. Raw vectors
of `SL3LocalElementaryFactor` are rejected by method absence/type checks. The
constructor will first call `verify_sl3_local_realization`; then it checks that
the supplied selected variable matches the certificate and that the original
input matches the certificate target for the current bridge.

For ordinary replays, every Murthy record must have denominator one and
`SL3LocalElementaryFactorReplay.mode == :ordinary`. The adapter will materialize
the ordinary matrices and build a `QuillenLocalRealizationCertificate` using a
caller-supplied or default length-one `LocalCertificate`, denominator one,
coverage multiplier one, and an elementary correction only when the Murthy
ordinary product itself is a single elementary correction. If the product is not
a single Quillen elementary correction, the adapter still records the
materialized factors in the adapter handoff but does not fabricate a single
Quillen local-realization certificate.

For localized replays, the adapter stores denominator product, cleared product,
records, local-unit witness metadata, and the Murthy branch. It must not expose
ordinary Quillen factors. A helper that requires an ordinary Quillen local
certificate will reject localized handoffs with the message that #183 must
define the localized Quillen shape first.

## Public Surface

The adapter remains internal. Function names will use leading underscores and
will not be exported from `src/Suslin.jl`.

Planned internal entry points:

- `_murthy_quillen_local_adapter(certificate, original_input, selected_variable; kwargs...)`
- `_verify_murthy_quillen_local_adapter(adapter)::Bool`
- `_murthy_quillen_local_realization_certificate(adapter)`

## Data Flow

1. Murthy solver produces `SL3LocalRealizationCertificate`.
2. Adapter verifies the Murthy certificate and extracts/rebuilds
   `SL3LocalElementaryFactorReplay`.
3. Adapter compares replay target, selected variable, stored factors, and
   denominator mode.
4. If ordinary, materialized factors and optional Quillen local certificate are
   stored.
5. If localized, replay metadata is stored and ordinary-only conversion throws.

## Error Handling

The adapter throws `ArgumentError` when:

- the certificate is not verified;
- the selected variable differs from the Murthy certificate;
- the original input differs from the current local target;
- the Murthy replay does not verify;
- a caller requests ordinary Quillen data for a localized replay;
- the caller supplies metadata that contradicts replayed denominator or witness
  data.

## Testing

Add expert tests in `test/expert/park_woodburn_quillen_route_adapter.jl` that:

- build a verified Murthy local certificate from the #206/#209/#210 fixture
  family;
- adapt an ordinary denominator-one Murthy certificate and verify the adapter
  replay, selected variable, local product, and witness metadata;
- adapt a localized Murthy certificate and verify that it stays a
  localized-replay handoff;
- reject a tampered Murthy local factor before constructing Quillen data;
- reject a mismatched selected variable;
- reject ordinary-only conversion for non-materializable localized records.

Preserve `test/expert/quillen_local_certificate.jl` and run it unchanged to
guard the #100 Quillen local certificate behavior.

## Out of Scope

- General Quillen patch assembly.
- Denominator cover discovery.
- Public `elementary_factorization` route changes.
- A complete #183 localized Quillen certificate shape.

