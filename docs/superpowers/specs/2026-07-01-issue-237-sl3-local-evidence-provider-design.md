# Issue 237 SL3 Local Evidence Provider Design

## Goal

Add an internal provider that starts from a verified `SL3LocalFormWitnessSelection`, asks the existing Murthy local solver for the matching local certificate, adapts that certificate with the existing Murthy-to-Quillen adapter, and returns evidence records aligned to the original `SL_3` driver context.

## Chosen Approach

Use a thin coordination layer in `src/algorithm/factorization.jl`. The provider will not compute Murthy branches, denominator covers, Quillen local-to-global patch assemblies, or transformed local-form composition itself. It will validate the #235 input context and #236 witness selection, construct `SL3LocalMurthyInputContext` through `sl3_local_murthy_input_context`, realize it through `realize_sl3_local_certificate`, adapt it through `_murthy_quillen_local_adapter`, and, only for ordinary materializable adapter mode, expose the verified Quillen local sequences that #219 consumers can use.

Alternatives considered:

- Put provider logic in `quillen_induction.jl`. This would keep all Quillen-facing records together, but it would make the Quillen layer know about #184 driver context records.
- Extend the Murthy solver to emit driver-context evidence directly. This would couple #182 to #184 and duplicate adapter responsibilities.
- Build the provider in `factorization.jl`. This matches existing `SL3RealizationInputContext` and witness-selection ownership and keeps the new code as context/provenance glue. This is the selected approach.

## Provider Shape

Add `SL3MurthyQuillenLocalEvidenceProvider` with these fields:

- original `SL3RealizationInputContext`
- verified `SL3LocalFormWitnessSelection`
- `SL3LocalMurthyInputContext`
- `SL3LocalRealizationCertificate`
- `MurthyQuillenLocalAdapter`
- optional ordinary Quillen local sequences
- local product, selected variable, denominator metadata, provider replay metadata, staged diagnostic, and verification

The constructor `_sl3_murthy_quillen_local_evidence_provider(selection; ...)` will require a verified witness selection whose support is `:supported` and whose #236 witness payload is replayable and bound to the #235 context. It will reject context/witness mismatches before returning Quillen-facing evidence. It will accept optional Murthy witness data (`witness`, `local_unit_witnesses`, `split_witness`, `bezout_witness`) and pass them to existing Murthy APIs unchanged.

## Data Flow

1. Verify `selection` with `_verify_sl3_local_form_witness_selection`.
2. Require `selection.context` to verify and require `selection.local_form_matrix == context.matrix` for this provider. This keeps the issue scoped to driver contexts whose selected #236 witness is already the local target for the #235 matrix.
3. Build the Murthy input context from `selection.entries` and `selection.selected_variable`.
4. Realize the Murthy certificate using `realize_sl3_local_certificate`.
5. Adapt the certificate with `_murthy_quillen_local_adapter`, carrying adapter witness metadata for #235 context, #236 witness, selected variable, and selection status; the provider replay metadata carries the provider source, local product, denominator product, and original matrix identity.
6. If the adapter mode is `:ordinary_quillen_factor_sequence`, expose the verified local sequences for the #219 consumption layer.
7. If the adapter mode is `:localized_replay_handoff`, return a verified staged provider record with a diagnostic explaining that localized denominator-cleared replay is not ordinary Quillen factor-sequence evidence.

## Verification

The provider verification will recompute all provider fields from the stored selection and optional inputs. It will check:

- context and selection still verify
- selected variable and selected-variable index match exactly
- local-form entries and local product match the context matrix
- Murthy input context, certificate, and adapter replay
- denominator product and denominator metadata match adapter/sequence replay
- provider metadata ties back to #235 context metadata and #236 witness metadata
- ordinary adapter records expose verified local sequences and localized records do not enter ordinary Quillen sequences

Tampering with the #236 witness record, Murthy local factors, selected variable, denominator/local-unit witness, or #235 context identity must make provider verification fail or provider construction throw before Quillen patch assembly can receive the evidence.

## Tests

Add `test/expert/park_woodburn_sl3_local_evidence_provider.jl` for the new provider. It will build a non-catalog multivariate `SL_3` driver context, select the #236 local-form witness, produce provider evidence, assert exact local product, selected variable, denominator metadata, Murthy provenance, #236 witness provenance, and #235 matrix metadata, then exercise negative controls.

Extend `test/expert/park_woodburn_quillen_route_adapter.jl` with a small provider-to-#219 assertion for the ordinary path and a localized diagnostic assertion for a Murthy fixture that cannot materialize ordinary sequence evidence.

Add the new expert file to `test/runtests.jl`.
