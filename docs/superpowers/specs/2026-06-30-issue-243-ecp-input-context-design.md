# Issue 243 ECP Input Context Design

## Context

Issue #243 adds an internal checked input boundary for ordinary polynomial ECP columns. The current reducer validates a column and immediately tries the supported exact routes. Parent #185 needs one reusable object that records the validated column, ring metadata, selected ECP variable hints, unimodularity witness, and the staged diagnostic boundary for unsupported-but-unimodular cases.

Issue #242 is already merged on `main` and provides the Park-Woodburn ECP mainline fixture catalog used here for supported and staged-positive examples.

## Approach Options

1. Add a minimal internal `ECPInputContext` with a replay verifier.
   This follows the existing certificate pattern: construct a provisional object, recompute verification metadata, store it, and make `verify_ecp_input_context` compare the stored verification with a fresh replay.

2. Extend `diagnose_unimodular_column_reduction` to return a richer diagnostic object.
   This would centralize some data but would mix validation context with the public diagnostic surface and would not naturally store selected-variable hints or witnesses.

3. Thread context fields through the existing reducer calls without a new object.
   This avoids a new type but keeps validation scattered, which is the problem #243 is meant to solve.

The chosen design is option 1. It is internal, scoped, and mirrors existing replay-verified ECP certificate code.

## Design

Add an unexported `ECPInputContext` struct in `src/algorithm/column_reduction.jl` with fields for:

- the coerced one-based column;
- the input ring and `_column_reduction_ring_profile(R)`;
- `variables = tuple(gens(R)...)` and normalized `variable_order`;
- `column_length`;
- a recomputed `_unimodular_witness`;
- `selected_variable_index` and `selected_variable` when supplied;
- `support_classification`, `staged_failure_reason`, and the full staged diagnostic from `diagnose_unimodular_column_reduction`;
- replay verification metadata.

Provide `ecp_input_context(v, R; variable_order = tuple(gens(R)...), selected_variable = nothing, unimodularity_witness = nothing)` plus an outer `ECPInputContext(v, R; ...)` constructor. These names are intentionally not exported.

The constructor rejects Laurent rings and unsupported preconditions before storing a context. It uses `_validated_unimodular_column`, `_column_reduction_ring_profile`, `_ecp_normalize_variable_order`, `_ecp_selected_variable_index`, `_unimodular_witness`, and `diagnose_unimodular_column_reduction`. A supplied witness is accepted only as a checked hint: it must be one-based, length-matched, coercible into `R`, and sum with the column to `one(R)`. The stored witness remains the canonical recomputed witness so the verifier can replay it exactly.

`verify_ecp_input_context(context)` recomputes every stored field from the stored column and ring:

- one-based indexing and length;
- coerced column entries;
- ring profile and generators;
- normalized variable order;
- selected-variable generator membership and variable-order membership;
- canonical unimodularity witness and witness identity;
- support classification, staged failure reason, and full diagnostic metadata.

The verifier returns `false` for tampered context fields and catches ordinary malformed objects without throwing, except for interrupts.

## Tests

Add `test/expert/ecp_input_context.jl` and register it in the expert group. The focused verification command is:

```bash
julia --project=. -e 'include("test/expert/ecp_input_context.jl")'
```

The test covers:

- one #242 supported ordinary multivariate column from `ECPMainlineFixtureCatalog`;
- one existing GF(2) staged ECP column from `ECPColumnFixtureCatalog`;
- one unsupported but unimodular ordinary #242 length-four boundary column;
- negative controls for non-unimodular input, length-two input, selected variable outside `gens(R)`, invalid supplied witness, and tampered stored witness.

## Out Of Scope

This change does not produce elementary factors, extract link witnesses, call #184 `SL_3` routes, alter Laurent behavior, export new public API names, or change public reducer behavior.

## Self-Review

- No placeholders remain.
- The chosen approach keeps the new API internal and follows existing replay-verifier patterns.
- The staged boundary is diagnostic-only and does not attempt unsupported reduction work.
- The tests map directly to the issue verification and negative controls.
