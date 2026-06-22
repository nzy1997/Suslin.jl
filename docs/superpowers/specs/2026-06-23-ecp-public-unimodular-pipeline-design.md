# ECP Public Unimodular Pipeline Design

## Goal

Wire the supported elementary-column-property stages into the public `reduce_unimodular_column(v, R)` path for fixture-backed exact polynomial columns, while preserving the existing factor-returning public API and legacy reducer behavior for unsupported or non-ECP inputs.

## Architecture

The public reducer will validate unimodularity first, then try a narrow staged ECP route for exact ordinary polynomial fixture families with deterministic in-tree witness data. If the staged route does not recognize the column, the reducer falls back to the existing `ecp_column_reduction_certificate(v, R)` path, preserving existing exact small-column, embedded-block, and Laurent behavior.

The staged route will be exposed through a non-exported expert API, `ecp_staged_column_reduction_certificate(v, R; ...)`, and a verifier, `verify_ecp_staged_column_reduction(certificate)`. This API will compose the existing #87 `ecp_link_witness`, #88 `ecp_link_step_certificate`, and #90 `ecp_induction_normality_certificate` verifiers instead of duplicating unchecked replay logic.

## Data Flow

1. `_validated_unimodular_column(v, R)` rejects malformed or non-unimodular input before any route attempts.
2. `_ecp_public_staged_reduction_certificate(column, R)` recognizes supported exact ordinary polynomial fixture families and supplies deterministic link-witness metadata.
3. `ecp_staged_column_reduction_certificate` builds a verified link-step certificate, obtains or verifies the lower-variable reduction, builds deterministic normality-witness data when one is not supplied, and calls `ecp_induction_normality_certificate`.
4. The staged certificate stores the route, final factors, final column, and replay summary. Its verifier checks the stored route by reusing the existing certificate verifiers and exact product replay.
5. `reduce_unimodular_column(v, R)` returns the staged certificate factors when this route applies, otherwise returns the existing certificate factors.

## Error Handling

Non-unimodular columns continue to throw `ArgumentError("v must be a unimodular column")` before any unsupported-route work. Unsupported ordinary or Laurent columns keep the existing staged `ArgumentError` behavior. Corrupted supplied link-witness or normality-witness data must fail through the reused replay constructors and throw `ArgumentError`; the public reducer must never return factors from unverified staged data.

## Scope

This design intentionally supports only the already merged fixture families from #87, #88, and #90. It does not implement Quillen local-to-global patching, general Park-Woodburn witness extraction, Laurent determinant correction, or factor minimization.

## Tests

Add `test/expert/elementary_column_property.jl` and register it in the expert group. The focused test will cover at least four GF(2)[x,y] hard-slice permutations that have no unit entries and do not pass the witness-unit shortcut. It will assert exact public factor replay, inspect the full staged route for the canonical case, check inverse-substitution factors over the original ring for a permuted case, and verify immediate non-unimodular rejection plus corrupt replay failure controls.

## Self-Review

No placeholders remain. The design is intentionally narrow and consistent with prior staged ECP PRs. Public API compatibility is preserved because `reduce_unimodular_column(v, R)` still returns only factors, and new inspection names remain expert/internal.
