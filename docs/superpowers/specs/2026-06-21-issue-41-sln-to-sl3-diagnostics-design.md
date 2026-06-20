# Issue 41 SL_n To SL_3 Diagnostics Design

## Goal

Add structured diagnostics for failed `SL_n -> local SL_3` reductions, focused on Laurent determinant-one cores from Issue #38 and reusable for later staged algorithm triage.

## Context

`reduce_sln_to_sl3(A; block_locations=nothing)` currently normalizes the input, walks default or supplied 3-coordinate blocks, and throws a staged `ArgumentError` if `realize_sl3_local` cannot solve a local obligation. The public error text intentionally hides the lower-level local solver reason. Issue #41 needs the hidden reason as data: block location, determinant status, local shape, and whether any alternative `3+3` coordinate partition works for the `6 x 6` Issue #38 core.

The Issue #39 fixture now provides the exact row-normalized and column-normalized determinant-one cores for the Issue #38 ToricBuilder `Q` block, plus a synthetic `40 x 40` Laurent block-local acceptance fixture that should remain accepted.

## Approaches Considered

1. Add a public diagnostic function and record types beside `reduce_sln_to_sl3`.

   This keeps existing failures stable while giving tests and triage code a structured API. The diagnostic path can reuse the same validation, normalization, block-location handling, and local solver helpers as the reduction path. This is the recommended approach because it is explicit, low risk, and does not force callers to catch a new exception type.

2. Replace staged `ArgumentError` with a custom exception carrying diagnostics.

   This attaches data directly to failure objects, but it changes the current error surface and risks breaking tests that only require stable staged failure text. It also makes diagnostics available only after a thrown failure, not as a proactive probe.

3. Parse local solver error strings after a failed reduction.

   This is small but brittle. It would not reliably report determinant status, partition search, or future local failure codes, and it would couple tests to user-facing prose.

## Design

Expose `diagnose_sln_to_sl3_reduction(A; block_locations=nothing, search_partitions=true)`. It returns an `SLNToSL3ReductionDiagnostic` with stable symbol fields:

- `status`: `:success` or `:failure`.
- `failure_code`: `nothing`, `:determinant_failure`, `:local_shape_failure`, `:local_solver_failure`, or `:reassembly_failure`.
- `determinant_status`: `:determinant_one`, `:determinant_not_one`, `:determinant_requires_correction`, or `:determinant_check_failed`.
- `determinant_classification`: the Laurent normalization classification when available.
- `block_diagnostics`: local records for examined non-identity obligations.
- `partition_search`: a named tuple reporting whether alternative `3+3` partitions were searched, how many were attempted, which partitions succeeded, and a stable status code.
- `message`: the underlying failure message when one exists.

Each `SL3LocalReductionDiagnostic` records:

- `block_location`: the exact 3-coordinate block.
- `status`: `:success` or `:failure`.
- `failure_code`: `nothing`, `:local_shape_failure`, or `:local_solver_failure`.
- `determinant_status`: determinant status for the local `3 x 3` target.
- `local_shape_reason`: `:embedded_2x2_with_trailing_identity` or `:not_embedded_2x2_with_trailing_identity`.
- `solver_status`: `:success`, `:not_attempted`, or `:failure`.
- `message`: local solver message when present.

The diagnostic path should leave `reduce_sln_to_sl3` errors compatible with current tests. The existing staged failure message must still include `failed to solve local SL_3 obligation on block [1, 2, 3]` for the Issue #38 row and column cores.

## Partition Search

When the normalized matrix is `6 x 6`, a local failure occurs, and `search_partitions=true`, diagnostics should try every unordered coordinate partition into two disjoint triples. A partition succeeds only if the existing reduction path succeeds for those two block locations, including exact reconstruction. For Issue #38 row and column cores, the search should attempt all 10 unordered partitions and report no successful partition.

For other sizes, or when the primary diagnostic path succeeds, partition search should report `searched=false`. This prevents an accidental combinatorial search on large fixtures such as the `40 x 40` acceptance case.

## Testing

Add `test/expert/sln_to_sl3_diagnostics.jl` and register it in the expert group. The test should:

- Probe the Issue #38 row-normalized and column-normalized cores from `test/fixtures/toricbuilder_issue38_cases.jl`.
- Assert global determinant status is `:determinant_one`.
- Assert the first failing block is `[1, 2, 3]`.
- Assert the local failure code is `:local_shape_failure` with reason `:not_embedded_2x2_with_trailing_identity`.
- Assert alternative `3+3` partitions were searched and no partition succeeded.
- Assert the existing staged reduction error text remains recognizable.
- Probe the synthetic `40 x 40` block-local Laurent acceptance fixture from `test/fixtures/laurent_large_acceptance_cases.jl` and assert diagnostics report success with no failure diagnostics and no partition search.

## Out Of Scope

Do not make the Issue #38 matrix factorize. Do not search general block decompositions beyond the `6 x 6` `3+3` diagnostic check. Do not change the local `SL_3` solver family.

## Self Review

- No placeholder requirements remain.
- The design keeps the current staged error text while adding a separate structured API.
- The partition search is bounded to `6 x 6` matrices and cannot affect the `40 x 40` negative control.
- The tests cover the Issue #38 row and column cores and the supported Laurent acceptance fixture.
