# Task 1-2 Report: ECP Public Unimodular Pipeline

## What Changed

- Added `test/expert/elementary_column_property.jl` as the focused acceptance test for the public ECP unimodular-column path and staged expert API.
- Registered the new expert test in `test/runtests.jl`.
- Updated `src/algorithm/column_reduction.jl` to:
  - add `ECPStagedColumnReductionCertificate`;
  - keep `reduce_unimodular_column(v, R)` factor-returning while preferring a verified staged ECP route when a deterministic public fixture witness is recognized;
  - add `ecp_staged_column_reduction_certificate(v, R; ...)`;
  - add `verify_ecp_staged_column_reduction(certificate)`;
  - add `_ecp_staged_column_reduction_replay_summary(certificate)`;
  - add `_ecp_public_staged_reduction_certificate(column, R)`;
  - add `_ecp_default_public_link_witness(column, R, selected_variable)` for the GF(2) fixture route;
  - add `_ecp_default_public_normality_witness(lower_factors, R)`;
  - extract `_ecp_column_reduction_certificate_validated(column, R)` so the public reducer validates once before staged-or-legacy dispatch.

## RED Output Summary

Command:

```bash
julia --project=. -e 'include("test/expert/elementary_column_property.jl")'
```

Observed RED:

- The test failed before any production edit.
- Failure reason: `UndefVarError: ecp_staged_column_reduction_certificate not defined in Suslin`.
- The failure occurred in the new acceptance test at the first staged API call.
- Summary reported by Julia: `12 passed, 0 failed, 1 errored`.

## GREEN Output Summary

Command:

```bash
julia --project=. -e 'include("test/expert/elementary_column_property.jl")'
```

Observed GREEN:

- Exit code: `0`
- Test summary: `27 passed, 27 total`
- The public reducer returned staged factors for the canonical GF(2) fixture and legacy factors remained available through `ecp_column_reduction_certificate`.

## Regression Command Results

1. `julia --project=. -e 'include("test/expert/unimodular_reduction_exact.jl")'`
   - Exit code: `0`
   - Summaries:
     - `exact unimodular reduction supports longer ordinary columns`: `8 passed`
     - `exact unimodular reduction supports Laurent-normalized columns`: `9 passed`
     - `exact unimodular reduction preserves old small cases`: `6 passed`
     - `exact unimodular reduction staged failures`: `11 passed`

2. `julia --project=. -e 'include("test/expert/ecp_variable_change_replay.jl")'`
   - Exit code: `0`
   - Summary: `ECP variable-change replay records`: `93 passed`

3. `julia --project=. -e 'include("test/expert/ecp_induction_normality.jl")'`
   - Exit code: `0`
   - Summary: `ECP induction and normality replay`: `27 passed`

## Files Changed

- `src/algorithm/column_reduction.jl`
- `test/runtests.jl`
- `test/expert/elementary_column_property.jl`
- `docs/superpowers/plans/2026-06-23-ecp-public-unimodular-pipeline-task-1-2-report.md`

## Self-Review

- The staged public route is intentionally narrow and fixture-gated, matching the plan instead of broadening support heuristically.
- The fallback path remains the existing certificate reducer, so non-fixture ordinary polynomial columns and Laurent consumers keep prior behavior.
- Replay verification for the staged certificate is stored and rechecked in the same style as the existing link-step and induction/normality certificates.
- The focused test exercises both success and tampered-witness failure cases, which is the right pressure for this API surface.

## Concerns

- `_ecp_default_public_link_witness` currently recognizes only the hard-coded GF(2) fixture route from the plan. That is intentional for this task, but expanding public staged coverage later will need a less ad hoc witness recognizer.

## Follow-up Fix: Staged Input Validation

### What Changed

- Added two focused regressions in `test/expert/elementary_column_property.jl`:
  - `lower_reduction = Any[]` must throw `ArgumentError` instead of surfacing a raw `BoundsError`.
  - `variable_order = ()` must throw a clearer staged `ArgumentError` mentioning `variable_order` or the need for at least one generator.
- Updated `ecp_staged_column_reduction_certificate` in `src/algorithm/column_reduction.jl` to:
  - reject empty `variable_order` before calling `first(normalized_order)`;
  - validate injected `lower_reduction` through `_ecp_verified_lower_reduction` before synthesizing the default normality witness;
  - pass the verified lower-column length into `_ecp_default_public_normality_witness` so the default witness no longer depends on `first(lower_factors)`.

### RED Output

Command:

```bash
julia --project=. -e 'include("test/expert/elementary_column_property.jl")'
```

Observed RED:

- Exit code: `1`
- Failure 1: `err isa ArgumentError` failed because the actual exception was `BoundsError(Any[], (1,))`.
- Failure 2: the `variable_order = ()` assertion failed because the error text did not mention `variable_order` or `at least one`.
- Summary: `28 passed, 2 failed, 0 errored`.

### GREEN Output

Command:

```bash
julia --project=. -e 'include("test/expert/elementary_column_property.jl")'
```

Observed GREEN:

- Exit code: `0`
- Summary: `public ECP unimodular-column pipeline | 30 passed, 30 total`

### Regression Output

1. `julia --project=. -e 'include("test/expert/ecp_induction_normality.jl")'`
   - Exit code: `0`
   - Summary: `ECP induction and normality replay | 27 passed, 27 total`

### Files Changed

- `src/algorithm/column_reduction.jl`
- `test/expert/elementary_column_property.jl`
- `docs/superpowers/plans/2026-06-23-ecp-public-unimodular-pipeline-task-1-2-report.md`

### Self-Review

- The constructor now uses the existing lower-reduction verifier instead of duplicating ad hoc factor extraction logic.
- The default normality witness is derived from verified lower-column shape, which removes the raw indexing hazard even if the verified factor sequence is empty.
- The change is scoped to staged expert input handling and does not alter the verified reduction route or factor order.

### Concerns

- I only reran the focused expert file and the requested induction/normality regression for this follow-up fix.

## Final Review Fix: Staged Replay Field Consistency

### What Changed

- Added focused regressions in `test/expert/elementary_column_property.jl` proving that top-level staged-certificate tampering now fails verification when any of these fields are changed independently of the nested verified certificates:
  - `monicity`
  - `lower_reduction`
  - `normality_witness`
- Strengthened the existing `lower_reduction = Any[]` regression to require an `ArgumentError` message mentioning `lower-variable reduction` or `v(0)`.
- Added a focused internal-helper regression showing that `_ecp_public_staged_reduction_certificate` no longer falls back to the legacy reducer after the deterministic public route has already been recognized; an injected staged-construction failure now propagates as an `ArgumentError`.
- Updated `src/algorithm/column_reduction.jl` so `_ecp_staged_column_reduction_replay_summary` now checks top-level staged fields against the nested verified route:
  - `original_column` must agree across the staged certificate, link-step certificate, and induction/normality certificate
  - `monicity` must equal the replayed monicity data derived from `link_step.link_witness`
  - `lower_reduction` must match the nested lower-reduction certificate or lower-variable factor sequence
  - `normality_witness` must match the nested induction/normality witness
- Updated `_ecp_public_staged_reduction_certificate` so it accepts optional internal dependency injection for focused testing, catches `ArgumentError` only while recognizing an unsupported deterministic route through `_ecp_default_public_link_witness`, and otherwise lets recognized-route staged-construction failures surface.

### RED Output

Command:

```bash
julia --project=. -e 'include("test/expert/elementary_column_property.jl")'
```

Observed RED:

- Exit code: `1`
- Three new tamper checks failed because `verify_ecp_staged_column_reduction` still returned `true` after top-level `monicity`, `lower_reduction`, and `normality_witness` were modified.
- The recognized-route propagation regression failed with a `MethodError` because `_ecp_public_staged_reduction_certificate` did not yet accept injected staged-construction inputs.
- Summary: `31 passed, 5 failed, 0 errored`

### GREEN Output

Command:

```bash
julia --project=. -e 'include("test/expert/elementary_column_property.jl")'
```

Observed GREEN:

- Exit code: `0`
- Summary: `public ECP unimodular-column pipeline | 36 passed, 36 total`

### Regression Output

1. `julia --project=. -e 'include("test/expert/ecp_induction_normality.jl")'`
   - Exit code: `0`
   - Summary: `ECP induction and normality replay | 27 passed, 27 total`

### Files Changed

- `src/algorithm/column_reduction.jl`
- `test/expert/elementary_column_property.jl`
- `docs/superpowers/plans/2026-06-23-ecp-public-unimodular-pipeline-task-1-2-report.md`

### Self-Review

- The staged replay verifier now treats the top-level staged fields as part of the trusted surface instead of assuming the nested certificates are the only state worth checking.
- The public helper change is intentionally narrow: unsupported-route recognition still falls back quietly, but once the deterministic route is recognized, downstream staged-construction failures are surfaced instead of being silently hidden by the legacy reducer.
- The only new internal surface is keyword-based dependency injection on the non-exported helper, which keeps the public API unchanged while making the supported-route failure mode testable.

### Concerns

- I reran the focused staged-pipeline expert file and the requested induction/normality regression only; I did not rerun the broader expert suite in this follow-up.
