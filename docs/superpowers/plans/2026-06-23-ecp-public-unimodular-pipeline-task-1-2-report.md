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
