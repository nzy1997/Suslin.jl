# Issue 147 Bounded Peel Evidence Design

## Goal

Record bounded, reproducible per-peel progress evidence for ToricBuilder
`case_007` and `case_009` without making either case part of the default report
or requiring either case to pass the Laurent `GL_n` certificate route.

The generated evidence must identify the furthest stage reached. When a bounded
run times out during `certificate_construction`, the row must also report the
active Laurent column-peel dimension, completed peel count, last completed peel
metadata, and last-column sparsity or term statistics.

## Context

Issue #131 tracks the ToricBuilder cache Q-block status. The current default
report exercises `case_001` through `case_006` and `case_010`; `case_007` and
`case_009` remain outside `DEFAULT_EXERCISED_CASE_IDS` because they are slow
optional evidence cases.

The local worktree already contains a draft implementation that writes
`peel.*` progress fields from the bounded worker and renders them in timeout
details. That draft proves the direction is useful, but it duplicates a large
part of the Laurent column-peel and Laurent `GL_n` certificate algorithm inside
`scripts/report_toricbuilder_cache_q_blocks.jl`. This design keeps the useful
progress schema while moving the instrumentation hook into the core algorithm
path so the report evidence follows the same route as normal certificate
construction.

## Approach Options

Recommended: add an observation-only progress callback to the private Laurent
column-peel path in `src/algorithm/laurent_column_peel.jl`, thread it through a
private Laurent `GL_n` certificate helper, and have the bounded report worker
use that helper to write `peel.*` progress. This keeps the report script thin
and prevents algorithm drift between `src/` and `scripts/`.

Alternative: keep the current script-local copy of the peel and certificate
logic. This is the shortest path from the draft diff to visible evidence, but
it creates a shadow implementation that can silently diverge from the real
certificate route.

Alternative: record only stage-level timings for `case_007` and `case_009`.
This is simple and avoids algorithm instrumentation, but it does not satisfy
issue #147 because `certificate_construction` timeouts would still not say
which peel dimension was active.

## Chosen Design

Use a core progress callback with no public behavior change.

`_factor_laurent_sl_column_peel(A; progress_callback = nothing)` should keep
returning the same `LaurentColumnPeelFactorization`. When the callback is
`nothing`, behavior and results are unchanged. When a callback is supplied, the
recursive peel path reports observation-only progress before starting a peel at
dimension `d` and after completing a peel step.

The callback record should include:

- `current_dimension`
- `completed_steps`
- `last_completed_dimension`
- `last_completed_elapsed_seconds`
- `last_completed_left_factors`
- `last_completed_right_factors`
- `last_column_nnz`
- `max_entry_terms`

The callback must not influence algorithm choices. If the callback throws, the
bounded worker may surface a `route_error`; losing the evidence channel should
not be disguised as successful algorithm evidence.

Add or extend a private Laurent `GL_n` certificate helper so the bounded worker
can pass the callback through certificate construction without changing the
public `laurent_gl_factorization_certificate(A)` API. The public function can
call the helper with `progress_callback = nothing`.

In `scripts/report_toricbuilder_cache_q_blocks.jl`, keep the bounded subprocess
model and the existing stage names. The worker callback should write the latest
progress to the worker progress file as `peel.*` keys. The parent process should
read these keys on timeout and append a concise `peel progress:` sentence to
Route Error Details and row evidence.

Do not add `case_007` or `case_009` to `DEFAULT_EXERCISED_CASE_IDS`.

## Components

- `src/algorithm/laurent_column_peel.jl`: owns the optional progress callback
  and all peel-step observation points.
- `src/algorithm/laurent_gl_certificate.jl`: owns the private certificate
  helper that threads the callback into the normalized determinant-one core
  factorization.
- `scripts/report_toricbuilder_cache_q_blocks.jl`: owns bounded worker
  orchestration, progress-file serialization/parsing, timeout row construction,
  markdown rendering, and command-line output paths.
- `test/internal/toricbuilder_cache_status_report.jl`: owns report-level parser,
  renderer, bounded worker, and validator coverage.
- `docs/audits/2026-06-27-qblock-case007-case009-peel-evidence.md`: records
  the actual `case_007` and `case_009` generated report paths and observed
  timings. A follow-up #131 comment may link to this artifact after review, but
  the local artifact is the required deliverable.

## Data Flow

1. The user runs the report script with explicit `--exercise=case_007` or
   `--exercise=case_009` plus `--timeout-seconds`.
2. The parent process starts a bounded worker and gives it a progress path.
3. The worker materializes the selected fixture and records
   determinant-classification and normalization timings as it does today.
4. During certificate construction, the worker calls the private certificate
   helper with a progress callback.
5. The core peel path emits observation records before and after peel steps.
6. The report callback writes the latest observation to the progress file as
   `peel.*` keys alongside the current stage timing data.
7. If the worker times out, the parent reads the progress file and turns the
   latest stage and peel fields into a structured timeout row.
8. The generated markdown includes the row in the case table, stage timing
   details, and Route Error Details.

## Error Handling And Evidence Rules

For determinant-classification or normalization timeouts, stage timing evidence
is sufficient because the peel path has not started.

For `certificate_construction` timeouts, issue #147 evidence is sufficient only
when Route Error Details include:

- `peel progress:`
- `current d=`
- `completed steps=`
- at least one last-column statistic, such as `last-column nnz=` or
  `max entry terms=`

A `certificate_construction` timeout without any `peel.*` keys can remain a
structured bounded timeout, but the issue #147 validator must reject it as
insufficient per-peel evidence.

For `case_007`, the expected useful outcome is a bounded
`certificate_construction` timeout with peel progress. For `case_009`, the
report should record the furthest reached stage; if it reaches
`certificate_construction`, the same peel-progress requirement applies.

## Tests

Add focused core coverage that exercises a small Laurent determinant-one case
with a callback and checks that progress records contain reasonable dimensions
and completed-step counts. Also check that running without a callback preserves
the existing factorization result and verification behavior.

Keep report-level synthetic coverage for `peel.*` parsing and rendering. A fake
bounded worker should write a progress file with
`current_stage=certificate_construction` and the peel fields, then time out.
The test should assert that Route Error Details and markdown include
`peel progress:`, `current d=`, `completed steps=`, and last-column statistics.

Add the issue #147 negative control: a fake bounded worker that times out during
`certificate_construction` but omits all `peel.*` keys. The validator must
reject that row as insufficient per-peel evidence even if the row is otherwise
a structured timeout.

Preserve existing default report tests that prove `case_007` and `case_009`
remain `not_exercised_in_default_report` in default runs.

## Verification

Required issue commands:

```bash
julia --project=. scripts/report_toricbuilder_cache_q_blocks.jl --exercise=case_007 --timeout-seconds=180 --output=/tmp/qblock-case007-peel-progress.md
julia --project=. scripts/report_toricbuilder_cache_q_blocks.jl --exercise=case_009 --timeout-seconds=650 --output=/tmp/qblock-case009-progress.md
julia --project=. -e 'include("test/internal/toricbuilder_cache_status_report.jl")'
```

Expected `case_007` result: the generated report contains a structured
`case_007` row, not `not_exercised_in_default_report`, and includes determinant,
normalization, and certificate-construction stage timing. If it times out during
`certificate_construction`, Route Error Details contain the required peel
progress fields.

Expected `case_009` result: the generated report contains a structured
`case_009` row, not `not_exercised_in_default_report`, and records the furthest
reached stage. If it reaches `certificate_construction`, Route Error Details
contain the required peel progress fields.

Run relevant core or expert tests for the callback-bearing Laurent column-peel
path. Run full package tests if the private helper threading touches shared
certificate behavior beyond the bounded report route.

## Out Of Scope

Do not require `case_007` or `case_009` to pass the Laurent `GL_n` certificate
route. Do not add either case to the default report. Do not optimize determinant
classification, normalization, or certificate-construction performance. Do not
broaden public Laurent `GL_n` API claims.

## Automatic Decisions

The visual companion was skipped because this is a CLI/report/test design with
no visual design decision. The design intentionally replaces the current
script-local algorithm copy with core instrumentation because the evidence
should come from the real Laurent certificate route.
