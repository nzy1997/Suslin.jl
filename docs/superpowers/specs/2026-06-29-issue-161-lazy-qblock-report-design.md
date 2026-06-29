# Issue 161 Lazy Q-Block Report Design

## Context

The ToricBuilder cache Q-block status report currently exercises a conservative
default set through the eager Laurent GL certificate route. Large cases can spend
their bounded budget before useful peel evidence appears. Issue #160 added the
public lazy determinant certificate route:

```julia
laurent_gl_factorization_certificate(
    A;
    determinant_strategy = :lazy,
    correction_side = :row,
)
```

Issue #161 asks the report to expose that route when a caller explicitly
requests it, without claiming success unless exact certificate verification
succeeds.

## Approach Options

Recommended: add explicit report options
`--determinant-strategy=eager|lazy` and
`--correction-side=row|column`. Keep the default eager behavior and the existing
default exercised case list unchanged. Carry determinant strategy metadata in
every row, and for lazy exercised rows also carry correction side and determinant
source. This is compatible with existing generated reports while giving bounded
lazy runs structured evidence.

Alternative: switch the default report to lazy. This would make issue #161
visible without flags, but it is a behavior change for the audit file and could
move runtime boundaries for the default report.

Alternative: add only textual lazy evidence to the existing evidence string.
This is smaller, but it does not give tests or readers a structured way to tell
whether a timeout happened before or after deferred determinant work.

## Chosen Design

The CLI accepts:

```text
--determinant-strategy=eager
--determinant-strategy=lazy
--correction-side=row
--correction-side=column
```

`determinant_strategy = :eager` remains the default. `correction_side` is valid
only with `determinant_strategy = :lazy`; the parser rejects unsupported values
and eager correction-side requests with `ArgumentError`.

Report rows gain these stable fields:

- `determinant_strategy`: `:eager`, `:lazy`, or `:not_run`.
- `correction_side`: `:row`, `:column`, or `:not_run`.
- `determinant_source`: `:deferred_submatrix`, `:not_reached`, or `:not_run`.

Lazy exercised rows always include `determinant_strategy = :lazy` and the
resolved `correction_side`, including timeout and boundary rows. Lazy rows report
`determinant_source = :deferred_submatrix` only after the lazy certificate or
worker progress proves the deferred determinant route was reached. If a bounded
run stops before that, the row reports `:not_reached`; non-lazy and non-exercised
rows report `:not_run`.

The bounded worker progress file records the determinant strategy, correction
side, and determinant source alongside stage timing and peel progress. The eager
worker path keeps its current staged timing behavior. The lazy worker path calls
`laurent_gl_factorization_certificate(A; determinant_strategy = :lazy,
correction_side, progress_callback)` directly so it avoids eager determinant
classification and eager normalization stages. Those eager-only stages are
recorded as `:not_run`, while certificate construction and verification retain
their existing timing detail.

Markdown adds a compact lazy metadata table and enriches exercised evidence so a
reviewer can see determinant strategy, correction side, determinant source, and
stage timing without parsing error text.

## Validation Rules

Internal status-report tests enforce truthful pass claims:

- A row with `route_status = :gl_certificate_pass` must have `verified = true`.
- A lazy pass must include `determinant_strategy = :lazy`,
  `correction_side in (:row, :column)`, and
  `determinant_source = :deferred_submatrix`.
- A lazy pass with `determinant_source = :not_reached` is invalid.
- A timeout must never be accepted as a pass.
- Lazy timeout rows must still include determinant strategy, correction side,
  determinant source, and stage timing.

## Testing

Extend `test/internal/toricbuilder_cache_status_report.jl` with synthetic worker
rows for lazy success, lazy timeout before deferred determinant work, and lazy
timeout after deferred determinant work. Add negative controls for false pass
claims with `verified = false`, missing lazy metadata, and
`determinant_source = :not_reached`.

Add one small real lazy fixture check using `case_010` through the worker helper
with `determinant_strategy = :lazy` and `correction_side = :row`, because it is
already part of the conservative default exercised set and should not expand
`DEFAULT_EXERCISED_CASE_IDS`.

Required verification commands:

```bash
julia --project=. -e 'include("test/internal/toricbuilder_cache_status_report.jl")'
julia --project=. scripts/report_toricbuilder_cache_q_blocks.jl --exercise=case_009 --timeout-seconds=180 --determinant-strategy=lazy --output=/tmp/qblock-case009-lazy.md
julia --project=. -e 'using Pkg; Pkg.test()'
```

The case_009 report may pass, hit a staged boundary, or time out, but it must
include a structured `case_009` row, lazy determinant metadata, and stage timing
details.

## Scope

Do not add slow cases to `DEFAULT_EXERCISED_CASE_IDS`. Do not update
`docs/audits/2026-06-24-toricbuilder-cache-q-block-status.md` because the
default report behavior intentionally remains eager and conservative. Do not
commit `Manifest.toml`.

## Automatic Decisions

- Clarifying questions, design approval, and written spec review were resolved
  by the Standing Answer Policy because this is a non-interactive Agent Desk
  run.
- The visual companion was skipped because no visual question would clarify this
  CLI/reporting change.
- The explicit lazy strategy with eager default was selected because it matches
  issue #161 and avoids broad report behavior changes.
- `case_010` was selected as the small real lazy fixture because it is already
  exercised by the conservative default status-report tests.
