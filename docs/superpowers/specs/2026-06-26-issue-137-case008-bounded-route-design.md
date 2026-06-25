# Issue 137 Case 008 Bounded Route Design

## Goal

Classify ToricBuilder `case_008` through the bounded Q-block report path with a
stable route status, stage name, failure or timeout code, and timing data.

## Context

Issue #136 added the bounded worker mode for explicitly exercised slow Q-block
rows. The current branch starts after that work is merged. The default report
still leaves `case_008` as `not_exercised_in_default_report`, which should not
change.

Running the issue command before this change:

```bash
julia --project=. scripts/report_toricbuilder_cache_q_blocks.jl --exercise=case_008 --timeout-seconds=120 --output=/tmp/qblock-case008.md
```

produced a bounded row, but the row had `route_status == :route_error` because
the parent tried to deserialize the worker row from captured stdout and found a
`Float32` payload before the expected report row. The stage timing showed the
work reached determinant classification for about 100 seconds, so the report
had useful stage context but surfaced it as an infrastructure route error.

## Approach Options

Recommended: keep #136's bounded subprocess architecture and send the worker's
serialized report row through a dedicated result file instead of stdout. The
parent can still capture stdout and stderr as diagnostics, but package or
library stdout noise can no longer corrupt the row payload. Add focused
`case_008` coverage that accepts only `:gl_certificate_pass`,
`:certified_algorithm_boundary`, or structured `:timed_out`.

Alternative: keep stdout serialization and scan for the last valid serialized
object. This is fragile because Julia serialization is a stream format, and
arbitrary stdout bytes can make recovery ambiguous.

Alternative: classify any invalid worker stdout as
`:certified_algorithm_boundary`. This would satisfy the status shape but would
mislabel infrastructure corruption as an algorithm boundary.

## Chosen Design

Extend the bounded worker protocol with `--worker-result=<path>`. The parent
creates a temporary result path and passes it to the worker along with the
existing progress path. The worker serializes the final row to that result path.
For direct manual worker invocation, a missing result path keeps the existing
stdout serialization fallback.

The parent continues to capture stdout and stderr. If the worker exits
successfully, the parent deserializes the result file. If the result file is
missing or invalid, the parent still reports a route error because that is an
infrastructure failure. The `case_008` acceptance path should not hit this
fallback once the result-file protocol is used.

Add a status validator in `test/internal/toricbuilder_cache_status_report.jl`
for bounded route rows. For `case_008`, the validator accepts exactly:

- `route_status == :gl_certificate_pass` with `verified == true`;
- `route_status == :certified_algorithm_boundary` with a non-`not_run` stage
  status, stage elapsed time, and stable error details;
- `route_status == :timed_out` with a non-`not_run` timed-out stage, elapsed
  time, a stage name in `error_details`, and the configured timeout budget text.

The validator rejects raw `:route_error` rows and plain timeout text that lacks
a stage name or budget. This is the negative control required by the issue.

## Tests

Add focused tests for:

- the bounded worker reading a serialized row from `--worker-result` even when
  stdout contains unrelated bytes;
- `case_008` with a one-second budget producing a structured bounded row
  instead of `:route_error` or `:not_exercised_in_default_report`;
- the structured route validator rejecting a raw route error and an
  unstructured timeout row.

Keep routine tests bounded to one-second `case_008` and `case_007` probes. The
issue verification command remains the explicit 120-second `case_008` report
run.

Required verification:

```bash
julia --project=. scripts/report_toricbuilder_cache_q_blocks.jl --exercise=case_008 --timeout-seconds=120 --output=/tmp/qblock-case008.md
julia --project=. -e 'include("test/internal/toricbuilder_cache_status_report.jl")'
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Out Of Scope

Do not require `case_008` to pass the Laurent GL certificate route. Do not add a
new determinant or normalization algorithm. Do not add `case_008` to the default
unbounded exercised set.

## Automatic Decisions

- Clarifying questions were answered by the Standing Answer Policy because this
  is a non-interactive Agent Desk run.
- The visual companion was skipped because this CLI report feature has no
  visual design decision.
- The dedicated worker result-file approach was selected because it fixes the
  observed stdout payload corruption while preserving #136's subprocess budget
  model.
- The design was approved automatically under the Standing Answer Policy because
  it is conservative and directly targets the issue's verification command.
