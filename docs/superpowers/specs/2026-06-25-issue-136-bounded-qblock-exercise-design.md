# Issue 136 Bounded Q-Block Exercise Design

## Goal

Add an optional bounded exercise mode to
`scripts/report_toricbuilder_cache_q_blocks.jl` so slow ToricBuilder Q-block
fixtures can be probed explicitly with stable timeout or classification rows and
stage-level timing metadata.

## Context

The current report script exercises the default contract set:

```text
case_001, case_002, case_003, case_004, case_005, case_006, case_010
```

The larger recorded fixtures `case_007`, `case_008`, `case_009`, `case_011`,
and `case_012` stay out of the default report. Issue #131 recorded that a
manual probe of `case_007` and `case_008` was stopped after more than two
minutes inside determinant or characteristic-polynomial work. The current
script also accepts unknown `--exercise` case IDs silently because it only
checks whether catalog rows match the requested set.

## Approach Options

Recommended: keep default behavior unbounded and unchanged, but route
`--timeout-seconds=N` exercised rows through a single-case Julia subprocess.
The parent process polls the child and kills it when the budget expires. The
child writes stage progress before determinant classification, normalization,
certificate construction, and verification, so the parent can produce a stable
`:timed_out` row with partial stage timings if the child is terminated.

Alternative: add elapsed-time checks around the existing in-process stages. This
would record timings for successful stages but would not interrupt determinant
or characteristic-polynomial calls while they are running, so it does not meet
the hard-budget requirement.

Alternative: require users to wrap the script with an external `timeout`
command. That avoids code in this repository but is not portable on macOS and
would not let the report produce a clean `case_007` row after timeout.

## Chosen Design

Add `--timeout-seconds=<positive number>` to the report script. When omitted,
`build_report()` and the default CLI path continue to use the current in-process
exercise route and the same default exercised case set. When present, every
explicitly exercised row is evaluated by a child Julia process running the same
script in worker mode for exactly one case.

Validate every requested `--exercise` ID against the checked-in catalog before
any report file is written. Unknown IDs such as `case_999` throw a clear
`ArgumentError` listing the unknown IDs and exit nonzero.

The worker path performs these measured stages in order:

1. determinant classification via `classify_laurent_determinant(A)`;
2. normalization via `normalize_laurent_gl_matrix(A)`;
3. certificate construction via `laurent_gl_factorization_certificate(A)`;
4. verification via `verify_laurent_gl_factorization_certificate(certificate)`.

Each report row gains `stage_timings`, a named tuple keyed by
`:determinant_classification`, `:normalization`, `:certificate_construction`,
and `:verification`. Each stage records `status`, `elapsed_seconds`, and
`error_details`. Stages not reached record `:not_run` and `not_run`.

## Status Semantics

The bounded worker returns `:gl_certificate_pass` when certificate verification
succeeds, `:gl_certificate_fail` when verification returns false,
`:certified_algorithm_boundary` when one of the measured stages throws a staged
or unsupported algorithm boundary, and `:route_error` for unexpected worker
errors. If the parent kills the worker because the budget expires, the row
records `route_status == :timed_out` and marks the current stage as
`:timed_out`.

For bounded worker rows, `public_elementary_status` remains `:not_run`; the
bounded mode measures the certificate route named in the issue rather than
spending budget on the separate public `elementary_factorization` path.
Unbounded default rows keep the existing public status behavior.

## Markdown Report

Keep the existing case table columns intact and add a `Stage Timing Details`
section. The section has one row per report row and one column for each measured
stage, rendered as `status (seconds)` when a numeric elapsed time exists or just
the symbolic status otherwise. This keeps the primary status table stable while
making stage metadata visible in generated reports.

## Tests

Add focused status-report tests for:

- unknown `--exercise` IDs failing with `ArgumentError`;
- parsing and storing `--timeout-seconds`;
- bounded `case_007` reporting `:timed_out` within a one-second budget and
  exposing stage timing metadata;
- rendered markdown including the `Stage Timing Details` section;
- default `build_report()` still leaving `case_007` and the other slow cases as
  `:not_exercised_in_default_report`.

The issue verification command exercises the slow case explicitly:

```bash
julia --project=. scripts/report_toricbuilder_cache_q_blocks.jl --exercise=case_007 --timeout-seconds=1 --output=/tmp/qblock-case007-timeout.md
```

The default internal report test remains:

```bash
julia --project=. -e 'include("test/internal/toricbuilder_cache_status_report.jl")'
```

The package verification command remains:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Out Of Scope

Do not make any slow case pass. Do not add a general Laurent determinant or
characteristic-polynomial algorithm. Do not put slow cases into routine default
tests except through a bounded one-second worker test.

## Automatic Decisions

- Clarifying questions were answered by the Standing Answer Policy because this
  is a non-interactive Agent Desk run.
- The subprocess timeout approach was selected because the issue asks for a hard
  budget that can stop determinant or characteristic-polynomial work.
- The visual companion was skipped because no visual question applies to this
  CLI report feature.
- The design was approved automatically under the Standing Answer Policy because
  it is conservative, preserves the default exercised set, and directly matches
  the issue verification commands.
