# Issue 300 Case008 Post-D15 Bounded Report Design

## Goal

Make the explicit bounded `case_008` Q-block report acceptance check prove that
the old `current d=15` Laurent-column reducer boundary is no longer the active
outcome after issue #299.

## Context

Issue #137 added the bounded report path and a generic structured-row validator
for explicit `case_008` runs. Issue #299, merged by PR #309, added certified
support for the `case_008 d=15` Laurent row-preconditioning stage. The remaining
integration gate is the full bounded `case_008` report: it may pass, find a
later algorithm boundary below dimension 15, or time out after making progress,
but it must not keep accepting the old unsupported d15 reducer error as current
evidence.

No `AGENTS.md`, `CLAUDE.md`, `GEMINI.md`, or `CONVENTIONS.md` file is present
in this worktree. `gh issue view` for #300 and #299 is blocked by the sandbox
proxy, so this design uses the full issue #300 body supplied by Agent Desk plus
locally available PR #309 metadata.

## Approaches Considered

Recommended: keep the existing bounded report runner and add a dedicated
post-d15 acceptance predicate in `test/internal/toricbuilder_cache_status_report.jl`.
The predicate builds on the existing structured-row validator, accepts
`gl_certificate_pass`, accepts `certified_algorithm_boundary` only when the
details do not report the old d15 unsupported Laurent-column family and do not
claim `current d=15`, and accepts `timed_out` only when the timeout metadata
shows the peel moved below d15 or explicitly completed d15. This keeps the
change local to the status-report contract and negative controls.

Alternative: change the report generator to rewrite old d15 errors as a newer
boundary. This would hide real route evidence instead of validating it.

Alternative: add a separate `case_008` runner. This duplicates the bounded
worker infrastructure and conflicts with the issue's request to reuse the
existing report path.

Chosen approach: status validator hardening plus report evidence checks. The
default report remains unchanged; `case_008` stays explicit through
`--exercise=case_008`.

## Acceptance Rules

For explicit bounded `case_008`, first require the existing structured-row
rules:

- `:gl_certificate_pass` requires `verified == true` and `error_details ==
  "none"`.
- `:certified_algorithm_boundary` requires exactly one boundary stage with
  numeric elapsed time, stable error details, and evidence naming the stage.
- `:timed_out` requires exactly one timed-out stage, numeric elapsed time, the
  configured timeout text in both row and stage details, and the stage name in
  row details.

Then apply the post-d15 gate:

- pass rows are accepted;
- algorithm-boundary rows are rejected if route or stage details contain
  `unsupported exact unimodular column reduction for Laurent-normalized column
  of length 15`, or if the details pair `current d=15` with
  `unsupported_laurent_column_family`;
- algorithm-boundary rows are accepted only if details expose a current peel
  dimension below 15, a Laurent-normalized column length below 15, or completed
  d15 progress metadata;
- timeout rows at certificate construction are accepted only if progress text
  reports `current d=<n>` with `n < 15`, or reports `last completed d=15`;
- timeout rows without certificate-construction peel progress remain rejected
  for the post-d15 `case_008` gate, even though they can still be structurally
  valid bounded rows for other issues;
- timeout rows at stages before d15 reduction are rejected.

The old strings are intentionally matched in tests so that a future report
cannot accidentally treat the historical d15 unsupported boundary as accepted
post-#299 evidence.

## Tests

Add focused negative controls to
`test/internal/toricbuilder_cache_status_report.jl`:

- an algorithm-boundary row whose details contain the exact old d15 unsupported
  Laurent-column text must fail the post-d15 predicate;
- an algorithm-boundary row whose details say `current d=15` with
  `unsupported_laurent_column_family` must fail;
- an algorithm-boundary row whose details say `current d=14` with
  `unsupported_laurent_column_family` must pass;
- a timeout row whose progress says `current d=15` but has no completed-d15
  metadata must fail;
- a timeout row whose progress says `current d=14` must pass;
- a timeout row whose progress says `current d=15` and `last completed d=15`
  must pass;
- a generated explicit `case_008` row must satisfy the post-d15 predicate under
  the issue verification timeout.

Required verification:

```bash
julia --project=. scripts/report_toricbuilder_cache_q_blocks.jl --exercise=case_008 --timeout-seconds=180 --output=/tmp/qblock-case008-after-d15.md
julia --project=. -e 'include("test/internal/toricbuilder_cache_status_report.jl")'
julia --project=. test/runtests.jl all
julia --project=. -e 'using Pkg; Pkg.test()'
git diff --check
```

Inspect `/tmp/qblock-case008-after-d15.md` and confirm the `case_008` row is
structured and the route details do not contain the old d15 unsupported reducer
strings.

## Out Of Scope

Do not require `case_007`, `case_009`, `case_011`, or `case_012` to pass. Do
not add Steinberg factor-count optimization. Do not add `case_008` to
`DEFAULT_EXERCISED_CASE_IDS`.

## Automatic Decisions

- Visual companion skipped because the task is a CLI/report validator change
  with no visual design decision.
- Clarifying questions skipped because Agent Desk is non-interactive and the
  issue body gives a precise acceptance contract.
- Recommended approach selected: harden the existing bounded report validator,
  because it preserves the runner contract and directly encodes the post-d15
  evidence requirement.
- Design approval auto-approved under the Standing Answer Policy.
