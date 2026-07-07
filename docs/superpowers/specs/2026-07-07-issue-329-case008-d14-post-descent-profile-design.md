# Issue 329 Case008 D14 Post-Descent Profile Design

## Goal

Record a deterministic post-descent `case_008 d=14` Laurent profile/report by
replaying the certified entry-addition operation from the original d14 boundary
fixture. The report becomes stable test evidence for later Laurent link-witness
and endpoint-reduction work, but it must not claim reducibility or add public
reducer support.

## Context

No repository `AGENTS.md`, `CLAUDE.md`, or `GEMINI.md` file is present in this
worktree. Issue #329 is open with no comments. Relevant merged PR context:

- #330 promoted internal Laurent descent measure, replay, strict-decrease, and
  certificate-validation helpers.
- #331 exposed the certified single-step d14 Laurent descent diagnostic.

The worker branch starts at the #331 merge commit. The checked-in d14 fixture is
`test/fixtures/toricbuilder_case008_d14_column_boundary.jl`, and the existing
expert tests already define deterministic support, valuation, leading-monomial,
measure, and replay helpers for the same fixture.

## Approach Options

Recommended: add a focused expert-only test file
`test/expert/case008_d14_laurent_post_descent_profile.jl`. The file includes the
existing d14 descent-step certificate helpers, replays the certified operation
`target_index = 1`, `source_index = 2`, `coefficient = 1`, and
`exponent = (-1, 1)`, computes the post-descent profile from the replayed
column, and validates a supplied report by recomputing the before measure, after
measure, and all post-descent summaries.

Alternative: promote a reusable post-descent report helper into `src/`. This
would make later consumers easier to call, but the current issue only asks for a
fixture/report and doing this in production code risks implying a supported
Laurent-native reduction path.

Alternative: store the replayed post-descent column as after-data. This would be
straightforward, but it violates the requirement that the post-descent report be
derived from replay rather than trusted stored after-data.

Chosen approach: expert-only report derived from replay. It reuses the internal
Laurent descent helpers and existing expert summary helpers, keeps the after
column computed rather than stored, and leaves Laurent link witnesses and
endpoint reductions out of scope.

## Report Shape

`case008_d14_laurent_post_descent_profile_report(fixture)` returns a named tuple
with stable top-level fields:

- `case_id = "case_008"`;
- `dimension = 14`;
- `source_boundary = :case008_d14_original`;
- `ring_generators = ("u", "v")`;
- `operation_family = :entry_addition`;
- `operation`, carrying the certified operation and ring-generator metadata;
- `replay_status = :ok`;
- `before_measure`;
- `after_measure`;
- `measure_relation = :strict_decrease`;
- `post_descent_profile`;
- `post_descent_support_summary`;
- `post_descent_valuation_summary`;
- `post_descent_leading_monomial_summary`;
- `status = :post_descent_profile_report`.

The `post_descent_profile` and summary fields are computed from the replayed
column by `laurent_descent_step_profile(after_column, R; case_id)`. They are not
loaded from a stored after-column fixture.

## Locked Constants

The committed test asserts the exact post-descent measure constants produced by
replay:

- `whole_support_count = 7378`;
- `max_entry_terms = 3734`;
- `valuation_span = (97, 92)`;
- `leading_exponent = (49, -5)`;
- `leading_entry_index = 10`.

It also asserts exact post-descent entry term counts, whole-support bounds,
valuation ranges, full per-entry support-summary records, and the full
leading-monomial candidate summary. The validator then recomputes the complete
report from the original fixture and certified operation to reject stale or
hand-edited report data.

## Validation

`validate_case008_d14_laurent_post_descent_profile_report(report, fixture)`
returns `:ok` only when:

- the original d14 boundary fixture validates;
- report metadata matches `case_008`, dimension 14, source boundary, ring
  generators, status, replay status, operation family, and strict-decrease
  relation;
- the operation validates over the supplied ring;
- the replayed after column yields the recorded after measure;
- `strictly_decreases_laurent_measure(before_measure, after_measure)` is true;
- post-descent profile, support, valuation, and leading-monomial summaries all
  match recomputation from the replayed column.

Negative controls mutate the operation, operation ring generators, report ring
generators, after measure, full post-descent profile, support summary, valuation
summary, and leading-monomial summary. Each mutation must be rejected before the
report can be used as later Laurent link-witness input.

## Tests

Create `test/expert/case008_d14_laurent_post_descent_profile.jl` and register it
in the expert list in `test/runtests.jl` near the existing d14 descent-step
tests.

Focused verification:

```bash
julia --project=. -e 'include("test/expert/case008_d14_laurent_post_descent_profile.jl")'
```

Package verification:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Out Of Scope

Do not implement Laurent link witnesses, endpoint reductions,
normality/conjugation replay, determinant normalization, recursive peel
integration, or full `case_008` success. Do not require a local ToricBuilder
checkout. Do not add public API or exports.

## Automatic Decisions

- Visual companion skipped because this is a deterministic test-data/report
  task, not a visual design question.
- Clarifying questions skipped because Agent Desk is non-interactive and the
  issue body gives exact operation, fields, and verification commands.
- Recommended approach selected: expert-only replay-derived report, because it
  satisfies the issue without production API churn or unsupported reducer
  claims.
- Design approval auto-approved under the Standing Answer Policy.
