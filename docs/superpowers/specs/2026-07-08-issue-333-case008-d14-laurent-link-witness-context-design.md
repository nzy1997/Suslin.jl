# Issue 333 Case008 D14 Laurent Link-Witness Context Design

## Goal

Define one expert-only Laurent link-witness context for the replay-derived
post-descent `case_008 d=14` column. The context fixes the source markers,
measure state, pivot metadata, partner-index set, and witness schema that later
bounded search and certificate issues must use.

## Context

No repository `AGENTS.md`, `CLAUDE.md`, or `GEMINI.md` file is present in this
worktree. Issue #333 is open with no comments. Relevant dependency context:

- #327 promoted Laurent descent measure and replay helpers into internal code.
- #328 exposed the certified single-step d14 Laurent descent diagnostic and set
  `next_boundary = :laurent_link_witness`.
- #329 added
  `case008_d14_laurent_post_descent_profile_report()` and
  `validate_case008_d14_laurent_post_descent_profile_report(report)`.
- PR #332 merged #329 and verified the focused post-descent report test and full
  package test suite.

The local baseline command
`julia --project=. -e 'include("test/expert/case008_d14_laurent_post_descent_profile.jl")'`
passes after dependency instantiation. The new context should include that
source report file rather than copy the replayed post-descent column or profile.

## Approach Options

Recommended: add one focused expert test file,
`test/expert/case008_d14_laurent_link_witness_context.jl`. It includes the #329
post-descent profile report, validates the report with
`validate_case008_d14_laurent_post_descent_profile_report`, then constructs a
NamedTuple context from the validated report. This keeps the helper available to
later expert-only link-witness tests without exposing production API or implying
full Laurent ECP support.

Alternative: promote a context constructor into `src/algorithm`. This would
make future code easier to call, but it would add public or semi-public surface
before any link witness search/certificate exists and would blur this issue's
expert-only boundary.

Alternative: hard-code the expected context constants in a standalone fixture.
This would be simple, but it would violate the requirement to derive
`source_measure` from `report.after_measure` and pivot metadata from the
post-descent leading-monomial summary.

Chosen approach: expert-only NamedTuple context derived from a validated #329
report.

## Context Shape

`case008_d14_laurent_link_witness_context(report =
case008_d14_laurent_post_descent_profile_report())` returns a NamedTuple with
these stable fields:

- `case_id = "case_008"`;
- `dimension = 14`;
- `ring_generators = ("u", "v")`;
- `source_report_boundary = :case008_d14_original`;
- `source_boundary = :case008_d14_post_descent`;
- `boundary = :laurent_link_witness`;
- `source_report_status = :post_descent_profile_report`;
- `source_measure = report.after_measure`;
- `pivot_entry_index`;
- `pivot_leading_exponent`;
- `pivot_term_count`;
- `candidate_partner_indices`;
- `required_witness_fields`;
- `status = :link_witness_context`.

The helper first calls
`validate_case008_d14_laurent_post_descent_profile_report(report)` and rejects
any result other than `:ok`. It then derives the pivot fields from
`first(report.post_descent_leading_monomial_summary.candidates)`. Partner
indices are derived from `1:report.dimension` after excluding the pivot entry,
so they remain tied to the validated report's dimension and pivot while using
the stable ascending order expected by later search issues.

## Locked Constants

The committed test asserts:

- `source_report_boundary = :case008_d14_original`;
- `source_boundary = :case008_d14_post_descent`;
- `boundary = :laurent_link_witness`;
- `source_report_status = :post_descent_profile_report`;
- `source_measure.whole_support_count = 7378`;
- `source_measure.max_entry_terms = 3734`;
- `source_measure.valuation_span = (97, 92)`;
- `pivot_entry_index = 10`;
- `pivot_leading_exponent = (49, -5)`;
- `pivot_term_count = 3692`;
- `candidate_partner_indices = (1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 14)`;
- `required_witness_fields =
  (:family, :pivot_index, :partner_index, :coefficient, :exponent,
  :ring_generators)`;
- `status = :link_witness_context`.

## Validation

Add `validate_case008_d14_laurent_link_witness_context(context, report =
case008_d14_laurent_post_descent_profile_report())::Symbol` for later tests and
negative controls. It validates the source report, checks all required context
fields, confirms the required witness-field schema is complete and exact, and
compares the supplied context with recomputation from the same validated report.

Negative controls reject:

- a stale post-descent report;
- a report whose source validator returns a non-`:ok` result;
- `measure_relation != :strict_decrease`;
- swapped ring generators;
- a tampered first leading-monomial candidate pivot entry;
- a context whose required witness schema omits any required field.

## Tests

Create `test/expert/case008_d14_laurent_link_witness_context.jl` and register it
in the expert list in `test/runtests.jl` immediately after the #329
post-descent profile report.

Focused verification:

```bash
julia --project=. -e 'include("test/expert/case008_d14_laurent_link_witness_context.jl")'
```

Package verification:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Out Of Scope

Do not search for a Laurent link witness. Do not implement endpoint reductions,
normality/conjugation replay, determinant normalization, recursive peel
integration, full `case_008` support, public API, or production reducer changes.

## Automatic Decisions

- Visual companion skipped because this is a deterministic Julia test-helper
  change, not a visual design question.
- Clarifying questions skipped because Agent Desk is non-interactive and the
  issue body fixes the fields, constants, and verification command.
- Recommended approach selected: expert-only context helper derived from the
  validated #329 report.
- Design approval auto-approved under the Standing Answer Policy.
