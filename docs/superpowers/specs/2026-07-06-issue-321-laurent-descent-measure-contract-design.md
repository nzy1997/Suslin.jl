# Issue 321 Laurent Descent Measure Contract Design

## Goal

Define an expert-only, replayable Laurent descent measure contract for the
checked-in `case_008 d=14` boundary profile. The contract turns the existing
profile-only data from #317 into a deterministic measure tuple, exposes a
strict lexicographic decrease predicate, and rejects stale or tampered profile
summaries before accepting any measure.

## Context

No repository `AGENTS.md`, `CLAUDE.md`, `GEMINI.md`, or `CONVENTIONS.md` file is
present in this worktree. `gh issue view 321` was blocked by the sandbox proxy,
so the full issue body supplied by Agent Desk is the authoritative issue
context for this run. A `gh pr list` search for `321` returned only an unrelated
historical Dependabot PR, so there is no relevant open PR context to merge or
adapt.

Current `main` already contains the dependencies for this issue:

- #315 added the validated `case_008 d=14` Laurent boundary fixture.
- #316 exposed the staged `:laurent_native_ecp_boundary` diagnostic.
- #317 added `test/expert/case008_d14_laurent_descent_profile.jl`, including
  `case008_d14_laurent_descent_profile()` and
  `validate_laurent_descent_profile(profile, fixture)`.

The baseline package test command
`julia --project=. -e 'using Pkg; Pkg.test()'` passed before edits. The run
emitted existing Julia 1.12 world-age warnings in Quillen fixture tests.

## Approach Options

Recommended: create a focused expert test file
`test/expert/case008_d14_laurent_descent_measure_contract.jl`. It will include
the #317 profile helper, validate the supplied profile through
`validate_laurent_descent_profile`, derive the measure fields from validated
profile summaries, and keep the comparator local to expert tests. Register the
file in `test/runtests.jl` so the expert group covers the contract.

Alternative: append the measure helper to the existing profile test file. That
would reuse local helper scope, but it would mix two contracts with different
statuses and make the #317 profile-only test harder to read.

Alternative: expose a production measure type or reducer hook under `src/`.
That would imply public API or Laurent reducer support before the later Laurent
move issues exist.

Chosen approach: focused expert-only measure contract file. It keeps this PR
inside the issue scope, avoids production API churn, and makes future reducer
work consume a replayable contract without claiming `case_008 d=14` support.

## Measure Shape

`case008_d14_laurent_descent_measure(profile; fixture = boundary_fixture())`
returns a named tuple with stable fields:

- `case_id = "case_008"`;
- `dimension = 14`;
- `ring_generators = ("u", "v")`;
- `status = :measure_contract`;
- `order = :lexicographic_minimize`;
- `components = (:whole_support_count, :max_entry_terms, :valuation_span,
  :leading_exponent, :leading_entry_index)`;
- `whole_support_count = 7387`;
- `max_entry_terms = 3734`;
- `valuation_span = (97, 93)`;
- `leading_exponent = (49, -5)`;
- `leading_entry_index = 10`.

`valuation_span` is computed from the validated profile as
`(u.max - u.min, v.max - v.min)` in `ring_generators` order. The leading
component uses the first leading-monomial candidate recorded by #317, which is
already sorted by descending exponent and entry-index tie-breaker.

## Validation

The measure constructor must validate the input profile before deriving any
measure field:

1. Check the profile shape for the fields the measure consumes.
2. Check `status == :profile_only`, `case_id == "case_008"`,
   `dimension == 14`, and `ring_generators == ("u", "v")`.
3. Call `validate_laurent_descent_profile(profile, fixture)` and reject any
   result other than `:ok`.

Rejected inputs throw `ArgumentError` with the validation symbol included in
the message. Negative controls cover wrong status, swapped ring generators,
stale support summary, and tampered leading-monomial metadata.

## Comparator

`strictly_decreases_laurent_measure(before, after)::Bool` compares only
measure objects that declare:

- `status == :measure_contract`;
- `order == :lexicographic_minimize`;
- identical `components`.

The comparator builds the component value tuple in declared component order and
uses Julia tuple ordering. Nested tuple fields such as `valuation_span` and
`leading_exponent` therefore compare lexicographically in the same generator
order. `leading_entry_index` is the final tie-breaker. Equal measures do not
strictly decrease.

Tests use synthetic measure objects derived from the validated baseline measure
to prove the comparator has teeth:

- reducing `whole_support_count` is a strict decrease;
- an equal measure is not a strict decrease;
- increasing `max_entry_terms` or `valuation_span` is not a strict decrease.

## Tests

Create `test/expert/case008_d14_laurent_descent_measure_contract.jl`.

Positive checks assert the real d14 fixture measure records:

- `case_id == "case_008"`;
- `dimension == 14`;
- `ring_generators == ("u", "v")`;
- `status == :measure_contract`;
- `order == :lexicographic_minimize`;
- `components == (:whole_support_count, :max_entry_terms, :valuation_span,
  :leading_exponent, :leading_entry_index)`;
- baseline values `7387`, `3734`, `(97, 93)`, `(49, -5)`, and `10`.

Negative controls assert profile rejection for:

- `status != :profile_only`;
- swapped ring generators;
- stale support summary;
- tampered leading-monomial metadata.

Runner coverage asserts `test/runtests.jl` includes the new expert file.

## Verification

Run:

```bash
julia --project=. -e 'include("test/expert/case008_d14_laurent_descent_measure_contract.jl")'
julia --project=. test/runtests.jl expert
julia --project=. -e 'using Pkg; Pkg.test()'
git diff --check
```

Expected result: all commands exit 0.

## Out Of Scope

Do not implement Laurent link witnesses, endpoint reductions, recursive peel
integration, or production support for `case_008 d=14`. Do not make diagonal
monomial balancing or polynomialization part of this measure contract. Do not
change production reducer behavior under `src/`.

## Automatic Decisions

- Visual companion skipped because this is a deterministic test-only contract,
  not a visual design task.
- Clarifying questions skipped because Agent Desk is non-interactive and the
  issue body gives exact fields, values, negative controls, and verification
  commands.
- Recommended approach selected: a new expert-only measure contract file,
  because it keeps #317 profile helpers intact and avoids production API churn.
- Design approval auto-approved under the Standing Answer Policy.
- User review of this written spec auto-approved under the Standing Answer
  Policy so the non-interactive run can continue to implementation planning.
