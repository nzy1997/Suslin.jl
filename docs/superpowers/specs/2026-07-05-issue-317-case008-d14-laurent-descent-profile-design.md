# Issue 317 Case008 D14 Laurent Descent Profile Design

## Goal

Add an expert-only profile for the `case_008 d=14` Laurent boundary column.
The profile records Newton-support, valuation, leading-monomial, and term-count
data from the checked-in #315 fixture. It is acceptance evidence for future
Laurent-native ECP design only; it must not choose a final descent measure or
add reducer support.

## Context

No repository `AGENTS.md` is present in this worktree. Live GitHub context was
read through the GitHub connector because the sandbox blocks `gh` network
access. Issue #317 is open with no comments. Dependency issue #315 is closed by
PR #318 and provides `test/fixtures/toricbuilder_case008_d14_column_boundary.jl`.
Dependency issue #316 is closed by PR #319 and exposes the staged
`:laurent_native_ecp_boundary` diagnostic. Current `main` is already at the
merged #319 state.

The #315 fixture records a validated length-14 Laurent column over
`GF(2)[u^+/-1, v^+/-1]`, with 14 nonzero entries and maximum entry term count
3734. Issue #317 asks for deterministic profile data computed from that
fixture, not from a ToricBuilder checkout.

## Approach Options

Recommended: create a focused expert test
`test/expert/case008_d14_laurent_descent_profile.jl` with local test-only
helpers. The helpers compute the profile directly from the fixture, validate
the full profile by recomputing it, and include negative controls for fixture,
summary, term-count, ring-generator, and status tampering. Register the file in
the expert runner so it is covered by `test/runtests.jl expert`.

Alternative: put reusable helpers in the d14 fixture module. That would make
the fixture heavier and blur the line between boundary data and exploratory
profile analysis.

Alternative: add production diagnostic fields to `Suslin.diagnose_unimodular_column_reduction`.
That would expose an unsettled descent profile through production APIs and
risks implying Laurent-native reducer support.

Chosen approach: focused expert-only test helper. It keeps the evidence
offline and deterministic, follows existing case008 profile patterns, and
preserves production reducer behavior.

## Profile Shape

The profile is a named tuple with stable fields:

- `case_id = "case_008"`;
- `dimension = 14`;
- `ring_generators = ("u", "v")`;
- `nonzero_entries = 14`;
- `max_entry_terms = 3734`;
- `entry_term_counts`, one count per column entry;
- `valuation_ranges`, a named tuple keyed by `u` and `v`, each storing the
  minimum and maximum exponent over the whole column support;
- `newton_support_summary`, including per-entry support counts and bounds,
  whole-column support count, whole-column bounds, and generator order;
- `leading_monomial_candidates`, including candidate records and explicit
  ordering metadata;
- `candidate_measure_families =
  (:newton_support, :valuation, :leading_monomial)`;
- `status = :profile_only`.

Support extraction uses `exponents(entry)` from Oscar Laurent elements. Zero
entries have zero support count and `nothing` bounds, though the recorded
fixture should have no zero entries.

## Leading Candidate Ordering

The leading-monomial candidate helper records enough metadata to reproduce the
ordering:

- `generator_order = ("u", "v")`;
- `order = :lexicographic_descending`;
- `tie_breaker = :entry_index_ascending`;
- each candidate records `entry_index`, `leading_exponent`, and `term_count`.

For each entry, the leading exponent is the lexicographic maximum exponent
tuple in the generator order. The global candidate list is sorted by descending
leading exponent and then ascending entry index. This is only a candidate
ordering for future design; it is not the selected Laurent descent measure.

## Validation

`validate_laurent_descent_profile(profile, fixture)::Symbol` recomputes the
profile from the supplied fixture and compares the full named tuple. It returns
`:ok` for the recorded profile and stable rejection symbols for:

- invalid fixture metadata or fixture corruption;
- wrong case id, dimension, or ring-generator metadata;
- wrong nonzero or term-count summary;
- wrong support or valuation summary;
- wrong leading-candidate summary or ordering metadata;
- missing required candidate measure families;
- any status other than `:profile_only`.

This validation strategy makes stale or hand-edited summaries fail without
storing a second hard-coded copy of the large support data.

## Tests

Create `test/expert/case008_d14_laurent_descent_profile.jl`.

The positive test asserts:

- the d14 fixture validates through `ToricBuilderCase008D14ColumnBoundary`;
- the profile records `case_id == "case_008"`;
- `dimension == 14`;
- `ring_generators == ("u", "v")`;
- `nonzero_entries == 14`;
- `max_entry_terms == 3734`;
- candidate measure families include `:newton_support`, `:valuation`, and
  `:leading_monomial`;
- `status == :profile_only`;
- support, valuation, and leading-candidate summaries have the expected shape
  and are accepted by the validator.

Negative controls mutate:

- the fixture column;
- ring-generator metadata;
- term-count summary;
- support summary;
- leading-candidate ordering metadata;
- status from `:profile_only` to `:supported`.

Each control must be rejected by `validate_laurent_descent_profile`.

## Verification

Run:

```bash
julia --project=. -e 'include("test/expert/case008_d14_laurent_descent_profile.jl")'
julia --project=. test/runtests.jl expert
julia --project=. -e 'using Pkg; Pkg.test()'
git diff --check
```

Expected result: all commands exit 0.

## Out Of Scope

Do not implement Laurent ECP. Do not add Laurent link witnesses, endpoint
reductions, Laurent normality/conjugation replay, or recursive Laurent peel
integration. Do not choose the final Laurent descent measure. Do not make
diagonal monomial balancing or polynomialization the primary route.

## Automatic Decisions

- Visual companion skipped because this is a deterministic data-profile task.
- Clarifying questions skipped because Agent Desk is non-interactive and the
  issue body gives exact profile fields and verification commands.
- Recommended approach selected: focused expert-only profile file with local
  helpers and expert-runner registration, because this avoids production API
  churn and keeps the profile clearly test-only.
- Design approval auto-approved under the Standing Answer Policy.
- The spec date uses the issue-series date, 2026-07-05, to match the adjacent
  #315 and #316 Superpowers specs already in the repository.
