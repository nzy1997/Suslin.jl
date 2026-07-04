# Issue 298 Case008 D15 Preconditioning Search Design

## Goal

Add an expert-only bounded search for the `case_008` `d=15` full matrix that returns a real replay-verified elementary preconditioning candidate. The accepted candidate must transform target column `15` into a column accepted by `Suslin.diagnose_unimodular_column_reduction`.

## Context

Issue #296 added the offline `15 x 15` matrix fixture in `test/fixtures/toricbuilder_case008_d15_matrix_boundary.jl`. Issue #297 recorded that the current target column is Laurent-unimodular but unsupported: it has no unit entry, no Laurent unit-creation candidate, no witness unit from the existing witness solve, a normalized ordinary-polynomial column that is not unimodular, and no production row-preconditioning candidate.

The dependency diagnostics are intentionally bounded because full Laurent witness solving for arbitrary d15 column candidates is expensive. Local verification showed the existing d15 matrix fixture and d15 witness profile pass after `Pkg.instantiate()`, and exhaustive full diagnostics over naive right-column additions are too slow for a focused expert test.

No `AGENTS.md`, `CLAUDE.md`, or `GEMINI.md` file is present in this worktree. Follow the existing Julia expert-test style and keep exploratory search out of production APIs.

## Approaches Considered

Recommended: keep the search local to `test/expert/case008_d15_preconditioning_search.jl`. The helper first runs cheap, bounded right-column addition probes using coefficient candidates that begin with `one(R)`. Candidates whose cheap direct-unit and normalized-unimodularity checks fail are recorded as skipped instead of entering expensive witness solving. The helper then runs a bounded row-side synthesis pass: for a configured pivot row, solve a Laurent linear system over the other target-column entries for coefficients that make that pivot entry equal to `one(R)`, emit one elementary left-preconditioning step per nonzero coefficient, replay the steps against the full matrix, and verify the transformed column with the existing reducer diagnostic.

Alternative: run full reducer diagnostics for every right-column candidate. This is the most literal search, but it spends minutes in witness solving on known-unpromising d15 candidates and makes the expert test too slow.

Alternative: hard-code the found row-side coefficients. That would make the test faster but less reviewable; the coefficients are large Laurent polynomials, and recomputing them from the fixture is a better replay-verifiable handoff to #299.

Chosen approach: local expert-test helper with cheap right-side progress records and bounded row-side synthesis. It preserves the issue's production out-of-scope boundary, finds a real `:found` result, and keeps expensive work explicit.

## Search Interface

Add `case008_d15_preconditioning_search(fixture; kwargs...)` in the expert test file. It accepts:

- `max_depth::Integer`, default `1`, for right-column addition probes;
- `side_order`, default `(:right, :left)`;
- `operation_family::Symbol`, default `:column_addition`;
- `column_index`, default `fixture.current_peel_column_index`;
- `source_column_candidates`, default `1:14`;
- `coefficient_candidates`, default `(one(R),)`;
- `right_full_diagnostic_limit::Integer`, default `0`, so right-side candidates are cheap-probed and then skipped before witness solving;
- `row_synthesis_pivots`, default `(1,)`;
- `row_synthesis_max_steps::Integer`, default `14`.

The returned record includes:

- `status`;
- `bounds`;
- `attempt_count`;
- `steps`;
- `transformed_column`;
- `reducer_diagnostic`;
- `progress_summary`.

For `status == :found`, `steps` is a tuple of `Suslin.elementary_preconditioning_step` records. The final matrix is reconstructed with `Suslin.replay_elementary_preconditioning`, checked with `Suslin.verify_elementary_preconditioning`, and supplies `transformed_column`. The reducer diagnostic must report `status == :supported`.

For `status == :not_found`, `steps == ()`, `bounds` preserves the exact requested limits, and `progress_summary` records skipped right probes and exhausted row synthesis bounds.

## Determinism And Bounds

The search order is deterministic:

1. Validate the d15 matrix fixture.
2. Diagnose the unmodified target column once.
3. For right-side depth-1 probes, enumerate source columns in ascending order and coefficient candidates in the supplied order. Run cheap direct-unit and normalized-unimodularity checks. If both fail and `right_full_diagnostic_limit == 0`, record a skip with reason `:expensive_witness_diagnostic_skipped`.
4. If right probes do not find support, enumerate configured row-synthesis pivots. For each pivot, solve
   `sum(coeff[source] * column[source] for source != pivot) == one(R) + column[pivot]`.
5. Build left-side elementary steps targeting the pivot row for each nonzero coefficient, with sources in ascending order. Reject the candidate if it exceeds `row_synthesis_max_steps`.
6. Replay the steps from the original full matrix, verify the replay, diagnose the transformed target column, and return `:found` only when the diagnostic is supported.

The default row-synthesis bounds intentionally find pivot `1` with `14` nonzero left-side steps. The transformed target column has direct unit entry `1`, so the final reducer diagnostic is supported by the existing `:unit_entry` stage.

## Validation And Negative Controls

Add local result validator `_case008_d15_preconditioning_result_is_verified(original_matrix, result)` that rejects anything except a replay-verified `:found` result with nonempty steps and a supported transformed-column diagnostic.

Negative controls:

- Tamper the first found step's factor and assert `Suslin.verify_elementary_preconditioning` fails against the found final matrix.
- Construct a fake `:found` result with an unsupported reducer diagnostic and assert the validator rejects it.
- Run deterministic `:not_found` bounds with `source_column_candidates = ()` and `row_synthesis_pivots = ()`, and assert the exact bounds and skip/progress summary are reported instead of silently passing.

## Tests

Create `test/expert/case008_d15_preconditioning_search.jl`.

Focused verification:

```bash
julia --project=. -e 'include("test/expert/case008_d15_preconditioning_search.jl")'
```

Required Agent Desk verification:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

The expert file stays out of `test/runtests.jl` so the default suite remains bounded.

## Out Of Scope

Do not add a production reducer stage. Do not change peel ordering. Do not require the full `case_008` report to pass. Do not make exploratory search helpers public APIs.

## Automatic Decisions

- Visual companion: skipped because this is algebraic test-harness work with no visual design question.
- Clarifying questions: skipped because Agent Desk is non-interactive and the issue body defines the acceptance contract.
- Design approval: auto-approved under the Standing Answer Policy because there is no interactive user gate and the chosen approach is the conservative local expert-test implementation.
- Approach: selected local expert-test helper with cheap right probes plus bounded row synthesis because full right-side diagnostics are too slow and hard-coded coefficients are less maintainable.
- Execution approach: will use the recommended Superpowers option, subagent-driven development, when the implementation plan hands off execution.
