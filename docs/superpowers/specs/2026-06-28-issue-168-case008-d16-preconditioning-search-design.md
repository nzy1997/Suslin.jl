# Issue 168 Case008 D16 Preconditioning Search Design

## Goal

Add an expert-only bounded search for elementary preconditioning candidates
around the `case_008` `d=16` boundary. The search should produce a structured
record that either contains a replay-verified candidate sequence whose
transformed column is accepted by an existing reducer diagnostic, or reports
`:not_found` with the exact deterministic bounds that were searched.

## Context

Issue #166 added the full `16 x 16` Laurent matrix fixture and proved the
current last column is unsupported while alternative columns `8`, `9`, and
`10` are accepted by the current reducer. Issue #167 profiled the current last
column's Laurent witness and showed that the witness has no unit entry, so the
existing Laurent witness-unit stage cannot apply directly.

The repository already exposes side-aware elementary preconditioning helpers:
`elementary_preconditioning_step`, `replay_elementary_preconditioning`, and
`verify_elementary_preconditioning`. This issue should reuse those exact replay
and verification primitives rather than adding a production reducer path.

No `AGENTS.md`, `CLAUDE.md`, or `GEMINI.md` file is present in this worktree.
The local style is driven by the existing Julia tests and nearby Superpowers
specs.

## Approaches Considered

Recommended: keep the bounded search local to the requested expert test file,
`test/expert/case008_d16_preconditioning_search.jl`. The file can expose local
helpers for deterministic breadth-first search, result validation, synthetic
negative controls, and the issue-specific test assertions. This satisfies the
interface without making exploratory search semantics part of the package API.

Alternative: add a production internal helper under `src/algorithm`. That would
make the search easier to reuse later, but it prematurely stabilizes candidate
ordering, status fields, and reducer-support semantics for an exploratory issue.

Alternative: hard-code a known successful sequence if one is found. That would
be fast, but it would not satisfy the requested bounded search interface or the
`:not_found` reporting contract.

Chosen approach: local expert-test search helpers. The test file is allowed to
be expert-only and deterministic, and the issue explicitly says not to wire the
search into production reducer logic.

## Search Interface

Add local helper `case008_d16_preconditioning_search(fixture; kwargs...)` in the
expert test file. It accepts the matrix fixture and keyword bounds:

- `max_depth::Integer`
- `side::Symbol`
- `operation_family::Symbol`
- `coefficient_candidates`
- `column_index::Integer`
- `source_column_candidates`

The default issue search uses the full `d=16` matrix, searches `:right` column
operations, targets the current peel column `16`, uses source columns `1:15`,
and uses coefficient candidates `(one(R),)`. With these bounds it should find
the known supported alternative-column route by applying one right-side
elementary operation that adds column `8`, `9`, or `10` to column `16` over
`GF(2)[u^+/-1, v^+/-1]`.

The returned record must always include:

- `status`
- `bounds`
- `attempt_count`
- `steps`
- `transformed_column`
- `reducer_diagnostic`

For `status == :found`, `steps` is a tuple of the exact step records returned
by `elementary_preconditioning_step`, `transformed_column` is the target column
after replay, and `reducer_diagnostic.status == :supported`. The test rejects a
found result if replay verification fails or if the transformed-column
diagnostic is unsupported.

For `status == :not_found`, `steps == ()`, `attempt_count > 0`, and `bounds`
matches the exact requested search bounds. The transformed column and reducer
diagnostic are the original target-column diagnostic so callers still get a
useful failure record.

## Determinism And Bounds

The search is breadth-first by sequence depth. At each depth it enumerates
source columns in ascending order, then coefficient candidates in the supplied
order. Candidate steps are created with `elementary_preconditioning_step` and
verified through `verify_elementary_preconditioning` against their recorded
transformed matrix before the transformed target column is diagnosed.

Only `:right` column-addition search is needed for the issue fixture. The helper
should still make unsupported side or operation-family choices explicit by
throwing `ArgumentError`, so malformed bounds do not silently change semantics.

## Validation And Negative Controls

Add local validation helper `_preconditioning_result_is_verified(original,
result)` that returns `true` only when:

- `status == :found`
- `attempt_count > 0`
- `steps` is nonempty
- `replay_elementary_preconditioning(original, result.steps)` reconstructs the
  matrix that supplied `result.transformed_column`
- `verify_elementary_preconditioning(original, result.steps, final_matrix)` is
  true
- `result.reducer_diagnostic.status == :supported`

The d16 fixture test should assert the found candidate is fully replay-verified
and has a supported reducer diagnostic. A bounded `:not_found` test should use
an empty source-column list or depth `0` with nonzero attempts and exact bounds.

Negative controls:

- Tamper one found step's factor and prove `verify_elementary_preconditioning`
  fails against the found final matrix.
- Build a synthetic known candidate whose replay succeeds, then tamper the
  factor and prove verification fails. This keeps the negative control useful
  even if future bounds return `:not_found`.
- Construct a fake `:found` result with an unsupported transformed-column
  diagnostic and assert `_preconditioning_result_is_verified` rejects it.

## Tests

Create `test/expert/case008_d16_preconditioning_search.jl`.

The focused command is:

```bash
julia --project=. -e 'include("test/expert/case008_d16_preconditioning_search.jl")'
```

The test file remains expert-only and is not registered in `test/runtests.jl`,
matching the issue's exploratory scope and keeping default package tests focused
on public/internal coverage.

## Out Of Scope

Do not add a production reducer stage. Do not change peel ordering. Do not make
the preconditioning result an exported API. Do not rerun ToricBuilder cache
generation or rewrite existing fixtures.

## Automatic Decisions

- Visual companion: skipped because this is algebraic test-harness work.
- Clarifying questions: skipped because Agent Desk requires a non-interactive
  run and the issue body defines the result record and verification contract.
- Approach: selected local expert-test helpers in the requested file because it
  is the narrowest change and preserves the production out-of-scope boundary.
- Default search bounds: selected depth `1`, side `:right`, operation family
  `:column_addition`, target column `16`, source columns `1:15`, coefficient
  candidates `(one(R),)` because dependency issue #166 proved columns `8`, `9`,
  and `10` are already supported by the reducer.
