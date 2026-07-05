# Issue 315 Case008 D14 Laurent Boundary Design

## Goal

Add a compact offline ToricBuilder `case_008` fixture for the post-d15 `d=14`
Laurent boundary, plus a default internal validator that records the new active
certificate-construction layer without claiming full `case_008` success.

## Context

No `AGENTS.md`, `CLAUDE.md`, `GEMINI.md`, or `CONVENTIONS.md` file is present
in this worktree. `gh issue view` for #315 is blocked by the sandbox proxy, so
this design uses the full issue body supplied by Agent Desk plus local merged
context from the `case_008 d=15` fixture and post-d15 bounded report gate.

Issue #300 records that explicit bounded `case_008` evidence now passes
determinant classification and Laurent normalization, completes the old d15
layer, and reaches certificate construction at `current d=14`. This issue asks
to preserve that new boundary as offline fixture data. The fixture is a
boundary target, not a reducer improvement and not a promise that all of
`case_008` passes.

## Approach Options

Recommended: follow the existing `d=15` column-boundary fixture shape and store
the exact `d=14` last column as static Laurent expression strings. Generate the
strings once by replaying one certified peel step from the checked-in
`test/fixtures/toricbuilder_case008_d15_matrix_boundary.jl` fixture. The
validator then checks metadata, ring identity, exact column identity,
unimodularity, nonzero and term-count statistics, and post-d15 provenance
metadata.

Alternative: derive the `d=14` column at fixture load time from the d15 matrix
fixture. That would reduce new file size, but it would make the fixture depend
on a large matrix parse and current peel-step implementation during every
default internal test.

Alternative: run the live d14 Laurent reducer diagnostic in the validator. That
would be too expensive for a default fixture gate and drifts into the out-of-
scope d14 reducer work.

Chosen approach: static column-only fixture plus metadata validation. This
matches the local d15 template, keeps the test offline, and makes the post-d15
boundary reviewable without a local ToricBuilder checkout.

## Fixture Shape

Create `test/fixtures/toricbuilder_case008_d14_column_boundary.jl` with module
`ToricBuilderCase008D14ColumnBoundary`.

The module will define:

- source metadata for `case_008`, `case_008.jls`, source block
  `:column_transformation_upper_left_q_block`, source matrix dimensions
  `(30, 30)`, and source column-transformation dimensions `(60, 60)`;
- ring metadata for `GF(2)[u^+/-1, v^+/-1]`;
- `EXPECTED_PASSED_PEEL_DIMENSIONS ==
  (30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15)`;
- `FIRST_FAILING_PEEL_DIMENSION == 14`;
- `EXPECTED_LAST_COLUMN_NONZERO_COUNT == 14`;
- `EXPECTED_MAX_ENTRY_TERM_COUNT == 3734`;
- `EXPECTED_BOUNDARY_PROVENANCE` with certificate-construction status,
  `current_peel_dimension == 14`, `last_completed_peel_dimension == 15`,
  `failure_code == :unsupported_laurent_column_family`, and
  `old_d15_boundary_cleared == true`;
- `FAILING_COLUMN_ENTRIES`, containing exactly 14 static Laurent expression
  strings generated from the checked-in d15 matrix fixture;
- `boundary_fixture()`, returning a named tuple with parsed column entries;
- `validate_boundary_fixture(fixture)::Symbol`;
- `non_unimodular_negative_control(fixture)`, multiplying the column by
  `v + 1`.

The validator will not run a full d14 reducer diagnostic. It records the
current diagnostic boundary as provenance metadata because issue #315 is about
the explicit bounded report boundary, and d14 reducer construction is out of
scope.

## Validation

`validate_boundary_fixture` should return `:ok` for the recorded fixture and
stable rejection symbols for:

- missing metadata;
- wrong case id or source identity;
- wrong source dimensions;
- wrong passed peel dimensions;
- wrong current peel dimension;
- a fixture that still claims the old d15 current boundary;
- wrong column length;
- wrong ring;
- wrong column statistics;
- non-unimodular columns;
- unit-entry profiles other than zero unit entries;
- columns that do not match the stored column snapshot.

The old d15 negative control is explicit: a fixture with
`first_failing_peel_dimension = 15` and provenance
`failure_code = :unsupported_laurent_column_family` must validate as
`:old_d15_boundary`.

## Tests

Create `test/internal/toricbuilder_case008_d14_column_boundary.jl`.

The focused test should assert:

- `validate_boundary_fixture(fixture) == :ok`;
- `case_id == "case_008"`;
- source block and dimensions match the issue;
- the Laurent ring has generators `("u", "v")` over `GF(2)`;
- `first_failing_peel_dimension == 14`;
- passed dimensions include d15;
- `length(failing_column) == 14`;
- last-column nonzero count is `14`;
- maximum entry term count is `3734`;
- provenance records certificate construction at current d14 after d15;
- the old d15 unsupported boundary validates as `:old_d15_boundary`;
- corrupted exact-column data validates as `:wrong_column`;
- the non-unimodular control validates as `:not_unimodular`.

Register the internal test in `test/runtests.jl` next to the existing d15 and
d16 case_008 boundary validators. The required expert diagnostics command does
not need to run the d14 reducer; it should continue to pass unchanged.

## Verification

Run:

```bash
julia --project=. -e 'include("test/internal/toricbuilder_case008_d14_column_boundary.jl")'
julia --project=. -e 'include("test/expert/laurent_column_reduction_diagnostics.jl")'
julia --project=. test/runtests.jl
julia --project=. -e 'using Pkg; Pkg.test()'
git diff --check
```

Expected result: all commands exit 0.

## Out Of Scope

Do not implement a new d14 Laurent reducer. Do not claim full `case_008`
success. Do not require a local ToricBuilder checkout. Do not make diagonal
monomial balancing or polynomialization the primary Laurent algorithm.

## Automatic Decisions

- Visual companion skipped because this is a data-fixture and validator task.
- Clarifying questions skipped because Agent Desk is non-interactive and the
  issue body gives exact fixture metadata and verification commands.
- Recommended approach selected: static column-only fixture generated from the
  checked-in d15 matrix fixture, because it follows the d15 template and keeps
  default validation offline.
- The d14 reducer diagnostic is stored as bounded-report provenance metadata
  instead of executed by the validator, because the issue is a boundary fixture
  and d14 Laurent reducer work is out of scope.
- Design approval auto-approved under the Standing Answer Policy.
