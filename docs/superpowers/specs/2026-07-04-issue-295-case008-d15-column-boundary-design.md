# Issue 295 Case008 D15 Column Boundary Design

## Goal

Add a compact offline ToricBuilder `case_008` fixture for the current `d=15`
Laurent column-reduction boundary, plus a default internal validator test.

The fixture should let later reducer work load the failing column directly,
without rerunning the slow `30x30` normalization and peel path.

## Context

The repository already has `case_008` boundary fixtures for the earlier
`d=21` and `d=16` stages. The `d=16` column is now supported by the current
Laurent column reducer, and `test/fixtures/toricbuilder_case008_d16_matrix_boundary.jl`
stores a certified `16x16` matrix snapshot that can replay the next peel step.

Issue #295 asks for the next unsupported boundary:

- source case: `case_008`;
- source Q-block dimensions: `(30, 30)`;
- source column-transformation dimensions: `(60, 60)`;
- passed peel dimensions:
  `(30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16)`;
- first failing peel dimension: `15`;
- ring description: `GF(2)[u^+/-1, v^+/-1]`;
- stored data: only the failing length-15 column;
- current reducer diagnostic:
  `failure_code == :unsupported_laurent_column_family`.

## Approach Options

Recommended: replay one certified peel step from the checked-in `d=16` matrix
fixture, extract the next `d=15` last column, and store that column as static
Laurent expression strings in a new column-only fixture. This follows the issue
recommendation, keeps the default validator offline, and reuses the shape of
the existing `d=16` column boundary fixture.

Alternative: store the full derived `15x15` matrix. That would simplify future
column-choice analysis, but the issue explicitly keeps full matrix storage out
of scope.

Alternative: teach the reducer a new `d=15` stage now. That would move the
boundary rather than fixture it, and the issue explicitly defers reducer work.

## Chosen Design

Create `test/fixtures/toricbuilder_case008_d15_column_boundary.jl` with module
`ToricBuilderCase008D15ColumnBoundary`.

The module will define:

- source metadata for `case_008`, including source cache file, source block,
  source matrix dimensions, and source column-transformation dimensions;
- ring metadata for `GF(2)[u^+/-1, v^+/-1]`;
- `EXPECTED_PASSED_PEEL_DIMENSIONS`;
- `FIRST_FAILING_PEEL_DIMENSION`;
- `FAILING_COLUMN_ENTRIES`, containing exactly 15 Laurent expression strings;
- `boundary_fixture()`, returning a named tuple with parsed column entries;
- `validate_boundary_fixture(fixture)::Symbol`;
- `non_unimodular_negative_control(fixture)`, multiplying the column by `v + 1`.

The fixture will keep local Laurent-expression parsing helpers, matching the
`d=16` column fixture. It will not load the ToricBuilder cache catalog, and it
will not materialize a full `15x15` matrix during validation.

## Validation

`validate_boundary_fixture` should reject:

- missing metadata;
- wrong case id or source identity;
- wrong source dimensions;
- wrong passed peel dimensions;
- wrong failing peel dimension;
- wrong column length;
- wrong ring;
- non-unimodular columns;
- unit-entry profiles other than zero unit entries;
- columns whose diagnostic no longer matches the stored unsupported boundary;
- columns that do not match the stored column snapshot.

The diagnostic match is intentionally compact: the current expected behavior is
`status == :unsupported`, `failure_code == :unsupported_laurent_column_family`,
and `column_length == 15`.

## Tests

Create `test/internal/toricbuilder_case008_d15_column_boundary.jl`.

The focused test should assert the issue-required facts:

- `validate_boundary_fixture(fixture) == :ok`;
- `case_id == "case_008"`;
- `source_matrix_dimensions == (30, 30)`;
- `source_column_transformation_dimensions == (60, 60)`;
- `first_failing_peel_dimension == 15`;
- `passed_peel_dimensions ==
  (30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16)`;
- `length(failing_column) == 15`;
- `ring_description == "GF(2)[u^+/-1, v^+/-1]"`;
- `Suslin.is_unimodular_column(failing_column, ring)`;
- `count(is_unit, failing_column) == 0`;
- `Suslin.diagnose_unimodular_column_reduction(failing_column, ring).failure_code ==
  :unsupported_laurent_column_family`;
- `Suslin.reduce_unimodular_column(failing_column, ring)` throws while the
  boundary remains unsupported;
- the corrupted `v + 1` negative control validates as `:not_unimodular`;
- `Suslin.reduce_unimodular_column` throws for the corrupted column.

Register the internal test in `test/runtests.jl` so it runs in the default
public/internal suite.

## Verification

Run:

```bash
julia --project=. -e 'include("test/internal/toricbuilder_case008_d15_column_boundary.jl")'
julia --project=. test/runtests.jl
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected result: all commands exit 0.

## Out Of Scope

Do not implement a new reducer stage. Do not require the full `case_008` report
to pass. Do not store a full `15x15` matrix. Do not add the slow `30x30`
derivation path to default tests.

## Automatic Decisions

- Visual companion: skipped because this is a data-fixture and validator task.
- Clarifying questions: skipped because this is a non-interactive Agent Desk run
  and the issue body gives concrete acceptance criteria.
- Design approval: accepted automatically under the standing answer policy.
- Fixture scope: column-only, following the issue recommendation.
- Derivation source: one peel step from the checked-in `d=16` matrix fixture.
