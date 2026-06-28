# Issue 164 Case008 D16 Column Boundary Design

## Goal

Add a compact ToricBuilder `case_008` fixture for the current `d=16` Laurent
column-reduction boundary, plus a default internal validator test.

The fixture should let later reducer work load the failing column directly,
without rerunning the slow `30x30` normalization and peel path.

## Context

Issue #149 moved `case_008` past the old `d=21` Laurent witness-unit boundary.
The current bounded investigation now reaches a smaller boundary:

- source case: `case_008`;
- source Q-block dimensions: `(30, 30)`;
- successful peel dimensions: `(30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17)`;
- first failing peel dimension: `16`;
- failing input: only the last column of the `d=16` Laurent matrix;
- current reducer diagnostic: `failure_code == :unsupported_laurent_column_family`.

The issue explicitly keeps the full `d=16` matrix and any column-choice analysis
out of scope.

## Approach Options

Recommended: add a standalone column-only fixture file under `test/fixtures/`
and a default internal validator under `test/internal/`. The fixture records
source metadata, peel metadata, ring metadata, the exact length-16 column, and
the expected unsupported diagnostic. This matches the issue and avoids making
the default test suite pay for the slow derivation path.

Alternative: extend the existing `d=21` fixture to include the later `d=16`
boundary. This would mix a now-supported historical boundary with the current
unsupported boundary and would keep unnecessary matrix-heavy assumptions close
to the new test.

Alternative: store the full `d=16` matrix. That may be useful for a later
column-choice issue, but it is explicitly out of scope here.

## Chosen Design

Create `test/fixtures/toricbuilder_case008_d16_column_boundary.jl` with module
`ToricBuilderCase008D16ColumnBoundary`.

The module will define:

- source metadata for `case_008`, including source cache file, source block,
  source matrix dimensions, and source column-transformation dimensions;
- ring metadata for `GF(2)[u^+/-1, v^+/-1]`;
- `EXPECTED_PASSED_PEEL_DIMENSIONS`;
- `FIRST_FAILING_PEEL_DIMENSION`;
- `FAILING_COLUMN_ENTRIES`, containing exactly 16 Laurent expression strings;
- `boundary_fixture()`, returning a named tuple with parsed column entries;
- `validate_boundary_fixture(fixture)::Symbol`;
- `non_unimodular_negative_control(fixture)`, multiplying the column by `v + 1`.

The fixture will include local Laurent-expression parsing helpers so it does not
load or materialize the full ToricBuilder cache Q-block catalog.

## Validation

`validate_boundary_fixture` should reject:

- missing metadata;
- wrong case id;
- wrong source dimensions;
- wrong passed peel dimensions;
- wrong failing peel dimension;
- wrong column length;
- wrong ring;
- non-unimodular columns;
- columns whose diagnostic no longer matches the stored unsupported boundary;
- columns that do not match the stored column snapshot.

The non-unimodular negative control should return `:not_unimodular`, because it
multiplies every column entry by the nonunit `v + 1`.

## Tests

Create `test/internal/toricbuilder_case008_d16_column_boundary.jl`.

The focused test should assert the issue-required facts:

- `case_id == "case_008"`;
- `source_matrix_dimensions == (30, 30)`;
- `first_failing_peel_dimension == 16`;
- `passed_peel_dimensions == (30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17)`;
- `length(failing_column) == 16`;
- `Suslin.is_unimodular_column(failing_column, ring)`;
- `count(is_unit, failing_column) == 0`;
- `Suslin.diagnose_unimodular_column_reduction(failing_column, ring).failure_code == :unsupported_laurent_column_family`;
- the corrupted `v + 1` negative control validates as `:not_unimodular`.

Register the internal test in `test/runtests.jl` so it runs in the default
public/internal suite.

## Verification

Run:

```bash
julia --project=. -e 'include("test/internal/toricbuilder_case008_d16_column_boundary.jl")'
julia --project=. test/runtests.jl
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected result: all commands exit 0.

## Out Of Scope

Do not implement a new reducer stage. Do not require `case_008` to pass. Do not
store the full `d=16` matrix. Do not add the slow `30x30` derivation path to
the default tests.

## Automatic Decisions

- Visual companion: skipped because this is a data-fixture and validator task.
- Design approval: accepted automatically under the standing non-interactive
  answer policy.
- Fixture scope: column-only, following the issue recommendation.
