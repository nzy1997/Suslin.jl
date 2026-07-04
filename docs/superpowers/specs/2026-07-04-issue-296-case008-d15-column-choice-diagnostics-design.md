# Issue 296 Case008 D15 Column-Choice Diagnostics Design

## Goal

Add an offline ToricBuilder `case_008` `d=15` full matrix fixture and a focused
internal diagnostic that evaluates all 15 columns of that matrix as possible
peel columns.

The diagnostic should keep the current boundary explicit: column `15` is the
current peel column from issue #295, it remains unsupported by the current
Laurent reducer, and the report records whether any other column is already
supported or useful for a later preconditioning search.

## Context

Issue #295 added the compact `case_008` `d=15` last-column fixture. That fixture
is appropriate for default reducer regression coverage, but it cannot answer
whether the full `15 x 15` peel matrix contains a better column choice.

The repository already has the exact pattern needed for this issue:

- `test/fixtures/toricbuilder_case008_d16_matrix_boundary.jl` stores a sparse
  offline full-matrix snapshot.
- `test/internal/toricbuilder_case008_d16_column_choice.jl` validates that
  snapshot, evaluates every column, and keeps the report helper internal.
- `test/fixtures/toricbuilder_case008_d15_column_boundary.jl` stores the current
  issue #295 last-column boundary and expected unsupported reducer diagnostic.

## Approach Options

Recommended: create a separate
`test/fixtures/toricbuilder_case008_d15_matrix_boundary.jl` module plus a
focused internal diagnostic
`test/internal/toricbuilder_case008_d15_column_choice.jl`. Derive the matrix by
replaying one certified peel step from the existing `d=16` matrix fixture, store
only sparse string entries, validate the current last column against the issue
#295 compact fixture, and keep the diagnostic unregistered from default tests.

Alternative: extend the existing `d=15` column-boundary fixture with the full
matrix. This reduces the number of files but makes a compact default fixture
materialize much more data than its callers need.

Alternative: derive the `15 x 15` matrix at test runtime from the `d=16` matrix
fixture. This keeps less data in git but violates the issue's request for an
offline stored matrix snapshot and makes the diagnostic depend on replay logic
instead of validating an immutable snapshot.

## Chosen Design

Create `test/fixtures/toricbuilder_case008_d15_matrix_boundary.jl` with module
`ToricBuilderCase008D15MatrixBoundary`.

The module will include the existing issue #295 column-only fixture and reuse
its Laurent parsing helpers. It will define:

- source metadata for `case_008`;
- `EXPECTED_PASSED_PEEL_DIMENSIONS ==
  (30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16)`;
- `FIRST_FAILING_PEEL_DIMENSION == 15`;
- `CURRENT_PEEL_COLUMN_INDEX == 15`;
- `FAILING_INPUT_SPARSE_ENTRIES`, a stored sparse string snapshot of the full
  `15 x 15` matrix;
- `matrix_fixture()`, returning a named tuple with `failing_input_matrix`,
  `current_peel_column`, `current_peel_column_index`, `ring`, and metadata;
- `validate_matrix_fixture(fixture)::Symbol`;
- `corrupted_matrix_entry_negative_control(fixture)`, changing one matrix entry
  in a copied fixture so validation rejects the snapshot.

The sparse matrix snapshot will be generated offline by replaying exactly one
certified peel step from
`ToricBuilderCase008D16MatrixBoundary.matrix_fixture().failing_input_matrix`.
The generated column `15` must equal
`ToricBuilderCase008D15ColumnBoundary.boundary_fixture().failing_column`.

## Diagnostic

Create `test/internal/toricbuilder_case008_d15_column_choice.jl`.

The file will expose local helper
`case008_d15_column_choice_report(fixture =
ToricBuilderCase008D15MatrixBoundary.matrix_fixture())`. The helper must
validate the fixture first and throw `ArgumentError` for invalid fixtures.

The report returns:

```julia
(;
    case_id = "case_008",
    dimension = 15,
    current_peel_column_index = 15,
    candidates = (...),
)
```

Each candidate entry includes:

- `column_index`;
- `is_current_peel_column`;
- `is_unimodular`;
- `unit_entry_count`;
- `laurent_witness_outcome`;
- `laurent_witness_unit_index`;
- `normalized_precondition_status`;
- `normalized_failure_code`;
- `row_preconditioning_outcome`;
- `row_preconditioning_transformed_stage`;
- `status`;
- `failure_code`;
- `supported_by_current_reducer`.

The diagnostic will use `Suslin.diagnose_unimodular_column_reduction(column, R)`
for each column and extract the same stage-detail fields as the existing d16
column-choice diagnostic. It will not add public API.

## Validation

`validate_matrix_fixture` should reject:

- missing metadata;
- wrong case id or source identity;
- wrong source dimensions;
- wrong peel dimension or passed peel dimensions;
- wrong matrix size;
- wrong ring;
- a current peel column index other than `15`;
- a current peel column that differs from the matrix column at index `15`;
- a current peel column that differs from the issue #295 column-only fixture;
- a stored matrix snapshot that differs from `FAILING_INPUT_SPARSE_ENTRIES`.

The negative control should corrupt one copied matrix entry and validate as
`:wrong_snapshot` or `:wrong_column`. The report helper must throw before
producing a report for that corrupted fixture.

## Tests

The focused internal diagnostic test should assert:

- the matrix fixture validates as `:ok`;
- the matrix has dimensions `(15, 15)`;
- `current_peel_column_index == 15`;
- the current peel column equals both the matrix column `15` and the issue #295
  column-only fixture;
- exactly 15 candidate columns are checked;
- candidate indices are exactly `1:15`;
- exactly one candidate is the current peel column;
- each candidate exposes the required diagnostic fields;
- the current candidate is unimodular, has zero unit entries, has
  `status == :unsupported`, and has
  `failure_code == :unsupported_laurent_column_family`;
- supported candidates, if any, have `status == :supported` and no failure code;
- unsupported unimodular candidates have an explicit failure code;
- the corrupted copied fixture is rejected before report generation.

Do not register the new diagnostic in `test/runtests.jl`; the issue asks for a
focused command and the full matrix should stay out of default slow paths.

## Verification

Run:

```bash
julia --project=. -e 'include("test/internal/toricbuilder_case008_d15_column_choice.jl")'
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected result: both commands exit 0. The focused command validates the
fixture, checks exactly 15 candidate columns, identifies column `15` as the
current peel column, verifies that the current column still reports
`:unsupported_laurent_column_family`, and rejects corrupted matrix snapshots.

## Out Of Scope

Do not implement a new reducer stage. Do not run the full `30 x 30`
`case_008` report in default tests. Do not make the column-choice report a
public API. Do not register the focused full-matrix diagnostic in
`test/runtests.jl`.

## Automatic Decisions

- Visual companion: skipped because this is a data-fixture and diagnostic-test
  task.
- Clarifying questions: skipped because Agent Desk requires a non-interactive
  run and the issue body gives concrete acceptance criteria.
- Design approval: accepted automatically under the standing answer policy.
- Approach: selected a separate matrix fixture and focused diagnostic because
  this mirrors the existing d16 full-matrix diagnostic and keeps the compact
  issue #295 fixture unchanged.
- Runtime scope: the full-matrix diagnostic stays unregistered from default
  tests; only the issue verification command loads it.
