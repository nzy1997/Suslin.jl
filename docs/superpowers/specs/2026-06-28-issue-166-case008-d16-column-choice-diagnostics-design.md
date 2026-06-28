# Issue 166 Case008 D16 Column-Choice Diagnostics Design

## Goal

Add an offline ToricBuilder `case_008` `d=16` full matrix fixture and a focused
column-choice diagnostic that evaluates every possible peel column from that
matrix.

The diagnostic must make the current boundary explicit: it should identify the
current last column, report each column's unimodularity and reducer support
status, and prove whether any alternative column is currently accepted by
`diagnose_unimodular_column_reduction`.

## Context

Issue #164 added a compact `case_008` `d=16` column-only fixture for the current
last-column failure. Issue #165 added stage-level diagnostic metadata showing
that the last column is a Laurent unimodular column whose witness has no unit
entry and whose Laurent normalization is not an ordinary unimodular column.

That still leaves one question unanswered: whether the failure is caused by a
poor last-column choice. The new fixture must store the full `16 x 16` Laurent
matrix from the same peel path, not just the last column.

The full matrix should stay out of default slow paths. It should be loaded only
by a focused internal/offline diagnostic command.

## Approach Options

Recommended: create a separate fixture file
`test/fixtures/toricbuilder_case008_d16_matrix_boundary.jl` and a focused
diagnostic test `test/internal/toricbuilder_case008_d16_column_choice.jl`. The
matrix fixture stores sparse `(row, column, expression)` entries, materializes
the matrix on demand, validates it against the known last-column snapshot, and
keeps the existing default column fixture unchanged.

Alternative: extend `test/fixtures/toricbuilder_case008_d16_column_boundary.jl`
with the full matrix. This avoids one new fixture module but makes the existing
default internal boundary test materialize the full matrix whenever it calls
`boundary_fixture()`.

Alternative: derive the `d=16` matrix from the checked-in `d=21` fixture at test
runtime. That minimizes stored data but replays supported peel steps during the
diagnostic and fails the issue requirement to use a stored full `d=16` matrix
fixture.

## Chosen Design

Create `test/fixtures/toricbuilder_case008_d16_matrix_boundary.jl` with module
`ToricBuilderCase008D16MatrixBoundary`.

The module will include the existing column-only fixture and reuse its Laurent
parsing helpers. It will define:

- source metadata for `case_008`;
- `EXPECTED_PASSED_PEEL_DIMENSIONS`;
- `FIRST_FAILING_PEEL_DIMENSION == 16`;
- `CURRENT_PEEL_COLUMN_INDEX == 16`;
- `FAILING_INPUT_SPARSE_ENTRIES`, a stored sparse string snapshot of the full
  `16 x 16` matrix;
- `matrix_fixture()`, returning a named tuple with `failing_input_matrix`,
  `current_peel_column`, `current_peel_column_index`, `ring`, and metadata;
- `validate_matrix_fixture(fixture)::Symbol`;
- `corrupted_matrix_entry_negative_control(fixture)`, changing one stored matrix
  entry in a copied fixture so validation rejects it with an explicit code.

The matrix snapshot will be generated offline by replaying five supported peel
steps from the existing `d=21` fixture to dimensions `20, 19, 18, 17, 16`. The
generated `d=16` last column must equal the existing issue #164
`FAILING_COLUMN_ENTRIES` snapshot.

## Diagnostic

Create `test/internal/toricbuilder_case008_d16_column_choice.jl`.

The diagnostic will expose a helper such as `case008_d16_column_choice_report()`
that returns a named tuple:

```julia
(;
    case_id = "case_008",
    dimension = 16,
    current_peel_column_index = 16,
    candidates = (...),
)
```

Each candidate entry will include:

- `column_index`;
- `is_current_peel_column`;
- `is_unimodular`;
- `unit_entry_count`;
- `laurent_witness_outcome`;
- `laurent_witness_unit_index`;
- `normalized_precondition_status`;
- `normalized_failure_code`;
- `status`;
- `failure_code`;
- `supported_by_current_reducer`.

For non-unimodular columns, the diagnostic will still report the unit-entry
count and precondition failure code, with Laurent witness and normalization
fields set to `nothing` because reducer stages are not attempted.

For unimodular columns, the diagnostic will call
`Suslin.diagnose_unimodular_column_reduction(column, R)` and extract the
`:laurent_witness_unit` and `:laurent_normalization` stage details when present.

## Validation

`validate_matrix_fixture` should reject:

- missing metadata;
- wrong case id;
- wrong source dimensions;
- wrong peel dimension or passed peel dimensions;
- wrong matrix size;
- wrong ring;
- a current peel column index other than `16`;
- a current peel column that does not match the existing column-only fixture;
- a stored matrix snapshot that does not match `FAILING_INPUT_SPARSE_ENTRIES`.

The negative control should corrupt one copied matrix entry and validate as
`:wrong_snapshot` or `:wrong_column`. It must reject before any column-choice
report is produced.

## Tests

The focused diagnostic test should assert:

- the matrix fixture validates as `:ok`;
- the matrix has dimensions `(16, 16)`;
- the current peel column index is `16`;
- the current peel column equals the existing d16 column-only fixture;
- exactly 16 candidate columns are checked;
- candidate column indices are exactly `1:16`;
- each candidate has the required diagnostic fields;
- at least one candidate is the current last column;
- every supported candidate, if any, has `status == :supported`;
- if no candidate is supported, every candidate has either
  `is_unimodular == false` or `supported_by_current_reducer == false` with an
  explicit failure code;
- the corrupted copied fixture is rejected before report generation.

Do not register the new diagnostic in `test/runtests.jl`; the issue asks for a
focused command and the full matrix should stay out of default slow paths.

## Verification

Run:

```bash
julia --project=. -e 'include("test/internal/toricbuilder_case008_d16_column_choice.jl")'
julia --project=. -e 'include("test/internal/toricbuilder_case008_d16_column_boundary.jl")'
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected result: all commands exit 0. The first command checks exactly 16
candidate columns, identifies column `16` as the current peel column, and states
whether any alternative column is currently supported.

## Out Of Scope

Do not change peel ordering. Do not implement a new reducer stage. Do not make
`case_008` pass. Do not add ToricBuilder as a package dependency or rerun the
slow cache derivation in default tests.

## Automatic Decisions

- Visual companion: skipped because this is a data-fixture and diagnostic-test
  task.
- Clarifying questions: skipped because Agent Desk requires a non-interactive
  run and the issue body defines the observable interface.
- Approach: selected a separate matrix fixture so the existing default d16
  column-boundary fixture remains compact.
- Runtime scope: the all-column diagnostic stays focused and unregistered from
  `test/runtests.jl` because the issue asks to keep the full matrix out of
  default slow paths.
