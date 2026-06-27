# Issue 148 Case 008 D21 Laurent Column Boundary Design

## Goal

Extract the `case_008` Laurent-normalized column-reduction boundary at peel
dimension `d=21` into a small offline fixture. The fixture should give the next
algorithm repair a deterministic target without requiring a full `30x30`
`case_008` certificate run.

## Context

Issue #131 established that `case_008` materializes from
`test/fixtures/toricbuilder_cache_q_blocks.jl`, passes Laurent determinant
classification and normalization, then passes peel steps `d=30` through `d=22`.
The first unsupported Laurent column occurs at `d=21`, where
`diagnose_unimodular_column_reduction` reports an unsupported Laurent column
family. Issue #148 is test-support work only: it must preserve that boundary
as evidence, not repair the reducer.

The current project already has the pattern to mirror:

- `test/fixtures/toricbuilder_case010_column_boundary.jl`
- `test/internal/toricbuilder_case010_column_boundary.jl`

The `case_008` version should keep the same fixture-and-validator shape while
using a stored boundary snapshot so routine validation stays fast.

## Approach Options

Recommended: store the `d=21` boundary as a snapshot fixture with a strict local
validator. The fixture records the failing `21x21` input matrix, derives the
last column as the failing column, stores provenance metadata, and validates
the diagnostic without rerunning normalization and nine peel steps.

Alternative: store the snapshot and add a separate optional reconstruction
script that recomputes the boundary from `case_008` and compares it with the
snapshot. This strengthens provenance, but adds maintenance surface beyond the
focused issue.

Alternative: expose both a snapshot fixture and a live derivation function from
the fixture module. This is flexible, but it blurs the issue goal by tying the
offline fixture back to the slow certificate path.

The chosen design is the snapshot fixture with strict validation. If provenance
reconstruction becomes necessary, it should be a separate follow-up.

## Fixture Module

Add `test/fixtures/toricbuilder_case008_d21_column_boundary.jl` with a module
named `ToricBuilderCase008D21ColumnBoundary`.

The module will include `toricbuilder_cache_q_blocks.jl` for source metadata and
matrix materialization, but normal validation will not rerun the full
normalization or peel path. The module will expose:

- `boundary_fixture()`: returns a named tuple with `case_id`, `source_entry`,
  `original_matrix`, `source_matrix_dimensions`, `source_sparse_entry_count`,
  `normalization_provenance`, `passed_peel_dimensions`,
  `first_failing_peel_dimension`, `failing_input_matrix`, `failing_column`,
  `ring`, `ring_description`, and `expected_diagnostic`.
- `validate_boundary_fixture(fixture)`: returns `:ok` only for the expected
  `case_008` `d=21` boundary and returns stable rejection symbols for corrupted
  fixtures.
- `non_unimodular_negative_control(fixture = boundary_fixture())`: returns a
  fixture with the failing column multiplied by a nonunit such as `v + 1`, so
  validation must reject it as `:not_unimodular`.

The snapshot data should be stored as sparse coordinate entries over
`GF(2)[u^+/-1, v^+/-1]`, using the same Julia Laurent-expression string style
as `toricbuilder_cache_q_blocks.jl`. The fixture will materialize the `21x21`
matrix from those entries and set `failing_column` to its last column.

## Validation

`validate_boundary_fixture` should be strict but local. It will check:

- all required fields are present;
- `fixture.case_id == "case_008"`;
- `fixture.source_matrix_dimensions == (30, 30)`;
- `fixture.first_failing_peel_dimension == 21`;
- `fixture.passed_peel_dimensions == (30, 29, 28, 27, 26, 25, 24, 23, 22)`;
- `fixture.failing_input_matrix` is `21x21`;
- `length(fixture.failing_column) == 21`;
- the ring description and actual ring are `GF(2)[u^+/-1, v^+/-1]`;
- every failing-column entry belongs to the stored ring;
- `Suslin.is_unimodular_column(fixture.failing_column, fixture.ring)` is true;
- `Suslin.diagnose_unimodular_column_reduction(fixture.failing_column,
  fixture.ring)` returns `status == :unsupported` and
  `failure_code == :unsupported_laurent_column_family`.

Expected rejection symbols include `:missing_metadata`, `:wrong_case`,
`:wrong_peel_dimension`, `:wrong_column_length`, `:wrong_ring`,
`:not_unimodular`, and `:wrong_diagnostic`.

## Tests

Add `test/internal/toricbuilder_case008_d21_column_boundary.jl` and include it in
the `internal` group in `test/runtests.jl`.

The focused test will assert:

- the fixture file exists and loads;
- `fixture.case_id == "case_008"`;
- the original source matrix has size `30x30`;
- `fixture.first_failing_peel_dimension == 21`;
- `fixture.passed_peel_dimensions` records exactly `30` through `22`;
- the failing input matrix has size `21x21`;
- `length(fixture.failing_column) == 21`;
- the ring is `GF(2)[u^+/-1, v^+/-1]` with generators `("u", "v")`;
- the failing column is unimodular;
- `diagnose_unimodular_column_reduction` returns the unsupported Laurent-family
  diagnostic required by the issue;
- `validate_boundary_fixture(boundary_fixture()) == :ok`;
- the non-unimodular negative control is rejected as `:not_unimodular`;
- a wrong-ring or wrong-dimension fixture is rejected instead of accepted.

The issue verification command is:

```bash
julia --project=. -e 'include("test/internal/toricbuilder_case008_d21_column_boundary.jl")'
```

Repository-level verification remains:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Error Handling

Snapshot materialization should throw clear `ArgumentError`s for malformed
Laurent expressions, unexpected case-catalog shape, or invalid sparse-entry
coordinates. Validation itself should return stable rejection symbols for
expected negative controls instead of throwing. Interrupts should be rethrown.

Diagnostic validation should inspect the structured diagnostic fields rather
than matching only a message string.

## Out Of Scope

Do not fix the `case_008` Laurent column reducer. Do not require
`laurent_gl_factorization_certificate` to pass for the original `case_008`
matrix. Do not add `case_008` to the default ToricBuilder cache report. Do not
add a live reconstruction script in this issue.
