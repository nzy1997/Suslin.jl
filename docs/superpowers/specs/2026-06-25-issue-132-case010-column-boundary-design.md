# Issue 132 Case 010 Laurent Column Boundary Design

## Goal

Extract the exact `case_010` Laurent column-reduction boundary into offline test
support so the current algorithm limit is reproducible without rerunning the
full ToricBuilder cache report.

## Chosen Approach

Create a focused fixture module in
`test/fixtures/toricbuilder_case010_column_boundary.jl` and a focused internal
validator in `test/internal/toricbuilder_case010_column_boundary.jl`. The fixture
will materialize `case_010` from `test/fixtures/toricbuilder_cache_q_blocks.jl`,
call `normalize_laurent_gl_matrix`, then run the same Laurent column-peel step
sequence on the normalized matrix until a call to
`Suslin.reduce_unimodular_column(column, R)` throws the current unsupported
length-5 diagnostic.

This is narrower than adding public helpers or changing the status report. It
uses existing internal Suslin functions and fixture patterns, records enough
metadata for tests to explain the boundary, and leaves the algorithm behavior
unchanged.

## Data Shape

The fixture module will expose:

- `boundary_fixture()`: returns a named tuple containing `case_id`,
  `original_matrix`, `normalization`, `normalized_matrix`,
  `first_failing_peel_dimension`, `failing_column`, `ring`, and
  `expected_diagnostic`.
- `validate_boundary_fixture(fixture)`: returns `:ok` only when the fixture is
  the expected `case_010` boundary, and returns `:not_unimodular` before
  accepting a perturbed non-unimodular column.

Validation will check that the failing column has length `5`, lives over a
Laurent ring with variables `u` and `v`, is unimodular, and still makes
`Suslin.reduce_unimodular_column(column, R)` throw an `ArgumentError` containing
`unsupported exact unimodular column reduction`.

## Error Handling

The fixture extraction will rethrow unexpected failures. The validator will
return stable symbols for expected rejection cases, including `:not_unimodular`
for the negative control required by the issue. Unsupported-column acceptance is
only valid when the thrown diagnostic contains the current expected substring.

## Testing

The internal test will include the fixture module and assert:

- the source matrix comes from `case_010`;
- the normalized matrix has determinant one after Laurent normalization;
- the failing peel dimension is `5`;
- the failing column is over `GF(2)[u^+/-1, v^+/-1]`;
- `validate_boundary_fixture(boundary_fixture()) == :ok`;
- perturbing one failing-column entry makes the validator return
  `:not_unimodular`.

The verification command for the issue is:

```bash
julia --project=. -e 'include("test/internal/toricbuilder_case010_column_boundary.jl")'
```

The full repository verification remains:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Out Of Scope

Do not fix `case_010`. Do not change `reduce_unimodular_column`. Do not change
the ToricBuilder cache status report or its recorded `case_010` route status.
