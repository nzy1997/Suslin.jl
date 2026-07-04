# Issue 296 Case008 D15 Column-Choice Diagnostics Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an offline full `case_008` `d=15` matrix fixture and a focused internal diagnostic that evaluates all 15 columns as possible peel columns.

**Architecture:** Follow the existing d16 full-matrix fixture and column-choice diagnostic. The new fixture stores sparse Laurent expression strings, materializes the matrix on demand, validates the current column against the issue #295 d15 column-only fixture, and stays out of default test registration.

**Tech Stack:** Julia, Oscar.jl Laurent rings, existing Suslin private peel helper for one-off offline derivation, existing Suslin reducer diagnostics, `Test`.

## Global Constraints

- Repository: `nzy1997/Suslin.jl`.
- Issue: `#296 Add a case_008 d=15 full-matrix column-choice diagnostic`.
- Base branch: `main`.
- Worker branch: `agent/issue-296-add-a-case_008-d-15-full-matrix-column-choice-di-run-1`.
- Fixture module path must be `test/fixtures/toricbuilder_case008_d15_matrix_boundary.jl`.
- Internal diagnostic path must be `test/internal/toricbuilder_case008_d15_column_choice.jl`.
- `case_id == "case_008"`.
- `source_matrix_dimensions == (30, 30)`.
- `source_column_transformation_dimensions == (60, 60)`.
- `passed_peel_dimensions == (30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16)`.
- `first_failing_peel_dimension == 15`.
- `current_peel_column_index == 15`.
- `ring_description == "GF(2)[u^+/-1, v^+/-1]"`.
- The matrix fixture must store the exact offline `15 x 15` sparse snapshot derived by one certified peel step from `ToricBuilderCase008D16MatrixBoundary.matrix_fixture().failing_input_matrix`.
- The fixture validator must reject a corrupted stored matrix entry or current-column entry with `:wrong_snapshot` or `:wrong_column`.
- The report helper must validate the fixture and throw `ArgumentError` before producing a report for a corrupted fixture.
- The diagnostic must check exactly 15 candidate columns and expose the report fields requested in issue #296.
- The current column must report `failure_code == :unsupported_laurent_column_family`.
- Do not implement a reducer stage.
- Do not run the full `30 x 30` `case_008` report in default tests.
- Do not make the column-choice report public API.
- Do not register `test/internal/toricbuilder_case008_d15_column_choice.jl` in `test/runtests.jl`.
- The focused verification command is `julia --project=. -e 'include("test/internal/toricbuilder_case008_d15_column_choice.jl")'`.
- The full verification command is `julia --project=. -e 'using Pkg; Pkg.test()'`.

---

## File Structure

- Create `test/fixtures/toricbuilder_case008_d15_matrix_boundary.jl`: matrix fixture module, sparse string entries, materialization, validator, and corrupted-copy helper.
- Create `test/internal/toricbuilder_case008_d15_column_choice.jl`: focused diagnostic helper and tests. This file is intentionally not registered in `test/runtests.jl`.
- Keep `test/fixtures/toricbuilder_case008_d15_column_boundary.jl`, `test/runtests.jl`, and public API files unchanged.

---

### Task 1: Add The D15 Matrix Fixture And Column-Choice Diagnostic

**Files:**
- Create: `test/internal/toricbuilder_case008_d15_column_choice.jl`
- Create: `test/fixtures/toricbuilder_case008_d15_matrix_boundary.jl`
- Test: `test/fixtures/toricbuilder_case008_d15_column_boundary.jl`
- Test: `test/fixtures/toricbuilder_case008_d16_matrix_boundary.jl`

**Interfaces:**
- Consumes: `ToricBuilderCase008D16MatrixBoundary.matrix_fixture()` and `ToricBuilderCase008D16MatrixBoundary.validate_matrix_fixture(fixture)`.
- Consumes: `ToricBuilderCase008D15ColumnBoundary.boundary_fixture()` and `ToricBuilderCase008D15ColumnBoundary._parse_laurent_value(R, u, v, value)`.
- Produces: module `ToricBuilderCase008D15MatrixBoundary`.
- Produces: `matrix_fixture() -> NamedTuple`.
- Produces: `validate_matrix_fixture(fixture)::Symbol`.
- Produces: `corrupted_matrix_entry_negative_control(fixture = matrix_fixture())`.
- Produces: local helper `case008_d15_column_choice_report(fixture = ToricBuilderCase008D15MatrixBoundary.matrix_fixture()) -> NamedTuple`.

- [ ] **Step 1: Write the failing focused diagnostic test**

Create `test/internal/toricbuilder_case008_d15_column_choice.jl`:

```julia
using Test
using Suslin
using Oscar

const TORICBUILDER_CASE008_D15_MATRIX_BOUNDARY_PATH =
    joinpath(@__DIR__, "..", "fixtures", "toricbuilder_case008_d15_matrix_boundary.jl")
const TORICBUILDER_CASE008_D15_COLUMN_BOUNDARY_PATH =
    joinpath(@__DIR__, "..", "fixtures", "toricbuilder_case008_d15_column_boundary.jl")

function _case008_d15_stage_detail(diagnostic, stage::Symbol)
    hasproperty(diagnostic, :stage_details) || return nothing
    idx = findfirst(detail -> detail.stage == stage, diagnostic.stage_details)
    return idx === nothing ? nothing : diagnostic.stage_details[idx]
end

function _case008_d15_column(M, column_index::Int)
    return [M[row, column_index] for row in 1:nrows(M)]
end

function _case008_d15_detail_property(detail, field::Symbol)
    detail === nothing && return nothing
    return hasproperty(detail, field) ? getproperty(detail, field) : nothing
end

function _case008_d15_candidate_report(M, R, column_index::Int, current_index::Int)
    column = _case008_d15_column(M, column_index)
    unit_entry_count = count(is_unit, column)
    is_unimodular = Suslin.is_unimodular_column(column, R)
    diagnostic = Suslin.diagnose_unimodular_column_reduction(column, R)
    witness = _case008_d15_stage_detail(diagnostic, :laurent_witness_unit)
    normalization = _case008_d15_stage_detail(diagnostic, :laurent_normalization)
    row_preconditioning =
        _case008_d15_stage_detail(diagnostic, :laurent_elementary_row_preconditioning)
    return (;
        column_index,
        is_current_peel_column = column_index == current_index,
        is_unimodular,
        unit_entry_count,
        laurent_witness_outcome = _case008_d15_detail_property(witness, :outcome),
        laurent_witness_unit_index = _case008_d15_detail_property(witness, :witness_unit_index),
        normalized_precondition_status = _case008_d15_detail_property(normalization, :normalized_status),
        normalized_failure_code = _case008_d15_detail_property(normalization, :normalized_failure_code),
        row_preconditioning_outcome = _case008_d15_detail_property(row_preconditioning, :outcome),
        row_preconditioning_transformed_stage =
            _case008_d15_detail_property(row_preconditioning, :transformed_stage),
        status = diagnostic.status,
        failure_code = diagnostic.failure_code,
        supported_by_current_reducer = diagnostic.status == :supported,
    )
end

function case008_d15_column_choice_report(fixture = ToricBuilderCase008D15MatrixBoundary.matrix_fixture())
    validation = ToricBuilderCase008D15MatrixBoundary.validate_matrix_fixture(fixture)
    validation == :ok ||
        throw(ArgumentError("invalid case_008 d15 matrix fixture: $(validation)"))

    M = fixture.failing_input_matrix
    R = fixture.ring
    d = fixture.first_failing_peel_dimension
    candidates = tuple(
        (
            _case008_d15_candidate_report(
                M,
                R,
                column_index,
                fixture.current_peel_column_index,
            )
            for column_index in 1:d
        )...,
    )
    return (;
        case_id = fixture.case_id,
        dimension = d,
        current_peel_column_index = fixture.current_peel_column_index,
        candidates,
    )
end

@testset "ToricBuilder case_008 d=15 column-choice diagnostics" begin
    @test isfile(TORICBUILDER_CASE008_D15_MATRIX_BOUNDARY_PATH)
    @test isfile(TORICBUILDER_CASE008_D15_COLUMN_BOUNDARY_PATH)
    include(TORICBUILDER_CASE008_D15_MATRIX_BOUNDARY_PATH)

    fixture = ToricBuilderCase008D15MatrixBoundary.matrix_fixture()
    @test ToricBuilderCase008D15MatrixBoundary.validate_matrix_fixture(fixture) == :ok
    @test (nrows(fixture.failing_input_matrix), ncols(fixture.failing_input_matrix)) == (15, 15)
    @test fixture.current_peel_column_index == 15
    @test fixture.current_peel_column ==
          _case008_d15_column(fixture.failing_input_matrix, fixture.current_peel_column_index)

    column_fixture = ToricBuilderCase008D15ColumnBoundary.boundary_fixture()
    @test ToricBuilderCase008D15ColumnBoundary.validate_boundary_fixture(column_fixture) == :ok
    @test fixture.current_peel_column == column_fixture.failing_column

    report = case008_d15_column_choice_report(fixture)
    @test report.case_id == "case_008"
    @test report.dimension == 15
    @test report.current_peel_column_index == 15
    @test length(report.candidates) == 15
    @test Tuple(candidate.column_index for candidate in report.candidates) == Tuple(1:15)
    @test count(candidate -> candidate.is_current_peel_column, report.candidates) == 1

    required_fields = (
        :column_index,
        :is_current_peel_column,
        :is_unimodular,
        :unit_entry_count,
        :laurent_witness_outcome,
        :laurent_witness_unit_index,
        :normalized_precondition_status,
        :normalized_failure_code,
        :row_preconditioning_outcome,
        :row_preconditioning_transformed_stage,
        :status,
        :failure_code,
        :supported_by_current_reducer,
    )
    @test all(candidate -> all(field -> hasproperty(candidate, field), required_fields), report.candidates)

    current = only(filter(candidate -> candidate.is_current_peel_column, report.candidates))
    @test current.column_index == 15
    @test current.is_unimodular
    @test current.unit_entry_count == 0
    @test current.status == :unsupported
    @test current.failure_code == :unsupported_laurent_column_family
    @test current.supported_by_current_reducer == false

    supported = filter(candidate -> candidate.supported_by_current_reducer, report.candidates)
    supported_alternatives = filter(candidate -> candidate.column_index != 15, supported)
    @info "case_008 d15 column-choice diagnostic" current_peel_column_index = 15 supported_columns = Tuple(candidate.column_index for candidate in supported) supported_alternative_columns = Tuple(candidate.column_index for candidate in supported_alternatives)
    @test all(candidate -> candidate.status == :supported, supported)
    @test all(candidate -> candidate.failure_code === nothing, supported)
    @test all(
        candidate ->
            !candidate.is_unimodular ||
            candidate.supported_by_current_reducer ||
            candidate.failure_code !== nothing,
        report.candidates,
    )

    corrupted = ToricBuilderCase008D15MatrixBoundary.corrupted_matrix_entry_negative_control(fixture)
    rejection = ToricBuilderCase008D15MatrixBoundary.validate_matrix_fixture(corrupted)
    @test rejection in (:wrong_snapshot, :wrong_column)
    @test_throws ArgumentError case008_d15_column_choice_report(corrupted)

    wrong_column = copy(fixture.current_peel_column)
    wrong_column[1] += one(fixture.ring)
    wrong_column_fixture = merge(fixture, (; current_peel_column = wrong_column))
    @test ToricBuilderCase008D15MatrixBoundary.validate_matrix_fixture(wrong_column_fixture) == :wrong_column
    @test_throws ArgumentError case008_d15_column_choice_report(wrong_column_fixture)
end
```

- [ ] **Step 2: Run the focused test to verify RED**

Run:

```bash
julia --project=. -e 'include("test/internal/toricbuilder_case008_d15_column_choice.jl")'
```

Expected: FAIL because `test/fixtures/toricbuilder_case008_d15_matrix_boundary.jl` does not exist yet. The failure should be the first `@test isfile(...)` assertion, not a syntax error.

- [ ] **Step 3: Generate the stored sparse matrix entries**

Run this one-off command and keep the generated `/tmp/case008_d15_sparse_entries.jl` out of git:

```bash
julia --project=. -e '
using Oscar
using Suslin

include("test/fixtures/toricbuilder_case008_d16_matrix_boundary.jl")
include("test/fixtures/toricbuilder_case008_d15_column_boundary.jl")

function _quote_laurent_string(value)
    escaped = replace(string(value), "\\" => "\\\\", "\"" => "\\\"")
    return "\"" * escaped * "\""
end

d16_fixture = ToricBuilderCase008D16MatrixBoundary.matrix_fixture()
ToricBuilderCase008D16MatrixBoundary.validate_matrix_fixture(d16_fixture) == :ok ||
    error("invalid d16 matrix fixture")
step = Suslin._laurent_column_peel_step(d16_fixture.failing_input_matrix)
current = step.next_block
nrows(current) == 15 && ncols(current) == 15 ||
    error("expected one peel step to produce a 15 x 15 matrix, got $(nrows(current)) x $(ncols(current))")

column_fixture = ToricBuilderCase008D15ColumnBoundary.boundary_fixture()
last_column = [current[row, 15] for row in 1:15]
last_column == column_fixture.failing_column ||
    error("derived d15 matrix last column does not match the stored d15 column fixture")

R = base_ring(current)
open("/tmp/case008_d15_sparse_entries.jl", "w") do io
    println(io, "# Generated by replaying one certified peel step from the d=16 matrix boundary fixture.")
    println(io, "const FAILING_INPUT_SPARSE_ENTRIES = (")
    for row in 1:15, col in 1:15
        entry = current[row, col]
        entry == zero(R) && continue
        println(io, "    (", row, ", ", col, ", ", _quote_laurent_string(entry), "),")
    end
    println(io, ")")
end
println("wrote /tmp/case008_d15_sparse_entries.jl")
'
```

Expected: command exits 0 and prints `wrote /tmp/case008_d15_sparse_entries.jl`.

- [ ] **Step 4: Add the matrix fixture module**

Create `test/fixtures/toricbuilder_case008_d15_matrix_boundary.jl` by following `test/fixtures/toricbuilder_case008_d16_matrix_boundary.jl` with these exact d15 changes:

```julia
module ToricBuilderCase008D15MatrixBoundary

using Oscar
using Suslin

include(joinpath(@__DIR__, "toricbuilder_case008_d15_column_boundary.jl"))

const CASE_ID = ToricBuilderCase008D15ColumnBoundary.CASE_ID
const SOURCE_CACHE_FILE = ToricBuilderCase008D15ColumnBoundary.SOURCE_CACHE_FILE
const SOURCE_BLOCK = ToricBuilderCase008D15ColumnBoundary.SOURCE_BLOCK
const SOURCE_MATRIX_DIMENSIONS = ToricBuilderCase008D15ColumnBoundary.SOURCE_MATRIX_DIMENSIONS
const SOURCE_COLUMN_TRANSFORMATION_DIMENSIONS =
    ToricBuilderCase008D15ColumnBoundary.SOURCE_COLUMN_TRANSFORMATION_DIMENSIONS
const EXPECTED_RING_DESCRIPTION =
    ToricBuilderCase008D15ColumnBoundary.EXPECTED_RING_DESCRIPTION
const EXPECTED_PASSED_PEEL_DIMENSIONS =
    ToricBuilderCase008D15ColumnBoundary.EXPECTED_PASSED_PEEL_DIMENSIONS
const FIRST_FAILING_PEEL_DIMENSION =
    ToricBuilderCase008D15ColumnBoundary.FIRST_FAILING_PEEL_DIMENSION
const CURRENT_PEEL_COLUMN_INDEX = FIRST_FAILING_PEEL_DIMENSION
const REQUIRED_MATRIX_FIELDS = (
    :case_id,
    :source_case,
    :source_cache_file,
    :source_block,
    :source_matrix_dimensions,
    :source_column_transformation_dimensions,
    :passed_peel_dimensions,
    :first_failing_peel_dimension,
    :failing_input_matrix,
    :current_peel_column,
    :current_peel_column_index,
    :ring,
    :ring_description,
)
```

Immediately after these constants, paste the complete generated `FAILING_INPUT_SPARSE_ENTRIES` tuple from `/tmp/case008_d15_sparse_entries.jl`.

After the tuple, copy the helper structure from the d16 matrix fixture:

```julia
function _parse_laurent_value(R, u, v, value::AbstractString)
    return ToricBuilderCase008D15ColumnBoundary._parse_laurent_value(R, u, v, value)
end

function _snapshot_matrix()
    R, (u, v) = Suslin.suslin_laurent_polynomial_ring(GF(2), ["u", "v"])
    M = zero_matrix(R, FIRST_FAILING_PEEL_DIMENSION, FIRST_FAILING_PEEL_DIMENSION)
    for (row, col, value) in FAILING_INPUT_SPARSE_ENTRIES
        M[row, col] = _parse_laurent_value(R, u, v, value)
    end
    return M
end

function _matrix_size(M)
    return (nrows(M), ncols(M))
end

function _column(M, column_index::Int)
    return [M[row, column_index] for row in 1:nrows(M)]
end

function _has_required_matrix_fields(fixture)::Bool
    return all(field -> hasproperty(fixture, field), REQUIRED_MATRIX_FIELDS)
end

function _column_fixture()
    return ToricBuilderCase008D15ColumnBoundary.boundary_fixture()
end

function matrix_fixture()
    column_fixture = _column_fixture()
    failing_input_matrix = _snapshot_matrix()
    R = base_ring(failing_input_matrix)
    return (;
        case_id = CASE_ID,
        source_case = CASE_ID,
        source_cache_file = SOURCE_CACHE_FILE,
        source_block = SOURCE_BLOCK,
        source_matrix_dimensions = SOURCE_MATRIX_DIMENSIONS,
        source_column_transformation_dimensions = SOURCE_COLUMN_TRANSFORMATION_DIMENSIONS,
        passed_peel_dimensions = EXPECTED_PASSED_PEEL_DIMENSIONS,
        first_failing_peel_dimension = FIRST_FAILING_PEEL_DIMENSION,
        failing_input_matrix,
        current_peel_column = _column(failing_input_matrix, CURRENT_PEEL_COLUMN_INDEX),
        current_peel_column_index = CURRENT_PEEL_COLUMN_INDEX,
        ring = R,
        ring_description = column_fixture.ring_description,
    )
end

function validate_matrix_fixture(fixture)::Symbol
    _has_required_matrix_fields(fixture) || return :missing_metadata
    fixture.case_id == CASE_ID || return :wrong_case
    fixture.source_case == CASE_ID || return :wrong_case
    fixture.source_cache_file == SOURCE_CACHE_FILE || return :wrong_source
    fixture.source_block == SOURCE_BLOCK || return :wrong_source
    fixture.source_matrix_dimensions == SOURCE_MATRIX_DIMENSIONS || return :wrong_source_matrix_dimensions
    fixture.source_column_transformation_dimensions == SOURCE_COLUMN_TRANSFORMATION_DIMENSIONS ||
        return :wrong_source_matrix_dimensions
    fixture.first_failing_peel_dimension == FIRST_FAILING_PEEL_DIMENSION ||
        return :wrong_peel_dimension
    fixture.passed_peel_dimensions == EXPECTED_PASSED_PEEL_DIMENSIONS ||
        return :wrong_passed_peel_dimensions
    fixture.current_peel_column_index == CURRENT_PEEL_COLUMN_INDEX || return :wrong_column
    _matrix_size(fixture.failing_input_matrix) ==
        (FIRST_FAILING_PEEL_DIMENSION, FIRST_FAILING_PEEL_DIMENSION) ||
        return :wrong_failing_input_dimension
    length(fixture.current_peel_column) == FIRST_FAILING_PEEL_DIMENSION ||
        return :wrong_column_length
    fixture.ring_description == EXPECTED_RING_DESCRIPTION || return :wrong_ring

    R = fixture.ring
    ToricBuilderCase008D15ColumnBoundary._is_expected_uv_laurent_ring(R) ||
        return :wrong_ring
    base_ring(fixture.failing_input_matrix) == R || return :wrong_ring
    fixture.current_peel_column == _column(fixture.failing_input_matrix, CURRENT_PEEL_COLUMN_INDEX) ||
        return :wrong_column

    column_fixture = _column_fixture()
    fixture.current_peel_column == column_fixture.failing_column || return :wrong_column
    fixture.failing_input_matrix == _snapshot_matrix() || return :wrong_snapshot
    return :ok
end

function corrupted_matrix_entry_negative_control(fixture = matrix_fixture())
    corrupted_matrix = copy(fixture.failing_input_matrix)
    corrupted_matrix[1, 1] += one(fixture.ring)
    return merge(fixture, (; failing_input_matrix = corrupted_matrix))
end

end
```

- [ ] **Step 5: Run the focused test to verify GREEN**

Run:

```bash
julia --project=. -e 'include("test/internal/toricbuilder_case008_d15_column_choice.jl")'
```

Expected: PASS. The output may include one `@info` line listing supported columns, and the testset must have zero failures.

- [ ] **Step 6: Run the existing dependent focused tests**

Run:

```bash
julia --project=. -e 'include("test/internal/toricbuilder_case008_d15_column_boundary.jl")'
julia --project=. -e 'include("test/internal/toricbuilder_case008_d16_column_choice.jl")'
```

Expected: both commands exit 0.

- [ ] **Step 7: Commit**

```bash
git add test/fixtures/toricbuilder_case008_d15_matrix_boundary.jl test/internal/toricbuilder_case008_d15_column_choice.jl
git commit -m "test: add case008 d15 matrix column-choice diagnostic"
```

## Plan Self-Review

- Spec coverage: the task creates the required full-matrix fixture, internal report helper, focused diagnostic test, matrix corruption negative control, current-column corruption negative control, and verification commands.
- Placeholder scan: no placeholder markers or deferred implementation steps remain; the generated sparse entries are produced by an exact command and pasted into the fixture.
- Type consistency: fixture and diagnostic names match the d16 matrix pattern with `D15`, `15`, and issue #295 d15 boundary constants.
