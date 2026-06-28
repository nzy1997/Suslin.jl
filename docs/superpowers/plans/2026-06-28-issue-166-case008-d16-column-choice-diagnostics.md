# Issue 166 Case008 D16 Column-Choice Diagnostics Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an offline full `case_008` `d=16` matrix fixture and a focused diagnostic that evaluates all 16 columns as possible peel columns.

**Architecture:** Keep the existing default d16 column fixture compact. Add a separate matrix fixture that stores sparse Laurent string entries and materializes the matrix only for the focused diagnostic test. The diagnostic helper validates the fixture before reporting and uses the existing `diagnose_unimodular_column_reduction` stage details.

**Tech Stack:** Julia, Oscar.jl Laurent rings, existing Suslin internal diagnostics, `Test`.

## Global Constraints

- Repository: `nzy1997/Suslin.jl`.
- Issue: `#166 Add full case008 d16 matrix column-choice diagnostics`.
- Base branch: `main`.
- Worker branch: `agent/issue-166-add-full-case008-d16-matrix-column-choice-diagno-run-1`.
- Do not change peel ordering.
- Do not implement a new reducer stage.
- Do not add the full matrix to default slow paths.
- Do not register `test/internal/toricbuilder_case008_d16_column_choice.jl` in `test/runtests.jl`.
- The focused verification command is `julia --project=. -e 'include("test/internal/toricbuilder_case008_d16_column_choice.jl")'`.
- The diagnostic must check exactly 16 candidate columns.
- The diagnostic must identify current last column index `16`.
- If no candidate is supported, every candidate must be asserted as either non-unimodular or unsupported with an explicit failure code.
- Negative control: corrupt one stored matrix entry in a copied fixture and require validator rejection before producing a column-choice report.

---

## File Structure

- Create `test/fixtures/toricbuilder_case008_d16_matrix_boundary.jl`: matrix-only fixture module, stored sparse string entries, materialization, validation, and corrupted-copy helper.
- Create `test/internal/toricbuilder_case008_d16_column_choice.jl`: focused internal diagnostic and tests. This file is intentionally not registered in `test/runtests.jl`.
- Keep `test/fixtures/toricbuilder_case008_d16_column_boundary.jl` unchanged except for using it from the new matrix fixture.

---

### Task 1: Add The Failing Column-Choice Diagnostic Test

**Files:**
- Create: `test/internal/toricbuilder_case008_d16_column_choice.jl`

**Interfaces:**
- Consumes: planned `ToricBuilderCase008D16MatrixBoundary.matrix_fixture()`, `validate_matrix_fixture(fixture)`, and `corrupted_matrix_entry_negative_control(fixture)`.
- Produces: `case008_d16_column_choice_report(fixture = ToricBuilderCase008D16MatrixBoundary.matrix_fixture())` returning a named tuple with `case_id`, `dimension`, `current_peel_column_index`, and `candidates`.

- [ ] **Step 1: Write the failing test and diagnostic helper**

Create `test/internal/toricbuilder_case008_d16_column_choice.jl` with this content:

```julia
using Test
using Suslin
using Oscar

const TORICBUILDER_CASE008_D16_MATRIX_BOUNDARY_PATH =
    joinpath(@__DIR__, "..", "fixtures", "toricbuilder_case008_d16_matrix_boundary.jl")

function _case008_d16_stage_detail(diagnostic, stage::Symbol)
    hasproperty(diagnostic, :stage_details) || return nothing
    idx = findfirst(detail -> detail.stage == stage, diagnostic.stage_details)
    return idx === nothing ? nothing : diagnostic.stage_details[idx]
end

function _case008_d16_column(M, column_index::Int)
    return [M[row, column_index] for row in 1:nrows(M)]
end

function _case008_d16_detail_property(detail, field::Symbol)
    detail === nothing && return nothing
    return hasproperty(detail, field) ? getproperty(detail, field) : nothing
end

function _case008_d16_candidate_report(M, R, column_index::Int, current_index::Int)
    column = _case008_d16_column(M, column_index)
    unit_entry_count = count(is_unit, column)
    is_unimodular = Suslin.is_unimodular_column(column, R)
    diagnostic = Suslin.diagnose_unimodular_column_reduction(column, R)
    witness = _case008_d16_stage_detail(diagnostic, :laurent_witness_unit)
    normalization = _case008_d16_stage_detail(diagnostic, :laurent_normalization)
    return (;
        column_index,
        is_current_peel_column = column_index == current_index,
        is_unimodular,
        unit_entry_count,
        laurent_witness_outcome = _case008_d16_detail_property(witness, :outcome),
        laurent_witness_unit_index = _case008_d16_detail_property(witness, :witness_unit_index),
        normalized_precondition_status = _case008_d16_detail_property(normalization, :normalized_status),
        normalized_failure_code = _case008_d16_detail_property(normalization, :normalized_failure_code),
        status = diagnostic.status,
        failure_code = diagnostic.failure_code,
        supported_by_current_reducer = diagnostic.status == :supported,
    )
end

function case008_d16_column_choice_report(fixture = ToricBuilderCase008D16MatrixBoundary.matrix_fixture())
    validation = ToricBuilderCase008D16MatrixBoundary.validate_matrix_fixture(fixture)
    validation == :ok ||
        throw(ArgumentError("invalid case_008 d16 matrix fixture: $(validation)"))

    M = fixture.failing_input_matrix
    R = fixture.ring
    d = fixture.first_failing_peel_dimension
    candidates = tuple(
        (
            _case008_d16_candidate_report(
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

@testset "ToricBuilder case_008 d=16 column-choice diagnostics" begin
    @test isfile(TORICBUILDER_CASE008_D16_MATRIX_BOUNDARY_PATH)
    include(TORICBUILDER_CASE008_D16_MATRIX_BOUNDARY_PATH)

    fixture = ToricBuilderCase008D16MatrixBoundary.matrix_fixture()
    @test ToricBuilderCase008D16MatrixBoundary.validate_matrix_fixture(fixture) == :ok
    @test (nrows(fixture.failing_input_matrix), ncols(fixture.failing_input_matrix)) == (16, 16)
    @test fixture.current_peel_column_index == 16
    @test fixture.current_peel_column ==
          _case008_d16_column(fixture.failing_input_matrix, fixture.current_peel_column_index)

    report = case008_d16_column_choice_report(fixture)
    @test report.case_id == "case_008"
    @test report.dimension == 16
    @test report.current_peel_column_index == 16
    @test length(report.candidates) == 16
    @test Tuple(candidate.column_index for candidate in report.candidates) == Tuple(1:16)
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
        :status,
        :failure_code,
        :supported_by_current_reducer,
    )
    @test all(candidate -> all(field -> hasproperty(candidate, field), required_fields), report.candidates)

    current = only(filter(candidate -> candidate.is_current_peel_column, report.candidates))
    @test current.column_index == 16
    @test current.is_unimodular
    @test current.unit_entry_count == 0
    @test current.status == :unsupported
    @test current.failure_code == :unsupported_laurent_column_family
    @test current.laurent_witness_outcome == :witness_without_unit
    @test current.normalized_precondition_status == :precondition_failed
    @test current.normalized_failure_code == :not_unimodular
    @test !current.supported_by_current_reducer

    supported = filter(candidate -> candidate.supported_by_current_reducer, report.candidates)
    supported_alternatives = filter(candidate -> candidate.column_index != 16, supported)
    @info "case_008 d16 column-choice diagnostic" current_peel_column_index = 16 supported_columns = Tuple(candidate.column_index for candidate in supported) supported_alternative_columns = Tuple(candidate.column_index for candidate in supported_alternatives)
    @test all(candidate -> candidate.status == :supported, supported)

    if isempty(supported)
        @test all(
            candidate ->
                !candidate.is_unimodular ||
                (!candidate.supported_by_current_reducer && candidate.failure_code !== nothing),
            report.candidates,
        )
    else
        @test all(candidate -> candidate.failure_code === nothing, supported)
    end

    corrupted = ToricBuilderCase008D16MatrixBoundary.corrupted_matrix_entry_negative_control(fixture)
    rejection = ToricBuilderCase008D16MatrixBoundary.validate_matrix_fixture(corrupted)
    @test rejection in (:wrong_snapshot, :wrong_column)
    @test_throws ArgumentError case008_d16_column_choice_report(corrupted)
end
```

- [ ] **Step 2: Run the focused test to verify RED**

Run:

```bash
julia --project=. -e 'include("test/internal/toricbuilder_case008_d16_column_choice.jl")'
```

Expected: FAIL because `test/fixtures/toricbuilder_case008_d16_matrix_boundary.jl` does not exist yet. The failure should be the `@test isfile(...)` assertion, not a syntax error.

- [ ] **Step 3: Commit the red test**

```bash
git add test/internal/toricbuilder_case008_d16_column_choice.jl
git commit -m "test: add case008 d16 column-choice diagnostic"
```

---

### Task 2: Add The Full D16 Matrix Fixture And Validator

**Files:**
- Create: `test/fixtures/toricbuilder_case008_d16_matrix_boundary.jl`
- Test: `test/internal/toricbuilder_case008_d16_column_choice.jl`

**Interfaces:**
- Consumes: `ToricBuilderCase008D16ColumnBoundary.boundary_fixture()` from `test/fixtures/toricbuilder_case008_d16_column_boundary.jl`.
- Produces: `ToricBuilderCase008D16MatrixBoundary.matrix_fixture()`, `validate_matrix_fixture(fixture)::Symbol`, and `corrupted_matrix_entry_negative_control(fixture)`.

- [ ] **Step 1: Generate the stored sparse matrix entries**

Run this one-off command and keep the generated `/tmp/case008_d16_sparse_entries.jl` out of git:

```bash
julia --project=. -e '
using Oscar
using Suslin

include("test/fixtures/toricbuilder_case008_d21_column_boundary.jl")
include("test/fixtures/toricbuilder_case008_d16_column_boundary.jl")

function _quote_laurent_string(value)
    escaped = replace(string(value), "\\" => "\\\\", "\"" => "\\\"")
    return "\"" * escaped * "\""
end

fixture = ToricBuilderCase008D21ColumnBoundary.boundary_fixture()
current = fixture.failing_input_matrix
passed = Int[]
while nrows(current) > 16
    push!(passed, nrows(current))
    step = Suslin._laurent_column_peel_step(current)
    global current = step.next_block
end

column_fixture = ToricBuilderCase008D16ColumnBoundary.boundary_fixture()
last_column = [current[row, 16] for row in 1:16]
last_column == column_fixture.failing_column ||
    error("derived d16 matrix last column does not match the stored d16 column fixture")
Tuple(passed) == (21, 20, 19, 18, 17) ||
    error("unexpected d16 peel path: $(Tuple(passed))")

R = base_ring(current)
open("/tmp/case008_d16_sparse_entries.jl", "w") do io
    println(io, "# Generated by replaying case_008 peel dimensions $(Tuple(passed)) from the d=21 fixture.")
    println(io, "const FAILING_INPUT_SPARSE_ENTRIES = (")
    for row in 1:16, col in 1:16
        entry = current[row, col]
        entry == zero(R) && continue
        println(io, "    (", row, ", ", col, ", ", _quote_laurent_string(entry), "),")
    end
    println(io, ")")
end

println("wrote /tmp/case008_d16_sparse_entries.jl")
println("nnz=", count(idx -> current[idx[1], idx[2]] != zero(R), Iterators.product(1:16, 1:16)))
'
```

Expected: exits 0, prints `nnz=256`, and writes 256 sparse entries.

- [ ] **Step 2: Add the matrix fixture module**

Create `test/fixtures/toricbuilder_case008_d16_matrix_boundary.jl` with this structure, inserting the generated `FAILING_INPUT_SPARSE_ENTRIES` constant from `/tmp/case008_d16_sparse_entries.jl` where indicated:

```julia
module ToricBuilderCase008D16MatrixBoundary

using Oscar
using Suslin

include(joinpath(@__DIR__, "toricbuilder_case008_d16_column_boundary.jl"))

const CASE_ID = ToricBuilderCase008D16ColumnBoundary.CASE_ID
const SOURCE_CACHE_FILE = ToricBuilderCase008D16ColumnBoundary.SOURCE_CACHE_FILE
const SOURCE_BLOCK = ToricBuilderCase008D16ColumnBoundary.SOURCE_BLOCK
const SOURCE_MATRIX_DIMENSIONS = ToricBuilderCase008D16ColumnBoundary.SOURCE_MATRIX_DIMENSIONS
const SOURCE_COLUMN_TRANSFORMATION_DIMENSIONS =
    ToricBuilderCase008D16ColumnBoundary.SOURCE_COLUMN_TRANSFORMATION_DIMENSIONS
const EXPECTED_RING_DESCRIPTION =
    ToricBuilderCase008D16ColumnBoundary.EXPECTED_RING_DESCRIPTION
const EXPECTED_PASSED_PEEL_DIMENSIONS =
    ToricBuilderCase008D16ColumnBoundary.EXPECTED_PASSED_PEEL_DIMENSIONS
const FIRST_FAILING_PEEL_DIMENSION =
    ToricBuilderCase008D16ColumnBoundary.FIRST_FAILING_PEEL_DIMENSION
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

# Paste generated FAILING_INPUT_SPARSE_ENTRIES here.

function _parse_laurent_value(R, u, v, value::AbstractString)
    return ToricBuilderCase008D16ColumnBoundary._parse_laurent_value(R, u, v, value)
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
    return ToricBuilderCase008D16ColumnBoundary.boundary_fixture()
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
    ToricBuilderCase008D16ColumnBoundary._is_expected_uv_laurent_ring(R) ||
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

- [ ] **Step 3: Run the focused test to verify GREEN**

Run:

```bash
julia --project=. -e 'include("test/internal/toricbuilder_case008_d16_column_choice.jl")'
```

Expected: PASS. The test output includes an `@info` line with `current_peel_column_index = 16`, `supported_columns`, and `supported_alternative_columns`.

- [ ] **Step 4: Run the existing d16 boundary test**

Run:

```bash
julia --project=. -e 'include("test/internal/toricbuilder_case008_d16_column_boundary.jl")'
```

Expected: PASS, proving the existing compact fixture remains valid.

- [ ] **Step 5: Commit the matrix fixture**

```bash
git add test/fixtures/toricbuilder_case008_d16_matrix_boundary.jl
git commit -m "test: add case008 d16 matrix fixture"
```

---

### Task 3: Final Verification

**Files:**
- No code edits expected.

**Interfaces:**
- Consumes: focused diagnostic and fixture from Tasks 1 and 2.
- Produces: final verification evidence for the PR.

- [ ] **Step 1: Run issue verification**

Run:

```bash
julia --project=. -e 'include("test/internal/toricbuilder_case008_d16_column_choice.jl")'
```

Expected: PASS and exactly 16 candidate columns checked.

- [ ] **Step 2: Run required full package verification**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: PASS.

- [ ] **Step 3: Inspect diff scope**

Run:

```bash
git status --short
git diff --stat origin/main..HEAD
```

Expected: only the issue #166 design/plan docs, the new focused diagnostic test, and the new matrix fixture are in scope.

---

## Self-Review

- Spec coverage: Task 1 covers the 16-column diagnostic, current column identity, explicit no-supported-candidate proof, and negative-control report rejection. Task 2 covers the stored matrix fixture, validator, copied-entry corruption, and same-peel-path derivation. Task 3 covers required verification.
- Placeholder scan: no unfinished markers or unspecified public API additions remain.
- Type consistency: the plan consistently uses `matrix_fixture()`, `validate_matrix_fixture(fixture)::Symbol`, `corrupted_matrix_entry_negative_control(fixture)`, and `case008_d16_column_choice_report(...)`.
- Scope: the new column-choice diagnostic is focused and intentionally not registered in `test/runtests.jl`, preserving default slow paths.

## Execution Handoff

Plan complete and saved to `docs/superpowers/plans/2026-06-28-issue-166-case008-d16-column-choice-diagnostics.md`. Two execution options:

1. **Subagent-Driven (recommended)** - dispatch a fresh subagent per task, review between tasks, fast iteration.
2. **Inline Execution** - execute tasks in this session using executing-plans, batch execution with checkpoints.

Automatic choice: Subagent-Driven, because it is marked recommended and subagent tools are available in this environment.
