# Issue 164 Case008 D16 Column Boundary Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a compact default-tested `case_008` `d=16` failing-column fixture.

**Architecture:** Add one standalone fixture module that parses stored Laurent expression strings into a length-16 column and validates metadata, unimodularity, and the current unsupported diagnostic. Add one internal test that pins the fixture contract and register it in the default internal suite.

**Tech Stack:** Julia, Oscar, Suslin, Test stdlib, existing ToricBuilder fixture/test conventions.

## Global Constraints

- Do not implement a new reducer stage.
- Do not require `case_008` to pass.
- Do not store the full `d=16` matrix.
- Do not add the slow `30x30` derivation path to default tests.
- The fixture must be column-only for issue #164.
- The focused test must assert `case_id == "case_008"`.
- The focused test must assert `source_matrix_dimensions == (30, 30)`.
- The focused test must assert `first_failing_peel_dimension == 16`.
- The focused test must assert `passed_peel_dimensions == (30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17)`.
- The focused test must assert `length(failing_column) == 16`.
- The focused test must assert `is_unimodular_column(failing_column, ring) == true`.
- The focused test must assert `count(is_unit, failing_column) == 0`.
- The focused test must assert `diagnose_unimodular_column_reduction(failing_column, ring).failure_code == :unsupported_laurent_column_family`.
- The negative control must corrupt the fixture by multiplying the failing column by a nonunit such as `v + 1`, and `validate_boundary_fixture(corrupted)` must explicitly reject it rather than pass.

---

## File Structure

- `test/internal/toricbuilder_case008_d16_column_boundary.jl`: new focused internal validator test.
- `test/fixtures/toricbuilder_case008_d16_column_boundary.jl`: new column-only fixture module and validator helpers.
- `test/runtests.jl`: register the focused internal test in the default internal group.

---

### Task 1: Add The Failing Internal Contract Test

**Files:**
- Create: `test/internal/toricbuilder_case008_d16_column_boundary.jl`

**Interfaces:**
- Consumes: `ToricBuilderCase008D16ColumnBoundary.boundary_fixture()`, `ToricBuilderCase008D16ColumnBoundary.validate_boundary_fixture(fixture)`, `ToricBuilderCase008D16ColumnBoundary.non_unimodular_negative_control(fixture)`.
- Produces: a failing test that pins the issue #164 fixture contract.

- [ ] **Step 1: Write the failing test**

Create `test/internal/toricbuilder_case008_d16_column_boundary.jl`:

```julia
using Test
using Suslin
using Oscar

const TORICBUILDER_CASE008_D16_COLUMN_BOUNDARY_PATH =
    joinpath(@__DIR__, "..", "fixtures", "toricbuilder_case008_d16_column_boundary.jl")

@testset "ToricBuilder case_008 d=16 Laurent column boundary" begin
    @test isfile(TORICBUILDER_CASE008_D16_COLUMN_BOUNDARY_PATH)

    include(TORICBUILDER_CASE008_D16_COLUMN_BOUNDARY_PATH)
    fixture = ToricBuilderCase008D16ColumnBoundary.boundary_fixture()

    @test fixture.case_id == "case_008"
    @test fixture.source_matrix_dimensions == (30, 30)
    @test fixture.source_column_transformation_dimensions == (60, 60)
    @test fixture.first_failing_peel_dimension == 16
    @test fixture.passed_peel_dimensions == (30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17)
    @test length(fixture.failing_column) == 16
    @test fixture.ring_description == "GF(2)[u^+/-1, v^+/-1]"
    @test Tuple(string.(gens(fixture.ring))) == ("u", "v")
    @test Suslin.is_unimodular_column(fixture.failing_column, fixture.ring)
    @test count(is_unit, fixture.failing_column) == 0

    diagnostic = Suslin.diagnose_unimodular_column_reduction(
        fixture.failing_column,
        fixture.ring,
    )
    @test diagnostic.status == :unsupported
    @test diagnostic.failure_code == :unsupported_laurent_column_family
    @test diagnostic.column_length == 16
    @test diagnostic.ring_profile.kind == :laurent_polynomial
    @test diagnostic.ring_profile.generators == ("u", "v")

    @test ToricBuilderCase008D16ColumnBoundary.validate_boundary_fixture(fixture) == :ok

    corrupted = ToricBuilderCase008D16ColumnBoundary.non_unimodular_negative_control(fixture)
    @test !Suslin.is_unimodular_column(corrupted.failing_column, corrupted.ring)
    @test ToricBuilderCase008D16ColumnBoundary.validate_boundary_fixture(corrupted) == :not_unimodular

    wrong_dimension = merge(fixture, (; first_failing_peel_dimension = 15))
    @test ToricBuilderCase008D16ColumnBoundary.validate_boundary_fixture(wrong_dimension) == :wrong_peel_dimension

    wrong_passed = merge(fixture, (; passed_peel_dimensions = (30, 29)))
    @test ToricBuilderCase008D16ColumnBoundary.validate_boundary_fixture(wrong_passed) == :wrong_passed_peel_dimensions

    wrong_ring, (wrong_u, wrong_v) = Suslin.suslin_laurent_polynomial_ring(GF(3), ["u", "v"])
    wrong_ring_fixture = merge(
        fixture,
        (;
            ring = wrong_ring,
            failing_column = [idx == 1 ? wrong_u : idx == 2 ? wrong_v : zero(wrong_ring) for idx in 1:16],
        ),
    )
    @test ToricBuilderCase008D16ColumnBoundary.validate_boundary_fixture(wrong_ring_fixture) == :wrong_ring
end
```

- [ ] **Step 2: Run the focused test to verify RED**

Run:

```bash
julia --project=. -e 'include("test/internal/toricbuilder_case008_d16_column_boundary.jl")'
```

Expected: FAIL because `test/fixtures/toricbuilder_case008_d16_column_boundary.jl` does not exist yet.

---

### Task 2: Add The Column-Only Fixture Module

**Files:**
- Create: `test/fixtures/toricbuilder_case008_d16_column_boundary.jl`
- Test: `test/internal/toricbuilder_case008_d16_column_boundary.jl`

**Interfaces:**
- Consumes: the issue #164 d=16 column values derived from the current slow peel path.
- Produces: `boundary_fixture()`, `validate_boundary_fixture(fixture)::Symbol`, and `non_unimodular_negative_control(fixture)`.

- [ ] **Step 1: Materialize the generated column constants**

Use the already-derived constants from `/tmp/case008_d16_column_data.jl` if present. If it is absent, regenerate it with a one-off local derivation command; do not add the derivation command to any default test.

The generated data must define:

```julia
const EXPECTED_PASSED_PEEL_DIMENSIONS = (30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17)
const FIRST_FAILING_PEEL_DIMENSION = 16
const FAILING_COLUMN_ENTRIES = (
    # exactly 16 Laurent expression strings
)
```

- [ ] **Step 2: Create the fixture module**

Create `test/fixtures/toricbuilder_case008_d16_column_boundary.jl` with:

```julia
module ToricBuilderCase008D16ColumnBoundary

using Oscar
using Suslin

const CASE_ID = "case_008"
const SOURCE_CACHE_FILE = "case_008.jls"
const SOURCE_BLOCK = :column_transformation_upper_left_q_block
const SOURCE_MATRIX_DIMENSIONS = (30, 30)
const SOURCE_COLUMN_TRANSFORMATION_DIMENSIONS = (60, 60)
const EXPECTED_RING_DESCRIPTION = "GF(2)[u^+/-1, v^+/-1]"
const EXPECTED_DIAGNOSTIC = (;
    status = :unsupported,
    failure_code = :unsupported_laurent_column_family,
)
const REQUIRED_BOUNDARY_FIELDS = (
    :case_id,
    :source_case,
    :source_cache_file,
    :source_block,
    :source_matrix_dimensions,
    :source_column_transformation_dimensions,
    :passed_peel_dimensions,
    :first_failing_peel_dimension,
    :failing_column,
    :ring,
    :ring_description,
    :expected_diagnostic,
)

# The generated data file supplies EXPECTED_PASSED_PEEL_DIMENSIONS,
# FIRST_FAILING_PEEL_DIMENSION, and FAILING_COLUMN_ENTRIES verbatim before
# these helper functions.

function _integer_exponent(ex)
    ex isa Integer && return Int(ex)
    if ex isa Expr && ex.head == :call && ex.args[1] == :- && length(ex.args) == 2 && ex.args[2] isa Integer
        return -Int(ex.args[2])
    end
    throw(ArgumentError("unsupported Laurent exponent expression: $(ex)"))
end

function _eval_laurent_expr(ex, R, variables)
    if ex isa Integer
        return ex == 0 ? zero(R) : ex == 1 ? one(R) : R(ex)
    elseif ex isa Symbol
        haskey(variables, ex) || throw(ArgumentError("unknown Laurent variable: $(ex)"))
        return variables[ex]
    elseif ex isa Expr && ex.head == :call
        op = ex.args[1]
        if op == :+
            value = zero(R)
            for arg in ex.args[2:end]
                value += _eval_laurent_expr(arg, R, variables)
            end
            return value
        elseif op == :*
            value = one(R)
            for arg in ex.args[2:end]
                value *= _eval_laurent_expr(arg, R, variables)
            end
            return value
        elseif op == :^ && length(ex.args) == 3
            return _eval_laurent_expr(ex.args[2], R, variables)^_integer_exponent(ex.args[3])
        elseif op == :- && length(ex.args) == 2
            return -_eval_laurent_expr(ex.args[2], R, variables)
        end
    end
    throw(ArgumentError("unsupported Laurent expression: $(ex)"))
end

function _parse_laurent_value(R, u, v, value::AbstractString)
    return _eval_laurent_expr(Meta.parse(value), R, Dict(:u => u, :v => v))
end

function _ring_generator_names(R)
    return Tuple(string.(gens(R)))
end

function _is_expected_uv_laurent_ring(R)::Bool
    try
        return Suslin._is_laurent_polynomial_ring(R) &&
            base_ring(R) == GF(2) &&
            _ring_generator_names(R) == ("u", "v")
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end

function _column_entries_are_over_ring(column, R)::Bool
    try
        return all(entry -> R(entry) == entry, column)
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end

function _failing_column(R, u, v)
    return [_parse_laurent_value(R, u, v, value) for value in FAILING_COLUMN_ENTRIES]
end

function boundary_fixture()
    R, (u, v) = Suslin.suslin_laurent_polynomial_ring(GF(2), ["u", "v"])
    return (;
        case_id = CASE_ID,
        source_case = CASE_ID,
        source_cache_file = SOURCE_CACHE_FILE,
        source_block = SOURCE_BLOCK,
        source_matrix_dimensions = SOURCE_MATRIX_DIMENSIONS,
        source_column_transformation_dimensions = SOURCE_COLUMN_TRANSFORMATION_DIMENSIONS,
        passed_peel_dimensions = EXPECTED_PASSED_PEEL_DIMENSIONS,
        first_failing_peel_dimension = FIRST_FAILING_PEEL_DIMENSION,
        failing_column = _failing_column(R, u, v),
        ring = R,
        ring_description = EXPECTED_RING_DESCRIPTION,
        expected_diagnostic = EXPECTED_DIAGNOSTIC,
    )
end

function _has_required_boundary_fields(fixture)::Bool
    return all(field -> hasproperty(fixture, field), REQUIRED_BOUNDARY_FIELDS)
end

function _diagnostic_matches_expected(column, R, expected)::Bool
    diagnostic = Suslin.diagnose_unimodular_column_reduction(column, R)
    return diagnostic.status == expected.status &&
        diagnostic.failure_code == expected.failure_code &&
        diagnostic.column_length == FIRST_FAILING_PEEL_DIMENSION
end

function validate_boundary_fixture(fixture)::Symbol
    _has_required_boundary_fields(fixture) || return :missing_metadata
    fixture.case_id == CASE_ID || return :wrong_case
    fixture.source_case == CASE_ID || return :wrong_case
    fixture.source_matrix_dimensions == SOURCE_MATRIX_DIMENSIONS || return :wrong_source_matrix_dimensions
    fixture.source_column_transformation_dimensions == SOURCE_COLUMN_TRANSFORMATION_DIMENSIONS ||
        return :wrong_source_matrix_dimensions
    fixture.first_failing_peel_dimension == FIRST_FAILING_PEEL_DIMENSION || return :wrong_peel_dimension
    fixture.passed_peel_dimensions == EXPECTED_PASSED_PEEL_DIMENSIONS || return :wrong_passed_peel_dimensions
    length(fixture.failing_column) == FIRST_FAILING_PEEL_DIMENSION || return :wrong_column_length
    fixture.ring_description == EXPECTED_RING_DESCRIPTION || return :wrong_ring

    R = fixture.ring
    _is_expected_uv_laurent_ring(R) || return :wrong_ring
    _column_entries_are_over_ring(fixture.failing_column, R) || return :wrong_ring
    Suslin.is_unimodular_column(fixture.failing_column, R) || return :not_unimodular
    count(is_unit, fixture.failing_column) == 0 || return :wrong_unit_profile
    _diagnostic_matches_expected(fixture.failing_column, R, fixture.expected_diagnostic) ||
        return :wrong_diagnostic

    u, v = gens(R)
    fixture.failing_column == _failing_column(R, u, v) || return :wrong_column

    return :ok
end

function non_unimodular_negative_control(fixture = boundary_fixture())
    _, v = gens(fixture.ring)
    nonunit = v + one(fixture.ring)
    return merge(fixture, (; failing_column = [nonunit * entry for entry in fixture.failing_column]))
end

end
```

- [ ] **Step 3: Run the focused test to verify GREEN**

Run:

```bash
julia --project=. -e 'include("test/internal/toricbuilder_case008_d16_column_boundary.jl")'
```

Expected: PASS.

---

### Task 3: Register The Internal Test

**Files:**
- Modify: `test/runtests.jl`
- Test: `test/runtests.jl`

**Interfaces:**
- Consumes: `test/internal/toricbuilder_case008_d16_column_boundary.jl`.
- Produces: default internal-suite coverage for the new fixture.

- [ ] **Step 1: Register the new internal test**

In `test/runtests.jl`, add:

```julia
"internal/toricbuilder_case008_d16_column_boundary.jl",
```

after:

```julia
"internal/toricbuilder_case008_d21_column_boundary.jl",
```

- [ ] **Step 2: Run default tests**

Run:

```bash
julia --project=. test/runtests.jl
```

Expected: PASS.

- [ ] **Step 3: Run package tests**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: PASS.
