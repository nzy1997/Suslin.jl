# Issue 132 Case 010 Laurent Column Boundary Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a focused offline fixture and validator for the exact normalized length-5 Laurent column that currently blocks `case_010`.

**Architecture:** Keep all new behavior in test support. A fixture module materializes `case_010`, normalizes it, follows the existing Laurent column-peel step path until the first unsupported `reduce_unimodular_column` boundary, and exposes validator helpers for the internal test.

**Tech Stack:** Julia, Oscar, Suslin internal test fixtures, Test stdlib.

## Global Constraints

- Source fixture must be `test/fixtures/toricbuilder_cache_q_blocks.jl`.
- New fixture path must be `test/fixtures/toricbuilder_case010_column_boundary.jl`.
- New issue verification command must be `julia --project=. -e 'include("test/internal/toricbuilder_case010_column_boundary.jl")'`.
- Full verification command must be `julia --project=. -e 'using Pkg; Pkg.test()'`.
- The extracted boundary column must have length `5`.
- The extracted boundary column must be over `GF(2)[u^+/-1, v^+/-1]`.
- The extracted boundary column must be unimodular.
- `Suslin.reduce_unimodular_column(column, R)` must still throw an `ArgumentError` containing `unsupported exact unimodular column reduction`.
- The validator must return `:not_unimodular` for a non-unimodular negative-control column.
- Do not fix `case_010`.
- Do not change ToricBuilder status-report statuses.
- Do not add public Suslin APIs for this issue.
- A one-entry zero perturbation of the extracted length-5 column remains unimodular in the current algebra; do not claim it is a non-unimodular negative control.

---

### Task 1: Add Case 010 Boundary Fixture And Internal Test

**Files:**
- Modify: `Project.toml`
- Create: `test/fixtures/toricbuilder_case010_column_boundary.jl`
- Create: `test/internal/toricbuilder_case010_column_boundary.jl`
- Modify: `test/runtests.jl`

**Interfaces:**
- Consumes: `ToricBuilderCacheQBlocks.catalog()`, `ToricBuilderCacheQBlocks.materialize_matrix(entry)`, `Suslin.normalize_laurent_gl_matrix(A)`, `Suslin._laurent_column_peel_step(current)`, `Suslin.is_unimodular_column(column, R)`, `Suslin.reduce_unimodular_column(column, R)`.
- Produces: `ToricBuilderCase010ColumnBoundary.boundary_fixture()`, `ToricBuilderCase010ColumnBoundary.validate_boundary_fixture(fixture)`, `ToricBuilderCase010ColumnBoundary.non_unimodular_negative_control(fixture)`, and `ToricBuilderCase010ColumnBoundary.single_entry_zero_perturbations(fixture)`.

- [ ] **Step 1: Write the failing internal test**

Create `test/internal/toricbuilder_case010_column_boundary.jl` with this content:

```julia
using Test
using Suslin
using Oscar

const TORICBUILDER_CASE010_COLUMN_BOUNDARY_PATH =
    joinpath(@__DIR__, "..", "fixtures", "toricbuilder_case010_column_boundary.jl")

@testset "ToricBuilder case_010 Laurent column boundary" begin
    @test isfile(TORICBUILDER_CASE010_COLUMN_BOUNDARY_PATH)

    include(TORICBUILDER_CASE010_COLUMN_BOUNDARY_PATH)
    fixture = ToricBuilderCase010ColumnBoundary.boundary_fixture()

    @test fixture.case_id == "case_010"
    @test (nrows(fixture.original_matrix), ncols(fixture.original_matrix)) == (6, 6)
    @test Suslin.verify_laurent_gl_normalization(fixture.original_matrix, fixture.normalization)
    @test det(fixture.normalized_matrix) == one(fixture.ring)
    @test fixture.first_failing_peel_dimension == 5
    @test length(fixture.failing_column) == 5
    @test fixture.ring_description == "GF(2)[u^+/-1, v^+/-1]"
    @test Tuple(string.(gens(fixture.ring))) == ("u", "v")
    @test Suslin.is_unimodular_column(fixture.failing_column, fixture.ring)
    @test ToricBuilderCase010ColumnBoundary.validate_boundary_fixture(fixture) == :ok

    err = try
        Suslin.reduce_unimodular_column(fixture.failing_column, fixture.ring)
        nothing
    catch caught
        caught
    end
    @test err isa ArgumentError
    @test occursin(fixture.expected_diagnostic, sprint(showerror, err))

    negative = ToricBuilderCase010ColumnBoundary.non_unimodular_negative_control(fixture)
    @test !Suslin.is_unimodular_column(negative.failing_column, negative.ring)
    @test ToricBuilderCase010ColumnBoundary.validate_boundary_fixture(negative) == :not_unimodular

    for perturbed in ToricBuilderCase010ColumnBoundary.single_entry_zero_perturbations(fixture)
        @test Suslin.is_unimodular_column(perturbed, fixture.ring)
    end
end
```

- [ ] **Step 2: Run the test to verify RED**

Run:

```bash
julia --project=. -e 'include("test/internal/toricbuilder_case010_column_boundary.jl")'
```

Expected: FAIL because `test/fixtures/toricbuilder_case010_column_boundary.jl` does not exist and the `isfile` assertion is false.

- [ ] **Step 3: Add the fixture implementation**

Create `test/fixtures/toricbuilder_case010_column_boundary.jl` with this content:

```julia
module ToricBuilderCase010ColumnBoundary

using Oscar
using Suslin

include(joinpath(@__DIR__, "toricbuilder_cache_q_blocks.jl"))

const CASE_ID = "case_010"
const EXPECTED_RING_DESCRIPTION = "GF(2)[u^+/-1, v^+/-1]"
const EXPECTED_DIAGNOSTIC = "unsupported exact unimodular column reduction"
const REQUIRED_BOUNDARY_FIELDS = (
    :case_id,
    :original_matrix,
    :normalization,
    :normalized_matrix,
    :first_failing_peel_dimension,
    :failing_input_matrix,
    :failing_column,
    :ring,
    :ring_description,
    :expected_diagnostic,
)

function _case010_entry()
    matches = filter(entry -> entry.id == CASE_ID, ToricBuilderCacheQBlocks.catalog().cases)
    length(matches) == 1 ||
        throw(ArgumentError("expected exactly one ToricBuilder cache Q-block entry for $(CASE_ID)"))
    return only(matches)
end

function _ring_generator_names(R)
    return Tuple(string.(gens(R)))
end

function _is_expected_uv_laurent_ring(R)::Bool
    try
        return Suslin._is_laurent_polynomial_ring(R) && _ring_generator_names(R) == ("u", "v")
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

function _first_failing_boundary(normalized_matrix)
    current = normalized_matrix
    passed_dimensions = Int[]

    while nrows(current) > 2
        d = nrows(current)
        try
            step = Suslin._laurent_column_peel_step(current)
            push!(passed_dimensions, d)
            current = step.next_block
        catch err
            err isa InterruptException && rethrow()
            message = sprint(showerror, err)
            (err isa ArgumentError && occursin(EXPECTED_DIAGNOSTIC, message)) || rethrow()
            return (;
                first_failing_peel_dimension = d,
                failing_input_matrix = current,
                failing_column = [current[row, d] for row in 1:d],
                passed_peel_dimensions = Tuple(passed_dimensions),
                observed_diagnostic = message,
            )
        end
    end

    throw(ArgumentError("$(CASE_ID) Laurent column boundary did not fail before the final 2x2 block"))
end

function boundary_fixture()
    entry = _case010_entry()
    original_matrix = ToricBuilderCacheQBlocks.materialize_matrix(entry)
    normalization = Suslin.normalize_laurent_gl_matrix(original_matrix)
    boundary = _first_failing_boundary(normalization.normalized_matrix)

    return (;
        case_id = entry.id,
        source_entry = entry,
        original_matrix,
        normalization,
        normalized_matrix = normalization.normalized_matrix,
        ring = base_ring(original_matrix),
        ring_description = entry.ring.description,
        expected_diagnostic = EXPECTED_DIAGNOSTIC,
        boundary...,
    )
end

function _has_required_boundary_fields(fixture)::Bool
    return all(field -> hasproperty(fixture, field), REQUIRED_BOUNDARY_FIELDS)
end

function _column_reduction_diagnostic_status(column, R, expected_diagnostic)::Symbol
    try
        Suslin.reduce_unimodular_column(column, R)
        return :reduced
    catch err
        err isa InterruptException && rethrow()
        err isa ArgumentError || return :unexpected_error
        message = sprint(showerror, err)
        return occursin(expected_diagnostic, message) ? :expected_diagnostic : :unexpected_diagnostic
    end
end

function validate_boundary_fixture(fixture)::Symbol
    _has_required_boundary_fields(fixture) || return :missing_metadata
    fixture.case_id == CASE_ID || return :wrong_case
    fixture.first_failing_peel_dimension == 5 || return :wrong_peel_dimension
    length(fixture.failing_column) == 5 || return :wrong_column_length
    fixture.ring_description == EXPECTED_RING_DESCRIPTION || return :wrong_ring

    R = fixture.ring
    _is_expected_uv_laurent_ring(R) || return :wrong_ring
    _column_entries_are_over_ring(fixture.failing_column, R) || return :wrong_ring
    Suslin.is_unimodular_column(fixture.failing_column, R) || return :not_unimodular

    diagnostic_status = _column_reduction_diagnostic_status(
        fixture.failing_column,
        R,
        fixture.expected_diagnostic,
    )
    diagnostic_status == :expected_diagnostic || return diagnostic_status

    return :ok
end

function non_unimodular_negative_control(fixture = boundary_fixture())
    _, v = gens(fixture.ring)
    nonunit = v + one(fixture.ring)
    return merge(fixture, (; failing_column = [nonunit * entry for entry in fixture.failing_column]))
end

function single_entry_zero_perturbations(fixture = boundary_fixture())
    return [
        begin
            column = copy(fixture.failing_column)
            column[idx] = zero(fixture.ring)
            column
        end
        for idx in eachindex(fixture.failing_column)
    ]
end

end
```

- [ ] **Step 4: Register the internal test**

In `test/runtests.jl`, insert `"internal/toricbuilder_case010_column_boundary.jl",` immediately after `"internal/toricbuilder_cache_q_blocks.jl",`.

The surrounding block should become:

```julia
        "internal/toricbuilder_cache_q_blocks.jl",
        "internal/toricbuilder_case010_column_boundary.jl",
        "internal/toricbuilder_cache_status_report.jl",
```

- [ ] **Step 5: Run focused verification to verify GREEN**

Run:

```bash
julia --project=. -e 'include("test/internal/toricbuilder_case010_column_boundary.jl")'
```

Expected: PASS with a `ToricBuilder case_010 Laurent column boundary` test summary and no failures.

- [ ] **Step 6: Run internal group verification**

Run:

```bash
julia --project=. test/runtests.jl internal
```

Expected: PASS for the internal group, including the new `toricbuilder_case010_column_boundary.jl` file.

- [ ] **Step 7: Commit**

Run:

```bash
git add Project.toml test/fixtures/toricbuilder_case010_column_boundary.jl test/internal/toricbuilder_case010_column_boundary.jl test/runtests.jl
git commit -m "test: extract case010 Laurent column boundary"
```

Expected: commit succeeds and contains the stdlib test-target metadata fix, the new fixture, the new internal test, and the test registration.

---

## Self-Review

- Spec coverage: Task 1 covers the fixture, validator, diagnostic assertion, length-5 Laurent ring assertion, unimodularity assertion, negative-control rejection, and full-suite registration.
- Incomplete-marker scan: no open-ended implementation markers remain.
- Type consistency: the plan consistently uses `boundary_fixture()`, `validate_boundary_fixture(fixture)`, `non_unimodular_negative_control(fixture)`, and `single_entry_zero_perturbations(fixture)` from `ToricBuilderCase010ColumnBoundary`.
