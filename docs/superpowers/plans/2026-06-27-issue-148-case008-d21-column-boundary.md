# Issue 148 Case 008 D21 Laurent Column Boundary Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a fast offline fixture and validator for the `case_008` Laurent-normalized column-reduction boundary at peel dimension `d=21`.

**Architecture:** Keep the new behavior entirely in test support. A one-shot implementation helper derives the `d=21` boundary from the real `case_008` path and writes a snapshot fixture; routine tests load only the stored `21x21` snapshot and validate the unsupported diagnostic locally.

**Tech Stack:** Julia, Oscar, Suslin internals, `Test` stdlib, existing ToricBuilder sparse fixture catalog.

## Global Constraints

- Source case must be `case_008` from `test/fixtures/toricbuilder_cache_q_blocks.jl`.
- New fixture path must be `test/fixtures/toricbuilder_case008_d21_column_boundary.jl`.
- New validator path must be `test/internal/toricbuilder_case008_d21_column_boundary.jl`.
- The fixture module name must be `ToricBuilderCase008D21ColumnBoundary`.
- The fixture must use a stored `21x21` boundary snapshot for routine validation.
- Normal fixture validation must not rerun `normalize_laurent_gl_matrix` or the `30x30` peel path.
- `fixture.first_failing_peel_dimension == 21`.
- `fixture.passed_peel_dimensions == (30, 29, 28, 27, 26, 25, 24, 23, 22)`.
- `length(fixture.failing_column) == 21`.
- The ring must be `GF(2)[u^+/-1, v^+/-1]` with generators `("u", "v")`.
- `Suslin.is_unimodular_column(fixture.failing_column, fixture.ring)` must be true.
- `Suslin.diagnose_unimodular_column_reduction(fixture.failing_column, fixture.ring)` must return `status == :unsupported` and `failure_code == :unsupported_laurent_column_family`.
- The validator must reject a non-unimodular negative control as `:not_unimodular`.
- Do not fix the `case_008` Laurent column reducer.
- Do not require `laurent_gl_factorization_certificate` to pass for the original `case_008` matrix.
- Do not add `case_008` to the default ToricBuilder cache report.
- Do not add a checked-in live reconstruction script.

---

### Task 1: Add Case 008 D21 Snapshot Fixture And Focused Validator

**Files:**
- Create: `test/internal/toricbuilder_case008_d21_column_boundary.jl`
- Create: `test/fixtures/toricbuilder_case008_d21_column_boundary.jl`

**Interfaces:**
- Consumes: `ToricBuilderCacheQBlocks.catalog()`, `ToricBuilderCacheQBlocks.materialize_matrix(entry)`, `ToricBuilderCacheQBlocks._sparse_laurent_matrix(R, u, v, rows, cols, entries)`, `Suslin.normalize_laurent_gl_matrix(A)`, `Suslin._laurent_column_peel_step(current)`, `Suslin.is_unimodular_column(column, R)`, `Suslin.diagnose_unimodular_column_reduction(column, R)`, `Suslin.suslin_laurent_polynomial_ring(base, names)`.
- Produces: `ToricBuilderCase008D21ColumnBoundary.boundary_fixture()::NamedTuple`, `ToricBuilderCase008D21ColumnBoundary.validate_boundary_fixture(fixture)::Symbol`, and `ToricBuilderCase008D21ColumnBoundary.non_unimodular_negative_control(fixture = boundary_fixture())::NamedTuple`.

- [ ] **Step 1: Write the failing focused validator test**

Create `test/internal/toricbuilder_case008_d21_column_boundary.jl` with this content:

```julia
using Test
using Suslin
using Oscar

const TORICBUILDER_CASE008_D21_COLUMN_BOUNDARY_PATH =
    joinpath(@__DIR__, "..", "fixtures", "toricbuilder_case008_d21_column_boundary.jl")

@testset "ToricBuilder case_008 d=21 Laurent column boundary" begin
    @test isfile(TORICBUILDER_CASE008_D21_COLUMN_BOUNDARY_PATH)

    include(TORICBUILDER_CASE008_D21_COLUMN_BOUNDARY_PATH)
    fixture = ToricBuilderCase008D21ColumnBoundary.boundary_fixture()

    @test fixture.case_id == "case_008"
    @test fixture.source_matrix_dimensions == (30, 30)
    @test (nrows(fixture.original_matrix), ncols(fixture.original_matrix)) == (30, 30)
    @test fixture.first_failing_peel_dimension == 21
    @test fixture.passed_peel_dimensions == (30, 29, 28, 27, 26, 25, 24, 23, 22)
    @test (nrows(fixture.failing_input_matrix), ncols(fixture.failing_input_matrix)) == (21, 21)
    @test length(fixture.failing_column) == 21
    @test fixture.failing_column == [fixture.failing_input_matrix[row, 21] for row in 1:21]
    @test fixture.ring_description == "GF(2)[u^+/-1, v^+/-1]"
    @test Tuple(string.(gens(fixture.ring))) == ("u", "v")
    @test Suslin.is_unimodular_column(fixture.failing_column, fixture.ring)

    diagnostic = Suslin.diagnose_unimodular_column_reduction(
        fixture.failing_column,
        fixture.ring,
    )
    @test diagnostic.status == :unsupported
    @test diagnostic.failure_code == :unsupported_laurent_column_family
    @test diagnostic.column_length == 21
    @test diagnostic.ring_profile.kind == :laurent_polynomial
    @test diagnostic.ring_profile.generators == ("u", "v")

    @test ToricBuilderCase008D21ColumnBoundary.validate_boundary_fixture(fixture) == :ok

    negative = ToricBuilderCase008D21ColumnBoundary.non_unimodular_negative_control(fixture)
    @test !Suslin.is_unimodular_column(negative.failing_column, negative.ring)
    @test ToricBuilderCase008D21ColumnBoundary.validate_boundary_fixture(negative) == :not_unimodular

    wrong_dimension = merge(fixture, (; first_failing_peel_dimension = 20))
    @test ToricBuilderCase008D21ColumnBoundary.validate_boundary_fixture(wrong_dimension) == :wrong_peel_dimension

    wrong_ring, (wrong_u, wrong_v) = Suslin.suslin_laurent_polynomial_ring(GF(3), ["u", "v"])
    wrong_ring_fixture = merge(
        fixture,
        (;
            ring = wrong_ring,
            failing_column = [
                idx == 1 ? wrong_u :
                idx == 2 ? wrong_v :
                idx == 3 ? one(wrong_ring) :
                zero(wrong_ring)
                for idx in 1:21
            ],
        ),
    )
    @test ToricBuilderCase008D21ColumnBoundary.validate_boundary_fixture(wrong_ring_fixture) == :wrong_ring
end
```

- [ ] **Step 2: Run the focused test to verify RED**

Run:

```bash
julia --project=. -e 'include("test/internal/toricbuilder_case008_d21_column_boundary.jl")'
```

Expected: FAIL because `test/fixtures/toricbuilder_case008_d21_column_boundary.jl` does not exist and the `isfile` assertion is false.

- [ ] **Step 3: Write the one-shot snapshot generator**

Create `.tmp_issue148_snapshot_generator.jl` with this content:

```julia
using Oscar
using Suslin

include(joinpath(@__DIR__, "test", "fixtures", "toricbuilder_cache_q_blocks.jl"))

const CASE_ID = "case_008"
const TARGET_DIMENSION = 21
const PASSED_DIMENSIONS = (30, 29, 28, 27, 26, 25, 24, 23, 22)
const FIXTURE_PATH = joinpath(@__DIR__, "test", "fixtures", "toricbuilder_case008_d21_column_boundary.jl")

function case008_entry()
    matches = filter(entry -> entry.id == CASE_ID, ToricBuilderCacheQBlocks.catalog().cases)
    length(matches) == 1 ||
        throw(ArgumentError("expected exactly one ToricBuilder cache Q-block entry for " * CASE_ID))
    return only(matches)
end

function matrix_sparse_entries(M)
    entries = Tuple{Int, Int, String}[]
    for row in 1:nrows(M), col in 1:ncols(M)
        value = M[row, col]
        iszero(value) && continue
        push!(entries, (row, col, string(value)))
    end
    return entries
end

function entries_literal(entries)
    lines = String["const FAILING_INPUT_SPARSE_ENTRIES = ("]
    for (row, col, value) in entries
        push!(lines, "    (" * string(row) * ", " * string(col) * ", " * repr(value) * "),")
    end
    push!(lines, ")")
    return join(lines, "\n")
end

function derive_boundary_snapshot()
    entry = case008_entry()
    original_matrix = ToricBuilderCacheQBlocks.materialize_matrix(entry)
    normalization = Suslin.normalize_laurent_gl_matrix(original_matrix)
    current = normalization.normalized_matrix
    passed_dimensions = Int[]

    while nrows(current) > TARGET_DIMENSION
        push!(passed_dimensions, nrows(current))
        current = Suslin._laurent_column_peel_step(current).next_block
    end

    Tuple(passed_dimensions) == PASSED_DIMENSIONS ||
        error("unexpected passed peel dimensions: " * repr(Tuple(passed_dimensions)))
    nrows(current) == TARGET_DIMENSION ||
        error("expected boundary dimension " * string(TARGET_DIMENSION) * ", got " * string(nrows(current)))

    R = base_ring(current)
    column = [current[row, TARGET_DIMENSION] for row in 1:TARGET_DIMENSION]
    Suslin.is_unimodular_column(column, R) ||
        error("derived boundary column is not unimodular")
    diagnostic = Suslin.diagnose_unimodular_column_reduction(column, R)
    diagnostic.status == :unsupported ||
        error("expected unsupported diagnostic, got " * string(diagnostic.status))
    diagnostic.failure_code == :unsupported_laurent_column_family ||
        error("expected unsupported Laurent-family diagnostic, got " * string(diagnostic.failure_code))

    entries = matrix_sparse_entries(current)
    u, v = gens(R)
    round_trip = ToricBuilderCacheQBlocks._sparse_laurent_matrix(
        R,
        u,
        v,
        TARGET_DIMENSION,
        TARGET_DIMENSION,
        entries,
    )
    round_trip == current ||
        error("sparse string snapshot did not round-trip to the derived boundary matrix")

    return entries
end

function fixture_template(entries_text)
    return """
module ToricBuilderCase008D21ColumnBoundary

using Oscar
using Suslin

include(joinpath(@__DIR__, "toricbuilder_cache_q_blocks.jl"))

const CASE_ID = "case_008"
const EXPECTED_RING_DESCRIPTION = "GF(2)[u^+/-1, v^+/-1]"
const EXPECTED_PASSED_PEEL_DIMENSIONS = (30, 29, 28, 27, 26, 25, 24, 23, 22)
const FIRST_FAILING_PEEL_DIMENSION = 21
const REQUIRED_BOUNDARY_FIELDS = (
    :case_id,
    :source_entry,
    :original_matrix,
    :source_matrix_dimensions,
    :source_sparse_entry_count,
    :normalization_provenance,
    :passed_peel_dimensions,
    :first_failing_peel_dimension,
    :failing_input_matrix,
    :failing_column,
    :ring,
    :ring_description,
    :expected_diagnostic,
)

""" * entries_text * """

function _case008_entry()
    matches = filter(entry -> entry.id == CASE_ID, ToricBuilderCacheQBlocks.catalog().cases)
    length(matches) == 1 ||
        throw(ArgumentError("expected exactly one ToricBuilder cache Q-block entry for " * CASE_ID))
    return only(matches)
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

function _snapshot_matrix()
    R, (u, v) = Suslin.suslin_laurent_polynomial_ring(GF(2), ["u", "v"])
    return ToricBuilderCacheQBlocks._sparse_laurent_matrix(
        R,
        u,
        v,
        FIRST_FAILING_PEEL_DIMENSION,
        FIRST_FAILING_PEEL_DIMENSION,
        FAILING_INPUT_SPARSE_ENTRIES,
    )
end

function _matrix_size(M)
    return (nrows(M), ncols(M))
end

function _last_column(M)
    d = ncols(M)
    return [M[row, d] for row in 1:nrows(M)]
end

function _column_entries_are_over_ring(column, R)::Bool
    try
        return all(entry -> R(entry) == entry, column)
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end

function boundary_fixture()
    entry = _case008_entry()
    original_matrix = ToricBuilderCacheQBlocks.materialize_matrix(entry)
    failing_input_matrix = _snapshot_matrix()
    R = base_ring(failing_input_matrix)

    return (;
        case_id = entry.id,
        source_entry = entry,
        original_matrix,
        source_matrix_dimensions = entry.dimensions.matrix,
        source_sparse_entry_count = entry.sparse_entry_count,
        normalization_provenance = (;
            source_issue = 131,
            method = :normalize_laurent_gl_matrix,
            determinant_classification = :one,
            normalized_matrix_dimensions = entry.dimensions.matrix,
        ),
        passed_peel_dimensions = EXPECTED_PASSED_PEEL_DIMENSIONS,
        first_failing_peel_dimension = FIRST_FAILING_PEEL_DIMENSION,
        failing_input_matrix,
        failing_column = _last_column(failing_input_matrix),
        ring = R,
        ring_description = entry.ring.description,
        expected_diagnostic = (;
            status = :unsupported,
            failure_code = :unsupported_laurent_column_family,
        ),
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
    fixture.source_matrix_dimensions == (30, 30) || return :wrong_source_matrix_dimensions
    fixture.first_failing_peel_dimension == FIRST_FAILING_PEEL_DIMENSION || return :wrong_peel_dimension
    fixture.passed_peel_dimensions == EXPECTED_PASSED_PEEL_DIMENSIONS || return :wrong_passed_peel_dimensions
    _matrix_size(fixture.failing_input_matrix) == (FIRST_FAILING_PEEL_DIMENSION, FIRST_FAILING_PEEL_DIMENSION) ||
        return :wrong_failing_input_dimension
    length(fixture.failing_column) == FIRST_FAILING_PEEL_DIMENSION || return :wrong_column_length
    fixture.ring_description == EXPECTED_RING_DESCRIPTION || return :wrong_ring

    R = fixture.ring
    _is_expected_uv_laurent_ring(R) || return :wrong_ring
    base_ring(fixture.failing_input_matrix) == R || return :wrong_ring
    fixture.failing_column == _last_column(fixture.failing_input_matrix) || return :wrong_column
    _column_entries_are_over_ring(fixture.failing_column, R) || return :wrong_ring
    Suslin.is_unimodular_column(fixture.failing_column, R) || return :not_unimodular
    _diagnostic_matches_expected(fixture.failing_column, R, fixture.expected_diagnostic) ||
        return :wrong_diagnostic

    return :ok
end

function non_unimodular_negative_control(fixture = boundary_fixture())
    _, v = gens(fixture.ring)
    nonunit = v + one(fixture.ring)
    corrupted_matrix = copy(fixture.failing_input_matrix)
    d = FIRST_FAILING_PEEL_DIMENSION
    for row in 1:d
        corrupted_matrix[row, d] = nonunit * fixture.failing_column[row]
    end
    return merge(
        fixture,
        (;
            failing_input_matrix = corrupted_matrix,
            failing_column = _last_column(corrupted_matrix),
        ),
    )
end

end
"""
end

function main()
    entries = derive_boundary_snapshot()
    open(FIXTURE_PATH, "w") do io
        write(io, fixture_template(entries_literal(entries)))
    end
    println("wrote " * FIXTURE_PATH * " with " * string(length(entries)) * " sparse entries")
end

main()
```

- [ ] **Step 4: Run the generator to create the snapshot fixture**

Run:

```bash
julia --project=. .tmp_issue148_snapshot_generator.jl
```

Expected: PASS. The command prints a status line beginning with:

```text
wrote /Users/nzy/jcode/Suslin.jl/test/fixtures/toricbuilder_case008_d21_column_boundary.jl with
```

The line ends with ` sparse entries`, and the count after `with` and before ` sparse entries` must be a positive integer. The command may take roughly a minute because it performs the real `case_008` normalization and peel once.

- [ ] **Step 5: Remove the temporary generator**

Run:

```bash
rm -f .tmp_issue148_snapshot_generator.jl
```

Expected: PASS with no output. `git status --short .tmp_issue148_snapshot_generator.jl` prints no line.

- [ ] **Step 6: Run the focused test to verify GREEN**

Run:

```bash
julia --project=. -e 'include("test/internal/toricbuilder_case008_d21_column_boundary.jl")'
```

Expected: PASS. The testset `ToricBuilder case_008 d=21 Laurent column boundary` reports all assertions passing.

- [ ] **Step 7: Commit the focused fixture and validator**

Run:

```bash
git add test/fixtures/toricbuilder_case008_d21_column_boundary.jl test/internal/toricbuilder_case008_d21_column_boundary.jl
git commit -m "test: add case008 d21 boundary fixture"
```

Expected: PASS. The commit includes exactly the two new files from this task.

### Task 2: Register Internal Test And Run Repository Verification

**Files:**
- Modify: `test/runtests.jl`

**Interfaces:**
- Consumes: `ToricBuilderCase008D21ColumnBoundary.boundary_fixture()`, `ToricBuilderCase008D21ColumnBoundary.validate_boundary_fixture(fixture)`.
- Produces: the internal test group includes `internal/toricbuilder_case008_d21_column_boundary.jl`.

- [ ] **Step 1: Register the new internal test file**

Modify the `internal` array in `test/runtests.jl` so the ToricBuilder boundary entries read:

```julia
        "internal/toricbuilder_cache_q_blocks.jl",
        "internal/toricbuilder_case010_column_boundary.jl",
        "internal/toricbuilder_case008_d21_column_boundary.jl",
        "internal/toricbuilder_cache_status_report.jl",
```

- [ ] **Step 2: Run the focused issue verification command**

Run:

```bash
julia --project=. -e 'include("test/internal/toricbuilder_case008_d21_column_boundary.jl")'
```

Expected: PASS. The fixture loads without running the full `30x30` path, `validate_boundary_fixture(boundary_fixture()) == :ok`, and the negative control is rejected as `:not_unimodular`.

- [ ] **Step 3: Run the internal group**

Run:

```bash
julia --project=. test/runtests.jl internal
```

Expected: PASS. The internal group includes `internal/toricbuilder_case008_d21_column_boundary.jl` and all internal tests pass.

- [ ] **Step 4: Run full repository verification**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: PASS. Public and internal tests pass under the package test entry point.

- [ ] **Step 5: Check the diff for whitespace and accidental helper files**

Run:

```bash
git diff --check
git status --short
```

Expected: `git diff --check` prints no output. `git status --short` shows `test/runtests.jl` modified and no `.tmp_issue148_snapshot_generator.jl` file.

- [ ] **Step 6: Commit the test registration**

Run:

```bash
git add test/runtests.jl
git commit -m "test: register case008 d21 boundary validator"
```

Expected: PASS. The commit includes exactly `test/runtests.jl`.
