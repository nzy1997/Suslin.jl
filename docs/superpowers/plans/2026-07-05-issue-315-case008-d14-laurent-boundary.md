# Issue 315 Case008 D14 Laurent Boundary Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an offline `case_008 d=14` Laurent boundary column fixture and internal validator.

**Architecture:** Follow the existing d15 column-boundary fixture pattern. Store the generated d14 Laurent column as static expression strings, validate exact metadata and column statistics, and record the current d14 certificate-construction boundary as provenance metadata rather than running a new d14 reducer.

**Tech Stack:** Julia, Oscar, Suslin internal fixture helpers, Test stdlib.

## Global Constraints

- Fixture module path must be `test/fixtures/toricbuilder_case008_d14_column_boundary.jl`.
- Internal validator path must be `test/internal/toricbuilder_case008_d14_column_boundary.jl`.
- `case_id == "case_008"`.
- `source_block == :column_transformation_upper_left_q_block`.
- `source_matrix_dimensions == (30, 30)`.
- `source_column_transformation_dimensions == (60, 60)`.
- `passed_peel_dimensions == (30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15)`.
- `first_failing_peel_dimension == 14`.
- `ring_description == "GF(2)[u^+/-1, v^+/-1]"`.
- `last_column_nonzero_count == 14`.
- `max_entry_term_count == 3734`.
- Boundary provenance must record certificate construction at `current_peel_dimension == 14` after `last_completed_peel_dimension == 15`, with `failure_code == :unsupported_laurent_column_family` and `old_d15_boundary_cleared == true`.
- The fixture must be offline and must not require a local ToricBuilder checkout.
- Do not implement a new d14 Laurent reducer.
- Do not claim full `case_008` success.
- Do not make diagonal monomial balancing or polynomialization the primary Laurent algorithm.
- The focused verification commands are `julia --project=. -e 'include("test/internal/toricbuilder_case008_d14_column_boundary.jl")'` and `julia --project=. -e 'include("test/expert/laurent_column_reduction_diagnostics.jl")'`.
- The package verification command is `julia --project=. -e 'using Pkg; Pkg.test()'`.

---

## File Structure

- Create: `test/fixtures/toricbuilder_case008_d14_column_boundary.jl`
  - Owns d14 boundary metadata, static column entries, parsing helpers, validation, and negative-control construction.
- Create: `test/internal/toricbuilder_case008_d14_column_boundary.jl`
  - Owns focused assertions for the issue contract and negative controls.
- Modify: `test/runtests.jl`
  - Registers the new internal validator in the default internal suite.
- Modify: `test/expert/laurent_column_reduction_diagnostics.jl`
  - Only if needed for a lightweight metadata assertion. Do not run the d14 reducer diagnostic from this file.

### Task 1: Add Case008 D14 Boundary Fixture And Validator

**Files:**
- Create: `test/fixtures/toricbuilder_case008_d14_column_boundary.jl`
- Create: `test/internal/toricbuilder_case008_d14_column_boundary.jl`
- Modify: `test/runtests.jl`

**Interfaces:**
- Consumes: generated static entry tuple at `/tmp/case008_d14_column_entries.jl`; if missing, regenerate it with:

```bash
julia --project=. -e 'using Suslin, Oscar; include("test/fixtures/toricbuilder_case008_d15_matrix_boundary.jl"); f=ToricBuilderCase008D15MatrixBoundary.matrix_fixture(); @assert ToricBuilderCase008D15MatrixBoundary.validate_matrix_fixture(f) == :ok; step=Suslin._laurent_column_peel_step(f.failing_input_matrix); B=step.next_block; d=nrows(B); col=[B[row,d] for row in 1:d]; stats=Suslin._laurent_column_peel_column_stats(col); @assert d == 14; @assert stats.last_column_nnz == 14; @assert stats.max_entry_terms == 3734; open("/tmp/case008_d14_column_entries.jl", "w") do io; println(io, "const FAILING_COLUMN_ENTRIES = ("); for entry in col; println(io, "    ", repr(string(entry)), ","); end; println(io, ")"); end; println("wrote /tmp/case008_d14_column_entries.jl")'
```

- Produces: module `ToricBuilderCase008D14ColumnBoundary`.
- Produces: `boundary_fixture() -> NamedTuple`.
- Produces: `validate_boundary_fixture(fixture)::Symbol`.
- Produces: `non_unimodular_negative_control(fixture = boundary_fixture())`.

- [ ] **Step 1: Write the failing internal validator test**

Create `test/internal/toricbuilder_case008_d14_column_boundary.jl`:

```julia
using Test
using Suslin
using Oscar

const TORICBUILDER_CASE008_D14_COLUMN_BOUNDARY_PATH =
    joinpath(@__DIR__, "..", "fixtures", "toricbuilder_case008_d14_column_boundary.jl")

@testset "ToricBuilder case_008 d=14 Laurent column boundary" begin
    @test isfile(TORICBUILDER_CASE008_D14_COLUMN_BOUNDARY_PATH)

    include(TORICBUILDER_CASE008_D14_COLUMN_BOUNDARY_PATH)
    fixture = ToricBuilderCase008D14ColumnBoundary.boundary_fixture()

    @test fixture.case_id == "case_008"
    @test fixture.source_block == :column_transformation_upper_left_q_block
    @test fixture.source_matrix_dimensions == (30, 30)
    @test fixture.source_column_transformation_dimensions == (60, 60)
    @test fixture.first_failing_peel_dimension == 14
    @test fixture.passed_peel_dimensions == (30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15)
    @test length(fixture.failing_column) == 14
    @test fixture.ring_description == "GF(2)[u^+/-1, v^+/-1]"
    @test Tuple(string.(gens(fixture.ring))) == ("u", "v")
    @test base_ring(fixture.ring) == GF(2)
    @test Suslin.is_unimodular_column(fixture.failing_column, fixture.ring)
    @test count(!iszero, fixture.failing_column) == 14
    @test count(is_unit, fixture.failing_column) == 0
    @test ToricBuilderCase008D14ColumnBoundary.column_statistics(fixture.failing_column) ==
          (; last_column_nonzero_count = 14, max_entry_term_count = 3734)

    provenance = fixture.boundary_provenance
    @test provenance.source == :explicit_bounded_case008_report
    @test provenance.stage == :certificate_construction
    @test provenance.current_peel_dimension == 14
    @test provenance.last_completed_peel_dimension == 15
    @test provenance.failure_code == :unsupported_laurent_column_family
    @test provenance.old_d15_boundary_cleared == true
    @test ToricBuilderCase008D14ColumnBoundary.validate_boundary_fixture(fixture) == :ok

    corrupted = ToricBuilderCase008D14ColumnBoundary.corrupted_column_negative_control(fixture)
    @test ToricBuilderCase008D14ColumnBoundary.validate_boundary_fixture(corrupted) == :wrong_column

    non_unimodular = ToricBuilderCase008D14ColumnBoundary.non_unimodular_negative_control(fixture)
    @test !Suslin.is_unimodular_column(non_unimodular.failing_column, non_unimodular.ring)
    @test ToricBuilderCase008D14ColumnBoundary.validate_boundary_fixture(non_unimodular) == :not_unimodular

    old_d15 = merge(
        fixture,
        (;
            first_failing_peel_dimension = 15,
            boundary_provenance = merge(
                fixture.boundary_provenance,
                (;
                    current_peel_dimension = 15,
                    last_completed_peel_dimension = 16,
                    old_d15_boundary_cleared = false,
                ),
            ),
        ),
    )
    @test ToricBuilderCase008D14ColumnBoundary.validate_boundary_fixture(old_d15) == :old_d15_boundary

    wrong_passed = merge(fixture, (; passed_peel_dimensions = (30, 29)))
    @test ToricBuilderCase008D14ColumnBoundary.validate_boundary_fixture(wrong_passed) == :wrong_passed_peel_dimensions

    wrong_ring, (wrong_u, wrong_v) = Suslin.suslin_laurent_polynomial_ring(GF(3), ["u", "v"])
    wrong_ring_fixture = merge(
        fixture,
        (;
            ring = wrong_ring,
            failing_column = [idx == 1 ? wrong_u : idx == 2 ? wrong_v : zero(wrong_ring) for idx in 1:14],
        ),
    )
    @test ToricBuilderCase008D14ColumnBoundary.validate_boundary_fixture(wrong_ring_fixture) == :wrong_ring
end
```

Modify `test/runtests.jl` to include:

```julia
"internal/toricbuilder_case008_d14_column_boundary.jl",
```

near the existing d15 and d16 case_008 boundary validators.

- [ ] **Step 2: Run the focused test to verify RED**

Run:

```bash
julia --project=. -e 'include("test/internal/toricbuilder_case008_d14_column_boundary.jl")'
```

Expected: FAIL because `test/fixtures/toricbuilder_case008_d14_column_boundary.jl` does not exist yet, and the first `@test isfile(...)` is false.

- [ ] **Step 3: Add the fixture module**

Create `test/fixtures/toricbuilder_case008_d14_column_boundary.jl` by following `test/fixtures/toricbuilder_case008_d15_column_boundary.jl` with these exact changes:

```julia
module ToricBuilderCase008D14ColumnBoundary

using Oscar
using Suslin

const CASE_ID = "case_008"
const SOURCE_CACHE_FILE = "case_008.jls"
const SOURCE_BLOCK = :column_transformation_upper_left_q_block
const SOURCE_MATRIX_DIMENSIONS = (30, 30)
const SOURCE_COLUMN_TRANSFORMATION_DIMENSIONS = (60, 60)
const EXPECTED_RING_DESCRIPTION = "GF(2)[u^+/-1, v^+/-1]"
const EXPECTED_PASSED_PEEL_DIMENSIONS = (30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15)
const FIRST_FAILING_PEEL_DIMENSION = 14
const EXPECTED_LAST_COLUMN_NONZERO_COUNT = 14
const EXPECTED_MAX_ENTRY_TERM_COUNT = 3734
const EXPECTED_BOUNDARY_PROVENANCE = (;
    source = :explicit_bounded_case008_report,
    source_case = CASE_ID,
    stage = :certificate_construction,
    route_status = :certified_algorithm_boundary,
    current_peel_dimension = 14,
    last_completed_peel_dimension = 15,
    failure_code = :unsupported_laurent_column_family,
    old_d15_boundary_cleared = true,
)
```

Define `REQUIRED_BOUNDARY_FIELDS` to include:

```julia
(
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
    :boundary_provenance,
    :last_column_nonzero_count,
    :max_entry_term_count,
)
```

Immediately after those constants, paste the complete `FAILING_COLUMN_ENTRIES` tuple from `/tmp/case008_d14_column_entries.jl`.

Copy the parsing and ring validation helpers from the d15 column fixture. Add:

```julia
function _term_count(entry)::Int
    iszero(entry) && return 0
    return length(collect(coefficients(entry)))
end

function column_statistics(column)
    return (;
        last_column_nonzero_count = count(!iszero, column),
        max_entry_term_count = maximum(_term_count, column; init = 0),
    )
end
```

`boundary_fixture()` must return all required fields, including:

```julia
boundary_provenance = EXPECTED_BOUNDARY_PROVENANCE,
last_column_nonzero_count = EXPECTED_LAST_COLUMN_NONZERO_COUNT,
max_entry_term_count = EXPECTED_MAX_ENTRY_TERM_COUNT,
```

`validate_boundary_fixture(fixture)::Symbol` must check metadata in this order:

```julia
_has_required_boundary_fields(fixture) || return :missing_metadata
fixture.case_id == CASE_ID || return :wrong_case
fixture.source_case == CASE_ID || return :wrong_case
fixture.source_cache_file == SOURCE_CACHE_FILE || return :wrong_source
fixture.source_block == SOURCE_BLOCK || return :wrong_source
fixture.source_matrix_dimensions == SOURCE_MATRIX_DIMENSIONS || return :wrong_source_matrix_dimensions
fixture.source_column_transformation_dimensions == SOURCE_COLUMN_TRANSFORMATION_DIMENSIONS ||
    return :wrong_source_matrix_dimensions
_claims_old_d15_boundary(fixture) && return :old_d15_boundary
fixture.first_failing_peel_dimension == FIRST_FAILING_PEEL_DIMENSION || return :wrong_peel_dimension
fixture.passed_peel_dimensions == EXPECTED_PASSED_PEEL_DIMENSIONS || return :wrong_passed_peel_dimensions
length(fixture.failing_column) == FIRST_FAILING_PEEL_DIMENSION || return :wrong_column_length
fixture.ring_description == EXPECTED_RING_DESCRIPTION || return :wrong_ring
fixture.boundary_provenance == EXPECTED_BOUNDARY_PROVENANCE || return :wrong_boundary_provenance
fixture.last_column_nonzero_count == EXPECTED_LAST_COLUMN_NONZERO_COUNT || return :wrong_column_statistics
fixture.max_entry_term_count == EXPECTED_MAX_ENTRY_TERM_COUNT || return :wrong_column_statistics
```

Then check the ring and column:

```julia
R = fixture.ring
_is_expected_uv_laurent_ring(R) || return :wrong_ring
_column_entries_are_over_ring(fixture.failing_column, R) || return :wrong_ring
stats = column_statistics(fixture.failing_column)
stats.last_column_nonzero_count == EXPECTED_LAST_COLUMN_NONZERO_COUNT || return :wrong_column_statistics
stats.max_entry_term_count == EXPECTED_MAX_ENTRY_TERM_COUNT || return :wrong_column_statistics
Suslin.is_unimodular_column(fixture.failing_column, R) || return :not_unimodular
count(is_unit, fixture.failing_column) == 0 || return :wrong_unit_profile
u, v = gens(R)
fixture.failing_column == _failing_column(R, u, v) || return :wrong_column
return :ok
```

Define `_claims_old_d15_boundary(fixture)::Bool` so it returns true when either:

```julia
fixture.first_failing_peel_dimension == 15
```

or the provenance has both `current_peel_dimension == 15` and
`failure_code == :unsupported_laurent_column_family`, or has
`old_d15_boundary_cleared != true`. The valid d14 boundary still records
`:unsupported_laurent_column_family` as the current certificate-construction
failure code, so that symbol alone must not be treated as an old-d15 claim.

Add:

```julia
function corrupted_column_negative_control(fixture = boundary_fixture())
    corrupted = copy(fixture.failing_column)
    corrupted[1] += one(fixture.ring)
    return merge(fixture, (; failing_column = corrupted))
end

function non_unimodular_negative_control(fixture = boundary_fixture())
    _, v = gens(fixture.ring)
    nonunit = v + one(fixture.ring)
    return merge(fixture, (; failing_column = [nonunit * entry for entry in fixture.failing_column]))
end
```

- [ ] **Step 4: Run the focused test to verify GREEN**

Run:

```bash
julia --project=. -e 'include("test/internal/toricbuilder_case008_d14_column_boundary.jl")'
```

Expected: PASS, with all tests in the testset succeeding.

- [ ] **Step 5: Run the issue expert diagnostics command**

Run:

```bash
julia --project=. -e 'include("test/expert/laurent_column_reduction_diagnostics.jl")'
```

Expected: PASS. Do not add an expensive d14 reducer diagnostic if this already passes.

- [ ] **Step 6: Run default internal registration**

Run:

```bash
julia --project=. test/runtests.jl
```

Expected: PASS.

- [ ] **Step 7: Commit task changes**

Run:

```bash
git add docs/superpowers/specs/2026-07-05-issue-315-case008-d14-laurent-boundary-design.md docs/superpowers/plans/2026-07-05-issue-315-case008-d14-laurent-boundary.md test/fixtures/toricbuilder_case008_d14_column_boundary.jl test/internal/toricbuilder_case008_d14_column_boundary.jl test/runtests.jl
git commit -m "test: add case008 d14 Laurent boundary fixture"
```

## Self-Review

- Spec coverage: the task creates the required fixture, validator, negative controls, and test registration while preserving the expert diagnostics command.
- Placeholder scan: no deferred implementation steps remain; static entries are supplied by `/tmp/case008_d14_column_entries.jl`.
- Type consistency: fixture interfaces match the names used by the internal validator.
