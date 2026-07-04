# Issue 295 Case008 D15 Column Boundary Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a compact offline `case_008` `d=15` Laurent column boundary fixture and default internal validator.

**Architecture:** Follow the existing `d=16` column-boundary fixture shape. The implementation stores only static Laurent expression strings in the fixture and validates the parsed length-15 column against metadata, unimodularity, zero unit entries, and the unsupported Laurent reducer diagnostic.

**Tech Stack:** Julia, Oscar, Suslin internal diagnostics, Test stdlib.

## Global Constraints

- Fixture module path must be `test/fixtures/toricbuilder_case008_d15_column_boundary.jl`.
- Internal validator path must be `test/internal/toricbuilder_case008_d15_column_boundary.jl`.
- `case_id == "case_008"`.
- `source_matrix_dimensions == (30, 30)`.
- `source_column_transformation_dimensions == (60, 60)`.
- `passed_peel_dimensions == (30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16)`.
- `first_failing_peel_dimension == 15`.
- `ring_description == "GF(2)[u^+/-1, v^+/-1]"`.
- The fixture must be column-only; do not store a full `15x15` matrix.
- Do not implement a new reducer stage.
- The negative control must multiply every entry by `v + 1`, validate as `:not_unimodular`, and make `Suslin.reduce_unimodular_column` throw.
- The focused verification command must be `julia --project=. -e 'include("test/internal/toricbuilder_case008_d15_column_boundary.jl")'`.
- The full verification command must be `julia --project=. -e 'using Pkg; Pkg.test()'`.

---

### Task 1: Add Case008 D15 Column Boundary Fixture

**Files:**
- Create: `test/fixtures/toricbuilder_case008_d15_column_boundary.jl`
- Create: `test/internal/toricbuilder_case008_d15_column_boundary.jl`
- Modify: `test/runtests.jl`

**Interfaces:**
- Consumes: `Suslin.suslin_laurent_polynomial_ring(GF(2), ["u", "v"])`, `Suslin.is_unimodular_column`, `Suslin.diagnose_unimodular_column_reduction`, `Suslin.reduce_unimodular_column`.
- Consumes: generated static entry tuple at `/tmp/case008_d15_column_entries.jl`; if missing, regenerate it with:

```bash
julia --project=. -e 'using Oscar, Suslin; include("test/fixtures/toricbuilder_case008_d16_matrix_boundary.jl"); fixture = ToricBuilderCase008D16MatrixBoundary.matrix_fixture(); @assert ToricBuilderCase008D16MatrixBoundary.validate_matrix_fixture(fixture) == :ok; step = Suslin._laurent_column_peel_step(fixture.failing_input_matrix); B = step.next_block; d = nrows(B); col = [B[row, d] for row in 1:d]; open("/tmp/case008_d15_column_entries.jl", "w") do io; println(io, "const FAILING_COLUMN_ENTRIES = ("); for entry in col; println(io, "    ", repr(string(entry)), ","); end; println(io, ")"); end; println("wrote /tmp/case008_d15_column_entries.jl with ", length(col), " entries")'
```

- Produces: module `ToricBuilderCase008D15ColumnBoundary`.
- Produces: `boundary_fixture() -> NamedTuple`.
- Produces: `validate_boundary_fixture(fixture)::Symbol`.
- Produces: `non_unimodular_negative_control(fixture = boundary_fixture())`.

- [ ] **Step 1: Write the failing internal validator test**

Create `test/internal/toricbuilder_case008_d15_column_boundary.jl`:

```julia
using Test
using Suslin
using Oscar

const TORICBUILDER_CASE008_D15_COLUMN_BOUNDARY_PATH =
    joinpath(@__DIR__, "..", "fixtures", "toricbuilder_case008_d15_column_boundary.jl")

@testset "ToricBuilder case_008 d=15 Laurent column boundary" begin
    @test isfile(TORICBUILDER_CASE008_D15_COLUMN_BOUNDARY_PATH)

    include(TORICBUILDER_CASE008_D15_COLUMN_BOUNDARY_PATH)
    fixture = ToricBuilderCase008D15ColumnBoundary.boundary_fixture()

    @test fixture.case_id == "case_008"
    @test fixture.source_matrix_dimensions == (30, 30)
    @test fixture.source_column_transformation_dimensions == (60, 60)
    @test fixture.first_failing_peel_dimension == 15
    @test fixture.passed_peel_dimensions == (30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16)
    @test length(fixture.failing_column) == 15
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
    @test diagnostic.column_length == 15
    @test diagnostic.ring_profile.kind == :laurent_polynomial
    @test diagnostic.ring_profile.generators == ("u", "v")
    @test ToricBuilderCase008D15ColumnBoundary.validate_boundary_fixture(fixture) == :ok
    @test_throws ArgumentError Suslin.reduce_unimodular_column(
        fixture.failing_column,
        fixture.ring,
    )

    corrupted = ToricBuilderCase008D15ColumnBoundary.non_unimodular_negative_control(fixture)
    @test !Suslin.is_unimodular_column(corrupted.failing_column, corrupted.ring)
    @test ToricBuilderCase008D15ColumnBoundary.validate_boundary_fixture(corrupted) == :not_unimodular
    @test_throws ArgumentError Suslin.reduce_unimodular_column(
        corrupted.failing_column,
        corrupted.ring,
    )

    wrong_dimension = merge(fixture, (; first_failing_peel_dimension = 16))
    @test ToricBuilderCase008D15ColumnBoundary.validate_boundary_fixture(wrong_dimension) == :wrong_peel_dimension

    wrong_passed = merge(fixture, (; passed_peel_dimensions = (30, 29)))
    @test ToricBuilderCase008D15ColumnBoundary.validate_boundary_fixture(wrong_passed) == :wrong_passed_peel_dimensions

    wrong_ring, (wrong_u, wrong_v) = Suslin.suslin_laurent_polynomial_ring(GF(3), ["u", "v"])
    wrong_ring_fixture = merge(
        fixture,
        (;
            ring = wrong_ring,
            failing_column = [idx == 1 ? wrong_u : idx == 2 ? wrong_v : zero(wrong_ring) for idx in 1:15],
        ),
    )
    @test ToricBuilderCase008D15ColumnBoundary.validate_boundary_fixture(wrong_ring_fixture) == :wrong_ring
end
```

Modify `test/runtests.jl` to include the new internal validator near the existing `case_008` boundary validators:

```julia
include("internal/toricbuilder_case008_d15_column_boundary.jl")
```

- [ ] **Step 2: Run the focused test to verify RED**

Run:

```bash
julia --project=. -e 'include("test/internal/toricbuilder_case008_d15_column_boundary.jl")'
```

Expected: FAIL because `test/fixtures/toricbuilder_case008_d15_column_boundary.jl` does not exist yet, and the first `@test isfile(...)` is false.

- [ ] **Step 3: Add the fixture module**

Create `test/fixtures/toricbuilder_case008_d15_column_boundary.jl` by following `test/fixtures/toricbuilder_case008_d16_column_boundary.jl` with these exact changes:

```julia
module ToricBuilderCase008D15ColumnBoundary

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

# Generated by replaying one certified peel step from the d=16 matrix boundary fixture.
const EXPECTED_PASSED_PEEL_DIMENSIONS = (30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16)
const FIRST_FAILING_PEEL_DIMENSION = 15
```

Immediately after those constants, paste the complete `FAILING_COLUMN_ENTRIES` tuple from `/tmp/case008_d15_column_entries.jl`.

Then copy the parsing, ring validation, `boundary_fixture`, `validate_boundary_fixture`, and `non_unimodular_negative_control` helpers from the `d=16` column fixture. Update all module names and dimension constants to `D15`/`15`, and simplify `_diagnostic_matches_expected` so it only checks:

```julia
diagnostic.status == expected.status &&
    diagnostic.failure_code == expected.failure_code &&
    diagnostic.column_length == FIRST_FAILING_PEEL_DIMENSION
```

- [ ] **Step 4: Run the focused test to verify GREEN**

Run:

```bash
julia --project=. -e 'include("test/internal/toricbuilder_case008_d15_column_boundary.jl")'
```

Expected: PASS, with all tests in the testset succeeding.

- [ ] **Step 5: Run the default test entrypoint**

Run:

```bash
julia --project=. test/runtests.jl
```

Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add test/fixtures/toricbuilder_case008_d15_column_boundary.jl test/internal/toricbuilder_case008_d15_column_boundary.jl test/runtests.jl
git commit -m "test: add case008 d15 Laurent column boundary"
```

## Plan Self-Review

- Spec coverage: the task creates the required fixture, validator, negative control, and test registration.
- Placeholder scan: no placeholder markers or deferred implementation steps remain.
- Type consistency: fixture interfaces match the existing `d=16` boundary fixture and the internal validator consumes those exact names.
