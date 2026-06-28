# case008 d16 Laurent Witness Profile Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add the expert witness-profile test requested by issue #167 for the d=16 `case_008` Laurent witness family.

**Architecture:** Keep the profile helpers local to one expert test file. The test computes exact algebraic profile fields from the fixture column and the raw Laurent witness, then exercises a small synthetic Laurent column that should pass through the existing witness-unit reduction route.

**Tech Stack:** Julia, Test stdlib, Suslin.jl internals, Oscar Laurent polynomial APIs (`length`, `exponents`, `gcd`, `is_unit`).

## Global Constraints

- Repository has no `AGENTS.md`; follow existing Julia test style in `test/expert` and `test/internal`.
- Do not implement preconditioning search or a new reducer stage.
- Do not add a public API or production diagnostic function for this exploratory profile.
- Add the expert test at `test/expert/case008_d16_laurent_witness_profile.jl`.
- The d=16 test must assert `witness_length == 16`, `witness_unit_entry_count == 0`, and `column ⋅ witness == 1`.
- The profile must record every witness entry failing `is_unit` as a concrete reason the existing witness-unit stage cannot apply.
- The same test file must include a small synthetic Laurent column with a known unit witness; its profile must report `witness_unit_entry_count >= 1`, and `Suslin._reduce_via_laurent_witness_unit_certificate(column, R)` must reduce it successfully.
- Keep the d=16 witness solve single-pass in the d=16 profile test because local probing measured the solve at about 8.4 seconds.
- Verification command from the issue: `julia --project=. -e 'include("test/expert/case008_d16_laurent_witness_profile.jl")'`.
- Required final verification command from Agent Desk: `julia --project=. -e 'using Pkg; Pkg.test()'`.

---

## File Structure

- Create `test/expert/case008_d16_laurent_witness_profile.jl`: local profile helpers, exact d=16 profile assertions, and synthetic witness-unit route control.
- Do not modify `src/`; this issue is a profiling test only.

---

### Task 1: Add Expert Laurent Witness Profile Test

**Files:**
- Create: `test/expert/case008_d16_laurent_witness_profile.jl`

**Interfaces:**
- Consumes: `ToricBuilderCase008D16ColumnBoundary.boundary_fixture()`, `Suslin._laurent_unimodular_witness(column, R)`, `Suslin._reduce_via_laurent_witness_unit_certificate(column, R)`, `Suslin._apply_reduction_factors(factors, column, R)`, and `Suslin._target_reduced_column(R, n)`.
- Produces: a local test helper `_laurent_witness_profile(column, witness, R)::NamedTuple` with stable fields used only by this expert test.

- [ ] **Step 1: Write the failing expert test shell**

Create `test/expert/case008_d16_laurent_witness_profile.jl` with this initial content. This intentionally references `_laurent_witness_profile` before it exists so the first run proves the test exercises the missing profile behavior.

```julia
using Test
using Suslin
using Oscar

include(joinpath(@__DIR__, "..", "fixtures", "toricbuilder_case008_d16_column_boundary.jl"))

@testset "case_008 d=16 Laurent witness profile" begin
    fixture = ToricBuilderCase008D16ColumnBoundary.boundary_fixture()
    R = fixture.ring
    column = fixture.failing_column
    witness = Suslin._laurent_unimodular_witness(column, R)

    @test witness !== nothing
    profile = _laurent_witness_profile(column, witness, R)

    @test profile.witness_length == 16
    @test profile.witness_unit_entry_count == 0
    @test profile.column_dot_witness == one(R)
end
```

- [ ] **Step 2: Run the expert test and verify RED**

Run:

```bash
julia --project=. -e 'include("test/expert/case008_d16_laurent_witness_profile.jl")'
```

Expected: FAIL with `UndefVarError: _laurent_witness_profile not defined`.

- [ ] **Step 3: Implement the full profile test**

Replace the entire file with:

```julia
using Test
using Suslin
using Oscar

include(joinpath(@__DIR__, "..", "fixtures", "toricbuilder_case008_d16_column_boundary.jl"))

const CASE008_D16_EXPECTED_WITNESS_TERM_COUNTS = (
    45,
    55,
    42,
    51,
    64,
    50,
    52,
    66,
    84,
    56,
    48,
    54,
    46,
    92,
    47,
    17,
)

const CASE008_D16_EXPECTED_WITNESS_SUPPORT_BOUNDS = (
    (; min_exponents = (-6, -6), max_exponents = (5, 4)),
    (; min_exponents = (-6, -6), max_exponents = (5, 4)),
    (; min_exponents = (-6, -6), max_exponents = (5, 4)),
    (; min_exponents = (-6, -6), max_exponents = (5, 4)),
    (; min_exponents = (-7, -6), max_exponents = (6, 4)),
    (; min_exponents = (-6, -6), max_exponents = (5, 4)),
    (; min_exponents = (-6, -6), max_exponents = (5, 4)),
    (; min_exponents = (-6, -7), max_exponents = (6, 4)),
    (; min_exponents = (-8, -8), max_exponents = (7, 5)),
    (; min_exponents = (-6, -6), max_exponents = (5, 4)),
    (; min_exponents = (-6, -6), max_exponents = (5, 3)),
    (; min_exponents = (-6, -6), max_exponents = (5, 4)),
    (; min_exponents = (-6, -6), max_exponents = (5, 3)),
    (; min_exponents = (-8, -8), max_exponents = (7, 4)),
    (; min_exponents = (-6, -6), max_exponents = (5, 4)),
    (; min_exponents = (-2, -3), max_exponents = (3, 2)),
)

function _laurent_column_dot(column::AbstractVector, witness::AbstractVector, R)
    length(column) == length(witness) ||
        throw(ArgumentError("column and witness lengths must match"))

    value = zero(R)
    for idx in eachindex(column)
        value += column[idx] * witness[idx]
    end
    return value
end

function _laurent_support_bounds(entry)
    entry_exponents = collect(exponents(entry))
    isempty(entry_exponents) && return nothing

    dimension = length(first(entry_exponents))
    return (;
        min_exponents = ntuple(idx -> minimum(exponent[idx] for exponent in entry_exponents), dimension),
        max_exponents = ntuple(idx -> maximum(exponent[idx] for exponent in entry_exponents), dimension),
    )
end

function _gcd_is_unit_or_nothing(values)
    state = iterate(values)
    state === nothing && return nothing

    current, next_state = state
    try
        while true
            state = iterate(values, next_state)
            state === nothing && return is_unit(current)
            value, next_state = state
            current = gcd(current, value)
        end
    catch err
        err isa MethodError || err isa ErrorException || err isa ArgumentError || rethrow()
        return nothing
    end
end

function _column_witness_entry_gcd_units(column::AbstractVector, witness::AbstractVector)
    return Tuple(_gcd_is_unit_or_nothing((column[idx], witness[idx])) for idx in eachindex(column))
end

function _elementary_pair_unit_attempts(witness::AbstractVector)
    attempts = NamedTuple[]
    for left in 1:(length(witness) - 1), right in (left + 1):length(witness)
        plus = witness[left] + witness[right]
        minus = witness[left] - witness[right]
        push!(
            attempts,
            (;
                indices = (left, right),
                operation = :+,
                is_unit = is_unit(plus),
                term_count = length(plus),
                support_bounds = _laurent_support_bounds(plus),
            ),
        )
        push!(
            attempts,
            (;
                indices = (left, right),
                operation = :-,
                is_unit = is_unit(minus),
                term_count = length(minus),
                support_bounds = _laurent_support_bounds(minus),
            ),
        )
    end
    return Tuple(attempts)
end

function _laurent_witness_profile(column::AbstractVector, witness::AbstractVector, R)
    length(column) == length(witness) ||
        throw(ArgumentError("column and witness lengths must match"))

    witness_entry_is_unit = Tuple(is_unit(entry) for entry in witness)
    witness_unit_indices = Tuple(idx for idx in eachindex(witness) if witness_entry_is_unit[idx])
    witness_unit_entry_count = length(witness_unit_indices)
    column_dot_witness = _laurent_column_dot(column, witness, R)

    return (;
        witness_length = length(witness),
        witness_unit_entry_count,
        witness_unit_indices,
        witness_entry_term_counts = Tuple(length(entry) for entry in witness),
        witness_entry_support_bounds = Tuple(_laurent_support_bounds(entry) for entry in witness),
        witness_entry_is_unit,
        witness_gcd_is_unit = _gcd_is_unit_or_nothing(witness),
        column_witness_entry_gcd_is_unit = _column_witness_entry_gcd_units(column, witness),
        column_dot_witness,
        column_dot_witness_is_one = column_dot_witness == one(R),
        existing_witness_unit_stage_applicable = witness_unit_entry_count > 0,
        elementary_pair_unit_attempts = _elementary_pair_unit_attempts(witness),
        unit_obstruction_reason = witness_unit_entry_count == 0 ? :no_witness_unit_entry : :none,
    )
end

@testset "case_008 d=16 Laurent witness profile" begin
    fixture = ToricBuilderCase008D16ColumnBoundary.boundary_fixture()
    R = fixture.ring
    column = fixture.failing_column
    witness = Suslin._laurent_unimodular_witness(column, R)

    @test witness !== nothing
    profile = _laurent_witness_profile(column, witness, R)

    @test profile.witness_length == 16
    @test profile.witness_unit_entry_count == 0
    @test profile.witness_unit_indices == ()
    @test profile.witness_entry_is_unit == ntuple(_ -> false, 16)
    @test profile.witness_entry_term_counts == CASE008_D16_EXPECTED_WITNESS_TERM_COUNTS
    @test profile.witness_entry_support_bounds == CASE008_D16_EXPECTED_WITNESS_SUPPORT_BOUNDS
    @test profile.witness_gcd_is_unit === true
    @test profile.column_witness_entry_gcd_is_unit == ntuple(_ -> true, 16)
    @test profile.column_dot_witness == one(R)
    @test profile.column_dot_witness_is_one
    @test !profile.existing_witness_unit_stage_applicable
    @test profile.unit_obstruction_reason == :no_witness_unit_entry
    @test length(profile.elementary_pair_unit_attempts) == 240
    @test all(attempt -> !attempt.is_unit, profile.elementary_pair_unit_attempts)
    @test all(
        attempt -> attempt.support_bounds === nothing ||
            length(attempt.support_bounds.min_exponents) == 2 &&
            length(attempt.support_bounds.max_exponents) == 2,
        profile.elementary_pair_unit_attempts,
    )
end

@testset "synthetic Laurent unit-witness profile control" begin
    R, (u, v) = Suslin.suslin_laurent_polynomial_ring(GF(2), ["u", "v"])
    column = [u + v, u + v + one(R), zero(R)]
    known_witness = [one(R), one(R), zero(R)]

    @test count(is_unit, column) == 0
    @test _laurent_column_dot(column, known_witness, R) == one(R)

    known_profile = _laurent_witness_profile(column, known_witness, R)
    @test known_profile.witness_length == 3
    @test known_profile.witness_unit_entry_count >= 1
    @test known_profile.witness_unit_indices == (1, 2)
    @test known_profile.existing_witness_unit_stage_applicable
    @test known_profile.unit_obstruction_reason == :none

    solver_witness = Suslin._laurent_unimodular_witness(column, R)
    @test solver_witness !== nothing
    solver_profile = _laurent_witness_profile(column, solver_witness, R)
    @test solver_profile.witness_unit_entry_count >= 1
    @test solver_profile.column_dot_witness == one(R)

    certificate = Suslin._reduce_via_laurent_witness_unit_certificate(column, R)
    @test certificate !== nothing
    @test certificate.stage.kind == :witness_unit
    @test certificate.stage.pivot_index in solver_profile.witness_unit_indices

    target = Suslin._target_reduced_column(R, length(column))
    @test Suslin._apply_reduction_factors(certificate.factors, column, R) == target
    @test certificate.stage.output_column == target
end
```

- [ ] **Step 4: Run the expert test and verify GREEN**

Run:

```bash
julia --project=. -e 'include("test/expert/case008_d16_laurent_witness_profile.jl")'
```

Expected: PASS, exit 0.

- [ ] **Step 5: Run the default package test entry point**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: PASS, exit 0.

- [ ] **Step 6: Commit the implementation**

Run:

```bash
git add test/expert/case008_d16_laurent_witness_profile.jl
git commit -m "test: profile case008 d16 Laurent witness"
```

Expected: commit created with only the expert test file staged for this task.

---

## Self-Review

- Spec coverage: Task 1 creates the requested expert file, records length/unit/term/support/gcd/combination profile fields, checks `column ⋅ witness == 1`, records the no-unit obstruction, and includes the synthetic unit-witness reduction route control.
- Red-flag scan: the plan contains no incomplete implementation instructions.
- Type consistency: `_laurent_witness_profile(column, witness, R)` returns the stable field names listed in the design spec and used by the tests.
