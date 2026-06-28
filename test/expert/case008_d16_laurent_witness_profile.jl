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
    witness_nonunit_obstructions = Tuple(
        (; index = idx, is_unit = false, reason = :non_unit_witness_entry)
        for idx in eachindex(witness_entry_is_unit)
        if !witness_entry_is_unit[idx]
    )
    column_dot_witness = _laurent_column_dot(column, witness, R)

    return (;
        witness_length = length(witness),
        witness_unit_entry_count,
        witness_unit_indices,
        witness_nonunit_obstructions,
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
    @test length(profile.witness_nonunit_obstructions) == 16
    @test Tuple(record.index for record in profile.witness_nonunit_obstructions) == ntuple(identity, 16)
    @test all(record -> record.is_unit === false, profile.witness_nonunit_obstructions)
    @test all(record -> record.reason === :non_unit_witness_entry, profile.witness_nonunit_obstructions)
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
    @test known_profile.witness_nonunit_obstructions == (
        (; index = 3, is_unit = false, reason = :non_unit_witness_entry),
    )
    @test known_profile.witness_unit_entry_count >= 1
    @test known_profile.witness_unit_indices == (1, 2)
    @test known_profile.existing_witness_unit_stage_applicable
    @test known_profile.unit_obstruction_reason == :none

    solver_witness = Suslin._laurent_unimodular_witness(column, R)
    @test solver_witness !== nothing
    solver_profile = _laurent_witness_profile(column, solver_witness, R)
    @test length(solver_profile.witness_nonunit_obstructions) < solver_profile.witness_length
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
