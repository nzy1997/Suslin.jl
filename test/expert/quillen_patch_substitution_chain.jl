using Test
using Suslin
using Oscar

const QUILLEN_CHAIN_CATALOG_PATH =
    joinpath(@__DIR__, "..", "fixtures", "quillen_patch_cases.jl")
if !isdefined(Main, :QuillenPatchFixtureCatalog)
    include(QUILLEN_CHAIN_CATALOG_PATH)
end

function chain_product(factors, R, n)
    product = identity_matrix(R, n)
    for factor in factors
        product *= factor
    end
    return product
end

function chain_rebuild_step(step; kwargs...)
    fields = merge((
        step_index = step.step_index,
        selected_variable = step.selected_variable,
        raw_denominator = step.raw_denominator,
        exponent = step.exponent,
        powered_denominator = step.powered_denominator,
        coverage_multiplier = step.coverage_multiplier,
        sign_convention = step.sign_convention,
        previous_coefficient = step.previous_coefficient,
        next_coefficient = step.next_coefficient,
        previous_matrix = step.previous_matrix,
        next_matrix = step.next_matrix,
        bracket_target = step.bracket_target,
        replay_metadata = step.replay_metadata,
    ), kwargs)
    return Suslin.QuillenPatchSubstitutionStep(
        fields.step_index,
        fields.selected_variable,
        fields.raw_denominator,
        fields.exponent,
        fields.powered_denominator,
        fields.coverage_multiplier,
        fields.sign_convention,
        fields.previous_coefficient,
        fields.next_coefficient,
        fields.previous_matrix,
        fields.next_matrix,
        fields.bracket_target,
        fields.replay_metadata,
    )
end

function chain_rebuild(chain; kwargs...)
    fields = merge((
        original_matrix = chain.original_matrix,
        ring = chain.ring,
        size = chain.size,
        selected_variable = chain.selected_variable,
        sign_convention = chain.sign_convention,
        solver_result = chain.solver_result,
        cumulative_coefficients = chain.cumulative_coefficients,
        intermediate_matrices = chain.intermediate_matrices,
        steps = chain.steps,
        bracket_matrices = chain.bracket_matrices,
        base_term = chain.base_term,
        metadata = chain.metadata,
        replay_metadata = chain.replay_metadata,
        verification = chain.verification,
    ), kwargs)
    return Suslin.QuillenPatchSubstitutionChain(
        fields.original_matrix,
        fields.ring,
        fields.size,
        fields.selected_variable,
        fields.sign_convention,
        fields.solver_result,
        fields.cumulative_coefficients,
        fields.intermediate_matrices,
        fields.steps,
        fields.bracket_matrices,
        fields.base_term,
        fields.metadata,
        fields.replay_metadata,
        fields.verification,
    )
end

@testset "Park-Woodburn substitution chain" begin
    entries = Main.QuillenPatchFixtureCatalog.cases_by_id()
    entry = entries["quillen-two-open-cover-qq"]
    R = entry.ring.object
    X = entry.substitution_variable
    raw = [data.denominator for data in entry.denominator_data]
    solver = Suslin.solve_quillen_denominator_cover(R, raw; max_exponent = 2)

    chain = Suslin.quillen_patch_substitution_chain(
        entry.base_matrix,
        X,
        solver;
        metadata = (; fixture_id = entry.id, consumer_issue_id = "#217"),
    )

    @test chain isa Suslin.QuillenPatchSubstitutionChain
    @test Suslin.verify_quillen_denominator_cover_solver_result(solver)
    @test Suslin.verify_quillen_patch_substitution_chain(chain)
    replay = Suslin.replay_quillen_patch_substitution_chain(chain)
    @test replay.overall_ok
    @test replay.final_coefficient_ok
    @test replay.base_term_ok
    @test replay.telescope_ok

    r = raw[1]
    expected_coefficients = [one(R), one(R) - r, zero(R)]
    @test chain.cumulative_coefficients == expected_coefficients
    @test chain.sign_convention == :park_woodburn_minus
    @test chain.selected_variable == X
    @test chain.solver_result === solver
    @test chain.metadata == (; fixture_id = entry.id, consumer_issue_id = "#217")
    @test chain.replay_metadata.sign_convention == :park_woodburn_minus
    @test chain.replay_metadata.denominator_count == length(raw)
    @test chain.replay_metadata.coverage_sum == one(R)

    expected_matrices = [
        Suslin._quillen_substitute_matrix_scaled_variable(entry.base_matrix, X, coefficient)
        for coefficient in expected_coefficients
    ]
    @test chain.intermediate_matrices == expected_matrices
    @test chain.base_term ==
          Suslin._quillen_substitute_matrix_scaled_variable(entry.base_matrix, X, zero(R))
    @test chain.base_term == last(chain.intermediate_matrices)

    @test length(chain.steps) == length(raw)
    @test length(chain.bracket_matrices) == length(raw)
    for (idx, step) in enumerate(chain.steps)
        @test step.step_index == idx
        @test step.selected_variable == X
        @test step.raw_denominator == solver.raw_denominators[idx]
        @test step.exponent == solver.exponent
        @test step.powered_denominator == solver.powered_denominators[idx]
        @test step.coverage_multiplier == solver.coverage_multipliers[idx]
        @test step.sign_convention == :park_woodburn_minus
        @test step.previous_coefficient == chain.cumulative_coefficients[idx]
        @test step.next_coefficient == chain.cumulative_coefficients[idx + 1]
        @test step.previous_matrix == chain.intermediate_matrices[idx]
        @test step.next_matrix == chain.intermediate_matrices[idx + 1]
        @test step.bracket_target == chain.bracket_matrices[idx]
        @test step.previous_matrix * step.bracket_target == step.next_matrix
        @test step.replay_metadata.step_index == idx
        @test step.replay_metadata.exponent == solver.exponent
        @test step.replay_metadata.sign_convention == :park_woodburn_minus
    end

    telescope = entry.base_matrix
    for bracket in chain.bracket_matrices
        telescope *= bracket
    end
    @test telescope == chain.base_term
    @test chain_product(chain.bracket_matrices, R, chain.size) ==
          inv(entry.base_matrix) * chain.base_term

    corrupted_exponent_steps = copy(chain.steps)
    corrupted_exponent_steps[1] = chain_rebuild_step(
        chain.steps[1];
        exponent = chain.steps[1].exponent + 1,
    )
    @test !Suslin.verify_quillen_patch_substitution_chain(
        chain_rebuild(chain; steps = corrupted_exponent_steps),
    )

    corrupted_multiplier_steps = copy(chain.steps)
    corrupted_multiplier_steps[1] = chain_rebuild_step(
        chain.steps[1];
        coverage_multiplier = chain.steps[1].coverage_multiplier + one(R),
    )
    @test !Suslin.verify_quillen_patch_substitution_chain(
        chain_rebuild(chain; steps = corrupted_multiplier_steps),
    )

    @test !Suslin.verify_quillen_patch_substitution_chain(
        chain_rebuild(chain; sign_convention = :park_woodburn_plus),
    )

    corrupted_intermediates = copy(chain.intermediate_matrices)
    corrupted_intermediates[2] =
        corrupted_intermediates[2] * elementary_matrix(chain.size, 1, 2, one(R), R)
    @test !Suslin.verify_quillen_patch_substitution_chain(
        chain_rebuild(chain; intermediate_matrices = corrupted_intermediates),
    )

    wrong_variable = collect(gens(R))[2]
    @test wrong_variable != X
    @test !Suslin.verify_quillen_patch_substitution_chain(
        chain_rebuild(chain; selected_variable = wrong_variable),
    )
end
