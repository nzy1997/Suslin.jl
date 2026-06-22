using Test
using Suslin
using Oscar

include(joinpath(@__DIR__, "..", "fixtures", "ecp_column_cases.jl"))

const VARIABLE_CHANGE_REPLAY_CASE_IDS = (
    "ecp-variable-change-monic-gf2",
    "ecp-variable-change-permuted-gf2",
)

function _vc_column(entry)
    return [getproperty(entry.entries, name) for name in entry.column_order]
end

function _vc_target_column(R, n::Int)
    target = zero_matrix(R, n, 1)
    target[n, 1] = one(R)
    return target
end

function _vc_factor_product(factors, R, n::Int)
    product = identity_matrix(R, n)
    for factor in factors
        product *= factor
    end
    return product
end

function _vc_apply_factors(factors, column, R)
    return _vc_factor_product(factors, R, length(column)) * matrix(R, length(column), 1, collect(column))
end

function _vc_stage(cert)
    stages = [stage for stage in cert.stages if stage.kind == :monicity_normalization]
    @test length(stages) == 1
    return only(stages)
end

function _vc_map_values(substitution_map)
    return tuple((entry.value for entry in substitution_map)...)
end

function _vc_map_variables(substitution_map)
    return tuple((entry.variable for entry in substitution_map)...)
end

function _replace_tuple_entry(values, idx::Int, value)
    return tuple((i == idx ? value : values[i] for i in eachindex(values))...)
end

function _replace_substitution_value(substitution_map, idx::Int, value)
    replacement = collect(substitution_map)
    replacement[idx] = merge(replacement[idx], (; value = value))
    return tuple(replacement...)
end

function _tamper_stage_field(cert, field::Symbol, value)
    stages = collect(cert.stages)
    stage_idx = findfirst(stage -> stage.kind == :monicity_normalization, stages)
    stage_idx === nothing && error("certificate has no variable-change stage")
    stages[stage_idx] = merge(stages[stage_idx], NamedTuple{(field,)}((value,)))
    return Suslin.ECPColumnReductionCertificate(
        cert.original_column,
        cert.ring,
        tuple(stages...),
        cert.factors,
        cert.final_column,
        cert.verification,
    )
end

function _tamper_inverse_substitution_map(cert)
    stages = collect(cert.stages)
    stage_idx = findfirst(stage -> stage.kind == :monicity_normalization, stages)
    stage_idx === nothing && error("certificate has no variable-change stage")
    stage = stages[stage_idx]
    inverse_substitution = _replace_substitution_value(
        stage.inverse_substitution,
        stage.source_variable_index,
        stage.inverse_substitution[stage.source_variable_index].value + one(cert.ring),
    )
    return _tamper_stage_field(cert, :inverse_substitution, inverse_substitution)
end

function _tamper_selected_monic_index(cert)
    stages = collect(cert.stages)
    stage_idx = findfirst(stage -> stage.kind == :monicity_normalization, stages)
    stage_idx === nothing && error("certificate has no variable-change stage")
    stage = stages[stage_idx]
    bad_index = stage.selected_monic_index == length(stage.transformed_column) ? 1 : stage.selected_monic_index + 1
    return _tamper_stage_field(cert, :selected_monic_index, bad_index)
end

function _assert_variable_change_stage(entry)
    column = _vc_column(entry)
    R = entry.ring.object
    cert = Suslin.ecp_column_reduction_certificate(column, R)
    stage = _vc_stage(cert)
    ring_gens = tuple(gens(R)...)
    target = _vc_target_column(R, length(column))
    required_fields = (
        :variable_order,
        :source_variable_index,
        :source_variable,
        :target_variable_index,
        :target_variable,
        :shift_polynomial,
        :forward_substitution,
        :inverse_substitution,
        :transformed_column,
        :selected_monic_index,
        :selected_monic_entry,
        :first_coordinate_strategy,
        :first_coordinate_move_factors,
        :first_coordinate_column,
        :variable_change_verification,
    )
    missing_fields = filter(field -> !hasproperty(stage, field), required_fields)

    @test !any(is_unit, column)
    @test Suslin._reduce_supported_unimodular_column_certificate(column, R) === nothing
    @test isempty(missing_fields)
    !isempty(missing_fields) && return
    @test Suslin.verify_ecp_column_reduction(cert)

    @test stage.variable_order == ring_gens
    @test Tuple(Symbol(string(gen)) for gen in stage.variable_order) == entry.variable_order
    @test stage.source_variable_index == stage.variable_index
    @test stage.source_variable == ring_gens[stage.source_variable_index]
    @test stage.target_variable_index == stage.last_variable_index
    @test stage.target_variable == ring_gens[end]
    @test stage.shift_polynomial == stage.shift_sign * stage.target_variable^stage.shift_power
    @test _vc_map_variables(stage.forward_substitution) == ring_gens
    @test _vc_map_variables(stage.inverse_substitution) == ring_gens
    @test _vc_map_values(stage.forward_substitution) == stage.forward_values
    @test _vc_map_values(stage.inverse_substitution) == stage.inverse_values

    recomputed_transformed = tuple((
        R(evaluate(entry, collect(stage.forward_values))) for entry in column
    )...)
    @test stage.transformed_column == recomputed_transformed
    @test 1 <= stage.selected_monic_index <= length(stage.transformed_column)
    @test stage.selected_monic_entry == stage.transformed_column[stage.selected_monic_index]
    @test Suslin._is_monic_in_last_variable(stage.selected_monic_entry, R)

    expected_strategy = stage.selected_monic_index == 1 ? :already_first : :not_moved
    @test stage.first_coordinate_strategy == expected_strategy
    @test isempty(stage.first_coordinate_move_factors)
    @test stage.first_coordinate_column == stage.transformed_column

    @test _vc_apply_factors(stage.transformed_factors, collect(stage.transformed_column), R) == target
    @test all(factor -> base_ring(factor) == R, stage.inverse_substituted_factors)
    @test _vc_apply_factors(stage.inverse_substituted_factors, column, R) == target
    @test stage.variable_change_verification.selected_monic_ok == true
    @test stage.variable_change_verification.transformed_reduction_ok == true
    @test stage.variable_change_verification.original_reduction_ok == true

    tamper_cases = (
        (:variable_order, _replace_tuple_entry(ring_gens, 1, ring_gens[end])),
        (:source_variable_index, stage.source_variable_index == 1 ? 2 : 1),
        (:source_variable, ring_gens[end]),
        (:target_variable_index, stage.target_variable_index == 1 ? 2 : 1),
        (:target_variable, ring_gens[1]),
        (:shift_polynomial, zero(parent(stage.shift_polynomial))),
        (:forward_substitution, _replace_substitution_value(
            stage.forward_substitution,
            stage.source_variable_index,
            stage.forward_substitution[stage.source_variable_index].value + one(R),
        )),
        (:inverse_substitution, _replace_substitution_value(
            stage.inverse_substitution,
            stage.source_variable_index,
            stage.inverse_substitution[stage.source_variable_index].value + one(R),
        )),
        (:transformed_column, _replace_tuple_entry(stage.transformed_column, 1, target[end, 1])),
        (:selected_monic_index, stage.selected_monic_index == length(stage.transformed_column) ? 1 : stage.selected_monic_index + 1),
        (:selected_monic_entry, zero(R)),
        (:first_coordinate_strategy, stage.first_coordinate_strategy == :already_first ? :not_moved : :already_first),
        (:first_coordinate_move_factors, (identity_matrix(R, length(column)),)),
        (:first_coordinate_column, _replace_tuple_entry(stage.first_coordinate_column, 1, target[end, 1])),
        (:variable_change_verification, merge(stage.variable_change_verification, (; selected_monic_ok = false))),
    )

    for (field, bad_value) in tamper_cases
        @test !Suslin.verify_ecp_column_reduction(_tamper_stage_field(cert, field, bad_value))
    end

    @test !Suslin.verify_ecp_column_reduction(_tamper_inverse_substitution_map(cert))
    @test !Suslin.verify_ecp_column_reduction(_tamper_selected_monic_index(cert))
end

@testset "ECP variable-change replay records" begin
    cases = ECPColumnFixtureCatalog.cases_by_id()
    @test Set(VARIABLE_CHANGE_REPLAY_CASE_IDS) ⊆ Set(keys(cases))
    for id in VARIABLE_CHANGE_REPLAY_CASE_IDS
        _assert_variable_change_stage(cases[id])
    end
end
