using Test
using Suslin
using Oscar

include(joinpath(@__DIR__, "..", "fixtures", "ecp_column_cases.jl"))
include(joinpath(@__DIR__, "..", "fixtures", "ecp_mainline_cases.jl"))

function _mn_column(entry)
    return [getproperty(entry.entries, name) for name in entry.column_order]
end

function _mn_mainline_column(entry)
    return [getproperty(entry.column_entries, name) for name in entry.column_order]
end

function _mn_factor_product(factors, R, n::Int)
    product = identity_matrix(R, n)
    for factor in factors
        product *= factor
    end
    return product
end

function _mn_apply_factors(factors, column, R)
    return _mn_factor_product(factors, R, length(column)) * matrix(R, length(column), 1, collect(column))
end

function _mn_map_values(substitution_map)
    return tuple((entry.value for entry in substitution_map)...)
end

function _mn_replace_tuple_entry(values, idx::Int, value)
    return tuple((i == idx ? value : values[i] for i in eachindex(values))...)
end

function _mn_replace_record_field(record, field::Symbol, value)
    fields = fieldnames(typeof(record))
    idx = findfirst(==(field), fields)
    idx === nothing && error("unknown ECP monicity normalization field $(field)")
    values = [getfield(record, name) for name in fields]
    values[idx] = value
    return typeof(record)(values...)
end

function _mn_assert_inverse_substitution_maps(record, R)
    gens_R = tuple(gens(R)...)
    forward_values = _mn_map_values(record.forward_substitution)
    inverse_values = _mn_map_values(record.inverse_substitution)
    @test tuple((evaluate(evaluate(gen, forward_values), inverse_values) for gen in gens_R)...) == gens_R
    @test tuple((evaluate(evaluate(gen, inverse_values), forward_values) for gen in gens_R)...) == gens_R
end

function _mn_assert_replays(record, column, R)
    target = zero_matrix(R, length(column), 1)
    target[length(column), 1] = one(R)

    @test Suslin.verify_ecp_monicity_normalization(record)
    if record.selected_monic_index == 1
        @test isempty(record.coordinate_move_factors)
    else
        @test !isempty(record.coordinate_move_factors)
    end
    @test _mn_apply_factors(record.factors, column, R) == target
    @test record.normalized_column[1] == record.selected_monic_entry
    _mn_assert_inverse_substitution_maps(record, R)
end

function _mn_context_from_mainline(entry)
    column = _mn_mainline_column(entry)
    R = entry.ring.object
    return Suslin.ecp_input_context(
        column,
        R;
        variable_order = entry.ring.generators,
        selected_variable = entry.selected_variable.generator,
        unimodularity_witness = entry.unimodularity.coefficients,
    )
end

function _mn_context_from_column(entry; selected_variable)
    column = _mn_column(entry)
    R = entry.ring.object
    return Suslin.ecp_input_context(
        column,
        R;
        variable_order = entry.ring.generators,
        selected_variable,
    )
end

@testset "ECP monicity normalization" begin
    mainline_cases = ECPMainlineFixtureCatalog.cases_by_id()
    column_cases = ECPColumnFixtureCatalog.cases_by_id()

    already_first_entry = mainline_cases["ecp-mainline-gf2-hard-slice"]
    already_first_ctx = _mn_context_from_mainline(already_first_entry)
    already_first_record = Suslin.ecp_monicity_normalization(
        already_first_ctx;
        selected_variable = gens(already_first_ctx.ring)[1],
        selected_monic_index = 1,
        max_shift_power = 2,
    )
    @test already_first_record.selected_variable == gens(already_first_ctx.ring)[1]
    @test already_first_record.selected_monic_index == 1
    @test already_first_record.normalized_column[1] == already_first_record.selected_monic_entry
    _mn_assert_replays(already_first_record, already_first_ctx.column, already_first_ctx.ring)
    @test !Suslin.verify_ecp_monicity_normalization(
        _mn_replace_record_field(
            already_first_record,
            :inverse_substitution,
            _mn_replace_tuple_entry(
                already_first_record.inverse_substitution,
                1,
                already_first_record.inverse_substitution[1].value + one(already_first_ctx.ring),
            ),
        ),
    )

    bounded_entry = column_cases["ecp-variable-change-monic-gf2"]
    bounded_ctx = _mn_context_from_column(
        bounded_entry;
        selected_variable = gens(bounded_entry.ring.object)[2],
    )
    bounded_record = Suslin.ecp_monicity_normalization(
        bounded_ctx;
        selected_variable = gens(bounded_ctx.ring)[2],
        selected_monic_index = 1,
        max_shift_power = 2,
    )
    @test bounded_record.selected_variable == gens(bounded_ctx.ring)[2]
    @test bounded_record.selected_monic_index == 1
    @test bounded_record.normalized_column[1] == bounded_record.selected_monic_entry
    _mn_assert_replays(bounded_record, bounded_ctx.column, bounded_ctx.ring)
    @test !Suslin.verify_ecp_monicity_normalization(
        _mn_replace_record_field(bounded_record, :selected_variable, gens(bounded_ctx.ring)[1]),
    )
    @test !Suslin.verify_ecp_monicity_normalization(
        _mn_replace_record_field(bounded_record, :selected_monic_entry, zero(bounded_ctx.ring)),
    )

    permuted_entry = column_cases["ecp-variable-change-permuted-gf2"]
    permuted_ctx = _mn_context_from_column(
        permuted_entry;
        selected_variable = gens(permuted_entry.ring.object)[2],
    )
    permuted_record = Suslin.ecp_monicity_normalization(
        permuted_ctx;
        selected_variable = gens(permuted_ctx.ring)[2],
    )
    @test permuted_record.selected_monic_index != 1
    @test !isempty(permuted_record.coordinate_move_factors)
    @test permuted_record.normalized_column[1] == permuted_record.selected_monic_entry
    _mn_assert_replays(permuted_record, permuted_ctx.column, permuted_ctx.ring)
    @test !Suslin.verify_ecp_monicity_normalization(
        _mn_replace_record_field(permuted_record, :coordinate_move_factors, empty(permuted_record.coordinate_move_factors)),
    )
    @test !Suslin.verify_ecp_monicity_normalization(
        _mn_replace_record_field(permuted_record, :inverse_substitution, _mn_replace_tuple_entry(
            permuted_record.inverse_substitution,
            1,
            permuted_record.inverse_substitution[1].value + one(permuted_ctx.ring),
        )),
    )
end
