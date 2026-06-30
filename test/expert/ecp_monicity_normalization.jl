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

function _mn_column_matrix(column, R)
    return matrix(R, length(column), 1, collect(column))
end

function _mn_shift_inverse_column(entry)
    column = _mn_column(entry)
    R = entry.ring.object
    x, y = gens(R)
    inverse_values = (x + y^2, y)
    return [evaluate(column_entry, collect(inverse_values)) for column_entry in column]
end

function _mn_factor_product(factors, R, n::Int)
    product = identity_matrix(R, n)
    for factor in factors
        product *= factor
    end
    return product
end

function _mn_apply_factors(factors, column, R)
    return _mn_factor_product(factors, R, length(column)) * _mn_column_matrix(column, R)
end

function _mn_map_values(substitution_map)
    return tuple((entry.value for entry in substitution_map)...)
end

function _mn_assert_square_factor_family(factors, R, n::Int)
    @test all(factor -> base_ring(factor) == R, factors)
    @test all(factor -> nrows(factor) == n, factors)
    @test all(factor -> ncols(factor) == n, factors)
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
    @test tuple((evaluate(evaluate(gen, collect(forward_values)), collect(inverse_values)) for gen in gens_R)...) == gens_R
    @test tuple((evaluate(evaluate(gen, collect(inverse_values)), collect(forward_values)) for gen in gens_R)...) == gens_R
end

function _mn_assert_replays(record, column, R)
    target = zero_matrix(R, length(column), 1)
    target[length(column), 1] = one(R)
    normalized_column = _mn_column_matrix(record.normalized_column, R)
    inverse_values = _mn_map_values(record.inverse_substitution)
    inverse_normalized_column = _mn_column_matrix(
        [evaluate(entry, collect(inverse_values)) for entry in record.normalized_column],
        R,
    )

    @test Suslin.verify_ecp_monicity_normalization(record)
    _mn_assert_square_factor_family(record.coordinate_move_factors, R, length(column))
    _mn_assert_square_factor_family(record.inverse_substituted_coordinate_move_factors, R, length(column))
    _mn_assert_square_factor_family(record.inverse_substituted_reduction_factors, R, length(column))
    @test _mn_apply_factors(record.coordinate_move_factors, record.transformed_column, R) == normalized_column
    @test _mn_apply_factors(record.inverse_substituted_coordinate_move_factors, column, R) == inverse_normalized_column
    @test _mn_apply_factors(record.inverse_substituted_reduction_factors, inverse_normalized_column, R) == target
    @test _mn_apply_factors(record.factors, column, R) == target
    @test record.normalized_column[1] == record.selected_monic_entry
    @test record.verification.substitution_inverse_ok
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

function _mn_context_from_column(column, entry; selected_variable)
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
    @test !Suslin.verify_ecp_monicity_normalization(
        _mn_replace_record_field(already_first_record, :ring, nothing),
    )
    @test !Suslin.verify_ecp_monicity_normalization(
        _mn_replace_record_field(already_first_record, :shift_sign, nothing),
    )
    @test Suslin._ecp_monicity_normalization_full_reduction_certificate(
        [one(already_first_ctx.ring), zero(already_first_ctx.ring)],
        already_first_ctx.ring,
    ) === nothing
    @test !Suslin._ecp_substitution_maps_are_inverse(
        already_first_ctx.ring,
        (
            (; variable = gens(already_first_ctx.ring)[1]),
            (; variable = gens(already_first_ctx.ring)[2]),
        ),
        already_first_record.inverse_substitution,
    )

    bounded_entry = column_cases["ecp-variable-change-monic-gf2"]
    bounded_preimage_column = _mn_shift_inverse_column(bounded_entry)
    bounded_ctx = _mn_context_from_column(
        bounded_preimage_column,
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
        _mn_column(permuted_entry),
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

    @test_throws ArgumentError Suslin.ecp_monicity_normalization(
        already_first_ctx;
        selected_variable = gens(already_first_ctx.ring)[1],
        selected_monic_index = 1,
        max_shift_power = -1,
    )
    missed_hint = Suslin.ecp_monicity_normalization(
        already_first_ctx;
        selected_variable = gens(already_first_ctx.ring)[1],
        selected_monic_index = length(already_first_ctx.column) + 1,
        max_shift_power = 1,
    )
    @test missed_hint isa Suslin.ECPMonicityNormalizationFailure
    @test missed_hint.selected_monic_index_hint == length(already_first_ctx.column) + 1
    @test missed_hint.search_failure.attempted_candidates > 0

    R3, (z1, z2, z3) = Oscar.polynomial_ring(GF(2), ["z1", "z2", "z3"])
    partial_order_ctx = Suslin.ecp_input_context(
        [z1, z2, one(R3)],
        R3;
        variable_order = (z1, z2),
        selected_variable = z1,
    )
    partial_order_record = Suslin.ecp_monicity_normalization(
        partial_order_ctx;
        selected_variable = z1,
        selected_monic_index = 1,
    )
    @test !Suslin.verify_ecp_monicity_normalization(
        _mn_replace_record_field(
            _mn_replace_record_field(partial_order_record, :selected_variable_index, 3),
            :selected_variable,
            z3,
        ),
    )

    unsupported_entry = column_cases["ecp-unsupported-unimodular-gf2"]
    unsupported_ctx = _mn_context_from_column(
        _mn_column(unsupported_entry),
        unsupported_entry;
        selected_variable = gens(unsupported_entry.ring.object)[2],
    )
    exhausted = Suslin.ecp_monicity_normalization(
        unsupported_ctx;
        selected_variable = gens(unsupported_ctx.ring)[2],
        max_shift_power = 0,
    )
    @test exhausted isa Suslin.ECPMonicityNormalizationFailure
    @test exhausted.context === unsupported_ctx
    @test exhausted.search_failure isa Suslin.ECPMonicitySearchFailure
    @test exhausted.search_failure.kind == :monicity_search_exhausted
    @test exhausted.max_shift_power == 0
    @test exhausted.search_failure.attempted_candidates == 0
end
