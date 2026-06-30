using Test
using Suslin
using Oscar

include(joinpath(@__DIR__, "..", "fixtures", "ecp_column_cases.jl"))
include(joinpath(@__DIR__, "..", "fixtures", "ecp_mainline_cases.jl"))

function _input_context_fixture_column(entry)
    return [getproperty(entry.entries, name) for name in entry.column_order]
end

function _input_context_mainline_column(entry)
    return [getproperty(entry.column_entries, name) for name in entry.column_order]
end

function _input_context_witness_total(ctx)
    total = zero(ctx.ring)
    for idx in 1:ctx.column_length
        total += ctx.unimodularity_witness[idx] * ctx.column[idx]
    end
    return total
end

function _input_context_staged_failure_reason(diagnostic)
    return diagnostic.status == :unsupported ? diagnostic.failure_code : nothing
end

function _input_context_replace_field(ctx, field::Symbol, value)
    fields = fieldnames(typeof(ctx))
    idx = findfirst(==(field), fields)
    idx === nothing && error("unknown ECPInputContext field $(field)")
    values = [getfield(ctx, name) for name in fields]
    values[idx] = value
    return Suslin.ECPInputContext(values...)
end

function _assert_checked_context(
    ctx,
    column,
    R;
    variable_order = tuple(gens(R)...),
    selected_variable = nothing,
)
    diagnostic = Suslin.diagnose_unimodular_column_reduction(column, R)
    normalized_order = tuple(Suslin._ecp_normalize_variable_order(R, variable_order)...)

    @test ctx isa Suslin.ECPInputContext
    @test Suslin.verify_ecp_input_context(ctx)
    @test ctx.column == [Suslin._coerce_into_ring(R, column[idx], "column[$idx]") for idx in eachindex(column)]
    @test ctx.ring == R
    @test ctx.ring_profile == Suslin._column_reduction_ring_profile(R)
    @test ctx.variables == tuple(gens(R)...)
    @test ctx.variable_order == normalized_order
    @test ctx.column_length == length(column)
    @test _input_context_witness_total(ctx) == one(R)
    @test ctx.unimodularity_witness == Suslin._unimodular_witness(ctx.column, R)
    if selected_variable === nothing
        @test ctx.selected_variable_index === nothing
        @test ctx.selected_variable === nothing
    else
        selected_index = Suslin._ecp_selected_variable_index(R, selected_variable)
        @test ctx.selected_variable_index == selected_index
        @test ctx.selected_variable == gens(R)[selected_index]
        @test count(==(ctx.selected_variable), ctx.variable_order) == 1
    end
    @test ctx.support_classification == diagnostic.status
    @test ctx.staged_failure_reason == _input_context_staged_failure_reason(diagnostic)
    @test ctx.staged_diagnostic == diagnostic
    @test ctx.verification.overall_ok
    @test ctx.verification.one_based_indexing_ok
    @test ctx.verification.column_ok
    @test ctx.verification.ring_profile_ok
    @test ctx.verification.variables_ok
    @test ctx.verification.variable_order_ok
    @test ctx.verification.column_length_ok
    @test ctx.verification.unimodularity_witness_ok
    @test ctx.verification.selected_variable_ok
    @test ctx.verification.staged_diagnostic_ok
end

@testset "checked ECP input contexts" begin
    mainline_cases = ECPMainlineFixtureCatalog.cases_by_id()
    column_cases = ECPColumnFixtureCatalog.cases_by_id()

    supported_entry = mainline_cases["ecp-mainline-gf2-hard-slice"]
    supported_R = supported_entry.ring.object
    supported_column = _input_context_mainline_column(supported_entry)
    supported_ctx = Suslin.ecp_input_context(
        supported_column,
        supported_R;
        variable_order = supported_entry.ring.generators,
        selected_variable = supported_entry.selected_variable.generator,
        unimodularity_witness = supported_entry.unimodularity.coefficients,
    )
    _assert_checked_context(
        supported_ctx,
        supported_column,
        supported_R;
        variable_order = supported_entry.ring.generators,
        selected_variable = supported_entry.selected_variable.generator,
    )
    @test supported_ctx.support_classification == :supported
    @test supported_ctx.staged_failure_reason === nothing

    staged_entry = column_cases["ecp-unsupported-unimodular-gf2"]
    staged_R = staged_entry.ring.object
    staged_column = _input_context_fixture_column(staged_entry)
    staged_ctx = Suslin.ecp_input_context(
        staged_column,
        staged_R;
        variable_order = tuple(gens(staged_R)...),
        selected_variable = gens(staged_R)[1],
    )
    _assert_checked_context(
        staged_ctx,
        staged_column,
        staged_R;
        selected_variable = gens(staged_R)[1],
    )
    @test staged_ctx.support_classification == :unsupported
    @test staged_ctx.staged_failure_reason == :unsupported_polynomial_column_family

    unsupported_column = [zero(staged_R), gens(staged_R)[1]^2, gens(staged_R)[1] * gens(staged_R)[2] + one(staged_R), zero(staged_R)]
    @test Suslin.is_unimodular_column(unsupported_column, staged_R)
    unsupported_ctx = Suslin.ecp_input_context(
        unsupported_column,
        staged_R;
        selected_variable = gens(staged_R)[1],
    )
    _assert_checked_context(
        unsupported_ctx,
        unsupported_column,
        staged_R;
        selected_variable = gens(staged_R)[1],
    )
    @test unsupported_ctx.column_length == 4
    @test unsupported_ctx.support_classification == :unsupported
    @test unsupported_ctx.staged_failure_reason == :unsupported_polynomial_column_family

    length4_entry = mainline_cases["ecp-mainline-length4-coupled-qq"]
    length4_R = length4_entry.ring.object
    length4_column = _input_context_mainline_column(length4_entry)
    length4_ctx = Suslin.ECPInputContext(
        length4_column,
        length4_R;
        variable_order = length4_entry.ring.generators,
        selected_variable = length4_entry.selected_variable.generator,
        unimodularity_witness = length4_entry.unimodularity.coefficients,
    )
    _assert_checked_context(
        length4_ctx,
        length4_column,
        length4_R;
        variable_order = length4_entry.ring.generators,
        selected_variable = length4_entry.selected_variable.generator,
    )
    length4_diagnostic = Suslin.diagnose_unimodular_column_reduction(length4_column, length4_R)
    @test length4_ctx.support_classification == length4_diagnostic.status
    @test length4_ctx.staged_failure_reason == _input_context_staged_failure_reason(length4_diagnostic)

    tampered_witness = collect(length4_ctx.unimodularity_witness)
    tampered_witness[1] += one(length4_R)
    @test !Suslin.verify_ecp_input_context(
        _input_context_replace_field(length4_ctx, :unimodularity_witness, tampered_witness),
    )
end

@testset "ECP input context rejects invalid boundaries" begin
    column_cases = ECPColumnFixtureCatalog.cases_by_id()
    non_unimodular_entry = column_cases["ecp-non-unimodular-gf2"]
    non_unimodular_R = non_unimodular_entry.ring.object
    non_unimodular = _input_context_fixture_column(non_unimodular_entry)
    @test_throws ArgumentError Suslin.ecp_input_context(non_unimodular, non_unimodular_R)

    staged_entry = column_cases["ecp-unsupported-unimodular-gf2"]
    R = staged_entry.ring.object
    x, y = gens(R)
    staged_column = _input_context_fixture_column(staged_entry)
    @test_throws ArgumentError Suslin.ecp_input_context([one(R), x], R)
    @test_throws ArgumentError Suslin.ecp_input_context(staged_column, R; selected_variable = x + y)

    mainline_cases = ECPMainlineFixtureCatalog.cases_by_id()
    length4_entry = mainline_cases["ecp-mainline-length4-coupled-qq"]
    length4_R = length4_entry.ring.object
    length4_column = _input_context_mainline_column(length4_entry)
    bad_witness = collect(length4_entry.unimodularity.coefficients)
    bad_witness[1] = zero(length4_R)
    @test_throws ArgumentError Suslin.ecp_input_context(
        length4_column,
        length4_R;
        selected_variable = length4_entry.selected_variable.generator,
        unimodularity_witness = bad_witness,
    )
end
