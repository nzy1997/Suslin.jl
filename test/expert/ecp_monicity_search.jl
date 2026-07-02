using Test
using Suslin
using Oscar

include(joinpath(@__DIR__, "..", "fixtures", "ecp_column_cases.jl"))

function _search_column(entry)
    return [getproperty(entry.entries, name) for name in entry.column_order]
end

function _search_target_column(R, n::Int)
    target = zero_matrix(R, n, 1)
    target[n, 1] = one(R)
    return target
end

function _search_factor_product(factors, R, n::Int)
    product = identity_matrix(R, n)
    for factor in factors
        product *= factor
    end
    return product
end

function _search_apply_factors(factors, column, R)
    return _search_factor_product(factors, R, length(column)) * matrix(R, length(column), 1, collect(column))
end

function _search_stage(result)
    @test result isa Suslin.ECPMonicitySearchResult
    @test result.stage.kind == :monicity_normalization
    return result.stage
end

function _search_monicity_stage(cert)
    stages = [stage for stage in cert.stages if stage.kind == :monicity_normalization]
    @test length(stages) == 1
    return only(stages)
end

function _search_tamper_stage_field(cert, field::Symbol, value)
    stages = collect(cert.stages)
    stage_idx = findfirst(stage -> stage.kind == :monicity_normalization, stages)
    stage_idx === nothing && error("certificate has no monicity-normalization stage")
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

function _search_certificate_from_stage(column, R, stage)
    factors = stage.factors
    final_column = _search_apply_factors(factors, column, R)
    cert = Suslin.ECPColumnReductionCertificate(
        column,
        R,
        ((; kind = :validation, input_length = length(column), is_unimodular = true), stage),
        factors,
        final_column,
        nothing,
    )
    verification = Suslin._ecp_column_reduction_replay_summary(cert)
    return Suslin.ECPColumnReductionCertificate(
        column,
        R,
        cert.stages,
        factors,
        final_column,
        verification,
    )
end

function _assert_success_reduces(entry; variable_order = tuple(gens(entry.ring.object)...), max_shift_power = 3)
    column = _search_column(entry)
    R = entry.ring.object
    result = Suslin._deterministic_ecp_monicity_search(
        column,
        R;
        variable_order,
        max_shift_power,
    )
    stage = _search_stage(result)
    @test result.original_column == tuple(column...)
    @test result.ring == R
    @test result.variable_order == tuple(Suslin._ecp_normalize_variable_order(R, variable_order)...)
    @test result.max_shift_power == max_shift_power
    @test result.stage === stage
    @test stage.variable_order == result.variable_order
    @test stage.source_variable == result.source_variable
    @test stage.target_variable == result.target_variable
    @test stage.shift_power == result.shift_power
    @test stage.shift_sign == result.shift_sign
    @test stage.shift_polynomial == result.shift_polynomial
    @test stage.selected_monic_index == result.selected_monic_index
    @test stage.selected_monic_entry == result.selected_monic_entry
    @test Suslin._is_monic_in_variable(stage.selected_monic_entry, R, stage.target_variable_index)
    @test _search_apply_factors(result.factors, column, R) == _search_target_column(R, length(column))
    @test _search_apply_factors(stage.factors, column, R) == _search_target_column(R, length(column))
    cert = _search_certificate_from_stage(column, R, stage)
    @test Suslin.verify_ecp_column_reduction(cert)
    @test any(cert_stage -> cert_stage.kind == :monicity_normalization, cert.stages)
    public_cert = Suslin.ecp_column_reduction_certificate(column, R)
    @test Suslin.verify_ecp_column_reduction(public_cert)
    @test any(cert_stage -> cert_stage.kind == :ecp_pipeline, public_cert.stages)
    return result
end

function _target_x_fixture()
    R, (x, y) = Oscar.polynomial_ring(GF(2), ["x", "y"])
    return (;
        id = "ecp-variable-change-target-x-gf2",
        ring = (; object = R),
        variable_order = (:y, :x),
        entries = (
            a = y + x^2,
            b = x * y + y + one(R),
            c = y^2 + x * y + x + one(R),
        ),
        column_order = (:a, :b, :c),
    )
end

@testset "deterministic ECP monicity search" begin
    cases = ECPColumnFixtureCatalog.cases_by_id()

    old_bounded = _assert_success_reduces(cases["ecp-variable-change-monic-gf2"])
    @test old_bounded.source_variable == gens(old_bounded.ring)[1]
    @test old_bounded.target_variable == gens(old_bounded.ring)[2]
    @test old_bounded.shift_power == 2
    @test Suslin._ecp_first_monic_entry_index([zero(old_bounded.ring), old_bounded.selected_monic_entry], old_bounded.ring) == 2
    @test Suslin._leading_coefficient_in_last_variable(old_bounded.selected_monic_entry, old_bounded.ring) == one(old_bounded.ring)

    target_x_fixture = _target_x_fixture()
    target_x = _assert_success_reduces(
        target_x_fixture;
        variable_order = reverse(tuple(gens(target_x_fixture.ring.object)...)),
    )
    @test target_x.variable_order == reverse(tuple(gens(target_x.ring)...))
    @test target_x.source_variable == gens(target_x.ring)[2]
    @test target_x.target_variable == gens(target_x.ring)[1]
    @test target_x.shift_power == 2
    @test target_x.stage.target_variable_index == 1
    @test target_x.stage.source_variable_index == 2
    target_x_cert = _search_certificate_from_stage(_search_column(target_x_fixture), target_x.ring, target_x.stage)
    @test target_x_cert.verification.overall_ok == true
    @test Suslin.verify_ecp_column_reduction(target_x_cert)

    exhausted_entry = cases["ecp-unsupported-unimodular-gf2"]
    exhausted_column = _search_column(exhausted_entry)
    exhausted = Suslin._deterministic_ecp_monicity_search(
        exhausted_column,
        exhausted_entry.ring.object;
        variable_order = tuple(gens(exhausted_entry.ring.object)...),
        max_shift_power = 3,
    )
    @test exhausted isa Suslin.ECPMonicitySearchFailure
    @test exhausted.kind == :monicity_search_exhausted
    @test exhausted.variable_order == tuple(gens(exhausted_entry.ring.object)...)
    @test exhausted.max_shift_power == 3
    @test exhausted.source_variables == (gens(exhausted_entry.ring.object)[1],)
    @test exhausted.target_variable == gens(exhausted_entry.ring.object)[2]
    @test exhausted.shift_powers == (1, 2, 3)
    @test exhausted.shift_polynomials == (
        gens(exhausted_entry.ring.object)[2],
        gens(exhausted_entry.ring.object)[2],
        gens(exhausted_entry.ring.object)[2]^2,
        gens(exhausted_entry.ring.object)[2]^2,
        gens(exhausted_entry.ring.object)[2]^3,
        gens(exhausted_entry.ring.object)[2]^3,
    )
    @test exhausted.attempted_candidates == 6
    @test occursin("exhausted deterministic ECP monicity search", exhausted.message)

    low_bound = Suslin._deterministic_ecp_monicity_search(
        _search_column(cases["ecp-variable-change-monic-gf2"]),
        cases["ecp-variable-change-monic-gf2"].ring.object;
        variable_order = tuple(gens(cases["ecp-variable-change-monic-gf2"].ring.object)...),
        max_shift_power = 1,
    )
    @test low_bound isa Suslin.ECPMonicitySearchFailure
    @test low_bound.max_shift_power == 1

    missing_required_source = Suslin._deterministic_ecp_monicity_search(
        _search_column(cases["ecp-variable-change-monic-gf2"]),
        cases["ecp-variable-change-monic-gf2"].ring.object;
        variable_order = (last(gens(cases["ecp-variable-change-monic-gf2"].ring.object)),),
        max_shift_power = 3,
    )
    @test missing_required_source isa Suslin.ECPMonicitySearchFailure
    @test missing_required_source.variable_order == (last(gens(cases["ecp-variable-change-monic-gf2"].ring.object)),)
    @test missing_required_source.attempted_candidates == 0

    symbol_order = _assert_success_reduces(
        cases["ecp-variable-change-monic-gf2"];
        variable_order = (:x, :y),
    )
    @test symbol_order.variable_order == tuple(gens(symbol_order.ring)...)

    @test_throws ArgumentError Suslin._deterministic_ecp_monicity_search(
        _search_column(cases["ecp-variable-change-monic-gf2"]),
        cases["ecp-variable-change-monic-gf2"].ring.object;
        variable_order = (:x, :x),
    )
    @test_throws ArgumentError Suslin._deterministic_ecp_monicity_search(
        _search_column(cases["ecp-variable-change-monic-gf2"]),
        cases["ecp-variable-change-monic-gf2"].ring.object;
        variable_order = (:x, :z),
    )
end

@testset "monicity replay rejects malformed stage metadata" begin
    cases = ECPColumnFixtureCatalog.cases_by_id()
    entry = cases["ecp-variable-change-monic-gf2"]
    column = _search_column(entry)
    result = Suslin._deterministic_ecp_monicity_search(column, entry.ring.object)
    cert = _search_certificate_from_stage(column, entry.ring.object, result.stage)
    stage = _search_monicity_stage(cert)

    bad_certs = (
        _search_tamper_stage_field(cert, :source_variable_index, length(gens(cert.ring)) + 1),
        _search_tamper_stage_field(cert, :target_variable_index, 0),
        _search_tamper_stage_field(cert, :source_variable, gens(cert.ring)[end]),
        _search_tamper_stage_field(cert, :target_variable, gens(cert.ring)[1]),
        _search_tamper_stage_field(cert, :shift_power, -1),
        _search_tamper_stage_field(cert, :variable_order, nothing),
        _search_tamper_stage_field(cert, :variable_order, :x),
        _search_tamper_stage_field(cert, :variable_order, (stage.variable_order[1], stage.variable_order[1])),
    )

    for bad_cert in bad_certs
        @test Suslin._ecp_replay_stages(bad_cert).ok == false
        @test Suslin.verify_ecp_column_reduction(bad_cert) == false
    end
end
