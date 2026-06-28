using Test
using Suslin
using Oscar

include(joinpath(@__DIR__, "..", "fixtures", "toricbuilder_case008_d16_column_boundary.jl"))

function _case008_d16_reduction_product(factors, R, n::Int)
    product = identity_matrix(R, n)
    for factor in factors
        product *= factor
    end
    return product
end

function _case008_d16_apply_factors(factors, column, R)
    n = length(column)
    return _case008_d16_reduction_product(factors, R, n) * matrix(R, n, 1, collect(column))
end

function _case008_d16_target_column(R, n::Int)
    target = zero_matrix(R, n, 1)
    target[n, 1] = one(R)
    return target
end

function _case008_d16_tamper_first_factor(factors, R, n::Int)
    tampered = copy(factors)
    tampered[1] = identity_matrix(R, n)
    return tampered
end

function _case008_d16_tamper_certificate_first_factor(cert)
    tampered = _case008_d16_tamper_first_factor(
        cert.factors,
        cert.ring,
        length(cert.original_column),
    )
    return Suslin.ECPColumnReductionCertificate(
        cert.original_column,
        cert.ring,
        cert.stages,
        tampered,
        cert.final_column,
        cert.verification,
    )
end

function _case008_d16_tamper_stage_coefficient(cert)
    stages = collect(cert.stages)
    stage_idx = findfirst(stage -> stage.kind == :laurent_elementary_row_preconditioning, stages)
    stage_idx === nothing && error("missing row-preconditioning stage")
    stage = stages[stage_idx]
    stages[stage_idx] = merge(stage, (; coefficient = zero(cert.ring)))
    return Suslin.ECPColumnReductionCertificate(
        cert.original_column,
        cert.ring,
        tuple(stages...),
        cert.factors,
        cert.final_column,
        cert.verification,
    )
end

function _case008_d16_diagnostic_stage_detail(diagnostic, stage::Symbol)
    hasproperty(diagnostic, :stage_details) || return nothing
    idx = findfirst(detail -> detail.stage == stage, diagnostic.stage_details)
    return idx === nothing ? nothing : diagnostic.stage_details[idx]
end

@testset "case_008 d=16 Laurent column reduction" begin
    fixture = ToricBuilderCase008D16ColumnBoundary.boundary_fixture()
    R = fixture.ring
    column = fixture.failing_column
    target = _case008_d16_target_column(R, length(column))

    factors = Suslin.reduce_unimodular_column(column, R)
    @test _case008_d16_apply_factors(factors, column, R) == target

    tampered_factors = _case008_d16_tamper_first_factor(factors, R, length(column))
    @test _case008_d16_apply_factors(tampered_factors, column, R) != target

    certificate = Suslin.ecp_column_reduction_certificate(column, R)
    @test Suslin.verify_ecp_column_reduction(certificate)
    @test any(stage -> stage.kind == :laurent_elementary_row_preconditioning, certificate.stages)
    @test !any(stage -> stage.kind == :case008_special_case, certificate.stages)
    @test _case008_d16_apply_factors(
        certificate.factors,
        certificate.original_column,
        certificate.ring,
    ) == target
    @test !Suslin.verify_ecp_column_reduction(
        _case008_d16_tamper_certificate_first_factor(certificate),
    )
    @test !Suslin.verify_ecp_column_reduction(
        _case008_d16_tamper_stage_coefficient(certificate),
    )

    diagnostic = Suslin.diagnose_unimodular_column_reduction(column, R)
    @test diagnostic.status == :supported
    @test diagnostic.failure_code === nothing
    @test diagnostic.column_length == 16
    @test diagnostic.ring_profile.kind == :laurent_polynomial
    @test diagnostic.ring_profile.generators == ("u", "v")
    @test :laurent_elementary_row_preconditioning in diagnostic.attempted_stages
    @test !(:case008_special_case in diagnostic.attempted_stages)
    @test _case008_d16_diagnostic_stage_detail(diagnostic, :case008_special_case) === nothing

    negative = ToricBuilderCase008D16ColumnBoundary.non_unimodular_negative_control(fixture)
    negative_diagnostic = Suslin.diagnose_unimodular_column_reduction(
        negative.failing_column,
        negative.ring,
    )
    @test negative_diagnostic.status == :precondition_failed
    @test negative_diagnostic.failure_code == :not_unimodular
    @test isempty(negative_diagnostic.attempted_stages)
    @test_throws ArgumentError Suslin.reduce_unimodular_column(
        negative.failing_column,
        negative.ring,
    )
end
