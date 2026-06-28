using Test
using Suslin
using Oscar

include(joinpath(@__DIR__, "..", "fixtures", "toricbuilder_case008_d21_column_boundary.jl"))

function _case008_d21_reduction_product(factors, R, n::Int)
    product = identity_matrix(R, n)
    for factor in factors
        product *= factor
    end
    return product
end

function _case008_d21_apply_factors(factors, column, R)
    n = length(column)
    return _case008_d21_reduction_product(factors, R, n) * matrix(R, n, 1, collect(column))
end

function _case008_d21_target_column(R, n::Int)
    target = zero_matrix(R, n, 1)
    target[n, 1] = one(R)
    return target
end

function _case008_d21_tamper_first_factor(factors, R, n::Int)
    tampered = copy(factors)
    tampered[1] = identity_matrix(R, n)
    return tampered
end

function _case008_d21_tamper_certificate_first_factor(cert)
    tampered = _case008_d21_tamper_first_factor(
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

@testset "case_008 d=21 Laurent column reduction" begin
    fixture = ToricBuilderCase008D21ColumnBoundary.boundary_fixture()
    R = fixture.ring
    column = fixture.failing_column
    target = _case008_d21_target_column(R, length(column))

    factors = Suslin.reduce_unimodular_column(column, R)
    @test _case008_d21_apply_factors(factors, column, R) == target

    tampered_factors = _case008_d21_tamper_first_factor(factors, R, length(column))
    @test _case008_d21_apply_factors(tampered_factors, column, R) != target

    certificate = Suslin.ecp_column_reduction_certificate(column, R)
    @test Suslin.verify_ecp_column_reduction(certificate)
    @test any(stage -> stage.kind == :witness_unit, certificate.stages)
    @test _case008_d21_apply_factors(
        certificate.factors,
        certificate.original_column,
        certificate.ring,
    ) == target
    @test !Suslin.verify_ecp_column_reduction(
        _case008_d21_tamper_certificate_first_factor(certificate),
    )

    diagnostic = Suslin.diagnose_unimodular_column_reduction(column, R)
    @test diagnostic.status == :supported
    @test diagnostic.failure_code === nothing
    @test diagnostic.column_length == 21
    @test diagnostic.ring_profile.kind == :laurent_polynomial
    @test diagnostic.ring_profile.generators == ("u", "v")
    @test :laurent_witness_unit in diagnostic.attempted_stages

    unexpected_solver_failure = (A, B) -> throw(ArgumentError("unexpected Laurent solver failure"))
    @test_throws ArgumentError Suslin._laurent_unimodular_witness(
        column,
        R;
        solver = unexpected_solver_failure,
    )
end
