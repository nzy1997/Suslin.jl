using Test
using Suslin
using Oscar

const PARK_WOODBURN_SLN_PEEL_STEP_CATALOG_PATH =
    joinpath(@__DIR__, "..", "fixtures", "park_woodburn_sln_driver_cases.jl")

function _sln_step_product(factors, R, n::Int)
    product = identity_matrix(R, n)
    for factor in factors
        product *= factor
    end
    return product
end

function _sln_step_target_column(R, n::Int)
    target = zero_matrix(R, n, 1)
    target[n, 1] = one(R)
    return target
end

function _sln_step_replace(step; kwargs...)
    values = Dict{Symbol, Any}(name => getfield(step, name) for name in fieldnames(typeof(step)))
    for (name, value) in kwargs
        values[name] = value
    end
    return Suslin.PolynomialColumnPeelStep((values[name] for name in fieldnames(typeof(step)))...)
end

function _sln_step_tamper_ecp_certificate(cert)
    factors = copy(cert.factors)
    factors[1] = identity_matrix(cert.ring, length(cert.original_column))
    return Suslin.ECPColumnReductionCertificate(
        cert.original_column,
        cert.ring,
        cert.stages,
        factors,
        cert.final_column,
        cert.verification,
    )
end

function _sln_step_assert_rejected(step; kwargs...)
    @test !Suslin._polynomial_column_peel_step_verification(_sln_step_replace(step; kwargs...)).overall_ok
end

@testset "Park-Woodburn SLn ECP-backed polynomial peel step" begin
    if !isdefined(Main, :ParkWoodburnSLnDriverFixtureCatalog)
        include(PARK_WOODBURN_SLN_PEEL_STEP_CATALOG_PATH)
    end
    entries = ParkWoodburnSLnDriverFixtureCatalog.cases_by_id()
    entry = entries["sln-driver-sl4-gf2-ecp-mainline"]

    step = Suslin._polynomial_column_peel_step(entry.matrix)
    R = base_ring(step.input_matrix)
    left_product = _sln_step_product(step.left_factors, R, step.dimension)
    right_product = _sln_step_product(step.right_factors, R, step.dimension)
    recorded_column = matrix(R, step.dimension, 1, step.last_column)

    @test step.dimension == 4
    @test step.input_matrix == entry.matrix
    @test step.ecp_evidence isa Suslin.ECPColumnReductionCertificate
    @test step.left_certificate == step.ecp_evidence
    @test Suslin.verify_ecp_column_reduction(step.ecp_evidence)
    @test step.ecp_evidence.original_column == step.last_column
    @test step.ecp_evidence.factors == step.left_factors
    @test step.ecp_evidence.final_column == _sln_step_target_column(R, step.dimension)
    @test step.ecp_route_provenance.verifier == :verify_ecp_column_reduction
    @test step.ecp_route_provenance.status == :verified
    @test step.ecp_route_provenance.route == :embedded_three_block
    @test step.ecp_route_provenance.factor_count == length(step.left_factors)
    @test left_product * recorded_column == _sln_step_target_column(R, step.dimension)
    @test left_product * step.input_matrix == step.after_left_matrix
    @test step.after_left_matrix * right_product == step.peeled_matrix
    @test step.right_clearing_coefficients ==
          tuple((step.after_left_matrix[step.dimension, col] for col in 1:(step.dimension - 1))...)
    @test step.right_factors == [
        elementary_matrix(step.dimension, step.dimension, col, -step.right_clearing_coefficients[col], R)
        for col in 1:(step.dimension - 1)
        if step.right_clearing_coefficients[col] != zero(R)
    ]
    @test step.block_embedding_indices == collect(1:(step.dimension - 1))
    @test step.peeled_matrix == block_embedding(step.next_block, step.dimension, step.block_embedding_indices)
    @test det(step.next_block) == one(R)
    @test step.determinant_metadata.input_determinant == one(R)
    @test step.determinant_metadata.peeled_determinant == one(R)
    @test step.determinant_metadata.next_block_determinant == one(R)
    @test step.descent_metadata.input_dimension == step.dimension
    @test step.descent_metadata.next_dimension == step.dimension - 1
    @test step.descent_metadata.dimension_drop == 1
    @test Suslin._polynomial_column_peel_step_verification(step).overall_ok

    bad_column = copy(step.last_column)
    bad_column[1] += one(R)
    _sln_step_assert_rejected(step; last_column = bad_column)

    bad_left = copy(step.left_factors)
    bad_left[1] = identity_matrix(R, step.dimension)
    _sln_step_assert_rejected(step; left_factors = bad_left)

    bad_right_coefficients = collect(step.right_clearing_coefficients)
    bad_right_coefficients[1] += one(R)
    _sln_step_assert_rejected(step; right_clearing_coefficients = tuple(bad_right_coefficients...))

    bad_block = copy(step.next_block)
    bad_block[1, 1] += one(R)
    _sln_step_assert_rejected(step; next_block = bad_block)

    bad_provenance = merge(step.ecp_route_provenance, (; route = :tampered_route))
    _sln_step_assert_rejected(step; ecp_route_provenance = bad_provenance)

    bad_ecp = _sln_step_tamper_ecp_certificate(step.ecp_evidence)
    _sln_step_assert_rejected(step; ecp_evidence = bad_ecp, left_certificate = bad_ecp)
end
