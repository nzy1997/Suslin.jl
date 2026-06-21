using Test
using Suslin
using Oscar

include("../fixtures/toricbuilder_issue38_cases.jl")

function _issue59_product(factors, R, n::Int)
    product = identity_matrix(R, n)
    for factor in factors
        product *= factor
    end
    return product
end

function _issue59_tamper_correction(certificate)
    R = base_ring(certificate.original_matrix)
    bad_factor = copy(certificate.correction.factor)
    bad_factor[1, 1] = one(R)
    bad_correction = merge(certificate.correction, (; factor = bad_factor))
    return Suslin.LaurentGLFactorizationCertificate(
        certificate.original_matrix,
        certificate.determinant_profile,
        certificate.normalization,
        bad_correction,
        certificate.inverse_correction,
        certificate.normalized_core,
        certificate.core_factorization,
        certificate.core_factors,
        certificate.reconstructed_product,
        certificate.verification,
    )
end

function _issue59_tamper_core_factor(certificate)
    R = base_ring(certificate.original_matrix)
    n = nrows(certificate.original_matrix)
    bad_factors = copy(certificate.core_factors)
    bad_factors[1] = identity_matrix(R, n)
    return Suslin.LaurentGLFactorizationCertificate(
        certificate.original_matrix,
        certificate.determinant_profile,
        certificate.normalization,
        certificate.correction,
        certificate.inverse_correction,
        certificate.normalized_core,
        certificate.core_factorization,
        bad_factors,
        certificate.reconstructed_product,
        certificate.verification,
    )
end

function _issue59_tamper_reconstructed_product(certificate)
    R = base_ring(certificate.original_matrix)
    bad_product = copy(certificate.reconstructed_product)
    bad_product[1, 1] += one(R)
    return Suslin.LaurentGLFactorizationCertificate(
        certificate.original_matrix,
        certificate.determinant_profile,
        certificate.normalization,
        certificate.correction,
        certificate.inverse_correction,
        certificate.normalized_core,
        certificate.core_factorization,
        certificate.core_factors,
        bad_product,
        certificate.verification,
    )
end

@testset "Issue 38 Laurent GL certificate" begin
    entry = only(ToricBuilderIssue38Cases.catalog().cases)
    Q = entry.inputs.matrix
    R = entry.ring.object
    u, v = entry.ring.generators

    determinant_profile = classify_laurent_determinant(Q)
    @test determinant_profile.classification == :laurent_monomial_unit
    @test determinant_profile.determinant == u * v

    certificate = laurent_gl_factorization_certificate(Q)
    @test certificate.original_matrix == Q
    @test certificate.determinant_profile.classification == :laurent_monomial_unit
    @test certificate.determinant_profile.determinant == u * v
    @test certificate.correction == certificate.normalization.correction
    @test certificate.inverse_correction == certificate.correction.inverse_factor
    @test certificate.normalized_core == certificate.normalization.normalized_matrix
    @test det(certificate.normalized_core) == one(R)
    @test length(certificate.core_factors) > 0
    @test verify_factorization(certificate.normalized_core, certificate.core_factors)
    @test !verify_factorization(Q, certificate.core_factors)

    core_product = _issue59_product(certificate.core_factors, R, nrows(Q))
    @test core_product == certificate.normalized_core
    @test certificate.correction.factor * core_product == Q
    @test certificate.reconstructed_product == Q
    @test verify_laurent_gl_factorization_certificate(certificate)

    tampered_correction = _issue59_tamper_correction(certificate)
    @test !verify_laurent_gl_factorization_certificate(tampered_correction)

    tampered_factor = _issue59_tamper_core_factor(certificate)
    @test !verify_laurent_gl_factorization_certificate(tampered_factor)

    tampered_product = _issue59_tamper_reconstructed_product(certificate)
    @test !verify_laurent_gl_factorization_certificate(tampered_product)

    original_err = try
        elementary_factorization(Q)
        nothing
    catch err
        err
    end
    @test original_err isa ArgumentError
    @test occursin("Laurent GL_n normalization boundary", sprint(showerror, original_err))
end
