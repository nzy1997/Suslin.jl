using Test
using Suslin
using Oscar

include("../fixtures/laurent_lazy_determinant_cases.jl")

function _issue158_fixture(id::AbstractString)
    catalog = LaurentLazyDeterminantCases.catalog()
    return only(filter(entry -> entry.id == id, catalog.cases))
end

function _issue158_product(factors, R, n::Int)
    product = identity_matrix(R, n)
    for factor in factors
        product *= factor
    end
    return product
end

function _issue158_wrong_unrewritten_certificate(certificate)
    R = base_ring(certificate.original_matrix)
    n = nrows(certificate.original_matrix)
    deferred_n = nrows(certificate.normalized_deferred_core)
    peel_certificate = certificate.deferred_metadata.peel_certificate
    unrewritten_left = Suslin._inverse_elementary_sequence(peel_certificate.left_factors)
    embedded_core_factors = Suslin._embed_laurent_deferred_peel_factors(
        certificate.normalized_deferred_factors,
        R,
        n,
        deferred_n,
    )
    right_inverse = Suslin._inverse_elementary_sequence(peel_certificate.right_factors)
    wrong_factors = vcat(unrewritten_left, embedded_core_factors, right_inverse)
    wrong_product = _issue158_product(wrong_factors, R, n)

    return Suslin.LaurentLazyGLHoistCertificate(
        certificate.original_matrix,
        certificate.deferred_metadata,
        certificate.overall_determinant,
        certificate.correction,
        certificate.inverse_correction,
        certificate.normalized_deferred_core,
        certificate.normalized_deferred_factorization,
        certificate.normalized_deferred_factors,
        unrewritten_left,
        wrong_factors,
        wrong_product,
        certificate.correction.factor * wrong_product,
        certificate.verification,
    )
end

@testset "left elementary factors rewrite across Laurent diagonal correction" begin
    R, (x, y) = suslin_laurent_polynomial_ring(GF(2), ["x", "y"])
    diagonal = diagonal_matrix(R, [x, y, one(R)])
    coefficient = x + y
    factor = elementary_matrix(3, 2, 1, coefficient, R)

    rewritten = Suslin._rewrite_left_elementary_factor_across_diagonal(factor, diagonal)
    row, col, rewritten_coefficient = Suslin._elementary_factor_data(rewritten)

    @test (row, col) == (2, 1)
    @test rewritten_coefficient == y^-1 * coefficient * x
    @test factor * diagonal == diagonal * rewritten
end

@testset "lazy deferred Laurent correction hoists to original GL certificate" begin
    entry = _issue158_fixture("monomial-unit-row-column-cores")
    A = entry.inputs.matrix
    metadata = Suslin._factor_laurent_gl_lazy_determinant_peel(A)

    certificate = Suslin._laurent_gl_lazy_deferred_correction_certificate(metadata)

    @test certificate.original_matrix == A
    @test certificate.deferred_metadata == metadata
    @test certificate.overall_determinant == metadata.overall_determinant
    @test certificate.overall_determinant == entry.determinant_profile.expected_determinant
    @test certificate.correction.scope == :original_matrix
    @test certificate.correction.side == :left
    @test certificate.correction.factor * certificate.inverse_correction ==
        identity_matrix(base_ring(A), nrows(A))
    @test certificate.inverse_correction * certificate.correction.factor ==
        identity_matrix(base_ring(A), nrows(A))
    @test det(certificate.correction.factor) == certificate.overall_determinant
    @test certificate.normalized_deferred_core == metadata.normalized_deferred_core
    @test Suslin._verify_laurent_column_peel_replay(certificate.normalized_deferred_factorization)
    @test verify_factorization(
        certificate.normalized_deferred_core,
        certificate.normalized_deferred_factors,
    )

    expected_product = _issue158_product(
        certificate.elementary_factors,
        base_ring(A),
        nrows(A),
    )
    @test certificate.elementary_product == expected_product
    @test certificate.correction.factor * certificate.elementary_product == A
    @test certificate.reconstructed_product == A
    @test certificate.verification.overall_ok
    @test Suslin._verify_laurent_gl_lazy_deferred_correction_certificate(certificate)

    wrong = _issue158_wrong_unrewritten_certificate(certificate)
    @test !Suslin._verify_laurent_gl_lazy_deferred_correction_certificate(wrong)
    @test !Suslin._laurent_gl_lazy_deferred_correction_certificate_verification(wrong).rewritten_left_factors_ok
end
