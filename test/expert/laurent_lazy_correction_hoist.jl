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

function _issue158_noncommuting_deferred_metadata()
    R, (x, y) = suslin_laurent_polynomial_ring(GF(2), ["x", "y"])
    deferred = matrix(R, [
        x one(R);
        zero(R) one(R)
    ])
    target = block_embedding(deferred, 3, [1, 2])
    coefficient = y + one(R)
    left_factor = elementary_matrix(3, 1, 3, coefficient, R)
    input = elementary_matrix(3, 1, 3, -coefficient, R) * target
    right_factors = Suslin._expected_column_peel_right_factors(target, 3, R)
    step = Suslin.LaurentColumnPeelStep(
        3,
        input,
        [input[row, 3] for row in 1:3],
        typeof(left_factor)[left_factor],
        target,
        right_factors,
        target,
        deferred,
    )
    peel_certificate = Suslin.LaurentDeterminantDeferredPeelCertificate(
        input,
        Suslin.LaurentColumnPeelStep[step],
        deferred,
        :deferred_submatrix,
        nothing,
    )
    return Suslin._normalize_laurent_determinant_deferred_submatrix(peel_certificate)
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
    direct_certificate = Suslin._laurent_gl_lazy_deferred_correction_certificate(A)

    @test certificate.original_matrix == A
    @test certificate.deferred_metadata == metadata
    @test direct_certificate.original_matrix == A
    @test direct_certificate.reconstructed_product == A
    @test direct_certificate.verification.overall_ok
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

    left_inverse = Suslin._inverse_elementary_sequence(metadata.peel_certificate.left_factors)
    expected_rewritten_left = [
        Suslin._rewrite_left_elementary_factor_across_diagonal(factor, certificate.correction.factor)
        for factor in left_inverse
    ]
    @test certificate.rewritten_left_factors == expected_rewritten_left

    noncommuting_metadata = _issue158_noncommuting_deferred_metadata()
    noncommuting_certificate = Suslin._laurent_gl_lazy_deferred_correction_certificate(
        noncommuting_metadata,
    )
    @test Suslin._verify_laurent_gl_lazy_deferred_correction_certificate(noncommuting_certificate)
    wrong = _issue158_wrong_unrewritten_certificate(noncommuting_certificate)
    @test wrong.reconstructed_product != wrong.original_matrix
    @test !Suslin._verify_laurent_gl_lazy_deferred_correction_certificate(wrong)
    @test !Suslin._laurent_gl_lazy_deferred_correction_certificate_verification(wrong).rewritten_left_factors_ok

    outer_broken = Suslin._laurent_gl_lazy_deferred_correction_certificate_verification((
        original_matrix = nothing,
    ))
    @test !outer_broken.size_ok
    @test !outer_broken.overall_ok

    inner_broken = Suslin.LaurentLazyGLHoistCertificate(
        certificate.original_matrix,
        (;),
        certificate.overall_determinant,
        certificate.correction,
        certificate.inverse_correction,
        certificate.normalized_deferred_core,
        certificate.normalized_deferred_factorization,
        certificate.normalized_deferred_factors,
        certificate.rewritten_left_factors,
        certificate.elementary_factors,
        certificate.elementary_product,
        certificate.reconstructed_product,
        certificate.verification,
    )
    inner_verification = Suslin._laurent_gl_lazy_deferred_correction_certificate_verification(
        inner_broken,
    )
    @test inner_verification.size_ok
    @test !inner_verification.metadata_ok
    @test !inner_verification.overall_ok
end
