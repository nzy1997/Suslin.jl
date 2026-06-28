using Test
using Suslin
using Oscar

include("../fixtures/laurent_lazy_determinant_cases.jl")

function _issue159_fixture(id::AbstractString)
    catalog = LaurentLazyDeterminantCases.catalog()
    return only(filter(entry -> entry.id == id, catalog.cases))
end

function _issue159_product(factors, R, n::Int)
    product = identity_matrix(R, n)
    for factor in factors
        product *= factor
    end
    return product
end

function _issue159_rebuild(
    certificate;
    correction_side = certificate.correction_side,
    reconstruction_relation = certificate.reconstruction_relation,
    correction = certificate.correction,
    inverse_correction = certificate.inverse_correction,
    normalized_deferred_core = certificate.normalized_deferred_core,
    normalized_deferred_factorization = certificate.normalized_deferred_factorization,
    normalized_deferred_factors = certificate.normalized_deferred_factors,
    rewritten_left_factors = certificate.rewritten_left_factors,
    rewritten_right_factors = certificate.rewritten_right_factors,
    elementary_factors = certificate.elementary_factors,
    elementary_product = certificate.elementary_product,
    reconstructed_product = certificate.reconstructed_product,
    verification = certificate.verification,
)
    return Suslin.LaurentLazyGLHoistCertificate(
        certificate.original_matrix,
        certificate.deferred_metadata,
        certificate.overall_determinant,
        correction_side,
        reconstruction_relation,
        correction,
        inverse_correction,
        normalized_deferred_core,
        normalized_deferred_factorization,
        normalized_deferred_factors,
        rewritten_left_factors,
        rewritten_right_factors,
        elementary_factors,
        elementary_product,
        reconstructed_product,
        verification,
    )
end

@testset "lazy Laurent determinant correction side choices" begin
    entry = _issue159_fixture("issue-38-q-block-lazy-determinant")
    A = entry.inputs.matrix
    R = base_ring(A)
    n = nrows(A)

    row_certificate = Suslin._laurent_gl_lazy_deferred_correction_certificate(
        A;
        correction_side = :row,
    )
    column_certificate = Suslin._laurent_gl_lazy_deferred_correction_certificate(
        A;
        correction_side = :column,
    )

    @test row_certificate.correction_side == :row
    @test column_certificate.correction_side == :column
    @test row_certificate.correction.side == :left
    @test column_certificate.correction.side == :right
    @test row_certificate.reconstruction_relation == :left_correction_times_elementary_product
    @test column_certificate.reconstruction_relation == :elementary_product_times_right_correction
    @test row_certificate.overall_determinant == entry.determinant_profile.expected_determinant
    @test column_certificate.overall_determinant == entry.determinant_profile.expected_determinant
    @test row_certificate.overall_determinant == column_certificate.overall_determinant

    @test row_certificate.elementary_product ==
        _issue159_product(row_certificate.elementary_factors, R, n)
    @test column_certificate.elementary_product ==
        _issue159_product(column_certificate.elementary_factors, R, n)
    @test row_certificate.correction.factor * row_certificate.elementary_product == A
    @test column_certificate.elementary_product * column_certificate.correction.factor == A
    @test row_certificate.reconstructed_product == A
    @test column_certificate.reconstructed_product == A
    @test row_certificate.verification.overall_ok
    @test column_certificate.verification.overall_ok
    @test Suslin._verify_laurent_gl_lazy_deferred_correction_certificate(row_certificate)
    @test Suslin._verify_laurent_gl_lazy_deferred_correction_certificate(column_certificate)

    @test det(row_certificate.correction.factor) == row_certificate.overall_determinant
    @test det(column_certificate.correction.factor) == column_certificate.overall_determinant
    @test row_certificate.correction.factor * row_certificate.inverse_correction ==
        identity_matrix(R, n)
    @test column_certificate.correction.factor * column_certificate.inverse_correction ==
        identity_matrix(R, n)

    wrong_metadata = _issue159_rebuild(
        column_certificate;
        correction_side = :row,
        reconstruction_relation = :left_correction_times_elementary_product,
    )
    wrong_verification = Suslin._laurent_gl_lazy_deferred_correction_certificate_verification(
        wrong_metadata,
    )
    @test !wrong_verification.overall_ok
    @test !wrong_verification.correction_side_ok
    @test !wrong_verification.correction_ok
    @test !wrong_verification.rewritten_left_factors_ok
    @test !Suslin._verify_laurent_gl_lazy_deferred_correction_certificate(wrong_metadata)

    err = try
        Suslin._laurent_gl_lazy_deferred_correction_certificate(
            A;
            correction_side = :diagonal,
        )
        nothing
    catch caught
        caught
    end
    @test err isa ArgumentError
    @test occursin(":row", sprint(showerror, err))
    @test occursin(":column", sprint(showerror, err))
end
