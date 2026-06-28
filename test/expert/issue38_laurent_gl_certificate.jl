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

function _issue59_malformed_original_matrix_certificate()
    return Suslin.LaurentGLFactorizationCertificate(
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        Any[],
        nothing,
        nothing,
    )
end

function _issue59_malformed_normalization_certificate(certificate)
    return Suslin.LaurentGLFactorizationCertificate(
        certificate.original_matrix,
        certificate.determinant_profile,
        nothing,
        certificate.correction,
        certificate.inverse_correction,
        certificate.normalized_core,
        certificate.core_factorization,
        certificate.core_factors,
        certificate.reconstructed_product,
        certificate.verification,
    )
end

function _issue162_rebuild_lazy_certificate(
    certificate;
    overall_determinant = certificate.overall_determinant,
    determinant_source = certificate.determinant_source,
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
        overall_determinant,
        determinant_source,
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

function _issue162_tamper_lazy_correction(certificate)
    R = base_ring(certificate.original_matrix)
    bad_factor = copy(certificate.correction.factor)
    bad_factor[1, 1] = one(R)
    bad_correction = merge(certificate.correction, (; factor = bad_factor))
    return _issue162_rebuild_lazy_certificate(
        certificate;
        correction = bad_correction,
    )
end

function _issue162_tamper_lazy_correction_side(certificate)
    return _issue162_rebuild_lazy_certificate(
        certificate;
        correction_side = :column,
    )
end

function _issue162_tamper_lazy_hoisted_factor(certificate)
    R = base_ring(certificate.original_matrix)
    n = nrows(certificate.original_matrix)
    bad_factors = copy(certificate.elementary_factors)
    bad_factors[1] = identity_matrix(R, n)
    bad_product = _issue59_product(bad_factors, R, n)
    return _issue162_rebuild_lazy_certificate(
        certificate;
        elementary_factors = bad_factors,
        elementary_product = bad_product,
        reconstructed_product = certificate.correction.factor * bad_product,
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
    callback_records = Any[]
    callback_certificate = Suslin._laurent_gl_factorization_certificate(
        Q;
        progress_callback = record -> push!(callback_records, record),
    )
    @test verify_laurent_gl_factorization_certificate(callback_certificate)
    @test callback_certificate.reconstructed_product == certificate.reconstructed_product
    @test callback_certificate.normalized_core == certificate.normalized_core
    @test !isempty(callback_records)
    @test first(callback_records).current_dimension == nrows(callback_certificate.normalized_core)
    @test last(callback_records).current_dimension == 2

    normalized = normalize_laurent_gl_matrix(Q)
    from_normalization_records = Any[]
    from_normalization_certificate = Suslin._laurent_gl_factorization_certificate_from_normalization(
        Q,
        normalized;
        progress_callback = record -> push!(from_normalization_records, record),
    )
    @test verify_laurent_gl_factorization_certificate(from_normalization_certificate)
    @test from_normalization_certificate.normalization == normalized
    @test from_normalization_certificate.reconstructed_product == certificate.reconstructed_product
    @test !isempty(from_normalization_records)
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

    malformed_original = _issue59_malformed_original_matrix_certificate()
    @test !verify_laurent_gl_factorization_certificate(malformed_original)

    malformed_normalization = _issue59_malformed_normalization_certificate(certificate)
    @test !verify_laurent_gl_factorization_certificate(malformed_normalization)

    lazy_certificate = laurent_gl_factorization_certificate(
        Q;
        determinant_strategy = :lazy,
        correction_side = :row,
    )
    @test lazy_certificate isa LaurentLazyGLHoistCertificate
    @test lazy_certificate.original_matrix == Q
    @test lazy_certificate.overall_determinant == u * v
    @test lazy_certificate.determinant_source == :deferred_submatrix
    @test lazy_certificate.correction_side == :row
    @test lazy_certificate.reconstruction_relation ==
        :left_correction_times_elementary_product
    @test length(lazy_certificate.elementary_factors) > 0
    @test lazy_certificate.elementary_product ==
        _issue59_product(lazy_certificate.elementary_factors, R, nrows(Q))
    @test lazy_certificate.correction.factor * lazy_certificate.elementary_product == Q
    @test lazy_certificate.reconstructed_product == Q
    @test lazy_certificate.verification.overall_ok
    @test verify_laurent_gl_factorization_certificate(lazy_certificate)

    @test !verify_laurent_gl_factorization_certificate(
        _issue162_tamper_lazy_correction(lazy_certificate),
    )
    @test !verify_laurent_gl_factorization_certificate(
        _issue162_tamper_lazy_correction_side(lazy_certificate),
    )
    @test !verify_laurent_gl_factorization_certificate(
        _issue162_tamper_lazy_hoisted_factor(lazy_certificate),
    )

    original_err = try
        elementary_factorization(Q)
        nothing
    catch err
        err
    end
    @test original_err isa ArgumentError
    @test occursin("Laurent GL_n normalization boundary", sprint(showerror, original_err))
end
