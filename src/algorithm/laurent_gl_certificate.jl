struct LaurentGLFactorizationCertificate
    original_matrix
    determinant_profile
    normalization
    correction
    inverse_correction
    normalized_core
    core_factorization
    core_factors::Vector
    reconstructed_product
    verification

    function LaurentGLFactorizationCertificate(
        original_matrix,
        determinant_profile,
        normalization,
        correction,
        inverse_correction,
        normalized_core,
        core_factorization,
        core_factors::Vector,
        reconstructed_product,
        verification,
    )
        R = base_ring(original_matrix)
        product = correction.factor * _factor_product(core_factors, R, nrows(original_matrix))
        return new(
            original_matrix,
            determinant_profile,
            normalization,
            correction,
            inverse_correction,
            normalized_core,
            core_factorization,
            core_factors,
            product,
            verification,
        )
    end
end

function laurent_gl_factorization_certificate(A)
    normalization = normalize_laurent_gl_matrix(A)
    core = normalization.normalized_matrix
    core_factorization = _factor_laurent_sl_column_peel(core)
    core_factors = core_factorization.factors
    R = base_ring(A)
    reconstructed_product = normalization.correction.factor * _factor_product(core_factors, R, nrows(A))
    certificate = LaurentGLFactorizationCertificate(
        A,
        normalization.determinant_profile,
        normalization,
        normalization.correction,
        normalization.correction.inverse_factor,
        core,
        core_factorization,
        core_factors,
        reconstructed_product,
        nothing,
    )
    verification = _laurent_gl_factorization_certificate_verification(certificate)
    verification.overall_ok || error("internal Laurent GL factorization certificate verification failed")
    return LaurentGLFactorizationCertificate(
        certificate.original_matrix,
        certificate.determinant_profile,
        certificate.normalization,
        certificate.correction,
        certificate.inverse_correction,
        certificate.normalized_core,
        certificate.core_factorization,
        certificate.core_factors,
        certificate.reconstructed_product,
        verification,
    )
end

function verify_laurent_gl_factorization_certificate(certificate)::Bool
    return _laurent_gl_factorization_certificate_verification(certificate).overall_ok
end

function _laurent_gl_factorization_certificate_verification(certificate)
    size_ok = try
        nrows(certificate.original_matrix) == ncols(certificate.original_matrix)
    catch err
        err isa InterruptException && rethrow()
        false
    end

    normalization_ok = false
    determinant_profile_ok = false
    correction_ok = false
    inverse_correction_ok = false
    normalized_core_ok = false
    core_det_ok = false
    core_replay_ok = false
    core_factors_match_ok = false
    core_factors_ok = false
    reconstructed_product_ok = false

    if size_ok
        try
            A = certificate.original_matrix
            R = base_ring(A)
            n = nrows(A)
            recomputed = normalize_laurent_gl_matrix(A)
            normalization_ok = verify_laurent_gl_normalization(A, certificate.normalization)
            determinant_profile_ok =
                certificate.determinant_profile == recomputed.determinant_profile &&
                certificate.normalization.determinant_profile == recomputed.determinant_profile
            correction_ok =
                certificate.correction == recomputed.correction &&
                certificate.normalization.correction == recomputed.correction
            inverse_correction_ok =
                certificate.inverse_correction == recomputed.correction.inverse_factor &&
                certificate.correction.inverse_factor == recomputed.correction.inverse_factor
            normalized_core_ok =
                certificate.normalized_core == recomputed.normalized_matrix &&
                certificate.normalization.normalized_matrix == recomputed.normalized_matrix
            core_det_ok = det(certificate.normalized_core) == one(R)
            core_replay_ok =
                certificate.core_factorization.original_matrix == certificate.normalized_core &&
                _verify_laurent_column_peel_replay(certificate.core_factorization)
            core_factors_match_ok = _factor_sequences_equal(
                certificate.core_factors,
                certificate.core_factorization.factors,
            )
            core_factors_ok = verify_factorization(certificate.normalized_core, certificate.core_factors)
            rebuilt_core = _factor_product(certificate.core_factors, R, n)
            rebuilt = certificate.correction.factor * rebuilt_core
            reconstructed_product_ok =
                rebuilt_core == certificate.normalized_core &&
                rebuilt == A &&
                certificate.reconstructed_product == rebuilt
        catch err
            err isa InterruptException && rethrow()
        end
    end

    overall_ok = size_ok && normalization_ok && determinant_profile_ok &&
        correction_ok && inverse_correction_ok && normalized_core_ok &&
        core_det_ok && core_replay_ok && core_factors_match_ok &&
        core_factors_ok && reconstructed_product_ok

    return (
        overall_ok = overall_ok,
        size_ok = size_ok,
        normalization_ok = normalization_ok,
        determinant_profile_ok = determinant_profile_ok,
        correction_ok = correction_ok,
        inverse_correction_ok = inverse_correction_ok,
        normalized_core_ok = normalized_core_ok,
        core_det_ok = core_det_ok,
        core_replay_ok = core_replay_ok,
        core_factors_match_ok = core_factors_match_ok,
        core_factors_ok = core_factors_ok,
        reconstructed_product_ok = reconstructed_product_ok,
    )
end
