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
        return new(
            original_matrix,
            determinant_profile,
            normalization,
            correction,
            inverse_correction,
            normalized_core,
            core_factorization,
            core_factors,
            reconstructed_product,
            verification,
        )
    end
end

struct LaurentLazyGLHoistCertificate
    original_matrix
    deferred_metadata
    overall_determinant
    correction
    inverse_correction
    normalized_deferred_core
    normalized_deferred_factorization
    normalized_deferred_factors::Vector
    rewritten_left_factors::Vector
    elementary_factors::Vector
    elementary_product
    reconstructed_product
    verification
end

function _laurent_gl_factorization_certificate_from_normalization(
    A,
    normalization;
    progress_callback = nothing,
)
    core = normalization.normalized_matrix
    core_factorization = _factor_laurent_sl_column_peel(core; progress_callback)
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

function _laurent_diagonal_entries(diagonal_factor)
    n = _require_square_matrix(diagonal_factor, "diagonal factor")
    R = base_ring(diagonal_factor)
    entries = [diagonal_factor[i, i] for i in 1:n]
    for i in 1:n, j in 1:n
        i == j && continue
        diagonal_factor[i, j] == zero(R) ||
            throw(ArgumentError("diagonal factor must have zero off-diagonal entries"))
    end
    return entries
end

function _rewrite_left_elementary_factor_across_diagonal(factor, diagonal_factor)
    n = _require_square_matrix(diagonal_factor, "diagonal factor")
    nrows(factor) == n || throw(ArgumentError("elementary factor and diagonal factor sizes must match"))
    ncols(factor) == n || throw(ArgumentError("elementary factor and diagonal factor sizes must match"))
    _same_base_ring(base_ring(factor), base_ring(diagonal_factor)) ||
        throw(ArgumentError("elementary factor and diagonal factor must have the same base ring"))

    diagonal_entries = _laurent_diagonal_entries(diagonal_factor)
    row, col, coefficient = _elementary_factor_data(factor)
    rewritten_coefficient = inv(diagonal_entries[row]) * coefficient * diagonal_entries[col]
    return elementary_matrix(n, row, col, rewritten_coefficient, base_ring(factor))
end

function _rewrite_left_elementary_factors_across_diagonal(factors, diagonal_factor)
    R = base_ring(diagonal_factor)
    n = _require_square_matrix(diagonal_factor, "diagonal factor")
    rewritten = typeof(identity_matrix(R, n))[]
    diagonal_entries = _laurent_diagonal_entries(diagonal_factor)
    split_generator = any(entry != one(R) for entry in diagonal_entries) ? first(gens(R)) : nothing
    for factor in factors
        rewritten_factor = _rewrite_left_elementary_factor_across_diagonal(factor, diagonal_factor)
        if split_generator !== nothing && rewritten_factor == factor
            row, col, coefficient = _elementary_factor_data(rewritten_factor)
            leading_coefficient = coefficient * split_generator
            trailing_coefficient = coefficient - leading_coefficient
            if leading_coefficient != zero(R) && trailing_coefficient != zero(R)
                push!(rewritten, elementary_matrix(n, row, col, leading_coefficient, R))
                push!(rewritten, elementary_matrix(n, row, col, trailing_coefficient, R))
                continue
            end
        end
        push!(rewritten, rewritten_factor)
    end
    return rewritten
end

function _embed_laurent_deferred_correction(correction, original_dimension::Int)
    deferred_dimension = nrows(correction.factor)
    correction_side = correction.side
    correction_side == :left ||
        throw(ArgumentError("only left deferred Laurent corrections can be hoisted"))

    embedded_factor = block_embedding(
        correction.factor,
        original_dimension,
        collect(1:deferred_dimension),
    )
    embedded_inverse = block_embedding(
        correction.inverse_factor,
        original_dimension,
        collect(1:deferred_dimension),
    )
    return merge(correction, (;
        scope = :original_matrix,
        deferred_scope = correction.scope,
        factor = embedded_factor,
        inverse_factor = embedded_inverse,
        embedded_from_dimension = deferred_dimension,
    ))
end

function _laurent_lazy_hoist_elementary_factors(metadata, correction, normalized_deferred_factors)
    peel_certificate = metadata.peel_certificate
    R = base_ring(peel_certificate.original_matrix)
    original_dimension = nrows(peel_certificate.original_matrix)
    deferred_dimension = nrows(metadata.normalized_deferred_core)
    left_inverse = _inverse_elementary_sequence(peel_certificate.left_factors)
    rewritten_left = _rewrite_left_elementary_factors_across_diagonal(
        left_inverse,
        correction.factor,
    )
    embedded_core = _embed_laurent_deferred_peel_factors(
        normalized_deferred_factors,
        R,
        original_dimension,
        deferred_dimension,
    )
    right_inverse = _inverse_elementary_sequence(peel_certificate.right_factors)
    return (;
        rewritten_left,
        elementary_factors = vcat(rewritten_left, embedded_core, right_inverse),
    )
end

function _laurent_gl_lazy_deferred_correction_certificate(
    metadata::NamedTuple;
    progress_callback = nothing,
)
    _verify_laurent_determinant_deferred_submatrix_normalization(metadata) ||
        error("invalid deferred Laurent determinant normalization metadata")
    metadata.supported ||
        throw(ArgumentError("unsupported deferred Laurent determinant cannot be hoisted"))
    metadata.deferred_correction !== nothing ||
        throw(ArgumentError("supported deferred Laurent determinant metadata is missing a correction"))

    A = metadata.peel_certificate.original_matrix
    R = base_ring(A)
    n = nrows(A)
    correction = _embed_laurent_deferred_correction(metadata.deferred_correction, n)
    normalized_deferred_core = metadata.normalized_deferred_core
    normalized_deferred_factorization = _factor_laurent_sl_column_peel(
        normalized_deferred_core;
        progress_callback,
    )
    normalized_deferred_factors = normalized_deferred_factorization.factors
    assembly = _laurent_lazy_hoist_elementary_factors(
        metadata,
        correction,
        normalized_deferred_factors,
    )
    elementary_product = _factor_product(assembly.elementary_factors, R, n)
    reconstructed_product = correction.factor * elementary_product
    certificate = LaurentLazyGLHoistCertificate(
        A,
        metadata,
        metadata.overall_determinant,
        correction,
        correction.inverse_factor,
        normalized_deferred_core,
        normalized_deferred_factorization,
        normalized_deferred_factors,
        assembly.rewritten_left,
        assembly.elementary_factors,
        elementary_product,
        reconstructed_product,
        nothing,
    )
    verification = _laurent_gl_lazy_deferred_correction_certificate_verification(certificate)
    verification.overall_ok ||
        error("internal lazy Laurent GL correction hoist verification failed")
    return LaurentLazyGLHoistCertificate(
        certificate.original_matrix,
        certificate.deferred_metadata,
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
        verification,
    )
end

function _laurent_gl_lazy_deferred_correction_certificate(
    A;
    determinant_probe = classify_laurent_determinant,
    progress_callback = nothing,
)
    deferred_certificate = _laurent_determinant_deferred_peel_certificate(
        A;
        progress_callback,
    )
    metadata = _normalize_laurent_determinant_deferred_submatrix(
        deferred_certificate;
        determinant_probe,
    )
    return _laurent_gl_lazy_deferred_correction_certificate(
        metadata;
        progress_callback,
    )
end

function _laurent_gl_factorization_certificate(A; progress_callback = nothing)
    normalization = normalize_laurent_gl_matrix(A)
    return _laurent_gl_factorization_certificate_from_normalization(
        A,
        normalization;
        progress_callback,
    )
end

function laurent_gl_factorization_certificate(A)
    return _laurent_gl_factorization_certificate(A)
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

function _verify_laurent_gl_lazy_deferred_correction_certificate(certificate)::Bool
    return _laurent_gl_lazy_deferred_correction_certificate_verification(certificate).overall_ok
end

function _laurent_gl_lazy_deferred_correction_certificate_verification(certificate)
    size_ok = try
        nrows(certificate.original_matrix) == ncols(certificate.original_matrix)
    catch err
        err isa InterruptException && rethrow()
        false
    end

    metadata_ok = false
    correction_ok = false
    normalized_core_ok = false
    normalized_factorization_ok = false
    rewritten_left_factors_ok = false
    elementary_factors_ok = false
    elementary_product_ok = false
    reconstructed_product_ok = false

    if size_ok
        try
            A = certificate.original_matrix
            R = base_ring(A)
            n = nrows(A)
            metadata = certificate.deferred_metadata
            metadata_ok =
                metadata.peel_certificate.original_matrix == A &&
                metadata.supported &&
                metadata.determinant_source == :deferred_submatrix &&
                _verify_laurent_determinant_deferred_submatrix_normalization(metadata)

            expected_correction = _embed_laurent_deferred_correction(
                metadata.deferred_correction,
                n,
            )
            identity = identity_matrix(R, n)
            correction_ok =
                certificate.overall_determinant == metadata.overall_determinant &&
                certificate.correction == expected_correction &&
                certificate.inverse_correction == expected_correction.inverse_factor &&
                certificate.correction.factor * certificate.inverse_correction == identity &&
                certificate.inverse_correction * certificate.correction.factor == identity &&
                det(certificate.correction.factor) == certificate.overall_determinant

            normalized_core_ok =
                certificate.normalized_deferred_core == metadata.normalized_deferred_core &&
                det(certificate.normalized_deferred_core) ==
                    one(base_ring(certificate.normalized_deferred_core))

            normalized_factorization_ok =
                certificate.normalized_deferred_factorization.original_matrix ==
                    certificate.normalized_deferred_core &&
                _verify_laurent_column_peel_replay(certificate.normalized_deferred_factorization) &&
                _factor_sequences_equal(
                    certificate.normalized_deferred_factors,
                    certificate.normalized_deferred_factorization.factors,
                ) &&
                verify_factorization(
                    certificate.normalized_deferred_core,
                    certificate.normalized_deferred_factors,
                )

            expected_left_inverse = _inverse_elementary_sequence(
                metadata.peel_certificate.left_factors,
            )
            expected_rewritten_left = _rewrite_left_elementary_factors_across_diagonal(
                expected_left_inverse,
                certificate.correction.factor,
            )
            rewritten_left_factors_ok = _factor_sequences_equal(
                certificate.rewritten_left_factors,
                expected_rewritten_left,
            )

            assembly = _laurent_lazy_hoist_elementary_factors(
                metadata,
                certificate.correction,
                certificate.normalized_deferred_factors,
            )
            elementary_factors_ok = _factor_sequences_equal(
                certificate.elementary_factors,
                assembly.elementary_factors,
            )
            rebuilt_product = _factor_product(certificate.elementary_factors, R, n)
            elementary_product_ok = certificate.elementary_product == rebuilt_product
            reconstructed_product_ok =
                certificate.reconstructed_product ==
                    certificate.correction.factor * certificate.elementary_product &&
                certificate.reconstructed_product == A
        catch err
            err isa InterruptException && rethrow()
        end
    end

    overall_ok = size_ok && metadata_ok && correction_ok && normalized_core_ok &&
        normalized_factorization_ok && rewritten_left_factors_ok &&
        elementary_factors_ok && elementary_product_ok && reconstructed_product_ok

    return (;
        overall_ok,
        size_ok,
        metadata_ok,
        correction_ok,
        normalized_core_ok,
        normalized_factorization_ok,
        rewritten_left_factors_ok,
        elementary_factors_ok,
        elementary_product_ok,
        reconstructed_product_ok,
    )
end
