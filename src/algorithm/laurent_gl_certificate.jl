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
    determinant_source::Symbol
    correction_side::Symbol
    reconstruction_relation::Symbol
    correction
    inverse_correction
    normalized_deferred_core
    normalized_deferred_factorization
    normalized_deferred_factors::Vector
    rewritten_left_factors::Vector
    rewritten_right_factors::Vector
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
    for factor in factors
        push!(rewritten, _rewrite_left_elementary_factor_across_diagonal(factor, diagonal_factor))
    end
    return rewritten
end

function _rewrite_right_elementary_factor_across_diagonal(factor, diagonal_factor)
    n = _require_square_matrix(diagonal_factor, "diagonal factor")
    nrows(factor) == n || throw(ArgumentError("elementary factor and diagonal factor sizes must match"))
    ncols(factor) == n || throw(ArgumentError("elementary factor and diagonal factor sizes must match"))
    _same_base_ring(base_ring(factor), base_ring(diagonal_factor)) ||
        throw(ArgumentError("elementary factor and diagonal factor must have the same base ring"))

    diagonal_entries = _laurent_diagonal_entries(diagonal_factor)
    row, col, coefficient = _elementary_factor_data(factor)
    rewritten_coefficient = diagonal_entries[row] * coefficient * inv(diagonal_entries[col])
    return elementary_matrix(n, row, col, rewritten_coefficient, base_ring(factor))
end

function _rewrite_right_elementary_factors_across_diagonal(factors, diagonal_factor)
    R = base_ring(diagonal_factor)
    n = _require_square_matrix(diagonal_factor, "diagonal factor")
    rewritten = typeof(identity_matrix(R, n))[]
    for factor in factors
        push!(rewritten, _rewrite_right_elementary_factor_across_diagonal(factor, diagonal_factor))
    end
    return rewritten
end

function _laurent_gl_correction_side(correction_side)::Symbol
    correction_side isa Symbol ||
        throw(ArgumentError("correction_side must be :row or :column"))
    correction_side == :row && return :left
    correction_side == :column && return :right
    throw(ArgumentError("correction_side must be :row or :column"))
end

function _laurent_gl_external_correction_side(correction_side::Symbol)::Symbol
    correction_side == :left && return :row
    correction_side == :right && return :column
    throw(ArgumentError("unsupported internal Laurent correction side $(correction_side)"))
end

function _laurent_gl_reconstruction_relation(correction_side::Symbol)::Symbol
    correction_side == :left && return :left_correction_times_elementary_product
    correction_side == :right && return :elementary_product_times_right_correction
    throw(ArgumentError("unsupported internal Laurent correction side $(correction_side)"))
end

function _side_aware_identity_correction(R, n::Int, determinant, side::Symbol)
    correction = _identity_correction(R, n, determinant)
    return merge(correction, (; side))
end

function _right_diagonal_determinant_correction(R, n::Int, determinant)
    correction = _left_diagonal_determinant_correction(R, n, determinant)
    return merge(correction, (; side = :right, kind = :right_diagonal_determinant_correction))
end

function _laurent_lazy_deferred_correction(metadata, correction_side::Symbol)
    metadata.deferred_correction !== nothing ||
        throw(ArgumentError("supported deferred Laurent determinant metadata is missing a correction"))
    deferred = metadata.deferred_submatrix
    R = base_ring(deferred)
    n = nrows(deferred)
    seed_correction = metadata.deferred_correction
    determinant = seed_correction.determinant

    if metadata.determinant_classification == :one
        return merge(
            seed_correction,
            _side_aware_identity_correction(R, n, determinant, correction_side),
        )
    elseif metadata.determinant_classification == :laurent_monomial_unit
        correction = correction_side == :left ?
            _left_diagonal_determinant_correction(R, n, determinant) :
            _right_diagonal_determinant_correction(R, n, determinant)
        return merge(seed_correction, correction)
    end

    throw(ArgumentError("unsupported deferred Laurent determinant cannot be hoisted"))
end

function _laurent_lazy_normalized_deferred_core(metadata, correction)
    deferred = metadata.deferred_submatrix
    if metadata.determinant_classification == :one
        return deferred
    elseif correction.side == :left
        return correction.inverse_factor * deferred
    elseif correction.side == :right
        return deferred * correction.inverse_factor
    end

    throw(ArgumentError("unsupported internal Laurent correction side $(correction.side)"))
end

function _embed_laurent_deferred_correction(correction, original_dimension::Int)
    deferred_dimension = nrows(correction.factor)
    correction.side in (:left, :right) ||
        throw(ArgumentError("only left or right deferred Laurent corrections can be hoisted"))

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
    deferred_dimension = nrows(metadata.deferred_submatrix)
    left_inverse = _inverse_elementary_sequence(peel_certificate.left_factors)
    embedded_core = _embed_laurent_deferred_peel_factors(
        normalized_deferred_factors,
        R,
        original_dimension,
        deferred_dimension,
    )
    right_inverse = _inverse_elementary_sequence(peel_certificate.right_factors)
    empty_rewrites = typeof(identity_matrix(R, original_dimension))[]

    if correction.side == :left
        rewritten_left = _rewrite_left_elementary_factors_across_diagonal(
            left_inverse,
            correction.factor,
        )
        return (;
            rewritten_left,
            rewritten_right = empty_rewrites,
            elementary_factors = vcat(rewritten_left, embedded_core, right_inverse),
        )
    elseif correction.side == :right
        rewritten_right = _rewrite_right_elementary_factors_across_diagonal(
            right_inverse,
            correction.factor,
        )
        return (;
            rewritten_left = empty_rewrites,
            rewritten_right,
            elementary_factors = vcat(left_inverse, embedded_core, rewritten_right),
        )
    end

    throw(ArgumentError("unsupported internal Laurent correction side $(correction.side)"))
end

function _laurent_gl_lazy_deferred_correction_certificate(
    metadata::NamedTuple;
    correction_side = :row,
    progress_callback = nothing,
)
    _verify_laurent_determinant_deferred_submatrix_normalization(metadata) ||
        error("invalid deferred Laurent determinant normalization metadata")
    metadata.supported ||
        throw(ArgumentError("unsupported deferred Laurent determinant cannot be hoisted"))

    A = metadata.peel_certificate.original_matrix
    R = base_ring(A)
    n = nrows(A)
    internal_correction_side = _laurent_gl_correction_side(correction_side)
    deferred_correction = _laurent_lazy_deferred_correction(metadata, internal_correction_side)
    correction = _embed_laurent_deferred_correction(deferred_correction, n)
    normalized_deferred_core = _laurent_lazy_normalized_deferred_core(metadata, deferred_correction)
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
    reconstruction_relation = _laurent_gl_reconstruction_relation(internal_correction_side)
    reconstructed_product = internal_correction_side == :left ?
        correction.factor * elementary_product :
        elementary_product * correction.factor
    certificate = LaurentLazyGLHoistCertificate(
        A,
        metadata,
        metadata.overall_determinant,
        metadata.determinant_source,
        correction_side,
        reconstruction_relation,
        correction,
        correction.inverse_factor,
        normalized_deferred_core,
        normalized_deferred_factorization,
        normalized_deferred_factors,
        assembly.rewritten_left,
        assembly.rewritten_right,
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
        certificate.determinant_source,
        certificate.correction_side,
        certificate.reconstruction_relation,
        certificate.correction,
        certificate.inverse_correction,
        certificate.normalized_deferred_core,
        certificate.normalized_deferred_factorization,
        certificate.normalized_deferred_factors,
        certificate.rewritten_left_factors,
        certificate.rewritten_right_factors,
        certificate.elementary_factors,
        certificate.elementary_product,
        certificate.reconstructed_product,
        verification,
    )
end

function _laurent_gl_lazy_deferred_correction_certificate(
    A;
    correction_side = :row,
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
        correction_side,
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

function _laurent_gl_certificate_strategy(determinant_strategy)::Symbol
    determinant_strategy isa Symbol ||
        throw(ArgumentError("determinant_strategy must be :eager or :lazy"))
    determinant_strategy in (:eager, :lazy) && return determinant_strategy
    throw(ArgumentError("determinant_strategy must be :eager or :lazy"))
end

function laurent_gl_factorization_certificate(
    A;
    determinant_strategy = :eager,
    correction_side = nothing,
    progress_callback = nothing,
)
    strategy = _laurent_gl_certificate_strategy(determinant_strategy)
    if strategy == :eager
        correction_side === nothing ||
            throw(ArgumentError("correction_side is supported only with determinant_strategy = :lazy"))
        return _laurent_gl_factorization_certificate(
            A;
            progress_callback,
        )
    end

    side = correction_side === nothing ? :row : correction_side
    return _laurent_gl_lazy_deferred_correction_certificate(
        A;
        correction_side = side,
        progress_callback,
    )
end

function verify_laurent_gl_factorization_certificate(
    certificate::LaurentLazyGLHoistCertificate,
)::Bool
    return _laurent_gl_lazy_deferred_correction_certificate_verification(certificate).overall_ok
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
    determinant_source_ok = false
    correction_ok = false
    normalized_core_ok = false
    normalized_factorization_ok = false
    correction_side_ok = false
    reconstruction_relation_ok = false
    rewritten_left_factors_ok = false
    rewritten_right_factors_ok = false
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
            determinant_source_ok =
                certificate.determinant_source == :deferred_submatrix &&
                certificate.determinant_source == metadata.determinant_source

            expected_side = _laurent_gl_correction_side(certificate.correction_side)
            expected_relation = _laurent_gl_reconstruction_relation(expected_side)
            expected_deferred_correction = _laurent_lazy_deferred_correction(
                metadata,
                expected_side,
            )
            expected_correction = _embed_laurent_deferred_correction(
                expected_deferred_correction,
                n,
            )
            identity = identity_matrix(R, n)
            correction_side_ok =
                _laurent_gl_external_correction_side(expected_side) == certificate.correction_side &&
                certificate.correction.side == expected_side
            reconstruction_relation_ok =
                certificate.reconstruction_relation == expected_relation
            correction_ok =
                certificate.overall_determinant == metadata.overall_determinant &&
                certificate.correction == expected_correction &&
                certificate.inverse_correction == expected_correction.inverse_factor &&
                certificate.correction.factor * certificate.inverse_correction == identity &&
                certificate.inverse_correction * certificate.correction.factor == identity &&
                det(certificate.correction.factor) == certificate.overall_determinant

            expected_normalized_core = _laurent_lazy_normalized_deferred_core(
                metadata,
                expected_deferred_correction,
            )
            normalized_core_ok =
                certificate.normalized_deferred_core == expected_normalized_core &&
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
            expected_rewritten_left = expected_side == :left ?
                _rewrite_left_elementary_factors_across_diagonal(
                    expected_left_inverse,
                    certificate.correction.factor,
                ) :
                typeof(identity_matrix(R, n))[]
            rewritten_left_factors_ok = _factor_sequences_equal(
                certificate.rewritten_left_factors,
                expected_rewritten_left,
            )
            expected_right_inverse = _inverse_elementary_sequence(
                metadata.peel_certificate.right_factors,
            )
            expected_rewritten_right = expected_side == :right ?
                _rewrite_right_elementary_factors_across_diagonal(
                    expected_right_inverse,
                    certificate.correction.factor,
                ) :
                typeof(identity_matrix(R, n))[]
            rewritten_right_factors_ok = _factor_sequences_equal(
                certificate.rewritten_right_factors,
                expected_rewritten_right,
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
            expected_reconstructed_product = expected_side == :left ?
                certificate.correction.factor * certificate.elementary_product :
                certificate.elementary_product * certificate.correction.factor
            reconstructed_product_ok =
                certificate.reconstructed_product == expected_reconstructed_product &&
                certificate.reconstructed_product == A
        catch err
            err isa InterruptException && rethrow()
        end
    end

    overall_ok = size_ok && metadata_ok && determinant_source_ok &&
        correction_ok && normalized_core_ok &&
        normalized_factorization_ok && correction_side_ok &&
        reconstruction_relation_ok && rewritten_left_factors_ok &&
        rewritten_right_factors_ok &&
        elementary_factors_ok && elementary_product_ok && reconstructed_product_ok

    return (;
        overall_ok,
        size_ok,
        metadata_ok,
        determinant_source_ok,
        correction_ok,
        normalized_core_ok,
        normalized_factorization_ok,
        correction_side_ok,
        reconstruction_relation_ok,
        rewritten_left_factors_ok,
        rewritten_right_factors_ok,
        elementary_factors_ok,
        elementary_product_ok,
        reconstructed_product_ok,
    )
end
