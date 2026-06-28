struct LaurentColumnPeelStep
    dimension::Int
    input_matrix
    last_column::Vector
    left_factors::Vector
    after_left_matrix
    right_factors::Vector
    peeled_matrix
    next_block
end

struct LaurentDeterminantDeferredPeelCertificate
    original_matrix
    peel_steps::Vector{LaurentColumnPeelStep}
    deferred_submatrix
    determinant_source::Symbol
    left_factors::Vector
    right_factors::Vector
    left_product
    right_product
    target_matrix
    verification

    function LaurentDeterminantDeferredPeelCertificate(
        original_matrix,
        peel_steps::Vector{LaurentColumnPeelStep},
        deferred_submatrix,
        determinant_source::Symbol,
        verification,
    )
        replay = _laurent_determinant_deferred_peel_replay_data(
            original_matrix,
            peel_steps,
            deferred_submatrix,
        )
        return new(
            original_matrix,
            collect(peel_steps),
            deferred_submatrix,
            determinant_source,
            replay.left_factors,
            replay.right_factors,
            replay.left_product,
            replay.right_product,
            replay.target_matrix,
            verification,
        )
    end
end

struct LaurentColumnPeelFactorization
    original_matrix
    final_block
    final_local_target
    final_local_factors::Vector
    final_factors::Vector
    factors::Vector
    product
    peel_steps::Vector{LaurentColumnPeelStep}
    verification

    function LaurentColumnPeelFactorization(
        original_matrix,
        final_block,
        final_local_target,
        final_local_factors::Vector,
        final_factors::Vector,
        factors::Vector,
        product,
        peel_steps::Vector{LaurentColumnPeelStep},
        verification,
    )
        R = base_ring(original_matrix)
        replayed_factors = _replay_laurent_column_peel_factors(peel_steps, final_factors, R)
        replayed_product = _factor_product(replayed_factors, R, nrows(original_matrix))
        return new(
            original_matrix,
            final_block,
            final_local_target,
            final_local_factors,
            final_factors,
            replayed_factors,
            replayed_product,
            peel_steps,
            verification,
        )
    end
end

function _laurent_column_peel_factors(A)
    return _factor_laurent_sl_column_peel(A).factors
end

function _verify_laurent_column_peel_replay(certificate)::Bool
    return _laurent_column_peel_verification(certificate).overall_ok
end

function _validate_laurent_column_peel_input_shape_and_ring(A)
    nrows(A) == ncols(A) || throw(ArgumentError("Laurent column-peel factorization requires a square matrix"))
    nrows(A) >= 2 || throw(ArgumentError("Laurent column-peel factorization requires size at least 2"))
    R = base_ring(A)
    _is_laurent_polynomial_ring(R) || throw(ArgumentError("Laurent column-peel factorization requires a Laurent polynomial ring"))
    return nrows(A)
end

function _validate_laurent_column_peel_input_determinant_one(
    A;
    determinant_probe = classify_laurent_determinant,
)
    profile = determinant_probe(A)
    profile.classification == :one || throw(ArgumentError("Laurent column-peel factorization requires determinant-one input"))
    return nrows(A)
end

function _laurent_column_peel_empty_last_completed()
    return (;
        dimension = nothing,
        elapsed_seconds = nothing,
        left_factors = nothing,
        right_factors = nothing,
    )
end

function _laurent_column_peel_column_stats(column)
    nnz = count(!iszero, column)
    max_terms = 0
    for entry in column
        iszero(entry) && continue
        term_count = try
            length(collect(coefficients(entry)))
        catch err
            err isa MethodError || err isa ErrorException || rethrow()
            0
        end
        max_terms = max(max_terms, term_count)
    end
    return (; last_column_nnz = nnz, max_entry_terms = max_terms)
end

function _emit_laurent_column_peel_progress(progress_callback, current, completed_steps::Int, last_completed)
    progress_callback === nothing && return nothing
    d = nrows(current)
    last_column = [current[row, d] for row in 1:d]
    stats = _laurent_column_peel_column_stats(last_column)
    progress_callback((;
        current_dimension = d,
        completed_steps,
        last_completed_dimension = last_completed.dimension,
        last_completed_elapsed_seconds = last_completed.elapsed_seconds,
        last_completed_left_factors = last_completed.left_factors,
        last_completed_right_factors = last_completed.right_factors,
        last_column_nnz = stats.last_column_nnz,
        max_entry_terms = stats.max_entry_terms,
    ))
    return nothing
end

function _factor_laurent_sl_column_peel(
    A;
    progress_callback = nothing,
    determinant_probe = classify_laurent_determinant,
)
    _validate_laurent_column_peel_input_shape_and_ring(A)
    _emit_laurent_column_peel_progress(
        progress_callback,
        A,
        0,
        _laurent_column_peel_empty_last_completed(),
    )
    _validate_laurent_column_peel_input_determinant_one(A; determinant_probe)
    factors, steps, final_block, final_target, final_local, final_2x2 =
        _laurent_column_peel_recursive(
            A;
            progress_callback,
            completed_steps = 0,
            last_completed = _laurent_column_peel_empty_last_completed(),
            emit_current = false,
        )
    R = base_ring(A)
    product = _factor_product(factors, R, nrows(A))
    certificate = LaurentColumnPeelFactorization(
        A,
        final_block,
        final_target,
        final_local,
        final_2x2,
        factors,
        product,
        steps,
        nothing,
    )
    verification = _laurent_column_peel_verification(certificate)
    verification.overall_ok || error("internal Laurent column-peel verification failed")
    return LaurentColumnPeelFactorization(
        certificate.original_matrix,
        certificate.final_block,
        certificate.final_local_target,
        certificate.final_local_factors,
        certificate.final_factors,
        certificate.factors,
        certificate.product,
        certificate.peel_steps,
        verification,
    )
end

function _factor_laurent_gl_lazy_determinant_peel(
    A;
    determinant_probe = classify_laurent_determinant,
    progress_callback = nothing,
)
    _validate_laurent_column_peel_input_shape_and_ring(A)
    n = nrows(A)
    n >= 3 || throw(ArgumentError("Lazy Laurent determinant peel requires size at least 3 so one peel step can complete before determinant classification"))
    empty_last_completed = _laurent_column_peel_empty_last_completed()
    _emit_laurent_column_peel_progress(progress_callback, A, 0, empty_last_completed)

    step_started_at = time()
    step = _laurent_column_peel_step(A)
    step_elapsed = round(time() - step_started_at; digits = 3)
    completed_steps = 1
    last_completed = (;
        dimension = n,
        elapsed_seconds = step_elapsed,
        left_factors = length(step.left_factors),
        right_factors = length(step.right_factors),
    )
    _emit_laurent_column_peel_progress(progress_callback, step.next_block, completed_steps, last_completed)

    deferred_certificate = LaurentDeterminantDeferredPeelCertificate(
        A,
        LaurentColumnPeelStep[step],
        step.next_block,
        :deferred_submatrix,
        nothing,
    )
    deferred_metadata = _normalize_laurent_determinant_deferred_submatrix(
        deferred_certificate;
        determinant_probe,
    )
    deferred_metadata.determinant_classification == :one || return deferred_metadata

    next_factors, next_steps, final_block, final_target, final_local, final_2x2 =
        _laurent_column_peel_recursive(
            step.next_block;
            progress_callback,
            completed_steps,
            last_completed,
            emit_current = false,
        )
    R = base_ring(A)
    factors = vcat(
        _inverse_elementary_sequence(step.left_factors),
        _embed_upper_left_factors(next_factors, R, n),
        _inverse_elementary_sequence(step.right_factors),
    )
    product = _factor_product(factors, R, n)
    certificate = LaurentColumnPeelFactorization(
        A,
        final_block,
        final_target,
        final_local,
        final_2x2,
        factors,
        product,
        vcat(LaurentColumnPeelStep[step], next_steps),
        nothing,
    )
    verification = _laurent_column_peel_verification(certificate)
    verification.overall_ok || error("internal lazy Laurent column-peel verification failed")
    return LaurentColumnPeelFactorization(
        certificate.original_matrix,
        certificate.final_block,
        certificate.final_local_target,
        certificate.final_local_factors,
        certificate.final_factors,
        certificate.factors,
        certificate.product,
        certificate.peel_steps,
        verification,
    )
end

function _laurent_determinant_deferred_peel_certificate(
    A;
    min_steps::Int = 1,
    progress_callback = nothing,
)
    _validate_laurent_column_peel_input_shape_and_ring(A)
    min_steps >= 1 || throw(ArgumentError("determinant-deferred Laurent peel requires at least one peel step"))
    n = nrows(A)
    n - min_steps >= 2 || throw(ArgumentError("determinant-deferred Laurent peel must leave a deferred submatrix of size at least 2"))

    steps = LaurentColumnPeelStep[]
    current = A
    completed_steps = 0
    last_completed = _laurent_column_peel_empty_last_completed()
    _emit_laurent_column_peel_progress(progress_callback, current, completed_steps, last_completed)

    for _ in 1:min_steps
        step_started_at = time()
        step = _laurent_column_peel_step(current)
        step_elapsed = round(time() - step_started_at; digits = 3)
        push!(steps, step)
        completed_steps += 1
        last_completed = (;
            dimension = step.dimension,
            elapsed_seconds = step_elapsed,
            left_factors = length(step.left_factors),
            right_factors = length(step.right_factors),
        )
        current = step.next_block
        _emit_laurent_column_peel_progress(progress_callback, current, completed_steps, last_completed)
    end

    certificate = LaurentDeterminantDeferredPeelCertificate(
        A,
        steps,
        current,
        :deferred_submatrix,
        nothing,
    )
    verification = _laurent_determinant_deferred_peel_verification(certificate)
    verification.overall_ok || error("internal determinant-deferred Laurent peel verification failed")
    return LaurentDeterminantDeferredPeelCertificate(
        certificate.original_matrix,
        certificate.peel_steps,
        certificate.deferred_submatrix,
        certificate.determinant_source,
        verification,
    )
end

function _verify_laurent_determinant_deferred_peel_replay(certificate)::Bool
    return _laurent_determinant_deferred_peel_verification(certificate).overall_ok
end

function _is_supported_deferred_laurent_determinant_class(classification::Symbol)::Bool
    return classification == :one || classification == :laurent_monomial_unit
end

function _scoped_deferred_correction(correction)
    return merge(correction, (; scope = :deferred_submatrix))
end

function _deferred_laurent_determinant_boundary(certificate, determinant_profile)
    classification = determinant_profile.classification
    reason = classification == :non_unit ?
        :non_unit_deferred_determinant :
        :unsupported_deferred_unit_class
    deferred = certificate.deferred_submatrix
    return (;
        kind = :unsupported_deferred_laurent_determinant,
        determinant_source = certificate.determinant_source,
        overall_determinant = determinant_profile.determinant,
        determinant_classification = classification,
        deferred_submatrix_size = (nrows(deferred), ncols(deferred)),
        supported = false,
        reason,
    )
end

function _normalize_laurent_determinant_deferred_submatrix(
    certificate;
    determinant_probe = classify_laurent_determinant,
)
    _verify_laurent_determinant_deferred_peel_replay(certificate) ||
        error("invalid determinant-deferred Laurent peel certificate replay")

    deferred = certificate.deferred_submatrix
    R = base_ring(deferred)
    n = nrows(deferred)
    determinant_profile = determinant_probe(deferred)
    classification = determinant_profile.classification
    determinant = determinant_profile.determinant

    deferred_correction = nothing
    deferred_diagonal_correction = nothing
    normalized_deferred_core = nothing
    staged_boundary = nothing

    if classification == :one
        deferred_correction = _scoped_deferred_correction(
            _identity_correction(R, n, determinant),
        )
        normalized_deferred_core = deferred
    elseif classification == :laurent_monomial_unit
        deferred_correction = _scoped_deferred_correction(
            _left_diagonal_determinant_correction(R, n, determinant),
        )
        deferred_diagonal_correction = deferred_correction
        normalized_deferred_core = deferred_correction.inverse_factor * deferred
    else
        staged_boundary = _deferred_laurent_determinant_boundary(
            certificate,
            determinant_profile,
        )
    end

    metadata = (;
        peel_certificate = certificate,
        deferred_submatrix = deferred,
        determinant_source = certificate.determinant_source,
        determinant_profile,
        overall_determinant = determinant,
        determinant_classification = classification,
        supported = _is_supported_deferred_laurent_determinant_class(classification),
        deferred_correction,
        deferred_diagonal_correction,
        normalized_deferred_core,
        staged_boundary,
        verification = nothing,
    )
    verification = _laurent_determinant_deferred_submatrix_normalization_verification(metadata)
    verification.overall_ok ||
        error("internal deferred Laurent submatrix normalization metadata verification failed")
    return merge(metadata, (; verification))
end

function _verify_laurent_determinant_deferred_submatrix_normalization(metadata)::Bool
    return _laurent_determinant_deferred_submatrix_normalization_verification(metadata).overall_ok
end

function _laurent_determinant_deferred_submatrix_normalization_verification(metadata)
    certificate_ok = try
        _verify_laurent_determinant_deferred_peel_replay(metadata.peel_certificate)
    catch err
        err isa InterruptException && rethrow()
        false
    end

    determinant_ok = false
    correction_ok = false
    normalized_core_det_ok = false
    boundary_ok = false

    try
        deferred = metadata.peel_certificate.deferred_submatrix
        R = base_ring(deferred)
        n = nrows(deferred)
        profile = classify_laurent_determinant(deferred)
        classification = profile.classification
        determinant_ok =
            metadata.deferred_submatrix == deferred &&
            metadata.determinant_source == :deferred_submatrix &&
            metadata.overall_determinant == profile.determinant &&
            metadata.determinant_classification == classification &&
            metadata.determinant_profile == profile &&
            metadata.supported == _is_supported_deferred_laurent_determinant_class(classification)

        if metadata.supported
            correction = metadata.deferred_correction
            identity = identity_matrix(R, n)
            correction_ok =
                metadata.staged_boundary === nothing &&
                correction !== nothing &&
                correction.scope == :deferred_submatrix &&
                correction.determinant == metadata.overall_determinant &&
                nrows(correction.factor) == n &&
                ncols(correction.factor) == n &&
                nrows(correction.inverse_factor) == n &&
                ncols(correction.inverse_factor) == n &&
                correction.factor * correction.inverse_factor == identity &&
                correction.inverse_factor * correction.factor == identity

            if classification == :one
                correction_ok = correction_ok &&
                    correction.kind == :identity &&
                    metadata.deferred_diagonal_correction === nothing &&
                    metadata.normalized_deferred_core == deferred
            elseif classification == :laurent_monomial_unit
                correction_ok = correction_ok &&
                    correction.kind == :left_diagonal_determinant_correction &&
                    metadata.deferred_diagonal_correction == correction &&
                    correction.factor * metadata.normalized_deferred_core == deferred
            end

            normalized_core_det_ok =
                metadata.normalized_deferred_core !== nothing &&
                det(metadata.normalized_deferred_core) == one(R)
        else
            boundary = metadata.staged_boundary
            boundary_ok =
                metadata.deferred_correction === nothing &&
                metadata.deferred_diagonal_correction === nothing &&
                metadata.normalized_deferred_core === nothing &&
                boundary !== nothing &&
                boundary.kind == :unsupported_deferred_laurent_determinant &&
                boundary.determinant_source == :deferred_submatrix &&
                boundary.overall_determinant == metadata.overall_determinant &&
                boundary.determinant_classification == metadata.determinant_classification &&
                boundary.deferred_submatrix_size == (n, n) &&
                boundary.supported == false
        end
    catch err
        err isa InterruptException && rethrow()
    end

    supported_ok = metadata.supported ? (correction_ok && normalized_core_det_ok) : boundary_ok
    overall_ok = certificate_ok && determinant_ok && supported_ok
    return (;
        overall_ok,
        certificate_ok,
        determinant_ok,
        correction_ok,
        normalized_core_det_ok,
        boundary_ok,
    )
end

function _laurent_column_peel_recursive(
    current;
    progress_callback = nothing,
    completed_steps::Int = 0,
    last_completed = _laurent_column_peel_empty_last_completed(),
    emit_current::Bool = true,
)
    d = nrows(current)
    R = base_ring(current)
    emit_current &&
        _emit_laurent_column_peel_progress(progress_callback, current, completed_steps, last_completed)

    if d == 2
        final_block = current
        final_target = block_embedding(final_block, 3, [1, 2])
        X = _laurent_column_peel_generator(R)
        final_local = realize_sl3_local(final_target, X; check_monic = false)
        final_factors = _project_local_sl3_to_2x2(final_local, R)
        verify_factorization(final_block, final_factors) ||
            error("internal Laurent column-peel final 2x2 factorization failed exact verification")
        return final_factors, LaurentColumnPeelStep[], final_block, final_target, final_local, final_factors
    end

    step_started_at = time()
    step = _laurent_column_peel_step(current)
    step_elapsed = round(time() - step_started_at; digits = 3)
    next_completed_steps = completed_steps + 1
    next_last_completed = (;
        dimension = d,
        elapsed_seconds = step_elapsed,
        left_factors = length(step.left_factors),
        right_factors = length(step.right_factors),
    )
    next_factors, next_steps, final_block, final_target, final_local, final_2x2 =
        _laurent_column_peel_recursive(
            step.next_block;
            progress_callback,
            completed_steps = next_completed_steps,
            last_completed = next_last_completed,
            emit_current = true,
        )
    current_factors = vcat(
        _inverse_elementary_sequence(step.left_factors),
        _embed_upper_left_factors(next_factors, R, d),
        _inverse_elementary_sequence(step.right_factors),
    )
    return current_factors,
        vcat(LaurentColumnPeelStep[step], next_steps),
        final_block,
        final_target,
        final_local,
        final_2x2
end

function _laurent_column_peel_step(current)
    R = base_ring(current)
    d = nrows(current)
    last_column = [current[row, d] for row in 1:d]
    target_column = _column_peel_target_column(R, d)
    left_factors = matrix(R, d, 1, last_column) == target_column ?
        typeof(identity_matrix(R, d))[] :
        reduce_unimodular_column(last_column, R)
    left_product = _factor_product(left_factors, R, d)
    after_left = left_product * current
    recorded_column = matrix(R, d, 1, last_column)
    left_product * recorded_column == target_column ||
        throw(ArgumentError("Laurent column-peel left factors failed to send the last column to e_d"))
    after_left[:, d:d] == target_column ||
        throw(ArgumentError("Laurent column-peel left product failed to normalize the last column"))
    right_factors = _expected_column_peel_right_factors(after_left, d, R)
    right_product = _factor_product(right_factors, R, d)
    peeled = after_left * right_product
    next_block = matrix(R, [peeled[row, col] for row in 1:(d - 1), col in 1:(d - 1)])
    _is_valid_laurent_column_peel_step_data(
        d,
        current,
        last_column,
        left_factors,
        after_left,
        right_factors,
        peeled,
        next_block,
    ) || throw(ArgumentError("Laurent column-peel step failed exact replay"))
    return LaurentColumnPeelStep(d, current, last_column, left_factors, after_left, right_factors, peeled, next_block)
end

function _column_peel_target_column(R, d::Int)
    return matrix(R, d, 1, [row == d ? one(R) : zero(R) for row in 1:d])
end

function _expected_column_peel_right_factors(after_left, d::Int, R)
    right_factors = typeof(identity_matrix(R, d))[]
    for j in 1:(d - 1)
        coeff = -after_left[d, j]
        coeff == zero(R) && continue
        push!(right_factors, elementary_matrix(d, d, j, coeff, R))
    end
    return right_factors
end

function _is_valid_laurent_column_peel_step_data(d::Int, current, last_column, left_factors, after_left, right_factors, peeled, next_block)::Bool
    R = base_ring(current)
    actual_last_column = [current[row, d] for row in 1:d]
    last_column == actual_last_column || return false
    recorded_column = matrix(R, d, 1, last_column)
    target_column = _column_peel_target_column(R, d)
    left_product = try
        _factor_product(left_factors, R, d)
    catch err
        err isa InterruptException && rethrow()
        return false
    end
    left_product * current == after_left || return false
    left_product * recorded_column == target_column || return false
    after_left[:, d:d] == target_column || return false

    right_factors == _expected_column_peel_right_factors(after_left, d, R) || return false
    right_product = try
        _factor_product(right_factors, R, d)
    catch err
        err isa InterruptException && rethrow()
        return false
    end
    after_left * right_product == peeled || return false
    peeled[d, d] == one(R) || return false
    all(peeled[row, d] == zero(R) for row in 1:(d - 1)) || return false
    all(peeled[d, col] == zero(R) for col in 1:(d - 1)) || return false
    expected_next = matrix(R, [peeled[row, col] for row in 1:(d - 1), col in 1:(d - 1)])
    next_block == expected_next || return false
    return peeled == block_embedding(next_block, d, collect(1:(d - 1)))
end

function _inverse_elementary_sequence(factors)
    collected = collect(factors)
    isempty(collected) && return collected

    inverses = typeof(first(collected))[]
    for factor in Iterators.reverse(collected)
        row, col, coeff = _elementary_factor_data(factor)
        push!(inverses, elementary_matrix(nrows(factor), row, col, -coeff, base_ring(factor)))
    end
    return inverses
end

function _elementary_factor_data(factor)
    nrows(factor) == ncols(factor) || throw(ArgumentError("elementary factor must be square"))
    R = base_ring(factor)
    row = 0
    col = 0
    coeff = zero(R)

    for i in 1:nrows(factor), j in 1:ncols(factor)
        entry = factor[i, j]
        if i == j
            entry == one(R) || throw(ArgumentError("elementary factor diagonal must be one"))
            continue
        end
        entry == zero(R) && continue
        row == 0 || throw(ArgumentError("elementary factor must have a unique off-diagonal entry"))
        row = i
        col = j
        coeff = entry
    end

    row != 0 || throw(ArgumentError("elementary factor must have a nonzero off-diagonal entry"))
    return row, col, coeff
end

function _embed_upper_left_factors(factors, R, d::Int)
    collected = collect(factors)
    isempty(collected) && return typeof(identity_matrix(R, d))[]
    return [block_embedding(factor, d, collect(1:(d - 1))) for factor in collected]
end

function _embed_laurent_deferred_peel_factors(factors, R, original_dimension::Int, factor_dimension::Int)
    collected = collect(factors)
    if isempty(collected)
        return typeof(identity_matrix(R, original_dimension))[]
    end
    factor_dimension == original_dimension && return collected
    return [block_embedding(factor, original_dimension, collect(1:factor_dimension)) for factor in collected]
end

function _laurent_determinant_deferred_target(deferred_submatrix, original_dimension::Int)
    return block_embedding(deferred_submatrix, original_dimension, collect(1:nrows(deferred_submatrix)))
end

function _laurent_determinant_deferred_peel_replay_data(original_matrix, peel_steps, deferred_submatrix)
    R = base_ring(original_matrix)
    original_dimension = nrows(original_matrix)
    left_factors = typeof(identity_matrix(R, original_dimension))[]
    right_factors = typeof(identity_matrix(R, original_dimension))[]
    steps = collect(peel_steps)

    for step in Iterators.reverse(steps)
        append!(
            left_factors,
            _embed_laurent_deferred_peel_factors(
                step.left_factors,
                R,
                original_dimension,
                step.dimension,
            ),
        )
    end

    for step in steps
        append!(
            right_factors,
            _embed_laurent_deferred_peel_factors(
                step.right_factors,
                R,
                original_dimension,
                step.dimension,
            ),
        )
    end

    left_product = _factor_product(left_factors, R, original_dimension)
    right_product = _factor_product(right_factors, R, original_dimension)
    target_matrix = _laurent_determinant_deferred_target(deferred_submatrix, original_dimension)
    return (;
        left_factors,
        right_factors,
        left_product,
        right_product,
        target_matrix,
    )
end

function _laurent_column_peel_generator(R)
    ring_gens = collect(gens(R))
    isempty(ring_gens) && throw(ArgumentError("Laurent column-peel factorization requires a Laurent generator"))
    return ring_gens[1]
end

function _replay_inverse_or_original(factors)
    collected = collect(factors)
    try
        return _inverse_elementary_sequence(collected)
    catch err
        err isa InterruptException && rethrow()
        return collected
    end
end

function _replay_laurent_column_peel_factors(peel_steps, final_factors, R)
    replayed = collect(final_factors)

    for step in Iterators.reverse(collect(peel_steps))
        replayed = vcat(
            _replay_inverse_or_original(step.left_factors),
            _embed_upper_left_factors(replayed, R, step.dimension),
            _replay_inverse_or_original(step.right_factors),
        )
    end

    return replayed
end

function _factor_sequences_equal(left, right)::Bool
    length(left) == length(right) || return false
    return all(left[i] == right[i] for i in eachindex(left))
end

function _project_local_sl3_to_2x2(factors, R)
    collected = collect(factors)
    projected = typeof(matrix(R, 2, 2, [one(R), zero(R), zero(R), one(R)]))[]

    for factor in collected
        nrows(factor) == 3 || throw(ArgumentError("local SL_3 factor must be 3x3"))
        ncols(factor) == 3 || throw(ArgumentError("local SL_3 factor must be 3x3"))
        factor[1, 3] == zero(R) || throw(ArgumentError("local SL_3 factor touched the trailing identity column"))
        factor[2, 3] == zero(R) || throw(ArgumentError("local SL_3 factor touched the trailing identity column"))
        factor[3, 1] == zero(R) || throw(ArgumentError("local SL_3 factor touched the trailing identity row"))
        factor[3, 2] == zero(R) || throw(ArgumentError("local SL_3 factor touched the trailing identity row"))
        factor[3, 3] == one(R) || throw(ArgumentError("local SL_3 factor must preserve the trailing identity"))
        push!(projected, matrix(R, [
            factor[1, 1] factor[1, 2];
            factor[2, 1] factor[2, 2];
        ]))
    end

    return projected
end

function _laurent_column_peel_verification(certificate)
    R = base_ring(certificate.original_matrix)
    size_ok = nrows(certificate.original_matrix) == ncols(certificate.original_matrix) >= 2
    step_chain_ok = _laurent_column_peel_step_chain_ok(certificate.peel_steps, certificate.original_matrix)
    steps_ok = try
        all(step -> _is_valid_laurent_column_peel_step_data(
                step.dimension,
                step.input_matrix,
                step.last_column,
                step.left_factors,
                step.after_left_matrix,
                step.right_factors,
                step.peeled_matrix,
                step.next_block,
            ),
            certificate.peel_steps,
        )
    catch err
        err isa InterruptException && rethrow()
        false
    end
    final_metadata_ok = try
        certificate.final_local_target == block_embedding(certificate.final_block, 3, [1, 2]) &&
            certificate.final_factors == _project_local_sl3_to_2x2(certificate.final_local_factors, R)
    catch err
        err isa InterruptException && rethrow()
        false
    end
    final_local_ok = try
        verify_factorization(certificate.final_local_target, certificate.final_local_factors)
    catch err
        err isa InterruptException && rethrow()
        false
    end
    final_factors_ok = try
        verify_factorization(certificate.final_block, certificate.final_factors)
    catch err
        err isa InterruptException && rethrow()
        false
    end
    factor_sequence_ok = try
        replayed_factors = _replay_laurent_column_peel_factors(
            certificate.peel_steps,
            certificate.final_factors,
            R,
        )
        _factor_sequences_equal(certificate.factors, replayed_factors)
    catch err
        err isa InterruptException && rethrow()
        false
    end
    product_ok = try
        certificate.product == _factor_product(certificate.factors, R, nrows(certificate.original_matrix))
    catch err
        err isa InterruptException && rethrow()
        false
    end
    factors_ok = try
        verify_factorization(certificate.original_matrix, certificate.factors)
    catch err
        err isa InterruptException && rethrow()
        false
    end
    overall_ok = size_ok && step_chain_ok && steps_ok && final_metadata_ok &&
        final_local_ok && final_factors_ok && factor_sequence_ok && product_ok && factors_ok
    return (
        overall_ok = overall_ok,
        size_ok = size_ok,
        step_chain_ok = step_chain_ok,
        steps_ok = steps_ok,
        final_metadata_ok = final_metadata_ok,
        final_local_ok = final_local_ok,
        final_factors_ok = final_factors_ok,
        factor_sequence_ok = factor_sequence_ok,
        product_ok = product_ok,
        factors_ok = factors_ok,
    )
end

function _laurent_determinant_deferred_peel_verification(certificate)
    R = base_ring(certificate.original_matrix)
    n = nrows(certificate.original_matrix)
    size_ok = nrows(certificate.original_matrix) == ncols(certificate.original_matrix) >= 3
    determinant_source_ok = certificate.determinant_source == :deferred_submatrix
    step_chain_ok = _laurent_determinant_deferred_step_chain_ok(
        certificate.peel_steps,
        certificate.original_matrix,
        certificate.deferred_submatrix,
    )
    steps_ok = try
        all(step -> _is_valid_laurent_column_peel_step_data(
                step.dimension,
                step.input_matrix,
                step.last_column,
                step.left_factors,
                step.after_left_matrix,
                step.right_factors,
                step.peeled_matrix,
                step.next_block,
            ),
            certificate.peel_steps,
        )
    catch err
        err isa InterruptException && rethrow()
        false
    end
    deferred_shape_ok = nrows(certificate.deferred_submatrix) == ncols(certificate.deferred_submatrix) &&
        nrows(certificate.deferred_submatrix) >= 2 &&
        nrows(certificate.deferred_submatrix) < n
    replay = try
        _laurent_determinant_deferred_peel_replay_data(
            certificate.original_matrix,
            certificate.peel_steps,
            certificate.deferred_submatrix,
        )
    catch err
        err isa InterruptException && rethrow()
        nothing
    end
    replay_metadata_ok = replay !== nothing &&
        _factor_sequences_equal(certificate.left_factors, replay.left_factors) &&
        _factor_sequences_equal(certificate.right_factors, replay.right_factors) &&
        certificate.left_product == replay.left_product &&
        certificate.right_product == replay.right_product &&
        certificate.target_matrix == replay.target_matrix
    target_ok = replay !== nothing &&
        certificate.target_matrix == _laurent_determinant_deferred_target(certificate.deferred_submatrix, n)
    relation_ok = try
        certificate.left_product * certificate.original_matrix * certificate.right_product ==
            certificate.target_matrix
    catch err
        err isa InterruptException && rethrow()
        false
    end
    overall_ok = size_ok && determinant_source_ok && step_chain_ok && steps_ok &&
        deferred_shape_ok && replay_metadata_ok && target_ok && relation_ok
    return (;
        overall_ok,
        size_ok,
        determinant_source_ok,
        step_chain_ok,
        steps_ok,
        deferred_shape_ok,
        replay_metadata_ok,
        target_ok,
        relation_ok,
    )
end

function _laurent_column_peel_step_chain_ok(steps, original_matrix)::Bool
    collected = collect(steps)
    current = original_matrix

    for step in collected
        step.input_matrix == current || return false
        step.dimension == nrows(current) || return false
        step.dimension == ncols(current) || return false
        current = step.next_block
    end

    return nrows(current) == 2 && ncols(current) == 2
end

function _laurent_determinant_deferred_step_chain_ok(steps, original_matrix, deferred_submatrix)::Bool
    collected = collect(steps)
    isempty(collected) && return false
    current = original_matrix

    for step in collected
        step.input_matrix == current || return false
        step.dimension == nrows(current) || return false
        step.dimension == ncols(current) || return false
        current = step.next_block
    end

    return current == deferred_submatrix
end
