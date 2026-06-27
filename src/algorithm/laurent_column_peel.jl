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

function _validate_laurent_column_peel_input(A)
    nrows(A) == ncols(A) || throw(ArgumentError("Laurent column-peel factorization requires a square matrix"))
    nrows(A) >= 2 || throw(ArgumentError("Laurent column-peel factorization requires size at least 2"))
    R = base_ring(A)
    _is_laurent_polynomial_ring(R) || throw(ArgumentError("Laurent column-peel factorization requires a Laurent polynomial ring"))
    profile = classify_laurent_determinant(A)
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

function _factor_laurent_sl_column_peel(A; progress_callback = nothing)
    _validate_laurent_column_peel_input(A)
    factors, steps, final_block, final_target, final_local, final_2x2 =
        _laurent_column_peel_recursive(
            A;
            progress_callback,
            completed_steps = 0,
            last_completed = _laurent_column_peel_empty_last_completed(),
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

function _laurent_column_peel_recursive(
    current;
    progress_callback = nothing,
    completed_steps::Int = 0,
    last_completed = _laurent_column_peel_empty_last_completed(),
)
    d = nrows(current)
    R = base_ring(current)
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
    left_factors = reduce_unimodular_column(last_column, R)
    left_product = _factor_product(left_factors, R, d)
    after_left = left_product * current
    recorded_column = matrix(R, d, 1, last_column)
    left_product * recorded_column == _column_peel_target_column(R, d) ||
        throw(ArgumentError("Laurent column-peel left factors failed to send the last column to e_d"))
    after_left[:, d:d] == _column_peel_target_column(R, d) ||
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
