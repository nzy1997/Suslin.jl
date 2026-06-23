struct PolynomialColumnPeelStep
    dimension::Int
    input_matrix
    last_column::Vector
    left_factors::Vector
    after_left_matrix
    right_factors::Vector
    peeled_matrix
    next_block
end

struct PolynomialColumnPeelCertificate
    original_matrix
    peel_steps::Vector{PolynomialColumnPeelStep}
    final_block
    final_certificate
    final_factors::Vector
    factors::Vector
    product
    verification
end

function _polynomial_column_peel_certificate(A; final_route=nothing)
    _validate_polynomial_column_peel_input(A)
    _validate_polynomial_column_peel_final_route(final_route)

    peel_steps, final_block, final_certificate, final_factors =
        _polynomial_column_peel_recursive(A; final_route=final_route)
    isempty(peel_steps) &&
        throw(ArgumentError("polynomial column-peel certificate requires at least one real peel step"))
    R = base_ring(A)
    factors = _replay_polynomial_column_peel_factors(peel_steps, final_factors, R)
    product = _factor_product(factors, R, nrows(A))

    certificate = PolynomialColumnPeelCertificate(
        A,
        peel_steps,
        final_block,
        final_certificate,
        final_factors,
        factors,
        product,
        nothing,
    )
    verification = _polynomial_column_peel_core_verification(certificate)
    verification.overall_core_ok || error("internal polynomial column-peel verification failed")
    return PolynomialColumnPeelCertificate(
        certificate.original_matrix,
        certificate.peel_steps,
        certificate.final_block,
        certificate.final_certificate,
        certificate.final_factors,
        certificate.factors,
        certificate.product,
        verification,
    )
end

function _verify_polynomial_column_peel_certificate(cert)::Bool
    try
        return _polynomial_column_peel_verification(cert).overall_ok
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end

function _polynomial_column_peel_verification(cert)
    core = _polynomial_column_peel_core_verification(cert)
    stored_verification_ok = cert.verification == core
    return merge(core, (;
        stored_verification_ok,
        overall_ok = core.overall_core_ok && stored_verification_ok,
    ))
end

function _validate_polynomial_column_peel_input(A)
    nrows(A) == ncols(A) ||
        throw(ArgumentError("polynomial column-peel certificate requires a square matrix"))
    nrows(A) >= 3 ||
        throw(ArgumentError("polynomial column-peel certificate requires size at least 3"))
    R = base_ring(A)
    _is_laurent_polynomial_ring(R) &&
        throw(ArgumentError("polynomial column-peel certificate does not support Laurent polynomial rings"))
    _factorization_ring_profile(R) == :polynomial ||
        throw(ArgumentError("polynomial column-peel certificate requires an exact ordinary polynomial ring"))
    det(A) == one(R) ||
        throw(ArgumentError("polynomial column-peel certificate requires determinant-one input"))
    return nrows(A)
end

function _validate_polynomial_column_peel_final_route(final_route)
    final_route === nothing && return nothing
    final_route isa Symbol ||
        throw(ArgumentError("polynomial column-peel final route must be a Symbol"))
    final_route in (:fast_local_sl3, :disjoint_local_blocks) ||
        throw(ArgumentError("polynomial column-peel final route must be :fast_local_sl3 or :disjoint_local_blocks"))
    return nothing
end

function _polynomial_column_peel_recursive(current; final_route=nothing)
    d = nrows(current)
    final_certificate = _polynomial_column_peel_try_final_route(current; final_route=final_route)
    if final_certificate !== nothing
        return PolynomialColumnPeelStep[],
            current,
            final_certificate,
            copy(final_certificate.factors)
    end

    d == 3 &&
        throw(ArgumentError("polynomial column-peel certificate requires a supported final route at size 3"))

    step = _polynomial_column_peel_step(current)
    next_steps, final_block, next_certificate, final_factors =
        _polynomial_column_peel_recursive(step.next_block; final_route=final_route)
    return vcat(PolynomialColumnPeelStep[step], next_steps),
        final_block,
        next_certificate,
        final_factors
end

function _polynomial_column_peel_try_final_route(current; final_route=nothing)
    candidate_routes =
        final_route === nothing ?
        (:fast_local_sl3, :disjoint_local_blocks) :
        (final_route,)

    for route in candidate_routes
        certificate = try
            _polynomial_factorization_route_certificate(
                current;
                route=route,
                allow_recursive_column_peel=false,
            )
        catch err
            err isa InterruptException && rethrow()
            err isa ArgumentError || rethrow()
            nothing
        end

        if certificate !== nothing &&
                certificate.status == :supported &&
                certificate.route in (:fast_local_sl3, :disjoint_local_blocks) &&
                certificate.matrix != identity_matrix(base_ring(current), nrows(current))
            return certificate
        end
    end

    return nothing
end

function _polynomial_column_peel_step(current)
    R = base_ring(current)
    d = nrows(current)
    last_column = [current[row, d] for row in 1:d]
    left_factors = reduce_unimodular_column(last_column, R)
    left_product = _factor_product(left_factors, R, d)
    recorded_column = matrix(R, d, 1, last_column)
    target_column = _column_peel_target_column(R, d)
    left_product * recorded_column == target_column ||
        throw(ArgumentError("polynomial column-peel left factors failed to send the last column to e_d"))
    after_left = left_product * current
    after_left[:, d:d] == target_column ||
        throw(ArgumentError("polynomial column-peel left product failed to normalize the last column"))
    right_factors = _expected_column_peel_right_factors(after_left, d, R)
    right_product = _factor_product(right_factors, R, d)
    peeled = after_left * right_product
    next_block = matrix(R, [peeled[row, col] for row in 1:(d - 1), col in 1:(d - 1)])
    _is_valid_polynomial_column_peel_step_data(
        d,
        current,
        last_column,
        left_factors,
        after_left,
        right_factors,
        peeled,
        next_block,
    ) || throw(ArgumentError("polynomial column-peel step failed exact replay"))
    return PolynomialColumnPeelStep(
        d,
        current,
        last_column,
        left_factors,
        after_left,
        right_factors,
        peeled,
        next_block,
    )
end

function _is_valid_polynomial_column_peel_step_data(
    d::Int,
    current,
    last_column,
    left_factors,
    after_left,
    right_factors,
    peeled,
    next_block,
)::Bool
    d >= 4 || return false
    nrows(current) == d && ncols(current) == d || return false
    R = base_ring(current)
    _is_laurent_polynomial_ring(R) && return false
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

function _replay_polynomial_column_peel_factors(peel_steps, final_factors, R)
    replayed = collect(final_factors)

    for step in Iterators.reverse(collect(peel_steps))
        replayed = vcat(
            _inverse_elementary_sequence(step.left_factors),
            _embed_upper_left_factors(replayed, R, step.dimension),
            _inverse_elementary_sequence(step.right_factors),
        )
    end

    return replayed
end

function _polynomial_column_peel_core_verification(cert)
    preconditions_ok = _polynomial_column_peel_preconditions_ok(cert)
    step_chain_ok = _polynomial_column_peel_step_chain_ok(
        cert.peel_steps,
        cert.original_matrix,
        cert.final_block,
    )
    steps_ok = try
        all(step -> _is_valid_polynomial_column_peel_step_data(
                step.dimension,
                step.input_matrix,
                step.last_column,
                step.left_factors,
                step.after_left_matrix,
                step.right_factors,
                step.peeled_matrix,
                step.next_block,
            ),
            cert.peel_steps,
        )
    catch err
        err isa InterruptException && rethrow()
        false
    end
    final_certificate_ok = _polynomial_column_peel_final_certificate_ok(cert)
    factor_sequence_ok = try
        replayed_factors = _replay_polynomial_column_peel_factors(
            cert.peel_steps,
            cert.final_factors,
            base_ring(cert.original_matrix),
        )
        _factor_sequences_equal(cert.factors, replayed_factors)
    catch err
        err isa InterruptException && rethrow()
        false
    end
    product_ok = try
        cert.product == _factor_product(cert.factors, base_ring(cert.original_matrix), nrows(cert.original_matrix))
    catch err
        err isa InterruptException && rethrow()
        false
    end
    factors_ok = try
        verify_factorization(cert.original_matrix, cert.factors)
    catch err
        err isa InterruptException && rethrow()
        false
    end
    overall_core_ok = preconditions_ok && step_chain_ok && steps_ok && final_certificate_ok &&
        factor_sequence_ok && product_ok && factors_ok
    return (
        overall_core_ok=overall_core_ok,
        preconditions_ok=preconditions_ok,
        step_chain_ok=step_chain_ok,
        steps_ok=steps_ok,
        final_certificate_ok=final_certificate_ok,
        factor_sequence_ok=factor_sequence_ok,
        product_ok=product_ok,
        factors_ok=factors_ok,
    )
end

function _polynomial_column_peel_preconditions_ok(cert)::Bool
    try
        A = cert.original_matrix
        nrows(A) == ncols(A) || return false
        nrows(A) >= 3 || return false
        isempty(cert.peel_steps) && return false
        R = base_ring(A)
        _is_laurent_polynomial_ring(R) && return false
        _factorization_ring_profile(R) == :polynomial || return false
        det(A) == one(R) || return false
        return true
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end

function _polynomial_column_peel_step_chain_ok(steps, original_matrix, final_block)::Bool
    collected = collect(steps)
    isempty(collected) && return false
    current = original_matrix

    for step in collected
        step.input_matrix == current || return false
        step.dimension == nrows(current) || return false
        step.dimension == ncols(current) || return false
        current = step.next_block
    end

    current == final_block || return false
    nrows(final_block) == ncols(final_block) || return false
    return nrows(final_block) >= 3
end

function _polynomial_column_peel_final_certificate_ok(cert)::Bool
    try
        final_certificate = cert.final_certificate
        final_certificate isa PolynomialFactorizationRouteCertificate || return false
        final_certificate.status == :supported || return false
        final_certificate.route in (:fast_local_sl3, :disjoint_local_blocks) || return false
        final_certificate.matrix == cert.final_block || return false
        _verify_polynomial_factorization_route_certificate(final_certificate) || return false
        _factor_sequences_equal(cert.final_factors, final_certificate.factors) || return false
        cert.final_block == cert.original_matrix && return false
        return verify_factorization(cert.final_block, cert.final_factors)
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end
