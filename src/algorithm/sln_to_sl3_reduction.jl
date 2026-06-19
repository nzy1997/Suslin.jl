struct SL3LocalObligation
    block_location::Vector{Int}
    ring
    target_local_matrix
    required_assumptions::Vector{Symbol}
    embedded_target
    local_factors::Vector
    embedded_factors::Vector
    reassembly_data
end

struct SLNToSL3Reduction
    ring
    size::Int
    original_matrix
    normalized_matrix
    normalization
    obligations::Vector{SL3LocalObligation}
    factors::Vector
    product
    verification
end

function reduce_sln_to_sl3(A; block_locations=nothing)
    return _construct_sln_to_sl3_reduction(A, block_locations)
end

function verify_sln_to_sl3_reduction(reduction::SLNToSL3Reduction)::Bool
    return _sln_to_sl3_reduction_verification(reduction).overall_ok
end

function _construct_sln_to_sl3_reduction(A, block_locations)
    n = _validate_factorization_matrix(A)
    R = base_ring(A)
    ring_profile = _factorization_ring_profile(R)
    normalized_A, normalization = _normalize_factorization_input(A, ring_profile)
    normalized_R = base_ring(normalized_A)

    if ring_profile == :polynomial
        _require_polynomial_sl_determinant(normalized_A)
    end

    X = _reduction_generator(normalized_R, ring_profile)
    locations = _normalize_reduction_block_locations(n, block_locations)
    obligations = SL3LocalObligation[]
    for indices in locations
        _is_identity_local_block(normalized_A, normalized_R, indices) && continue
        push!(obligations, _build_sl3_local_obligation(normalized_A, normalized_R, indices, X))
    end

    factors = _compose_obligation_factors(obligations, normalized_R, n)
    product = _factor_product(factors, normalized_R, n)
    product == normalized_A || _throw_staged_sln_to_sl3_failure("embedded local SL_3 obligations did not exactly reconstruct the normalized matrix")

    reduction = SLNToSL3Reduction(
        normalized_R,
        n,
        A,
        normalized_A,
        normalization,
        obligations,
        factors,
        product,
        nothing,
    )

    verification = _sln_to_sl3_reduction_verification(reduction)
    verified_reduction = SLNToSL3Reduction(
        reduction.ring,
        reduction.size,
        reduction.original_matrix,
        reduction.normalized_matrix,
        reduction.normalization,
        reduction.obligations,
        reduction.factors,
        reduction.product,
        verification,
    )

    verification.overall_ok || error("internal SL_n to local SL_3 reduction verification failed")
    return verified_reduction
end

function _reduction_generator(R, ring_profile::Symbol)
    if ring_profile == :laurent
        _throw_staged_sln_to_sl3_failure("Laurent SL_n to local SL_3 reduction is not yet implemented")
    end

    ring_gens = collect(gens(R))
    length(ring_gens) == 1 || throw(ArgumentError("ordinary polynomial reduction currently requires a univariate base ring"))
    return ring_gens[1]
end

function _normalize_reduction_block_locations(n::Int, block_locations)
    locations = if block_locations === nothing
        [Int[first_idx, first_idx + 1, first_idx + 2] for first_idx in 1:3:(n - 2)]
    else
        [Int[index for index in location] for location in block_locations]
    end

    seen = Set{Int}()
    for location in locations
        length(location) == 3 || throw(ArgumentError("block locations must consist of 3 indices per local SL_3 obligation"))
        length(unique(location)) == 3 || throw(ArgumentError("block locations must contain distinct indices"))
        for index in location
            1 <= index <= n || throw(ArgumentError("block locations must lie within the matrix size"))
            index in seen && throw(ArgumentError("block locations must be pairwise disjoint"))
            push!(seen, index)
        end
    end

    return locations
end

function _build_sl3_local_obligation(A, R, indices::Vector{Int}, X)
    local_target = _principal_submatrix(A, indices)
    local_factors = try
        realize_sl3_local(local_target, X)
    catch err
        err isa InterruptException && rethrow()
        _throw_staged_sln_to_sl3_failure("failed to solve local SL_3 obligation on block $(indices)")
    end
    embedded_target = block_embedding(local_target, nrows(A), indices)
    embedded_factors = embed_factor_sequence(local_factors, nrows(A), indices)

    reassembly_data = _sl3_obligation_reassembly_data(
        local_factors,
        embedded_factors,
        R,
        nrows(A),
        local_target,
        embedded_target,
    )

    return SL3LocalObligation(
        copy(indices),
        R,
        local_target,
        Symbol[:univariate_base_ring, :determinant_one],
        embedded_target,
        local_factors,
        embedded_factors,
        reassembly_data,
    )
end

function _is_identity_local_block(A, R, indices::Vector{Int})::Bool
    return _principal_submatrix(A, indices) == identity_matrix(R, 3)
end

function _principal_submatrix(A, indices::Vector{Int})
    R = base_ring(A)
    return matrix(R, [A[row, col] for row in indices, col in indices])
end

function _compose_obligation_factors(obligations::AbstractVector, R, n::Int)
    isempty(obligations) && return typeof(identity_matrix(R, n))[]
    return compose_factor_sequences((obligation.embedded_factors for obligation in obligations)...)
end

function _same_factor_sequence(left, right)::Bool
    length(left) == length(right) || return false
    for idx in eachindex(left)
        left[idx] == right[idx] || return false
    end
    return true
end

function _sl3_obligation_reassembly_data(local_factors, embedded_factors, R, n::Int, local_target, embedded_target)
    local_product = _factor_product(local_factors, R, 3)
    embedded_product = _factor_product(embedded_factors, R, n)
    return (
        local_product = local_product,
        local_product_ok = local_product == local_target,
        embedded_product = embedded_product,
        embedded_product_ok = embedded_product == embedded_target,
    )
end

function _factor_product(factors, R, n::Int)
    product = identity_matrix(R, n)
    for factor in factors
        product *= factor
    end
    return product
end

function _sln_to_sl3_reduction_verification(reduction::SLNToSL3Reduction)
    obligation_locations_ok = try
        _normalize_reduction_block_locations(
            reduction.size,
            [obligation.block_location for obligation in reduction.obligations],
        )
        true
    catch err
        err isa InterruptException && rethrow()
        false
    end
    obligations_ok = try
        all(obligation -> _verify_sl3_local_obligation(
                obligation,
                reduction.normalized_matrix,
                reduction.ring,
                reduction.size,
            ),
            reduction.obligations,
        )
    catch err
        err isa InterruptException && rethrow()
        false
    end
    factors_ok = try
        verify_factorization(reduction.normalized_matrix, reduction.factors)
    catch err
        err isa InterruptException && rethrow()
        false
    end
    product_ok = try
        reduction.product == reduction.normalized_matrix
    catch err
        err isa InterruptException && rethrow()
        false
    end
    normalization_ok = try
        reduction.normalization === nothing || verify_laurent_gl_normalization(reduction.original_matrix, reduction.normalization)
    catch err
        err isa InterruptException && rethrow()
        false
    end
    original_reconstruction_ok = try
        reduction.normalization === nothing ?
            reduction.product == reduction.original_matrix :
            reduction.normalization.correction.factor * reduction.product == reduction.original_matrix
    catch err
        err isa InterruptException && rethrow()
        false
    end
    obligation_factors_ok = try
        _same_factor_sequence(
            _compose_obligation_factors(reduction.obligations, reduction.ring, reduction.size),
            reduction.factors,
        )
    catch err
        err isa InterruptException && rethrow()
        false
    end

    return (
        obligation_locations_ok = obligation_locations_ok,
        obligations_ok = obligations_ok,
        obligation_factors_ok = obligation_factors_ok,
        factors_ok = factors_ok,
        product_ok = product_ok,
        normalization_ok = normalization_ok,
        original_reconstruction_ok = original_reconstruction_ok,
        overall_ok = obligation_locations_ok && obligations_ok && obligation_factors_ok && factors_ok && product_ok && normalization_ok && original_reconstruction_ok,
    )
end

function _verify_sl3_local_obligation(obligation::SL3LocalObligation, normalized_matrix, R, n::Int)::Bool
    try
        _same_base_ring(obligation.ring, R) || return false
        obligation.required_assumptions == Symbol[:univariate_base_ring, :determinant_one] || return false

        target_local_matrix = _principal_submatrix(normalized_matrix, obligation.block_location)
        target_local_matrix == obligation.target_local_matrix || return false

        embedded_target = block_embedding(target_local_matrix, n, obligation.block_location)
        embedded_target == obligation.embedded_target || return false

        local_product = _factor_product(obligation.local_factors, obligation.ring, 3)
        embedded_product = _factor_product(obligation.embedded_factors, obligation.ring, n)
        reassembly_data = _sl3_obligation_reassembly_data(
            obligation.local_factors,
            obligation.embedded_factors,
            obligation.ring,
            n,
            obligation.target_local_matrix,
            obligation.embedded_target,
        )

        return local_product == obligation.target_local_matrix &&
               embedded_product == obligation.embedded_target &&
               reassembly_data == obligation.reassembly_data
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end

function _throw_staged_sln_to_sl3_failure(reason::AbstractString)
    throw(ArgumentError("staged SL_n to local SL_3 reduction failure: $(reason)"))
end
