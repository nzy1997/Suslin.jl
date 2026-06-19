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
    obligations = [
        _build_sl3_local_obligation(normalized_A, normalized_R, indices, X)
        for indices in locations
    ]

    factors = compose_factor_sequences((obligation.embedded_factors for obligation in obligations)...)
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

    local_product = _factor_product(local_factors, R, 3)
    embedded_product = _factor_product(embedded_factors, R, nrows(A))
    reassembly_data = (
        local_product = local_product,
        local_product_ok = local_product == local_target,
        embedded_product = embedded_product,
        embedded_product_ok = embedded_product == embedded_target,
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

function _principal_submatrix(A, indices::Vector{Int})
    R = base_ring(A)
    return matrix(R, [A[row, col] for row in indices, col in indices])
end

function _factor_product(factors, R, n::Int)
    product = identity_matrix(R, n)
    for factor in factors
        product *= factor
    end
    return product
end

function _sln_to_sl3_reduction_verification(reduction::SLNToSL3Reduction)
    obligations_ok = all(_verify_sl3_local_obligation, reduction.obligations)
    factors_ok = verify_factorization(reduction.normalized_matrix, reduction.factors)
    product_ok = reduction.product == reduction.normalized_matrix
    normalization_ok = reduction.normalization === nothing || verify_laurent_gl_normalization(reduction.original_matrix, reduction.normalization)
    original_reconstruction_ok = reduction.normalization === nothing ?
        reduction.product == reduction.original_matrix :
        reduction.normalization.correction.factor * reduction.product == reduction.original_matrix

    return (
        obligations_ok = obligations_ok,
        factors_ok = factors_ok,
        product_ok = product_ok,
        normalization_ok = normalization_ok,
        original_reconstruction_ok = original_reconstruction_ok,
        overall_ok = obligations_ok && factors_ok && product_ok && normalization_ok && original_reconstruction_ok,
    )
end

function _verify_sl3_local_obligation(obligation::SL3LocalObligation)::Bool
    return obligation.reassembly_data.local_product_ok && obligation.reassembly_data.embedded_product_ok
end

function _throw_staged_sln_to_sl3_failure(reason::AbstractString)
    throw(ArgumentError("staged SL_n to local SL_3 reduction failure: $(reason)"))
end
