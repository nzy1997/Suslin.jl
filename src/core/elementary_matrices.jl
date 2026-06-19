function elementary_matrix(n::Int, i::Int, j::Int, a, R)
    i == j && throw(ArgumentError("i and j must differ"))

    E = identity_matrix(R, n)
    coerced_a = _coerce_into_ring(R, a, "a")
    E[i, j] = coerced_a
    return E
end

function _require_square_matrix(M, label::AbstractString)
    nrows(M) == ncols(M) || throw(DimensionMismatch("$label must be square"))
    return nrows(M)
end

function _same_base_ring(left, right)::Bool
    if _is_laurent_polynomial_ring(left) || _is_laurent_polynomial_ring(right)
        return left === right
    end

    return left == right || left === right
end

function _embedding_indices(n::Int, block_size::Int, indices)
    block_size <= n || throw(DimensionMismatch("target size must be at least the block size"))
    length(indices) == block_size || throw(DimensionMismatch("number of indices must match the block size"))

    result = Int[]
    seen = Set{Int}()
    for index in indices
        index isa Integer || throw(ArgumentError("indices must be integers"))
        idx = Int(index)
        1 <= idx <= n || throw(ArgumentError("indices must be between 1 and the target size"))
        idx in seen && throw(ArgumentError("indices must be distinct"))
        push!(seen, idx)
        push!(result, idx)
    end
    return result
end

function block_embedding(block, n::Int, indices)
    block_size = _require_square_matrix(block, "block")
    target_indices = _embedding_indices(n, block_size, indices)

    R = base_ring(block)
    embedded = identity_matrix(R, n)
    for local_i in 1:block_size, local_j in 1:block_size
        embedded[target_indices[local_i], target_indices[local_j]] = block[local_i, local_j]
    end
    return embedded
end

function _factor_shape_and_ring(factor, label::AbstractString)
    factor_size = _require_square_matrix(factor, label)
    return factor_size, base_ring(factor)
end

function _require_matching_factor(factor, expected_size::Int, expected_ring, label::AbstractString)
    factor_size, factor_ring = _factor_shape_and_ring(factor, label)
    factor_size == expected_size || throw(ArgumentError("$label must have the same size as the first factor"))
    _same_base_ring(factor_ring, expected_ring) || throw(ArgumentError("$label must have the same base ring as the first factor"))
    return factor
end

function embed_factor_sequence(factors, n::Int, indices)
    collected = collect(factors)
    isempty(collected) && throw(ArgumentError("factor sequence must be nonempty"))

    factor_size, factor_ring = _factor_shape_and_ring(first(collected), "factor")
    _embedding_indices(n, factor_size, indices)
    for factor in Iterators.drop(collected, 1)
        _require_matching_factor(factor, factor_size, factor_ring, "factor")
    end

    return [block_embedding(factor, n, indices) for factor in collected]
end

function compose_factor_sequences(sequences...)
    isempty(sequences) && throw(ArgumentError("at least one factor sequence is required"))

    composed = nothing
    expected_size = nothing
    expected_ring = nothing

    for sequence in sequences
        for factor in sequence
            if composed === nothing
                expected_size, expected_ring = _factor_shape_and_ring(factor, "factor")
                composed = typeof(factor)[factor]
            else
                _require_matching_factor(factor, expected_size, expected_ring, "factor")
                push!(composed, factor)
            end
        end
    end

    composed === nothing && throw(ArgumentError("at least one factor is required"))
    return composed
end
