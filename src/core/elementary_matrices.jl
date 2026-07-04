function elementary_matrix(n::Int, i::Int, j::Int, a, R)
    i == j && throw(ArgumentError("i and j must differ"))

    E = identity_matrix(R, n)
    coerced_a = _coerce_into_ring(R, a, "a")
    E[i, j] = coerced_a
    return E
end

function _elementary_factor_term_exponents(value)
    try
        return collect(exponents(value))
    catch err
        err isa InterruptException && rethrow()
        err isa MethodError || err isa ErrorException || rethrow()
    end

    try
        return collect(AbstractAlgebra.exponent_vectors(value))
    catch err
        err isa InterruptException && rethrow()
        err isa MethodError || err isa ErrorException || rethrow()
        throw(ArgumentError("elementary factor analysis requires polynomial or Laurent polynomial entries"))
    end
end

function _elementary_factor_term_count(value)::Int
    iszero(value) && return 0
    return length(_elementary_factor_term_exponents(value))
end

function _elementary_factor_monomial_degree(raw_exponents)::Int
    return sum(abs, Int.(collect(raw_exponents)))
end

function _require_square_elementary_analysis_factor(factor)
    nrows(factor) == ncols(factor) || throw(ArgumentError("factor must be square"))
    return nrows(factor)
end

function max_elementary_factor_monomial_degree(factors)::Int
    max_degree = 0
    for factor in factors
        n = _require_square_elementary_analysis_factor(factor)
        for row in 1:n, col in 1:n
            row == col && continue
            value = factor[row, col]
            iszero(value) && continue
            for raw_exponents in _elementary_factor_term_exponents(value)
                max_degree = max(max_degree, _elementary_factor_monomial_degree(raw_exponents))
            end
        end
    end
    return max_degree
end

function total_elementary_factor_offdiagonal_monomials(factors)::Int
    total = 0
    for factor in factors
        n = _require_square_elementary_analysis_factor(factor)
        for row in 1:n, col in 1:n
            row == col && continue
            total += _elementary_factor_term_count(factor[row, col])
        end
    end
    return total
end

function _require_elementary_preconditioning_side(side)
    side isa Symbol || throw(ArgumentError("side must be :left or :right"))
    (side == :left || side == :right) || throw(ArgumentError("side must be :left or :right"))
    return side
end

function _require_preconditioning_index(index, limit::Int, label::AbstractString)
    index isa Integer || throw(ArgumentError("$label must be an integer"))
    idx = Int(index)
    1 <= idx <= limit || throw(ArgumentError("$label must be between 1 and $limit"))
    return idx
end

function _preconditioning_factor_size(A, side::Symbol)
    return side == :left ? nrows(A) : ncols(A)
end

function _preconditioning_factor_indices(side::Symbol, target::Int, source::Int)
    return side == :left ? (target, source) : (source, target)
end

function elementary_preconditioning_step(A, side, target, source, coefficient)
    checked_side = _require_elementary_preconditioning_side(side)
    R = base_ring(A)
    factor_size = _preconditioning_factor_size(A, checked_side)
    target_idx = _require_preconditioning_index(target, factor_size, "target")
    source_idx = _require_preconditioning_index(source, factor_size, "source")
    target_idx == source_idx && throw(ArgumentError("target and source must differ"))

    coerced_coefficient = _coerce_into_ring(R, coefficient, "coefficient")
    factor_row, factor_col = _preconditioning_factor_indices(checked_side, target_idx, source_idx)
    factor = elementary_matrix(factor_size, factor_row, factor_col, coerced_coefficient, R)
    transformed_matrix = checked_side == :left ? factor * A : A * factor

    return (;
        side = checked_side,
        target = target_idx,
        source = source_idx,
        coefficient = coerced_coefficient,
        factor,
        transformed_matrix,
    )
end

function _require_preconditioning_step_property(step, property::Symbol)
    property in propertynames(step) || throw(ArgumentError("preconditioning step must include $(property)"))
    return getproperty(step, property)
end

function _require_preconditioning_factor(current, side::Symbol, factor)
    expected_size = _preconditioning_factor_size(current, side)
    nrows(factor) == expected_size || throw(ArgumentError("preconditioning factor has wrong row count for $(side) side"))
    ncols(factor) == expected_size || throw(ArgumentError("preconditioning factor has wrong column count for $(side) side"))
    _same_base_ring(base_ring(factor), base_ring(current)) ||
        throw(ArgumentError("preconditioning factor must have the same base ring as the current matrix"))
    return factor
end

function _apply_preconditioning_factor(current, side, factor)
    checked_side = _require_elementary_preconditioning_side(side)
    checked_factor = _require_preconditioning_factor(current, checked_side, factor)
    return checked_side == :left ? checked_factor * current : current * checked_factor
end

function replay_elementary_preconditioning(A, steps)
    current = A
    for step in steps
        side = _require_preconditioning_step_property(step, :side)
        factor = _require_preconditioning_step_property(step, :factor)
        current = _apply_preconditioning_factor(current, side, factor)
    end
    return current
end

function verify_elementary_preconditioning(A, steps, expected)::Bool
    try
        replayed = replay_elementary_preconditioning(A, steps)
        nrows(replayed) == nrows(expected) || return false
        ncols(replayed) == ncols(expected) || return false
        _same_base_ring(base_ring(replayed), base_ring(expected)) || return false
        return replayed == expected
    catch err
        err isa InterruptException && rethrow()
        return false
    end
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

function _canonical_elementary_factor_record(factor)
    n = _require_square_matrix(factor, "elementary factor")
    R = base_ring(factor)
    row = 0
    col = 0
    coefficient = zero(R)

    for i in 1:n, j in 1:n
        entry = factor[i, j]
        if i == j
            entry == one(R) || throw(ArgumentError("elementary factor diagonal must be one"))
            continue
        end

        entry == zero(R) && continue
        row == 0 || throw(ArgumentError("elementary factor must have at most one nonzero off-diagonal entry"))
        row = i
        col = j
        coefficient = entry
    end

    row == 0 && return (; kind = :identity, n, ring = R)
    return (; kind = :elementary, n, ring = R, row, col, coefficient)
end

function _elementary_factor_record_matrix(record)
    record.kind == :identity && return identity_matrix(record.ring, record.n)
    record.kind == :elementary &&
        return elementary_matrix(record.n, record.row, record.col, record.coefficient, record.ring)
    throw(ArgumentError("unknown elementary factor record kind"))
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
