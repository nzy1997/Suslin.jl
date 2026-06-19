function reduce_unimodular_column(v::AbstractVector, R)
    column = _validated_unimodular_column(v, R)

    factors = _is_laurent_polynomial_ring(R) ?
        _reduce_laurent_unimodular_column(column, R) :
        _reduce_polynomial_unimodular_column_exact(column, R)

    factors !== nothing && return _checked_reduction_factors(factors, column, R, "public reducer")
    _throw_unsupported_unimodular_column_reduction(column, R)
end

function _validated_unimodular_column(v::AbstractVector, R)
    Base.require_one_based_indexing(v)
    n = length(v)
    n >= 3 || throw(ArgumentError("v must have length at least 3"))

    column = [_coerce_into_ring(R, v[idx], "v[$idx]") for idx in 1:n]
    is_unimodular_column(column, R) || throw(ArgumentError("v must be a unimodular column"))
    return column
end

function _factor_sequence_product(factors, R, n::Int)
    product = identity_matrix(R, n)
    for factor in factors
        nrows(factor) == n || throw(ArgumentError("reduction factor has wrong row count"))
        ncols(factor) == n || throw(ArgumentError("reduction factor has wrong column count"))
        _same_base_ring(base_ring(factor), R) || throw(ArgumentError("reduction factor has wrong base ring"))
        product *= factor
    end
    return product
end

function _apply_reduction_factors(factors, column::AbstractVector, R)
    n = length(column)
    return _factor_sequence_product(factors, R, n) * matrix(R, n, 1, collect(column))
end

function _target_reduced_column(R, n::Int)
    target = zero_matrix(R, n, 1)
    target[n, 1] = one(R)
    return target
end

function _checked_reduction_factors(factors, column::AbstractVector, R, stage::AbstractString)
    n = length(column)
    _apply_reduction_factors(factors, column, R) == _target_reduced_column(R, n) ||
        throw(ErrorException("internal error: $(stage) produced factors that do not reduce the column to e_n"))
    return factors
end

function _reduce_unit_witness_column(column::AbstractVector, pivot_idx::Int, R)
    n = length(column)
    u = column[pivot_idx]
    uinv = inv(u)
    elimination_factors = typeof(identity_matrix(R, n))[]

    for k in 1:n
        k == pivot_idx && continue
        coeff = -column[k] * uinv
        coeff == zero(R) && continue
        push!(elimination_factors, elementary_matrix(n, k, pivot_idx, coeff, R))
    end

    move_factors = typeof(identity_matrix(R, n))[]
    if pivot_idx != n
        push!(move_factors, elementary_matrix(n, pivot_idx, n, -one(R), R))
        push!(move_factors, elementary_matrix(n, n, pivot_idx, one(R), R))
    end

    normalization_factors = _unit_normalization_factors(n, u, uinv, R)
    return vcat(normalization_factors, move_factors, elimination_factors)
end

function _unimodular_witness(column::AbstractVector, R)
    witness_matrix = Oscar.coordinates(one(R), ideal(R, column))
    return [witness_matrix[1, idx] for idx in 1:ncols(witness_matrix)]
end

function _reduce_via_witness_unit(column::AbstractVector, witness::AbstractVector, pivot_idx::Int, R)
    n = length(column)
    u = witness[pivot_idx]
    uinv = inv(u)
    unit_creation_factors = typeof(identity_matrix(R, n))[]

    for k in 1:n
        k == pivot_idx && continue
        coeff = uinv * witness[k]
        coeff == zero(R) && continue
        push!(unit_creation_factors, elementary_matrix(n, pivot_idx, k, coeff, R))
    end

    created_column = copy(column)
    created_column[pivot_idx] = uinv
    reduction_factors = _reduce_unit_witness_column(created_column, pivot_idx, R)
    return vcat(reduction_factors, unit_creation_factors)
end

function _reduce_supported_unimodular_column(column::AbstractVector, R)
    unit_idx = findfirst(is_unit, column)
    unit_idx !== nothing && return _reduce_unit_witness_column(column, unit_idx, R)

    witness = _unimodular_witness(column, R)
    witness_unit_idx = findfirst(is_unit, witness)
    witness_unit_idx !== nothing && return _reduce_via_witness_unit(column, witness, witness_unit_idx, R)

    return nothing
end

function _reduce_polynomial_unimodular_column_exact(column::AbstractVector, R)
    factors = _reduce_exact_small_column(column, R)
    factors !== nothing && return factors

    block_factors = _reduce_via_supported_three_block(column, R)
    block_factors !== nothing && return block_factors

    return nothing
end

function _has_at_least_two_generators(R)::Bool
    try
        return ngens(R) >= 2
    catch err
        err isa MethodError || err isa ErrorException || rethrow()
        return false
    end
end

function _reduce_exact_small_column(column::AbstractVector, R)
    factors = _reduce_supported_unimodular_column(column, R)
    factors !== nothing && return _checked_reduction_factors(factors, column, R, "unit or witness reduction")

    _has_at_least_two_generators(R) || return nothing
    normalization_factors = _reduce_after_monicity_normalization(column, R)
    normalization_factors !== nothing &&
        return _checked_reduction_factors(normalization_factors, column, R, "monicity normalization reduction")

    return nothing
end

function _reduce_via_supported_three_block(column::AbstractVector, R)
    n = length(column)
    n > 3 || return nothing

    for i in 1:(n - 2), j in (i + 1):(n - 1), k in (j + 1):n
        indices = (i, j, k)
        subcolumn = [column[idx] for idx in indices]
        is_unimodular_column(subcolumn, R) || continue

        subfactors = _reduce_exact_small_column(subcolumn, R)
        subfactors === nothing && continue

        return _embedded_three_block_reduction(column, R, indices, subfactors)
    end

    return nothing
end

function _embedded_three_block_reduction(column::AbstractVector, R, indices, subfactors)
    n = length(column)
    pivot_idx = indices[end]
    embedded_factors = [block_embedding(factor, n, indices) for factor in subfactors]
    after_block = _apply_reduction_factors(embedded_factors, column, R)

    elimination_factors = typeof(identity_matrix(R, n))[]
    for row in 1:n
        row == pivot_idx && continue
        coeff = -after_block[row, 1]
        coeff == zero(R) && continue
        push!(elimination_factors, elementary_matrix(n, row, pivot_idx, coeff, R))
    end

    return _checked_reduction_factors(
        vcat(elimination_factors, embedded_factors),
        column,
        R,
        "embedded 3-entry reduction",
    )
end

function _reduce_after_monicity_normalization(column::AbstractVector, R)
    ring_gens = collect(gens(R))
    last_var = ring_gens[end]

    for var_idx in 1:(length(ring_gens) - 1), shift_power in 1:3, shift_sign in (one(R), -one(R))
        forward_values = copy(ring_gens)
        forward_values[var_idx] = ring_gens[var_idx] + shift_sign * last_var^shift_power
        transformed = [_coerce_into_ring(R, evaluate(entry, forward_values), "substituted column entry") for entry in column]
        any(entry -> _is_monic_in_last_variable(entry, R), transformed) || continue

        transformed_factors = _reduce_supported_unimodular_column(transformed, R)
        transformed_factors === nothing && continue

        inverse_values = copy(ring_gens)
        inverse_values[var_idx] = ring_gens[var_idx] - shift_sign * last_var^shift_power
        return [_substitute_matrix_entries(factor, inverse_values, R) for factor in transformed_factors]
    end

    return nothing
end

function _reduce_laurent_unimodular_column(column::AbstractVector, R)
    unit_idx = findfirst(is_unit, column)
    unit_idx !== nothing && return _reduce_unit_witness_column(column, unit_idx, R)

    normalization = normalize_laurent_object(column)
    poly_column = normalization.normalized_object
    P = normalization.metadata.polynomial_ring
    is_unimodular_column(poly_column, P) || return nothing

    poly_factors = _reduce_polynomial_unimodular_column_exact(poly_column, P)
    poly_factors === nothing && return nothing

    lifted_factors = [_lift_polynomial_reduction_factor(factor, R) for factor in poly_factors]
    shift = only(normalization.metadata.shift_monomials)
    inverse_shift = only(normalization.metadata.inverse_shift_monomials)
    normalization_factors = _unit_normalization_factors(length(column), inverse_shift, shift, R)

    return _checked_reduction_factors(
        vcat(normalization_factors, lifted_factors),
        column,
        R,
        "Laurent normalization reduction",
    )
end

function _lift_polynomial_reduction_factor(factor, R)
    entries = [
        _coerce_into_ring(R, factor[row, col], "lifted reduction factor entry")
        for col in 1:ncols(factor), row in 1:nrows(factor)
    ]
    return matrix(R, nrows(factor), ncols(factor), vec(entries))
end

function _is_monic_in_last_variable(p, R)
    iszero(p) && return false
    return _leading_coefficient_in_last_variable(p, R) == one(R)
end

function _leading_coefficient_in_last_variable(p, R)
    last_idx = ngens(R)
    target_degree = degree(p, last_idx)
    target_degree < 0 && return zero(R)

    ring_gens = collect(gens(R))
    total = zero(R)
    for (coeff, exponents) in zip(AbstractAlgebra.coefficients(p), AbstractAlgebra.exponent_vectors(p))
        exponents[last_idx] == target_degree || continue
        term = R(coeff)
        for idx in 1:(last_idx - 1)
            exponent = exponents[idx]
            exponent == 0 && continue
            term *= ring_gens[idx]^exponent
        end
        total += term
    end

    return total
end

function _substitute_matrix_entries(M, values::AbstractVector, R)
    entries = [_coerce_into_ring(R, evaluate(M[row, col], values), "substituted matrix entry") for col in 1:ncols(M), row in 1:nrows(M)]
    return matrix(R, nrows(M), ncols(M), vec(entries))
end

function _unit_normalization_factors(n::Int, u, uinv, R)
    factors = typeof(identity_matrix(R, n))[]
    u == one(R) && return factors

    i = n - 1
    j = n

    for (row, col, coeff) in (
        (i, j, u - one(R)),
        (j, i, one(R)),
        (i, j, uinv - one(R)),
        (j, i, -u),
    )
        push!(factors, elementary_matrix(n, row, col, coeff, R))
    end

    return factors
end

function _throw_unsupported_unimodular_column_reduction(column::AbstractVector, R)
    n = length(column)
    profile = _is_laurent_polynomial_ring(R) ? "Laurent-normalized" : "ordinary polynomial"
    throw(ArgumentError("unsupported exact unimodular column reduction for $(profile) column of length $(n): no supported unit, witness-unit, monicity-normalized, or 3-entry block reduction stage applies"))
end
