function reduce_unimodular_column(v::AbstractVector, R)
    Base.require_one_based_indexing(v)

    n = length(v)
    n >= 3 || throw(ArgumentError("v must have length at least 3"))

    column = [_coerce_into_ring(R, v[idx], "v[$idx]") for idx in 1:n]
    is_unimodular_column(column, R) || throw(ArgumentError("v must be a unimodular column"))
    factors = _reduce_supported_unimodular_column(column, R)
    factors !== nothing && return factors

    if ngens(R) >= 2
        normalization_factors = _reduce_after_monicity_normalization(column, R)
        normalization_factors !== nothing && return normalization_factors
    end

    throw(ArgumentError("reduce_unimodular_column currently supports only small unimodular columns reducible by direct unit witnesses or by a finite monicity-normalization search in the last variable"))
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
