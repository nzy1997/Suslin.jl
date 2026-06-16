function elementary_matrix(n::Int, i::Int, j::Int, a, R)
    i == j && throw(ArgumentError("i and j must differ"))

    E = identity_matrix(R, n)
    coerced_a = _coerce_into_ring(R, a, "a")
    E[i, j] = coerced_a
    return E
end
