function realize_cohn_type(n::Int, i::Int, j::Int, a, v::AbstractVector, R)
    n >= 3 || throw(ArgumentError("n must be at least 3"))
    1 <= i <= n || throw(ArgumentError("i must be between 1 and n"))
    1 <= j <= n || throw(ArgumentError("j must be between 1 and n"))
    i == j && throw(ArgumentError("i and j must differ"))
    Base.require_one_based_indexing(v)
    length(v) == 2 || throw(ArgumentError("v must contain exactly two entries"))

    t = findfirst(k -> k != i && k != j, 1:n)
    t === nothing && throw(ArgumentError("could not choose an auxiliary index distinct from i and j"))

    v1 = _coerce_into_ring(R, v[1], "v[1]")
    v2 = _coerce_into_ring(R, v[2], "v[2]")
    coerced_a = _coerce_into_ring(R, a, "a")

    return [
        elementary_matrix(n, i, t, -v1, R),
        elementary_matrix(n, j, t, -v2, R),
        elementary_matrix(n, t, i, -coerced_a * v2, R),
        elementary_matrix(n, t, j, coerced_a * v1, R),
        elementary_matrix(n, i, t, v1, R),
        elementary_matrix(n, j, t, v2, R),
        elementary_matrix(n, t, i, coerced_a * v2, R),
        elementary_matrix(n, t, j, -coerced_a * v1, R),
    ]
end
