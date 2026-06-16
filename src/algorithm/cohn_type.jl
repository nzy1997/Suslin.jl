function realize_cohn_type(n::Int, i::Int, j::Int, a, v::AbstractVector, R)
    n >= 3 || throw(ArgumentError("n must be at least 3"))
    1 <= i <= n || throw(ArgumentError("i must be between 1 and n"))
    1 <= j <= n || throw(ArgumentError("j must be between 1 and n"))
    i == j && throw(ArgumentError("i and j must differ"))
    Base.require_one_based_indexing(v)
    length(v) == n || throw(ArgumentError("v must contain exactly n entries"))

    t = findfirst(k -> k != i && k != j, 1:n)
    t === nothing && throw(ArgumentError("could not choose an auxiliary index distinct from i and j"))

    coerced_v = [_coerce_into_ring(R, v[idx], "v[$idx]") for idx in 1:n]
    vi = coerced_v[i]
    vj = coerced_v[j]
    coerced_a = _coerce_into_ring(R, a, "a")

    factors = [
        elementary_matrix(n, i, t, -vi, R),
        elementary_matrix(n, j, t, -vj, R),
        elementary_matrix(n, t, i, -coerced_a * vj, R),
        elementary_matrix(n, t, j, coerced_a * vi, R),
        elementary_matrix(n, i, t, vi, R),
        elementary_matrix(n, j, t, vj, R),
        elementary_matrix(n, t, i, coerced_a * vj, R),
        elementary_matrix(n, t, j, -coerced_a * vi, R),
    ]

    for l in 1:n
        (l == i || l == j) && continue
        coeff_li = coerced_a * coerced_v[l] * vj
        coeff_lj = -coerced_a * coerced_v[l] * vi
        coeff_li == zero(R) || push!(factors, elementary_matrix(n, l, i, coeff_li, R))
        coeff_lj == zero(R) || push!(factors, elementary_matrix(n, l, j, coeff_lj, R))
    end

    return factors
end
