function realize_conjugate_elementary(B, i::Int, j::Int, a)
    nrows(B) == ncols(B) || throw(ArgumentError("B must be square"))

    n = nrows(B)
    n >= 3 || throw(ArgumentError("B must have size at least 3"))
    1 <= i <= n || throw(ArgumentError("i must be between 1 and n"))
    1 <= j <= n || throw(ArgumentError("j must be between 1 and n"))
    i == j && throw(ArgumentError("i and j must differ"))

    R = base_ring(B)
    coerced_a = _coerce_into_ring(R, a, "a")
    Binv = inv(B)

    v = [B[row, i] for row in 1:n]
    w = [coerced_a * Binv[j, col] for col in 1:n]
    g = [Binv[i, col] for col in 1:n]

    _dot(v, w, R) == zero(R) || throw(ArgumentError("conjugated elementary matrix did not produce an orthogonal I + v w presentation"))
    _dot(g, v, R) == one(R) || throw(ArgumentError("could not extract a unimodular witness row for the selected column"))

    factors = typeof(identity_matrix(R, n))[]
    for p in 1:(n - 1), q in (p + 1):n
        coeff = w[p] * g[q] - w[q] * g[p]
        coeff == zero(R) && continue
        append!(factors, realize_cohn_type(n, p, q, coeff, v, R))
    end
    return factors
end

function _dot(v::AbstractVector, w::AbstractVector, R)
    length(v) == length(w) || throw(ArgumentError("vectors must have the same length"))

    total = zero(R)
    for idx in eachindex(v, w)
        total += v[idx] * w[idx]
    end
    return total
end
