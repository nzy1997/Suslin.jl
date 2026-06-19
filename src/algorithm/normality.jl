function realize_conjugate_elementary(B, i::Int, j::Int, a)
    nrows(B) == ncols(B) || throw(ArgumentError("B must be square"))

    n = nrows(B)
    n >= 3 || throw(ArgumentError("B must have size at least 3"))
    1 <= i <= n || throw(ArgumentError("i must be between 1 and n"))
    1 <= j <= n || throw(ArgumentError("j must be between 1 and n"))
    i == j && throw(ArgumentError("i and j must differ"))

    R = base_ring(B)
    coerced_a = _coerce_into_ring(R, a, "a")
    Binv = _inverse_matrix_over_base_ring(B)

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

function _inverse_matrix_over_base_ring(M)
    R = base_ring(M)
    _is_laurent_polynomial_ring(R) && return _laurent_inverse_matrix(M)

    return inv(M)
end

function _laurent_inverse_matrix(M)
    n = _require_square_matrix(M, "B")
    R = _require_laurent_polynomial_ring(base_ring(M); label="B base ring")
    determinant = det(M)
    determinant == zero(R) && throw(ArgumentError("B must be invertible over its Laurent polynomial ring"))

    determinant_inverse = try
        inv(determinant)
    catch err
        if err isa ArgumentError || err isa ErrorException || err isa MethodError
            throw(ArgumentError("B must have a Laurent-unit determinant"))
        end
        rethrow()
    end

    inverse = _adjugate_matrix(M)
    for row in 1:n, col in 1:n
        inverse[row, col] *= determinant_inverse
    end

    identity = identity_matrix(R, n)
    inverse * M == identity || throw(ArgumentError("failed to invert B over its Laurent polynomial ring"))
    M * inverse == identity || throw(ArgumentError("failed to invert B over its Laurent polynomial ring"))
    return inverse
end

function _adjugate_matrix(M)
    n = _require_square_matrix(M, "matrix")
    R = base_ring(M)

    if n == 1
        adjugate = zero_matrix(R, 1, 1)
        adjugate[1, 1] = one(R)
        return adjugate
    end

    adjugate = zero_matrix(R, n, n)
    for row in 1:n, col in 1:n
        sign = isodd(row + col) ? -one(R) : one(R)
        adjugate[col, row] = sign * det(_matrix_without_row_col(M, row, col))
    end
    return adjugate
end

function _matrix_without_row_col(M, skipped_row::Int, skipped_col::Int)
    n = _require_square_matrix(M, "matrix")
    R = base_ring(M)
    minor = zero_matrix(R, n - 1, n - 1)
    minor_row = 1

    for row in 1:n
        row == skipped_row && continue
        minor_col = 1
        for col in 1:n
            col == skipped_col && continue
            minor[minor_row, minor_col] = M[row, col]
            minor_col += 1
        end
        minor_row += 1
    end

    return minor
end

function _dot(v::AbstractVector, w::AbstractVector, R)
    length(v) == length(w) || throw(ArgumentError("vectors must have the same length"))

    total = zero(R)
    for idx in eachindex(v, w)
        total += v[idx] * w[idx]
    end
    return total
end
