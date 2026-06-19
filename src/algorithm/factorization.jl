function elementary_factorization(A)
    laurent_input = _is_laurent_polynomial_ring(base_ring(A))
    if _is_laurent_polynomial_ring(base_ring(A))
        normalization = normalize_laurent_gl_matrix(A)
        A = normalization.normalized_matrix
    end

    nrows(A) == ncols(A) || throw(ArgumentError("A must be square"))
    nrows(A) == 3 || throw(ArgumentError("elementary_factorization currently supports only 3x3 matrices"))

    R = base_ring(A)
    ring_gens = collect(gens(R))
    if laurent_input
        throw(ArgumentError("Laurent matrices are normalized at the GL_n boundary, but the SL_n factorization core currently supports only polynomial rings"))
    end
    length(ring_gens) == 1 || throw(ArgumentError("elementary_factorization currently supports only univariate polynomial rings"))
    det(A) == one(R) || throw(ArgumentError("elementary_factorization currently supports only matrices in SL_3"))

    if A[1, 3] != zero(R) || A[2, 3] != zero(R) || A[3, 1] != zero(R) || A[3, 2] != zero(R) || A[3, 3] != one(R)
        throw(ArgumentError("elementary_factorization currently supports only the local SL_3 slice with zero third row and column except the (3,3) entry"))
    end

    return realize_sl3_local(A[1, 1], A[1, 2], A[2, 1], A[2, 2], ring_gens[1])
end

function verify_factorization(A, factors)::Bool
    nrows(A) == ncols(A) || throw(ArgumentError("A must be square"))

    R = base_ring(A)
    product = identity_matrix(R, nrows(A))

    for factor in factors
        nrows(factor) == nrows(A) || throw(ArgumentError("all factors must have the same size as A"))
        ncols(factor) == ncols(A) || throw(ArgumentError("all factors must have the same size as A"))
        _same_base_ring(base_ring(factor), R) || throw(ArgumentError("all factors must lie in the same base ring as A"))
        product *= factor
    end

    return product == A
end
