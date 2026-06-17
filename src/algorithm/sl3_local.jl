function realize_sl3_local(p, q, r, s, X; check_monic::Bool=true)
    R = parent(p)
    parent(q) === R || throw(ArgumentError("p, q, r, s, and X must lie in the same polynomial ring"))
    parent(r) === R || throw(ArgumentError("p, q, r, s, and X must lie in the same polynomial ring"))
    parent(s) === R || throw(ArgumentError("p, q, r, s, and X must lie in the same polynomial ring"))
    parent(X) === R || throw(ArgumentError("p, q, r, s, and X must lie in the same polynomial ring"))

    ring_gens = collect(gens(R))
    var_idx = findfirst(isequal(X), ring_gens)
    var_idx === nothing && throw(ArgumentError("X must be one of the polynomial ring generators"))

    if check_monic && !_is_monic_in_variable(p, var_idx, R)
        throw(ArgumentError("p must be monic in X"))
    end

    target = matrix(R, [
        p q zero(R);
        r s zero(R);
        zero(R) zero(R) one(R)
    ])
    det(target) == one(R) || throw(ArgumentError("constructed matrix must have determinant 1"))

    if s == one(R) && p == one(R) + q * r
        return [
            elementary_matrix(3, 1, 2, q, R),
            elementary_matrix(3, 2, 1, r, R),
        ]
    end

    if p == one(R) && s == one(R) + q * r
        return [
            elementary_matrix(3, 2, 1, r, R),
            elementary_matrix(3, 1, 2, q, R),
        ]
    end

    throw(ArgumentError("only the open Gauss-cell slices p = 1 + q*r, s = 1 or p = 1, s = 1 + q*r are supported"))
end

function _is_monic_in_variable(p, var_idx::Int, R)
    iszero(p) && return false

    target_degree = degree(p, var_idx)
    target_degree < 0 && return false

    ring_gens = collect(gens(R))
    total = zero(R)
    for (coeff, exponents) in zip(AbstractAlgebra.coefficients(p), AbstractAlgebra.exponent_vectors(p))
        exponents[var_idx] == target_degree || continue
        term = R(coeff)
        for idx in 1:length(ring_gens)
            idx == var_idx && continue
            exponent = exponents[idx]
            exponent == 0 && continue
            term *= ring_gens[idx]^exponent
        end
        total += term
    end

    return total == one(R)
end
