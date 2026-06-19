function realize_sl3_local(A, X; check_monic::Bool=true)
    form = _recognize_sl3_local_matrix(A, X; check_monic)
    return _realize_sl3_local_form(form)
end

function realize_sl3_local(p, q, r, s, X; check_monic::Bool=true)
    form = _recognize_sl3_local_parameters(p, q, r, s, X; check_monic)
    return _realize_sl3_local_form(form)
end

function _recognize_sl3_local_matrix(A, X; check_monic::Bool=true)
    nrows(A) == 3 || throw(ArgumentError("local SL_3 special-form recognition failed: A must be 3x3"))
    ncols(A) == 3 || throw(ArgumentError("local SL_3 special-form recognition failed: A must be 3x3"))

    R = base_ring(A)
    parent(X) === R || throw(ArgumentError("A and X must lie in the same polynomial or Laurent polynomial ring"))
    if A[1, 3] != zero(R) || A[2, 3] != zero(R) ||
            A[3, 1] != zero(R) || A[3, 2] != zero(R) ||
            A[3, 3] != one(R)
        throw(ArgumentError("local SL_3 special-form recognition failed: A must be an embedded 2x2 block with trailing identity"))
    end

    return _recognize_sl3_local_parameters(A[1, 1], A[1, 2], A[2, 1], A[2, 2], X; check_monic)
end

function _recognize_sl3_local_parameters(p, q, r, s, X; check_monic::Bool=true)
    R = parent(p)
    parent(q) === R || throw(ArgumentError("p, q, r, s, and X must lie in the same polynomial or Laurent polynomial ring"))
    parent(r) === R || throw(ArgumentError("p, q, r, s, and X must lie in the same polynomial or Laurent polynomial ring"))
    parent(s) === R || throw(ArgumentError("p, q, r, s, and X must lie in the same polynomial or Laurent polynomial ring"))
    parent(X) === R || throw(ArgumentError("p, q, r, s, and X must lie in the same polynomial or Laurent polynomial ring"))

    ring_gens = collect(gens(R))
    var_idx = findfirst(isequal(X), ring_gens)
    var_idx === nothing && throw(ArgumentError("X must be one of the polynomial or Laurent polynomial ring generators"))

    if check_monic
        if _is_laurent_polynomial_ring(R)
            throw(ArgumentError("p monicity check is only supported for ordinary polynomial local SL_3 inputs; pass check_monic=false only when the caller has discharged the local monicity assumption"))
        end
        _is_monic_in_variable(p, var_idx, R) || throw(ArgumentError("p must be monic in X"))
    end

    target = matrix(R, [
        p q zero(R);
        r s zero(R);
        zero(R) zero(R) one(R)
    ])
    det(target) == one(R) || throw(ArgumentError("constructed matrix must have determinant 1"))

    if s == one(R) && p == one(R) + q * r
        return (; family = :open_s_one, R, p, q, r, s, target)
    end

    if p == one(R) && s == one(R) + q * r
        return (; family = :open_p_one, R, p, q, r, s, target)
    end

    s_inverse = _unit_inverse_or_nothing(s)
    if s_inverse !== nothing
        return (; family = :s_unit, R, p, q, r, s, target, pivot_inverse = s_inverse)
    end

    p_inverse = _unit_inverse_or_nothing(p)
    if p_inverse !== nothing
        return (; family = :p_unit, R, p, q, r, s, target, pivot_inverse = p_inverse)
    end

    _throw_staged_sl3_local_failure("supported families require one open unipotent slice or a unit diagonal pivot p or s")
end

function _realize_sl3_local_form(form)
    R = form.R

    factors = if form.family == :open_s_one
        [
            elementary_matrix(3, 1, 2, form.q, R),
            elementary_matrix(3, 2, 1, form.r, R),
        ]
    elseif form.family == :open_p_one
        [
            elementary_matrix(3, 2, 1, form.r, R),
            elementary_matrix(3, 1, 2, form.q, R),
        ]
    elseif form.family == :s_unit
        s_inverse = form.pivot_inverse
        vcat(
            [elementary_matrix(3, 1, 2, form.q * s_inverse, R)],
            _sl3_diagonal_unit_factors(s_inverse, R),
            [elementary_matrix(3, 2, 1, form.r * s_inverse, R)],
        )
    elseif form.family == :p_unit
        p_inverse = form.pivot_inverse
        vcat(
            [elementary_matrix(3, 2, 1, form.r * p_inverse, R)],
            _sl3_diagonal_unit_factors(form.p, R),
            [elementary_matrix(3, 1, 2, form.q * p_inverse, R)],
        )
    else
        _throw_staged_sl3_local_failure("unrecognized local solver family $(form.family)")
    end

    _verify_sl3_local_factorization(form.target, factors)
    return factors
end

function _sl3_diagonal_unit_factors(u, R)
    u_inverse = inv(u)
    return [
        elementary_matrix(3, 1, 2, u, R),
        elementary_matrix(3, 2, 1, -u_inverse, R),
        elementary_matrix(3, 1, 2, u, R),
        elementary_matrix(3, 1, 2, -one(R), R),
        elementary_matrix(3, 2, 1, one(R), R),
        elementary_matrix(3, 1, 2, -one(R), R),
    ]
end

function _unit_inverse_or_nothing(value)
    try
        is_unit(value) || return nothing
        return inv(value)
    catch err
        if err isa ArgumentError || err isa MethodError || err isa ErrorException
            return nothing
        end
        rethrow()
    end
end

function _verify_sl3_local_factorization(target, factors)
    R = base_ring(target)
    product = identity_matrix(R, 3)
    for factor in factors
        nrows(factor) == 3 || error("local SL_3 exact verification failed: factor has wrong row count")
        ncols(factor) == 3 || error("local SL_3 exact verification failed: factor has wrong column count")
        _same_base_ring(base_ring(factor), R) || error("local SL_3 exact verification failed: factor base ring mismatch")
        product *= factor
    end

    product == target || error("local SL_3 exact verification failed: factor product does not equal target")
    return nothing
end

function _throw_staged_sl3_local_failure(reason::AbstractString)
    throw(ArgumentError("staged local SL_3 solver failure: $(reason)"))
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
