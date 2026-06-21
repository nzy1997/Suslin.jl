struct SL3LocalRealizationCertificate
    target
    branch::Symbol
    factors::Vector
    selected_variable
    witness
end

function realize_sl3_local(A, X; check_monic::Bool=true)
    return realize_sl3_local_certificate(A, X; check_monic).factors
end

function realize_sl3_local(p, q, r, s, X; check_monic::Bool=true)
    return realize_sl3_local_certificate(p, q, r, s, X; check_monic).factors
end

function realize_sl3_local_certificate(A, X; check_monic::Bool=true)
    form = _recognize_sl3_local_matrix(A, X; check_monic)
    return _realize_sl3_local_certificate_form(form)
end

function realize_sl3_local_certificate(p, q, r, s, X; check_monic::Bool=true)
    form = _recognize_sl3_local_parameters(p, q, r, s, X; check_monic)
    return _realize_sl3_local_certificate_form(form)
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
        return (; family = :open_s_one, R, p, q, r, s, X, target)
    end

    if p == one(R) && s == one(R) + q * r
        return (; family = :open_p_one, R, p, q, r, s, X, target)
    end

    s_inverse = _unit_inverse_or_nothing(s)
    if s_inverse !== nothing
        return (; family = :s_unit, R, p, q, r, s, X, target, pivot_inverse = s_inverse)
    end

    p_inverse = _unit_inverse_or_nothing(p)
    if p_inverse !== nothing
        return (; family = :p_unit, R, p, q, r, s, X, target, pivot_inverse = p_inverse)
    end

    _throw_staged_sl3_local_failure("supported families require one open unipotent slice or a unit diagonal pivot p or s")
end

function _sl3_local_form_factors(form)
    R = form.R
    if form.family == :open_s_one
        return [
            elementary_matrix(3, 1, 2, form.q, R),
            elementary_matrix(3, 2, 1, form.r, R),
        ]
    elseif form.family == :open_p_one
        return [
            elementary_matrix(3, 2, 1, form.r, R),
            elementary_matrix(3, 1, 2, form.q, R),
        ]
    elseif form.family == :s_unit
        s_inverse = form.pivot_inverse
        return vcat(
            [elementary_matrix(3, 1, 2, form.q * s_inverse, R)],
            _sl3_diagonal_unit_factors(s_inverse, R),
            [elementary_matrix(3, 2, 1, form.r * s_inverse, R)],
        )
    elseif form.family == :p_unit
        p_inverse = form.pivot_inverse
        return vcat(
            [elementary_matrix(3, 2, 1, form.r * p_inverse, R)],
            _sl3_diagonal_unit_factors(form.p, R),
            [elementary_matrix(3, 1, 2, form.q * p_inverse, R)],
        )
    end

    _throw_staged_sl3_local_failure("unrecognized local solver family $(form.family)")
end

function _sl3_local_form_witness(form)
    if form.family in (:open_s_one, :open_p_one)
        return (; q = form.q, r = form.r)
    elseif form.family == :s_unit
        return (; pivot = form.s, pivot_inverse = form.pivot_inverse)
    elseif form.family == :p_unit
        return (; pivot = form.p, pivot_inverse = form.pivot_inverse)
    end

    _throw_staged_sl3_local_failure("unrecognized local solver family $(form.family)")
end

function _realize_sl3_local_certificate_form(form)
    factors = _sl3_local_form_factors(form)
    _verify_sl3_local_factorization(form.target, factors)
    certificate = SL3LocalRealizationCertificate(
        form.target,
        form.family,
        factors,
        form.X,
        _sl3_local_form_witness(form),
    )
    verify_sl3_local_realization(certificate) ||
        error("internal local SL_3 realization certificate verification failed")
    return certificate
end

function _realize_sl3_local_form(form)
    return _realize_sl3_local_certificate_form(form).factors
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
    is_unit(value) || return nothing
    return inv(value)
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

function verify_sl3_local_realization(certificate)::Bool
    try
        return _sl3_local_realization_verification(certificate).overall_ok
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end

function _sl3_local_realization_verification(certificate)
    target = certificate.target
    size_ok = nrows(target) == 3 && ncols(target) == 3
    R = size_ok ? base_ring(target) : nothing
    variable_ok =
        size_ok &&
        parent(certificate.selected_variable) === R &&
        certificate.selected_variable in collect(gens(R))
    shape_ok = size_ok && _sl3_local_target_entries(target) !== nothing
    determinant_ok = size_ok && det(target) == one(R)
    expected_factors =
        size_ok && variable_ok && shape_ok && determinant_ok ?
        _sl3_local_certificate_expected_factors(certificate) :
        nothing
    witness_ok = expected_factors !== nothing
    factors_match_ok = witness_ok && _sl3_local_factor_sequences_equal(certificate.factors, expected_factors)
    factors_ok = size_ok && verify_factorization(target, certificate.factors)
    overall_ok = size_ok && variable_ok && shape_ok && determinant_ok &&
        witness_ok && factors_match_ok && factors_ok
    return (;
        overall_ok,
        size_ok,
        variable_ok,
        shape_ok,
        determinant_ok,
        witness_ok,
        factors_match_ok,
        factors_ok,
    )
end

function _sl3_local_target_entries(target)
    R = base_ring(target)
    if target[1, 3] != zero(R) || target[2, 3] != zero(R) ||
            target[3, 1] != zero(R) || target[3, 2] != zero(R) ||
            target[3, 3] != one(R)
        return nothing
    end

    return (; p = target[1, 1], q = target[1, 2], r = target[2, 1], s = target[2, 2])
end

function _sl3_local_factor_sequences_equal(left, right)
    length(left) == length(right) || return false
    for idx in eachindex(left, right)
        left[idx] == right[idx] || return false
    end
    return true
end

function _sl3_local_certificate_expected_factors(certificate)
    _sl3_local_branch_witness_ok(certificate) || return nothing

    entries = _sl3_local_target_entries(certificate.target)
    R = base_ring(certificate.target)
    p = entries.p
    q = entries.q
    r = entries.r
    s = entries.s
    witness = certificate.witness

    if certificate.branch == :open_s_one
        return [
            elementary_matrix(3, 1, 2, witness.q, R),
            elementary_matrix(3, 2, 1, witness.r, R),
        ]
    elseif certificate.branch == :open_p_one
        return [
            elementary_matrix(3, 2, 1, witness.r, R),
            elementary_matrix(3, 1, 2, witness.q, R),
        ]
    elseif certificate.branch == :s_unit
        return vcat(
            [elementary_matrix(3, 1, 2, q * witness.pivot_inverse, R)],
            _sl3_diagonal_unit_factors(witness.pivot_inverse, R),
            [elementary_matrix(3, 2, 1, r * witness.pivot_inverse, R)],
        )
    elseif certificate.branch == :p_unit
        return vcat(
            [elementary_matrix(3, 2, 1, r * witness.pivot_inverse, R)],
            _sl3_diagonal_unit_factors(witness.pivot, R),
            [elementary_matrix(3, 1, 2, q * witness.pivot_inverse, R)],
        )
    end

    return nothing
end

function _sl3_local_branch_witness_ok(certificate)
    entries = _sl3_local_target_entries(certificate.target)
    entries === nothing && return false

    R = base_ring(certificate.target)
    p = entries.p
    q = entries.q
    r = entries.r
    s = entries.s
    witness = certificate.witness

    if certificate.branch == :open_s_one
        return witness.q == q && witness.r == r && s == one(R) && p == one(R) + q * r
    elseif certificate.branch == :open_p_one
        return witness.q == q && witness.r == r && p == one(R) && s == one(R) + q * r
    elseif certificate.branch == :s_unit
        return witness.pivot == s && witness.pivot * witness.pivot_inverse == one(R)
    elseif certificate.branch == :p_unit
        return witness.pivot == p && witness.pivot * witness.pivot_inverse == one(R)
    end

    return false
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
