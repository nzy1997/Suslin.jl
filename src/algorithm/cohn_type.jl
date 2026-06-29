struct CohnTypeRealizationCertificate
    n::Int
    i::Int
    j::Int
    a
    v::Vector
    ring
    auxiliary_index::Int
    target
    factors::Vector
    product
    verification
end

function _require_ordinary_polynomial_certificate_ring(R)
    _is_laurent_polynomial_ring(R) &&
        throw(ArgumentError("Cohn-type certificates require an ordinary polynomial ring"))
    try
        collect(gens(R))
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("Cohn-type certificates require an ordinary polynomial ring"))
    end
    return R
end

function _cohn_type_checked_data(n::Int, i::Int, j::Int, a, v::AbstractVector, R)
    n >= 3 || throw(ArgumentError("n must be at least 3"))
    1 <= i <= n || throw(ArgumentError("i must be between 1 and n"))
    1 <= j <= n || throw(ArgumentError("j must be between 1 and n"))
    i == j && throw(ArgumentError("i and j must differ"))
    Base.require_one_based_indexing(v)
    length(v) == n || throw(ArgumentError("v must contain exactly n entries"))

    t = findfirst(k -> k != i && k != j, 1:n)
    t === nothing && throw(ArgumentError("could not choose an auxiliary index distinct from i and j"))

    coerced_v = [_coerce_into_ring(R, v[idx], "v[$idx]") for idx in 1:n]
    coerced_a = _coerce_into_ring(R, a, "a")
    return (;
        n,
        i,
        j,
        a = coerced_a,
        v = coerced_v,
        ring = R,
        auxiliary_index = t,
    )
end

function _cohn_type_target_from_checked_data(data)
    target = identity_matrix(data.ring, data.n)
    vi = data.v[data.i]
    vj = data.v[data.j]
    for row in 1:data.n
        target[row, data.i] += data.a * data.v[row] * vj
        target[row, data.j] -= data.a * data.v[row] * vi
    end
    return target
end

function _cohn_type_factors_from_checked_data(data)
    n = data.n
    i = data.i
    j = data.j
    R = data.ring
    t = data.auxiliary_index
    vi = data.v[i]
    vj = data.v[j]
    a = data.a

    factors = [
        elementary_matrix(n, i, t, -vi, R),
        elementary_matrix(n, j, t, -vj, R),
        elementary_matrix(n, t, i, -a * vj, R),
        elementary_matrix(n, t, j, a * vi, R),
        elementary_matrix(n, i, t, vi, R),
        elementary_matrix(n, j, t, vj, R),
        elementary_matrix(n, t, i, a * vj, R),
        elementary_matrix(n, t, j, -a * vi, R),
    ]

    for l in 1:n
        (l == i || l == j) && continue
        coeff_li = a * data.v[l] * vj
        coeff_lj = -a * data.v[l] * vi
        coeff_li == zero(R) || push!(factors, elementary_matrix(n, l, i, coeff_li, R))
        coeff_lj == zero(R) || push!(factors, elementary_matrix(n, l, j, coeff_lj, R))
    end

    return factors
end

function _cohn_type_factor_product(factors, R, n::Int)
    product = identity_matrix(R, n)
    for factor in factors
        nrows(factor) == n || throw(ArgumentError("Cohn-type factor has wrong row count"))
        ncols(factor) == n || throw(ArgumentError("Cohn-type factor has wrong column count"))
        _same_base_ring(base_ring(factor), R) ||
            throw(ArgumentError("Cohn-type factor has wrong base ring"))
        product *= factor
    end
    return product
end

function realize_cohn_type(n::Int, i::Int, j::Int, a, v::AbstractVector, R)
    return _cohn_type_factors_from_checked_data(_cohn_type_checked_data(n, i, j, a, v, R))
end

function realize_cohn_type_certificate(n::Int, i::Int, j::Int, a, v::AbstractVector, R)
    _require_ordinary_polynomial_certificate_ring(R)
    data = _cohn_type_checked_data(n, i, j, a, v, R)
    target = _cohn_type_target_from_checked_data(data)
    factors = _cohn_type_factors_from_checked_data(data)
    product = _cohn_type_factor_product(factors, R, n)
    provisional = CohnTypeRealizationCertificate(
        data.n,
        data.i,
        data.j,
        data.a,
        data.v,
        data.ring,
        data.auxiliary_index,
        target,
        factors,
        product,
        nothing,
    )
    verification = _cohn_type_certificate_core_verification(provisional)
    verification.overall_core_ok ||
        error("internal Cohn-type certificate verification failed")
    return CohnTypeRealizationCertificate(
        provisional.n,
        provisional.i,
        provisional.j,
        provisional.a,
        provisional.v,
        provisional.ring,
        provisional.auxiliary_index,
        provisional.target,
        provisional.factors,
        provisional.product,
        verification,
    )
end

function _cohn_type_certificate_core_verification(cert)
    ring_ok = false
    checked_inputs_ok = false
    auxiliary_index_ok = false
    target_replay_ok = false
    factors_vector_ok = false
    product_replay_ok = false
    product_matches_stored_ok = false
    target_matches_product_ok = false
    factor_count = 0

    _require_ordinary_polynomial_certificate_ring(cert.ring)
    ring_ok = true

    data = _cohn_type_checked_data(cert.n, cert.i, cert.j, cert.a, cert.v, cert.ring)
    checked_inputs_ok = true
    auxiliary_index_ok = data.auxiliary_index == cert.auxiliary_index

    expected_target = _cohn_type_target_from_checked_data(data)
    target_replay_ok = cert.target == expected_target

    factors_vector_ok = cert.factors isa AbstractVector
    if factors_vector_ok
        factor_count = length(cert.factors)
        replayed_product = _cohn_type_factor_product(cert.factors, cert.ring, cert.n)
        product_replay_ok = true
        product_matches_stored_ok = replayed_product == cert.product
        target_matches_product_ok = target_replay_ok && replayed_product == cert.target
    end

    overall_core_ok =
        ring_ok &&
        checked_inputs_ok &&
        auxiliary_index_ok &&
        target_replay_ok &&
        factors_vector_ok &&
        product_replay_ok &&
        product_matches_stored_ok &&
        target_matches_product_ok

    return (;
        ring_ok,
        checked_inputs_ok,
        auxiliary_index_ok,
        target_replay_ok,
        factors_vector_ok,
        product_replay_ok,
        product_matches_stored_ok,
        target_matches_product_ok,
        factor_count,
        overall_core_ok,
    )
end

function _cohn_type_certificate_verification(cert)
    core = _cohn_type_certificate_core_verification(cert)
    stored_verification_ok = cert.verification == core
    return merge(core, (;
        stored_verification_ok,
        overall_ok = core.overall_core_ok && stored_verification_ok,
    ))
end

function verify_cohn_type_certificate(cert)::Bool
    try
        return _cohn_type_certificate_verification(cert).overall_ok
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end
