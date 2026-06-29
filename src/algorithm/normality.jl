struct RankOneNormalityCertificate
    n::Int
    v::Vector
    w::Vector
    g::Vector
    ring
    orthogonality
    bezout
    cohn_coefficients::Vector
    child_certificates::Vector{CohnTypeRealizationCertificate}
    factors::Vector
    target
    product
    verification
end

function _rank_one_checked_data(v::AbstractVector, w::AbstractVector, g::AbstractVector, R)
    Base.require_one_based_indexing(v)
    Base.require_one_based_indexing(w)
    Base.require_one_based_indexing(g)
    n = length(v)
    n >= 3 || throw(ArgumentError("v, w, and g must have length at least 3"))
    length(w) == n || throw(ArgumentError("w must have the same length as v"))
    length(g) == n || throw(ArgumentError("g must have the same length as v"))
    _require_ordinary_polynomial_certificate_ring(R)
    coerced_v = [_coerce_into_ring(R, v[idx], "v[$idx]") for idx in 1:n]
    coerced_w = [_coerce_into_ring(R, w[idx], "w[$idx]") for idx in 1:n]
    coerced_g = [_coerce_into_ring(R, g[idx], "g[$idx]") for idx in 1:n]
    orthogonality = _dot(coerced_w, coerced_v, R)
    bezout = _dot(coerced_g, coerced_v, R)
    return (;
        n,
        v = coerced_v,
        w = coerced_w,
        g = coerced_g,
        ring = R,
        orthogonality,
        bezout,
    )
end

function _rank_one_target_from_checked_data(data)
    target = identity_matrix(data.ring, data.n)
    for row in 1:data.n, col in 1:data.n
        target[row, col] += data.v[row] * data.w[col]
    end
    return target
end

function _rank_one_cohn_coefficients_from_checked_data(data)
    return [
        (; i, j, a = data.w[i] * data.g[j] - data.w[j] * data.g[i])
        for i in 1:(data.n - 1) for j in (i + 1):data.n
    ]
end

function _rank_one_child_certificates_from_checked_data(data, coefficients)
    children = CohnTypeRealizationCertificate[]
    for entry in coefficients
        entry.a == zero(data.ring) && continue
        push!(children, realize_cohn_type_certificate(data.n, entry.i, entry.j, entry.a, data.v, data.ring))
    end
    return children
end

function _rank_one_factors_from_children(children)
    factors = Any[]
    for child in children
        append!(factors, child.factors)
    end
    return factors
end

function realize_rank_one_normality_certificate(v::AbstractVector, w::AbstractVector, g::AbstractVector, R)
    data = _rank_one_checked_data(v, w, g, R)
    data.orthogonality == zero(R) || throw(ArgumentError("rank-one inputs must satisfy w*v == 0"))
    data.bezout == one(R) || throw(ArgumentError("rank-one inputs must satisfy g*v == 1"))
    coefficients = _rank_one_cohn_coefficients_from_checked_data(data)
    children = _rank_one_child_certificates_from_checked_data(data, coefficients)
    factors = _rank_one_factors_from_children(children)
    target = _rank_one_target_from_checked_data(data)
    product = _cohn_type_factor_product(factors, R, data.n)
    provisional = RankOneNormalityCertificate(
        data.n,
        data.v,
        data.w,
        data.g,
        data.ring,
        data.orthogonality,
        data.bezout,
        coefficients,
        children,
        factors,
        target,
        product,
        nothing,
    )
    verification = _rank_one_certificate_core_verification(provisional)
    verification.overall_core_ok ||
        error("internal rank-one normality certificate verification failed")
    return RankOneNormalityCertificate(
        provisional.n,
        provisional.v,
        provisional.w,
        provisional.g,
        provisional.ring,
        provisional.orthogonality,
        provisional.bezout,
        provisional.cohn_coefficients,
        provisional.child_certificates,
        provisional.factors,
        provisional.target,
        provisional.product,
        verification,
    )
end

function _rank_one_certificate_core_verification(cert)
    ring_ok = false
    checked_inputs_ok = false
    orthogonality_ok = false
    bezout_ok = false
    coefficient_table_ok = false
    child_count = 0
    child_certificates_ok = false
    factor_sequence_ok = false
    factor_count = 0
    product_replay_ok = false
    product_matches_stored_ok = false
    target_replay_ok = false
    target_matches_product_ok = false

    _require_ordinary_polynomial_certificate_ring(cert.ring)
    ring_ok = true

    data = _rank_one_checked_data(cert.v, cert.w, cert.g, cert.ring)
    checked_inputs_ok = true
    orthogonality_ok = cert.orthogonality == zero(cert.ring) && cert.orthogonality == data.orthogonality
    bezout_ok = cert.bezout == one(cert.ring) && cert.bezout == data.bezout
    expected_coefficients = _rank_one_cohn_coefficients_from_checked_data(data)
    coefficient_table_ok = cert.cohn_coefficients == expected_coefficients
    expected_children = _rank_one_child_certificates_from_checked_data(data, expected_coefficients)
    child_count = length(cert.child_certificates)
    child_certificates_ok =
        child_count == length(expected_children) &&
        all(verify_cohn_type_certificate, cert.child_certificates) &&
        all(zip(cert.child_certificates, expected_children)) do (actual, expected)
            actual.n == expected.n &&
                actual.i == expected.i &&
                actual.j == expected.j &&
                actual.a == expected.a &&
                actual.v == expected.v &&
                _same_base_ring(actual.ring, expected.ring) &&
                actual.auxiliary_index == expected.auxiliary_index &&
                actual.target == expected.target &&
                actual.factors == expected.factors &&
                actual.product == expected.product &&
                actual.verification == expected.verification
        end
    expected_factors = _rank_one_factors_from_children(cert.child_certificates)
    factor_sequence_ok = cert.factors == expected_factors
    factor_count = length(cert.factors)
    replayed_product = _cohn_type_factor_product(cert.factors, cert.ring, cert.n)
    product_replay_ok = true
    product_matches_stored_ok = replayed_product == cert.product
    expected_target = _rank_one_target_from_checked_data(data)
    target_replay_ok = cert.target == expected_target
    target_matches_product_ok = target_replay_ok && replayed_product == cert.target
    overall_core_ok =
        ring_ok &&
        checked_inputs_ok &&
        orthogonality_ok &&
        bezout_ok &&
        coefficient_table_ok &&
        child_certificates_ok &&
        factor_sequence_ok &&
        product_replay_ok &&
        product_matches_stored_ok &&
        target_matches_product_ok
    return (;
        ring_ok,
        checked_inputs_ok,
        orthogonality_ok,
        bezout_ok,
        coefficient_table_ok,
        child_count,
        child_certificates_ok,
        factor_sequence_ok,
        factor_count,
        product_replay_ok,
        product_matches_stored_ok,
        target_replay_ok,
        target_matches_product_ok,
        overall_core_ok,
    )
end

function _rank_one_certificate_verification(cert)
    core = _rank_one_certificate_core_verification(cert)
    stored_verification_ok = cert.verification == core
    return merge(core, (;
        stored_verification_ok,
        overall_ok = core.overall_core_ok && stored_verification_ok,
    ))
end

function verify_rank_one_normality_certificate(cert)::Bool
    try
        return _rank_one_certificate_verification(cert).overall_ok
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end

struct ConjugatedElementaryNormalityCertificate
    n::Int
    A
    i::Int
    j::Int
    a
    ring
    determinant
    inverse_A
    elementary_matrix
    conjugation_convention::Symbol
    conjugation_target
    v::Vector
    w::Vector
    g::Vector
    rank_one_certificate::RankOneNormalityCertificate
    factors::Vector
    product
    verification
end

function _conjugate_elementary_checked_data(A, i::Int, j::Int, a)
    nrows(A) == ncols(A) || throw(ArgumentError("A must be square"))

    n = nrows(A)
    n >= 3 || throw(ArgumentError("A must have size at least 3"))
    1 <= i <= n || throw(ArgumentError("i must be between 1 and n"))
    1 <= j <= n || throw(ArgumentError("j must be between 1 and n"))
    i == j && throw(ArgumentError("i and j must differ"))

    R = base_ring(A)
    _require_ordinary_polynomial_certificate_ring(R)
    coerced_a = _coerce_into_ring(R, a, "a")
    determinant = det(A)
    determinant == one(R) || throw(ArgumentError("A must have determinant one"))
    inverse_A = inv(A)
    identity = identity_matrix(R, n)
    A * inverse_A == identity || throw(ArgumentError("A inverse replay failed"))
    inverse_A * A == identity || throw(ArgumentError("A inverse replay failed"))

    elementary = elementary_matrix(n, i, j, coerced_a, R)
    conjugation_target = A * elementary * inverse_A
    v = [A[row, i] for row in 1:n]
    w = [coerced_a * inverse_A[j, col] for col in 1:n]
    g = [inverse_A[i, col] for col in 1:n]
    rank_one_certificate = realize_rank_one_normality_certificate(v, w, g, R)
    rank_one_certificate.target == conjugation_target ||
        error("internal conjugated elementary target did not match rank-one target")

    return (;
        n,
        A,
        i,
        j,
        a = coerced_a,
        ring = R,
        determinant,
        inverse_A,
        elementary_matrix = elementary,
        conjugation_convention = :A_E_invA,
        conjugation_target,
        v,
        w,
        g,
        rank_one_certificate,
    )
end

function realize_conjugate_elementary_certificate(A, i::Int, j::Int, a)
    data = _conjugate_elementary_checked_data(A, i, j, a)
    factors = data.rank_one_certificate.factors
    product = _cohn_type_factor_product(factors, data.ring, data.n)
    provisional = ConjugatedElementaryNormalityCertificate(
        data.n,
        data.A,
        data.i,
        data.j,
        data.a,
        data.ring,
        data.determinant,
        data.inverse_A,
        data.elementary_matrix,
        data.conjugation_convention,
        data.conjugation_target,
        data.v,
        data.w,
        data.g,
        data.rank_one_certificate,
        factors,
        product,
        nothing,
    )
    verification = _conjugate_elementary_certificate_core_verification(provisional)
    verification.overall_core_ok ||
        error("internal conjugated elementary certificate verification failed")
    return ConjugatedElementaryNormalityCertificate(
        provisional.n,
        provisional.A,
        provisional.i,
        provisional.j,
        provisional.a,
        provisional.ring,
        provisional.determinant,
        provisional.inverse_A,
        provisional.elementary_matrix,
        provisional.conjugation_convention,
        provisional.conjugation_target,
        provisional.v,
        provisional.w,
        provisional.g,
        provisional.rank_one_certificate,
        provisional.factors,
        provisional.product,
        verification,
    )
end

function _conjugate_elementary_certificate_core_verification(cert)
    ring_ok = false
    checked_inputs_ok = false
    n_ok = false
    determinant_ok = false
    inverse_replay_ok = false
    elementary_matrix_ok = false
    convention_ok = false
    target_replay_ok = false
    vector_replay_ok = false
    rank_one_certificate_ok = false
    factor_sequence_ok = false
    factor_count = 0
    product_replay_ok = false
    product_matches_stored_ok = false
    target_matches_product_ok = false

    data = _conjugate_elementary_checked_data(cert.A, cert.i, cert.j, cert.a)
    checked_inputs_ok = true
    ring_ok = _same_base_ring(cert.ring, data.ring)
    n_ok = cert.n == data.n
    determinant_ok = cert.determinant == one(data.ring) && cert.determinant == data.determinant
    identity = identity_matrix(data.ring, data.n)
    inverse_replay_ok =
        cert.inverse_A == data.inverse_A &&
        cert.A * cert.inverse_A == identity &&
        cert.inverse_A * cert.A == identity
    elementary_matrix_ok = cert.elementary_matrix == data.elementary_matrix
    convention_ok =
        cert.conjugation_convention == :A_E_invA &&
        cert.conjugation_convention == data.conjugation_convention
    target_replay_ok = cert.conjugation_target == data.conjugation_target
    vector_replay_ok = cert.v == data.v && cert.w == data.w && cert.g == data.g
    rank_one_certificate_ok =
        verify_rank_one_normality_certificate(cert.rank_one_certificate) &&
        cert.rank_one_certificate.n == data.n &&
        cert.rank_one_certificate.v == data.v &&
        cert.rank_one_certificate.w == data.w &&
        cert.rank_one_certificate.g == data.g &&
        _same_base_ring(cert.rank_one_certificate.ring, data.ring) &&
        cert.rank_one_certificate.target == data.conjugation_target &&
        cert.rank_one_certificate.factors == data.rank_one_certificate.factors &&
        cert.rank_one_certificate.product == data.rank_one_certificate.product &&
        cert.rank_one_certificate.verification == data.rank_one_certificate.verification
    factor_sequence_ok = cert.factors == cert.rank_one_certificate.factors
    factor_count = length(cert.factors)
    replayed_product = _cohn_type_factor_product(cert.factors, data.ring, data.n)
    product_replay_ok = true
    product_matches_stored_ok = replayed_product == cert.product
    target_matches_product_ok = target_replay_ok && replayed_product == cert.conjugation_target

    overall_core_ok =
        ring_ok &&
        checked_inputs_ok &&
        n_ok &&
        determinant_ok &&
        inverse_replay_ok &&
        elementary_matrix_ok &&
        convention_ok &&
        target_replay_ok &&
        vector_replay_ok &&
        rank_one_certificate_ok &&
        factor_sequence_ok &&
        product_replay_ok &&
        product_matches_stored_ok &&
        target_matches_product_ok

    return (;
        ring_ok,
        checked_inputs_ok,
        n_ok,
        determinant_ok,
        inverse_replay_ok,
        elementary_matrix_ok,
        convention_ok,
        target_replay_ok,
        vector_replay_ok,
        rank_one_certificate_ok,
        factor_sequence_ok,
        factor_count,
        product_replay_ok,
        product_matches_stored_ok,
        target_matches_product_ok,
        overall_core_ok,
    )
end

function _conjugate_elementary_certificate_verification(cert)
    core = _conjugate_elementary_certificate_core_verification(cert)
    stored_verification_ok = cert.verification == core
    return merge(core, (;
        stored_verification_ok,
        overall_ok = core.overall_core_ok && stored_verification_ok,
    ))
end

function verify_conjugate_elementary_certificate(cert)::Bool
    try
        return _conjugate_elementary_certificate_verification(cert).overall_ok
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end

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

    is_unit(determinant) || throw(ArgumentError("B must have a Laurent-unit determinant"))
    determinant_inverse = inv(determinant)

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
