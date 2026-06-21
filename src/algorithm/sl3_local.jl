struct SL3LocalRealizationCertificate
    target
    branch::Symbol
    factors::Vector
    selected_variable
    witness
end

struct SL3LocalSplitLemmaReplay
    split_id::Symbol
    original_target
    first_child_target
    second_child_target
    prefix_factors::Vector
    middle_factors::Vector
    suffix_factors::Vector
    wrapper_factors::Vector
    reassembled_product
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

function sl3_local_split_lemma_replay(
    a,
    a_prime,
    b,
    c,
    d,
    c1,
    d1,
    c2,
    d2;
    split_id::Symbol = :murthy_split_lemma,
)
    R = parent(a)
    normalized = _sl3_local_split_lemma_normalize_inputs(R, a, a_prime, b, c, d, c1, d1, c2, d2)
    a, a_prime, b, c, d, c1, d1, c2, d2 = normalized

    a * a_prime * d - b * c == one(R) ||
        throw(ArgumentError("Murthy split lemma original determinant relation failed"))
    a * d1 - b * c1 == one(R) ||
        throw(ArgumentError("Murthy split lemma first child determinant relation failed"))
    a_prime * d2 - b * c2 == one(R) ||
        throw(ArgumentError("Murthy split lemma second child determinant relation failed"))

    original_target = _sl3_local_special_form_target(R, a * a_prime, b, c, d)
    first_child_target = _sl3_local_special_form_target(R, a, b, c1, d1)
    second_child_target = _sl3_local_special_form_target(R, a_prime, b, c2, d2)
    segments = _sl3_local_split_lemma_wrapper_segments(R, a, a_prime, b, c, d, c1, d1, c2, d2)
    reassembled_product = _sl3_local_split_lemma_reassembled_product(
        segments.prefix_factors,
        first_child_target,
        segments.middle_factors,
        second_child_target,
        segments.suffix_factors,
        R,
    )

    reassembled_product == original_target ||
        throw(ArgumentError("Murthy split lemma elementary wrapper replay failed exact reassembly"))

    replay = SL3LocalSplitLemmaReplay(
        split_id,
        original_target,
        first_child_target,
        second_child_target,
        segments.prefix_factors,
        segments.middle_factors,
        segments.suffix_factors,
        vcat(segments.prefix_factors, segments.middle_factors, segments.suffix_factors),
        reassembled_product,
        (; a, a_prime, b, c, d, c1, d1, c2, d2),
    )
    verify_sl3_local_split_lemma_replay(replay) ||
        error("internal Murthy split lemma replay verification failed")
    return replay
end

function sl3_local_split_lemma_certificate(
    replay::SL3LocalSplitLemmaReplay,
    first_child_certificate,
    second_child_certificate,
    X,
)
    verify_sl3_local_split_lemma_replay(replay) ||
        throw(ArgumentError("Murthy split lemma replay must verify before certificate construction"))
    verify_sl3_local_realization(first_child_certificate) ||
        throw(ArgumentError("first child certificate must verify"))
    verify_sl3_local_realization(second_child_certificate) ||
        throw(ArgumentError("second child certificate must verify"))
    first_child_certificate.target == replay.first_child_target ||
        throw(ArgumentError("first child certificate target does not match split replay"))
    second_child_certificate.target == replay.second_child_target ||
        throw(ArgumentError("second child certificate target does not match split replay"))

    R = base_ring(replay.original_target)
    parent(X) === R || throw(ArgumentError("split certificate variable must lie in the target ring"))
    first_child_certificate.selected_variable == X ||
        throw(ArgumentError("first child certificate variable does not match split certificate variable"))
    second_child_certificate.selected_variable == X ||
        throw(ArgumentError("second child certificate variable does not match split certificate variable"))

    factors = vcat(
        replay.prefix_factors,
        first_child_certificate.factors,
        replay.middle_factors,
        second_child_certificate.factors,
        replay.suffix_factors,
    )
    certificate = SL3LocalRealizationCertificate(
        replay.original_target,
        :murthy_split_lemma,
        factors,
        X,
        (; split = replay, first_child_certificate, second_child_certificate),
    )
    verify_sl3_local_realization(certificate) ||
        error("internal Murthy split lemma certificate verification failed")
    return certificate
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

function _sl3_local_split_lemma_normalize_inputs(R, values...)
    return tuple((_coerce_into_ring(R, value, "split lemma input") for value in values)...)
end

function _sl3_local_special_form_target(R, p, q, r, s)
    return matrix(R, [
        p q zero(R);
        r s zero(R);
        zero(R) zero(R) one(R)
    ])
end

function _sl3_local_split_lemma_wrapper_segments(R, a, a_prime, b, c, d, c1, d1, c2, d2)
    prefix_factors = [
        elementary_matrix(3, 2, 1, c * d1 * d2 - d * (c2 + a_prime * c1 * d2), R),
        elementary_matrix(3, 2, 3, d2 - one(R), R),
        elementary_matrix(3, 3, 2, one(R), R),
        elementary_matrix(3, 2, 3, -one(R), R),
    ]
    middle_factors = [
        elementary_matrix(3, 2, 3, one(R), R),
        elementary_matrix(3, 3, 2, -one(R), R),
        elementary_matrix(3, 2, 3, one(R), R),
    ]
    suffix_factors = [
        elementary_matrix(3, 2, 3, -one(R), R),
        elementary_matrix(3, 3, 2, one(R), R),
        elementary_matrix(3, 2, 3, a - one(R), R),
        elementary_matrix(3, 3, 1, -a_prime * c1, R),
        elementary_matrix(3, 3, 2, -d1, R),
    ]
    return (; prefix_factors, middle_factors, suffix_factors)
end

function _sl3_local_factor_product(factors, R)
    product = identity_matrix(R, 3)
    for factor in factors
        product *= factor
    end
    return product
end

function _sl3_local_split_lemma_reassembled_product(
    prefix_factors,
    first_child_target,
    middle_factors,
    second_child_target,
    suffix_factors,
    R,
)
    return _sl3_local_factor_product(prefix_factors, R) *
        first_child_target *
        _sl3_local_factor_product(middle_factors, R) *
        second_child_target *
        _sl3_local_factor_product(suffix_factors, R)
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

function verify_sl3_local_split_lemma_replay(replay)::Bool
    try
        return _sl3_local_split_lemma_verification(replay).overall_ok
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end

function _sl3_local_split_lemma_verification(replay)
    split_id_ok = replay.split_id isa Symbol
    original_size_ok = nrows(replay.original_target) == 3 && ncols(replay.original_target) == 3
    first_size_ok = nrows(replay.first_child_target) == 3 && ncols(replay.first_child_target) == 3
    second_size_ok = nrows(replay.second_child_target) == 3 && ncols(replay.second_child_target) == 3
    size_ok = original_size_ok && first_size_ok && second_size_ok
    R = size_ok ? base_ring(replay.original_target) : nothing
    target_ring_ok = size_ok &&
        _same_base_ring(base_ring(replay.first_child_target), R) &&
        _same_base_ring(base_ring(replay.second_child_target), R)
    witness = replay.witness
    witness_keys_ok = _sl3_local_split_lemma_witness_keys_ok(witness)
    witness_ring_ok = target_ring_ok && witness_keys_ok && _sl3_local_split_lemma_witness_ring_ok(witness, R)

    relations_ok = false
    targets_ok = false
    determinant_ok = false
    wrapper_segments_ok = false
    wrapper_factors_ok = false
    reassembled_product_ok = false
    if witness_ring_ok
        a = witness.a
        a_prime = witness.a_prime
        b = witness.b
        c = witness.c
        d = witness.d
        c1 = witness.c1
        d1 = witness.d1
        c2 = witness.c2
        d2 = witness.d2

        relations_ok =
            a * a_prime * d - b * c == one(R) &&
            a * d1 - b * c1 == one(R) &&
            a_prime * d2 - b * c2 == one(R)
        expected_original = _sl3_local_special_form_target(R, a * a_prime, b, c, d)
        expected_first = _sl3_local_special_form_target(R, a, b, c1, d1)
        expected_second = _sl3_local_special_form_target(R, a_prime, b, c2, d2)
        targets_ok =
            replay.original_target == expected_original &&
            replay.first_child_target == expected_first &&
            replay.second_child_target == expected_second
        determinant_ok =
            det(replay.original_target) == one(R) &&
            det(replay.first_child_target) == one(R) &&
            det(replay.second_child_target) == one(R)

        expected_segments = _sl3_local_split_lemma_wrapper_segments(
            R,
            a,
            a_prime,
            b,
            c,
            d,
            c1,
            d1,
            c2,
            d2,
        )
        expected_wrapper_factors = vcat(
            expected_segments.prefix_factors,
            expected_segments.middle_factors,
            expected_segments.suffix_factors,
        )
        wrapper_segments_ok =
            _sl3_local_factor_sequences_equal(replay.prefix_factors, expected_segments.prefix_factors) &&
            _sl3_local_factor_sequences_equal(replay.middle_factors, expected_segments.middle_factors) &&
            _sl3_local_factor_sequences_equal(replay.suffix_factors, expected_segments.suffix_factors)
        wrapper_factors_ok =
            _sl3_local_factor_sequences_equal(replay.wrapper_factors, expected_wrapper_factors)
        expected_reassembled = _sl3_local_split_lemma_reassembled_product(
            expected_segments.prefix_factors,
            expected_first,
            expected_segments.middle_factors,
            expected_second,
            expected_segments.suffix_factors,
            R,
        )
        reassembled_product_ok =
            replay.reassembled_product == expected_reassembled &&
            expected_reassembled == replay.original_target
    end

    overall_ok = split_id_ok && size_ok && target_ring_ok && witness_keys_ok &&
        witness_ring_ok && relations_ok && targets_ok && determinant_ok &&
        wrapper_segments_ok && wrapper_factors_ok && reassembled_product_ok
    return (;
        overall_ok,
        split_id_ok,
        size_ok,
        target_ring_ok,
        witness_keys_ok,
        witness_ring_ok,
        relations_ok,
        targets_ok,
        determinant_ok,
        wrapper_segments_ok,
        wrapper_factors_ok,
        reassembled_product_ok,
    )
end

function _sl3_local_split_lemma_witness_keys_ok(witness)
    witness isa NamedTuple || return false
    return keys(witness) == (:a, :a_prime, :b, :c, :d, :c1, :d1, :c2, :d2)
end

function _sl3_local_split_lemma_witness_ring_ok(witness, R)
    for key in keys(witness)
        parent(getproperty(witness, key)) === R || return false
    end
    return true
end

function _sl3_local_realization_verification(certificate)
    target = certificate.target
    size_ok = nrows(target) == 3 && ncols(target) == 3
    R = size_ok ? base_ring(target) : nothing
    # Replay checks selected_variable is a legal target-ring generator; branch
    # witness equations and exact factors carry the algebraic replay.
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
    elseif certificate.branch == :murthy_split_lemma
        return vcat(
            witness.split.prefix_factors,
            witness.first_child_certificate.factors,
            witness.split.middle_factors,
            witness.second_child_certificate.factors,
            witness.split.suffix_factors,
        )
    end
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
    _sl3_local_witness_keys_ok(certificate.branch, witness) || return false

    if certificate.branch == :open_s_one
        return witness.q == q && witness.r == r && s == one(R) && p == one(R) + q * r
    elseif certificate.branch == :open_p_one
        return witness.q == q && witness.r == r && p == one(R) && s == one(R) + q * r
    elseif certificate.branch == :s_unit
        return witness.pivot == s && witness.pivot * witness.pivot_inverse == one(R)
    elseif certificate.branch == :p_unit
        return witness.pivot == p && witness.pivot * witness.pivot_inverse == one(R)
    elseif certificate.branch == :murthy_split_lemma
        return _sl3_local_split_certificate_witness_ok(certificate)
    end
end

function _sl3_local_split_certificate_witness_ok(certificate)
    witness = certificate.witness
    split = witness.split
    first_child_certificate = witness.first_child_certificate
    second_child_certificate = witness.second_child_certificate

    verify_sl3_local_split_lemma_replay(split) || return false
    split.original_target == certificate.target || return false
    first_child_certificate.target == split.first_child_target || return false
    second_child_certificate.target == split.second_child_target || return false
    verify_sl3_local_realization(first_child_certificate) || return false
    verify_sl3_local_realization(second_child_certificate) || return false
    first_child_certificate.selected_variable == certificate.selected_variable || return false
    second_child_certificate.selected_variable == certificate.selected_variable || return false
    return true
end

function _sl3_local_witness_keys_ok(branch, witness)
    witness isa NamedTuple || return false
    if branch in (:open_s_one, :open_p_one)
        return keys(witness) == (:q, :r)
    elseif branch in (:s_unit, :p_unit)
        return keys(witness) == (:pivot, :pivot_inverse)
    elseif branch == :murthy_split_lemma
        return keys(witness) == (:split, :first_child_certificate, :second_child_certificate)
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
