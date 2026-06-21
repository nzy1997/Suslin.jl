struct SL3LocalRealizationCertificate
    target
    branch::Symbol
    factors::Vector
    selected_variable
    witness
end

struct SL3LocalQDegreeNormalization
    target
    quotient
    remainder
    normalized_target
    elementary_correction
    selected_variable
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

struct SL3LocalMurthyQUnitReduction
    target
    q0
    q0_inverse
    p0
    right_e21_coefficient
    eliminated_target
    elimination_factor
    inverse_elimination_factor
    p_prime
    split_certificate
    selected_variable
    degree_p::Int
    degree_p_prime::Int
    locality_witness
end

struct SL3LocalMurthyQ0NonunitReduction
    target
    p0
    q0
    p_prime
    q_prime
    resultant
    bezout_target
    child_link_target
    left_factor
    first_elementary_factor
    child_certificate
    selected_variable
    degree_p::Int
    degree_q::Int
    degree_p_prime::Int
    degree_q_prime::Int
    branch_unit
    branch_unit_inverse
    witness_source::Symbol
end

function Base.:(==)(left::SL3LocalQDegreeNormalization, right::SL3LocalQDegreeNormalization)
    return left.target == right.target &&
        left.quotient == right.quotient &&
        left.remainder == right.remainder &&
        left.normalized_target == right.normalized_target &&
        left.elementary_correction == right.elementary_correction &&
        left.selected_variable == right.selected_variable
end

function realize_sl3_local(A, X; check_monic::Bool=true, murthy_q0_nonunit_witness=nothing)
    return realize_sl3_local_certificate(
        A,
        X;
        check_monic,
        murthy_q0_nonunit_witness,
    ).factors
end

function realize_sl3_local(p, q, r, s, X; check_monic::Bool=true, murthy_q0_nonunit_witness=nothing)
    return realize_sl3_local_certificate(
        p,
        q,
        r,
        s,
        X;
        check_monic,
        murthy_q0_nonunit_witness,
    ).factors
end

function realize_sl3_local_certificate(A, X; check_monic::Bool=true, murthy_q0_nonunit_witness=nothing)
    form = _recognize_sl3_local_matrix(A, X; check_monic, murthy_q0_nonunit_witness)
    return _realize_sl3_local_certificate_form(form)
end

function realize_sl3_local_certificate(
    p,
    q,
    r,
    s,
    X;
    check_monic::Bool=true,
    murthy_q0_nonunit_witness=nothing,
)
    form = _recognize_sl3_local_parameters(
        p,
        q,
        r,
        s,
        X;
        check_monic,
        murthy_q0_nonunit_witness,
    )
    return _realize_sl3_local_certificate_form(form)
end

function sl3_local_q_degree_normalization(A, X; check_monic::Bool=true)
    nrows(A) == 3 || throw(ArgumentError("q-degree normalization requires a local SL_3 special-form matrix"))
    ncols(A) == 3 || throw(ArgumentError("q-degree normalization requires a local SL_3 special-form matrix"))
    entries = _sl3_local_target_entries(A)
    entries === nothing && throw(ArgumentError("q-degree normalization requires a local SL_3 special-form matrix"))
    return sl3_local_q_degree_normalization(entries.p, entries.q, entries.r, entries.s, X; check_monic)
end

function sl3_local_q_degree_normalization(p, q, r, s, X; check_monic::Bool=true)
    form = _recognize_sl3_local_q_degree_normalization_parameters(p, q, r, s, X; check_monic)
    quotient, remainder = _sl3_local_divrem_monic_in_variable(form.q, form.p, form.var_idx, form.R)
    normalized_target = _sl3_local_special_form_target(
        form.R,
        form.p,
        remainder,
        form.r,
        form.s - quotient * form.r,
    )
    record = SL3LocalQDegreeNormalization(
        form.target,
        quotient,
        remainder,
        normalized_target,
        elementary_matrix(3, 1, 2, quotient, form.R),
        form.X,
    )
    verify_sl3_local_q_degree_normalization(record) ||
        error("internal Murthy q-degree normalization verification failed")
    return record
end

function sl3_local_q_degree_normalization_certificate(record::SL3LocalQDegreeNormalization)
    verify_sl3_local_q_degree_normalization(record) ||
        throw(ArgumentError("Murthy q-degree normalization record must verify before certificate construction"))
    certificate = SL3LocalRealizationCertificate(
        record.target,
        :murthy_q_degree_normalization,
        [record.normalized_target, record.elementary_correction],
        record.selected_variable,
        (; normalization = record),
    )
    verify_sl3_local_realization(certificate) ||
        error("internal Murthy q-degree normalization certificate verification failed")
    return certificate
end

function sl3_local_q_degree_normalization_certificate(args...; check_monic::Bool=true)
    return sl3_local_q_degree_normalization_certificate(
        sl3_local_q_degree_normalization(args...; check_monic),
    )
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

function _recognize_sl3_local_matrix(A, X; check_monic::Bool=true, murthy_q0_nonunit_witness=nothing)
    nrows(A) == 3 || throw(ArgumentError("local SL_3 special-form recognition failed: A must be 3x3"))
    ncols(A) == 3 || throw(ArgumentError("local SL_3 special-form recognition failed: A must be 3x3"))

    R = base_ring(A)
    parent(X) === R || throw(ArgumentError("A and X must lie in the same polynomial or Laurent polynomial ring"))
    if A[1, 3] != zero(R) || A[2, 3] != zero(R) ||
            A[3, 1] != zero(R) || A[3, 2] != zero(R) ||
            A[3, 3] != one(R)
        throw(ArgumentError("local SL_3 special-form recognition failed: A must be an embedded 2x2 block with trailing identity"))
    end

    return _recognize_sl3_local_parameters(
        A[1, 1],
        A[1, 2],
        A[2, 1],
        A[2, 2],
        X;
        check_monic,
        murthy_q0_nonunit_witness,
    )
end

function _recognize_sl3_local_parameters(
    p,
    q,
    r,
    s,
    X;
    check_monic::Bool=true,
    murthy_q0_nonunit_witness=nothing,
)
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

    if _sl3_local_supports_murthy_q0_unit_branch(R, var_idx) &&
            _is_monic_in_variable(p, var_idx, R)
        degree_p = degree(p, var_idx)
        degree_q = degree(q, var_idx)
        if degree_q >= degree_p
            return (; family = :murthy_q0_unit, R, p, q, r, s, X, target, var_idx, murthy_q0_nonunit_witness)
        end

        q_inverse = _unit_inverse_or_nothing(q)
        p0 = _sl3_local_constant_coefficient(p, var_idx, R)
        # This is the terminal q-unit child after q(0)-unit elimination has
        # made the first entry divisible by X. Inputs with nonzero p(0) stay
        # on the recursive Murthy branch so the q(0)-unit step is replayed.
        if q_inverse !== nothing && p0 == zero(R)
            return (; family = :q_unit, R, p, q, r, s, X, target, pivot_inverse = q_inverse)
        end

        q0 = _sl3_local_constant_coefficient(q, var_idx, R)
        if _unit_inverse_or_nothing(q0) !== nothing
            return (; family = :murthy_q0_unit, R, p, q, r, s, X, target, var_idx, murthy_q0_nonunit_witness)
        end

        return (; family = :murthy_q0_nonunit_bezout_resultant, R, p, q, r, s, X, target, var_idx, murthy_q0_nonunit_witness)
    end

    q_inverse = _unit_inverse_or_nothing(q)
    if q_inverse !== nothing
        return (; family = :q_unit, R, p, q, r, s, X, target, pivot_inverse = q_inverse)
    end

    _throw_staged_sl3_local_failure("supported families require one open unipotent slice, a unit pivot p/s/q, or ordinary univariate monic p in the Murthy q(0)-unit or Bezout/resultant branches")
end

function _recognize_sl3_local_q_degree_normalization_parameters(p, q, r, s, X; check_monic::Bool=true)
    R = parent(p)
    parent(q) === R || throw(ArgumentError("p, q, r, s, and X must lie in the same polynomial ring"))
    parent(r) === R || throw(ArgumentError("p, q, r, s, and X must lie in the same polynomial ring"))
    parent(s) === R || throw(ArgumentError("p, q, r, s, and X must lie in the same polynomial ring"))
    parent(X) === R || throw(ArgumentError("p, q, r, s, and X must lie in the same polynomial ring"))
    _is_laurent_polynomial_ring(R) &&
        throw(ArgumentError("q-degree normalization is only supported for ordinary polynomial local SL_3 inputs"))

    ring_gens = collect(gens(R))
    var_idx = findfirst(isequal(X), ring_gens)
    var_idx === nothing && throw(ArgumentError("X must be one of the polynomial ring generators"))
    if check_monic
        _is_monic_in_variable(p, var_idx, R) || throw(ArgumentError("p must be monic in X"))
    end

    target = _sl3_local_special_form_target(R, p, q, r, s)
    det(target) == one(R) || throw(ArgumentError("constructed matrix must have determinant 1"))
    return (; R, p, q, r, s, X, var_idx, target)
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
    elseif form.family == :q_unit
        q_inverse = form.pivot_inverse
        return [
            elementary_matrix(3, 2, 1, (form.s - one(R)) * q_inverse, R),
            elementary_matrix(3, 1, 2, form.q, R),
            elementary_matrix(3, 2, 1, (form.p - one(R)) * q_inverse, R),
        ]
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
    elseif form.family == :q_unit
        return (; pivot = form.q, pivot_inverse = form.pivot_inverse)
    end

    _throw_staged_sl3_local_failure("unrecognized local solver family $(form.family)")
end

function _realize_sl3_local_certificate_form(form)
    if form.family == :murthy_q0_unit
        return _realize_sl3_local_murthy_q0_unit_certificate(form)
    elseif form.family == :murthy_q0_nonunit_bezout_resultant
        return _realize_sl3_local_murthy_q0_nonunit_certificate(form)
    end

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

function _realize_sl3_local_murthy_q0_unit_certificate(form)
    degree_q = degree(form.q, form.var_idx)
    degree_p = degree(form.p, form.var_idx)
    if degree_q >= degree_p
        normalization = sl3_local_q_degree_normalization(form.p, form.q, form.r, form.s, form.X)
        normalized_certificate = realize_sl3_local_certificate(
            normalization.normalized_target,
            form.X;
            murthy_q0_nonunit_witness = form.murthy_q0_nonunit_witness,
        )
        certificate = SL3LocalRealizationCertificate(
            form.target,
            :murthy_q0_unit,
            vcat(normalized_certificate.factors, [normalization.elementary_correction]),
            form.X,
            (; normalization, normalized_certificate, reduction = nothing),
        )
        verify_sl3_local_realization(certificate) ||
            error("internal Murthy q(0)-unit normalization certificate verification failed")
        return certificate
    end

    reduction = _sl3_local_murthy_q0_unit_reduction(form)
    certificate = SL3LocalRealizationCertificate(
        form.target,
        :murthy_q0_unit,
        vcat(reduction.split_certificate.factors, [reduction.inverse_elimination_factor]),
        form.X,
        (; normalization = nothing, normalized_certificate = nothing, reduction),
    )
    verify_sl3_local_realization(certificate) ||
        error("internal Murthy q(0)-unit certificate verification failed")
    return certificate
end

function _sl3_local_murthy_q0_unit_reduction(form)
    R = form.R
    degree_p = degree(form.p, form.var_idx)
    q0 = _sl3_local_constant_coefficient(form.q, form.var_idx, R)
    q0_inverse = _unit_inverse_or_nothing(q0)
    q0_inverse === nothing &&
        _throw_staged_sl3_local_failure("Murthy q(0)-nonunit Bezout/resultant branch is not implemented")

    p0 = _sl3_local_constant_coefficient(form.p, form.var_idx, R)
    right_e21_coefficient = -q0_inverse * p0
    elimination_factor = elementary_matrix(3, 2, 1, right_e21_coefficient, R)
    inverse_elimination_factor = elementary_matrix(3, 2, 1, -right_e21_coefficient, R)
    eliminated_target = form.target * elimination_factor
    eliminated_entries = _sl3_local_target_entries(eliminated_target)
    eliminated_entries === nothing &&
        error("internal Murthy q(0)-unit elimination left special form")

    eliminated_p = eliminated_entries.p
    p_prime = divexact(eliminated_p, form.X)
    eliminated_p == form.X * p_prime ||
        error("internal Murthy q(0)-unit elimination did not produce an X-divisible p entry")
    degree_p_prime = degree(p_prime, form.var_idx)
    degree_p_prime < degree_p ||
        error("internal Murthy q(0)-unit recursion guard failed to decrease degree")

    gcd_x_q, d1, minus_c1 = gcdx(form.X, form.q)
    gcd_x_q == one(R) ||
        error("internal Murthy q(0)-unit split witness gcd(X, q) was not 1")
    gcd_p_prime_q, d2, minus_c2 = gcdx(p_prime, form.q)
    gcd_p_prime_q == one(R) ||
        error("internal Murthy q(0)-unit split witness gcd(p_prime, q) was not 1")
    c1 = -minus_c1
    c2 = -minus_c2
    split_replay = sl3_local_split_lemma_replay(
        form.X,
        p_prime,
        form.q,
        eliminated_entries.r,
        eliminated_entries.s,
        c1,
        d1,
        c2,
        d2;
        split_id = :murthy_q0_unit_split,
    )
    first_child_certificate = realize_sl3_local_certificate(split_replay.first_child_target, form.X)
    second_child_certificate = realize_sl3_local_certificate(split_replay.second_child_target, form.X)
    split_certificate = sl3_local_split_lemma_certificate(
        split_replay,
        first_child_certificate,
        second_child_certificate,
        form.X,
    )

    reduction = SL3LocalMurthyQUnitReduction(
        form.target,
        q0,
        q0_inverse,
        p0,
        right_e21_coefficient,
        eliminated_target,
        elimination_factor,
        inverse_elimination_factor,
        p_prime,
        split_certificate,
        form.X,
        degree_p,
        degree_p_prime,
        (;
            eliminated_p,
            eliminated_r = eliminated_entries.r,
            eliminated_s = eliminated_entries.s,
        ),
    )
    verify_sl3_local_murthy_q_unit_reduction(reduction) ||
        error("internal Murthy q(0)-unit reduction verification failed")
    return reduction
end

function _realize_sl3_local_murthy_q0_nonunit_certificate(form)
    reduction = _sl3_local_murthy_q0_nonunit_reduction(form)
    certificate = SL3LocalRealizationCertificate(
        form.target,
        :murthy_q0_nonunit_bezout_resultant,
        vcat(
            [reduction.left_factor, reduction.first_elementary_factor],
            reduction.child_certificate.factors,
        ),
        form.X,
        (; reduction),
    )
    verify_sl3_local_realization(certificate) ||
        error("internal Murthy q(0)-nonunit Bezout/resultant certificate verification failed")
    return certificate
end

function _sl3_local_murthy_q0_nonunit_reduction(form)
    p_prime, q_prime, witness_source = _sl3_local_murthy_q0_nonunit_bezout_pair(form)
    resultant = p_prime * form.p - q_prime * form.q
    resultant == one(form.R) ||
        throw(ArgumentError("Murthy q(0)-nonunit Bezout equality p_prime*p - q_prime*q must equal 1"))

    degree_p = degree(form.p, form.var_idx)
    degree_q = degree(form.q, form.var_idx)
    degree_p_prime = degree(p_prime, form.var_idx)
    degree_q_prime = degree(q_prime, form.var_idx)
    degree_p_prime < degree_q ||
        throw(ArgumentError("Murthy q(0)-nonunit p_prime degree guard failed"))
    degree_q_prime < degree_p ||
        throw(ArgumentError("Murthy q(0)-nonunit q_prime degree guard failed"))

    branch_unit = _sl3_local_constant_coefficient(form.q + p_prime, form.var_idx, form.R)
    branch_unit_inverse = _unit_inverse_or_nothing(branch_unit)
    branch_unit_inverse === nothing &&
        _throw_staged_sl3_local_failure("Murthy q(0)-nonunit Bezout child q(0) is not a unit")

    bezout_target = _sl3_local_special_form_target(form.R, form.p, form.q, q_prime, p_prime)
    child_link_target = _sl3_local_special_form_target(
        form.R,
        form.p + q_prime,
        form.q + p_prime,
        q_prime,
        p_prime,
    )
    left_factor = elementary_matrix(3, 2, 1, form.r * p_prime - form.s * q_prime, form.R)
    first_elementary_factor = elementary_matrix(3, 1, 2, -one(form.R), form.R)
    child_entries = _sl3_local_target_entries(child_link_target)
    child_entries === nothing &&
        error("internal Murthy q(0)-nonunit child link target left special form")
    child_form = (;
        family = :murthy_q0_unit,
        R = form.R,
        p = child_entries.p,
        q = child_entries.q,
        r = child_entries.r,
        s = child_entries.s,
        X = form.X,
        target = child_link_target,
        var_idx = form.var_idx,
        murthy_q0_nonunit_witness = nothing,
    )
    child_certificate = _realize_sl3_local_murthy_q0_unit_certificate(child_form)

    form.target == left_factor * bezout_target ||
        error("internal Murthy q(0)-nonunit Bezout reduction first equality failed")
    bezout_target == first_elementary_factor * child_link_target ||
        error("internal Murthy q(0)-nonunit Bezout reduction q0-unit equality failed")
    verify_sl3_local_realization(child_certificate) ||
        error("internal Murthy q(0)-nonunit child certificate verification failed")

    reduction = SL3LocalMurthyQ0NonunitReduction(
        form.target,
        _sl3_local_constant_coefficient(form.p, form.var_idx, form.R),
        _sl3_local_constant_coefficient(form.q, form.var_idx, form.R),
        p_prime,
        q_prime,
        resultant,
        bezout_target,
        child_link_target,
        left_factor,
        first_elementary_factor,
        child_certificate,
        form.X,
        degree_p,
        degree_q,
        degree_p_prime,
        degree_q_prime,
        branch_unit,
        branch_unit_inverse,
        witness_source,
    )
    verify_sl3_local_murthy_q0_nonunit_reduction(reduction) ||
        error("internal Murthy q(0)-nonunit Bezout/resultant reduction verification failed")
    return reduction
end

function _sl3_local_murthy_q0_nonunit_bezout_pair(form)
    if form.murthy_q0_nonunit_witness !== nothing
        witness = form.murthy_q0_nonunit_witness
        hasproperty(witness, :p_prime) ||
            throw(ArgumentError("Murthy q(0)-nonunit witness must provide p_prime"))
        hasproperty(witness, :q_prime) ||
            throw(ArgumentError("Murthy q(0)-nonunit witness must provide q_prime"))
        return (
            _coerce_into_ring(form.R, witness.p_prime, "Murthy q(0)-nonunit p_prime witness"),
            _coerce_into_ring(form.R, witness.q_prime, "Murthy q(0)-nonunit q_prime witness"),
            :supplied_bezout_witness,
        )
    end

    g, a, b = gcdx(form.p, form.q)
    g_inverse = _unit_inverse_or_nothing(g)
    g_inverse === nothing &&
        _throw_staged_sl3_local_failure("Murthy q(0)-nonunit Bezout extraction did not produce a unit gcd")
    return (g_inverse * a, -g_inverse * b, :extracted_bezout_witness)
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

function verify_sl3_local_q_degree_normalization(record)::Bool
    try
        return _sl3_local_q_degree_normalization_verification(record).overall_ok
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end

function verify_sl3_local_murthy_q_unit_reduction(reduction)::Bool
    try
        return _sl3_local_murthy_q0_unit_reduction_verification(reduction).overall_ok
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end

function verify_sl3_local_murthy_q0_nonunit_reduction(reduction)::Bool
    try
        return _sl3_local_murthy_q0_nonunit_reduction_verification(reduction).overall_ok
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end

function _sl3_local_q_degree_normalization_verification(record)
    target = record.target
    normalized_target = record.normalized_target
    target_size_ok = nrows(target) == 3 && ncols(target) == 3
    normalized_size_ok = nrows(normalized_target) == 3 && ncols(normalized_target) == 3
    size_ok = target_size_ok && normalized_size_ok
    R = size_ok ? base_ring(target) : nothing
    normalized_ring_ok = size_ok && _same_base_ring(base_ring(normalized_target), R)
    ordinary_polynomial_ok = normalized_ring_ok && !_is_laurent_polynomial_ring(R)
    quotient_ring_ok = ordinary_polynomial_ok && parent(record.quotient) === R
    remainder_ring_ok = ordinary_polynomial_ok && parent(record.remainder) === R
    variable_ok = ordinary_polynomial_ok &&
        parent(record.selected_variable) === R &&
        record.selected_variable in collect(gens(R))
    correction_size_ok =
        ordinary_polynomial_ok &&
        nrows(record.elementary_correction) == 3 &&
        ncols(record.elementary_correction) == 3
    correction_ring_ok = correction_size_ok && _same_base_ring(base_ring(record.elementary_correction), R)
    entries = ordinary_polynomial_ok ? _sl3_local_target_entries(target) : nothing
    normalized_entries = ordinary_polynomial_ok ? _sl3_local_target_entries(normalized_target) : nothing
    shape_ok = entries !== nothing && normalized_entries !== nothing

    determinant_ok = false
    monic_ok = false
    division_ok = false
    degree_ok = false
    normalized_target_ok = false
    correction_ok = false
    replay_ok = false
    if shape_ok && quotient_ring_ok && remainder_ring_ok && variable_ok && correction_ring_ok
        var_idx = findfirst(isequal(record.selected_variable), collect(gens(R)))
        p = entries.p
        q = entries.q
        r = entries.r
        s = entries.s
        determinant_ok = det(target) == one(R) && det(normalized_target) == one(R)
        monic_ok = _is_monic_in_variable(p, var_idx, R)
        division_ok = q == record.quotient * p + record.remainder
        degree_ok = degree(record.remainder, var_idx) < degree(p, var_idx)
        normalized_target_ok =
            normalized_entries.p == p &&
            normalized_entries.q == record.remainder &&
            normalized_entries.r == r &&
            normalized_entries.s == s - record.quotient * r
        correction_ok =
            record.elementary_correction == elementary_matrix(3, 1, 2, record.quotient, R)
        replay_ok = normalized_target * record.elementary_correction == target
    end

    overall_ok = size_ok && normalized_ring_ok && ordinary_polynomial_ok &&
        quotient_ring_ok && remainder_ring_ok && variable_ok && correction_size_ok &&
        correction_ring_ok && shape_ok && determinant_ok && monic_ok &&
        division_ok && degree_ok && normalized_target_ok && correction_ok && replay_ok
    return (;
        overall_ok,
        size_ok,
        normalized_ring_ok,
        ordinary_polynomial_ok,
        quotient_ring_ok,
        remainder_ring_ok,
        variable_ok,
        correction_size_ok,
        correction_ring_ok,
        shape_ok,
        determinant_ok,
        monic_ok,
        division_ok,
        degree_ok,
        normalized_target_ok,
        correction_ok,
        replay_ok,
    )
end

function _sl3_local_murthy_q0_unit_reduction_verification(reduction)
    target = reduction.target
    eliminated_target = reduction.eliminated_target
    target_size_ok = nrows(target) == 3 && ncols(target) == 3
    eliminated_size_ok = nrows(eliminated_target) == 3 && ncols(eliminated_target) == 3
    size_ok = target_size_ok && eliminated_size_ok
    R = size_ok ? base_ring(target) : nothing
    eliminated_ring_ok = size_ok && _same_base_ring(base_ring(eliminated_target), R)
    ordinary_univariate_ok = eliminated_ring_ok &&
        _sl3_local_supports_murthy_q0_unit_branch(
            R,
            findfirst(isequal(reduction.selected_variable), collect(gens(R))) === nothing ?
                0 :
                findfirst(isequal(reduction.selected_variable), collect(gens(R))),
        )
    entries = ordinary_univariate_ok ? _sl3_local_target_entries(target) : nothing
    eliminated_entries = ordinary_univariate_ok ? _sl3_local_target_entries(eliminated_target) : nothing
    shape_ok = entries !== nothing && eliminated_entries !== nothing
    variable_ok = ordinary_univariate_ok &&
        parent(reduction.selected_variable) === R &&
        reduction.selected_variable in collect(gens(R))
    var_idx = variable_ok ? findfirst(isequal(reduction.selected_variable), collect(gens(R))) : nothing
    scalar_ring_ok = ordinary_univariate_ok &&
        parent(reduction.q0) === R &&
        parent(reduction.q0_inverse) === R &&
        parent(reduction.p0) === R &&
        parent(reduction.right_e21_coefficient) === R &&
        parent(reduction.p_prime) === R
    factor_size_ok = ordinary_univariate_ok &&
        nrows(reduction.elimination_factor) == 3 &&
        ncols(reduction.elimination_factor) == 3 &&
        nrows(reduction.inverse_elimination_factor) == 3 &&
        ncols(reduction.inverse_elimination_factor) == 3
    factor_ring_ok = factor_size_ok &&
        _same_base_ring(base_ring(reduction.elimination_factor), R) &&
        _same_base_ring(base_ring(reduction.inverse_elimination_factor), R)
    locality_witness = reduction.locality_witness
    locality_keys_ok = locality_witness isa NamedTuple &&
        keys(locality_witness) == (:eliminated_p, :eliminated_r, :eliminated_s)
    locality_ring_ok = scalar_ring_ok && locality_keys_ok &&
        parent(locality_witness.eliminated_p) === R &&
        parent(locality_witness.eliminated_r) === R &&
        parent(locality_witness.eliminated_s) === R

    determinant_ok = false
    constants_ok = false
    elimination_factor_ok = false
    eliminated_target_ok = false
    divisibility_ok = false
    degree_ok = false
    locality_witness_ok = false
    split_certificate_ok = false
    final_factors_ok = false
    if shape_ok && variable_ok && scalar_ring_ok && factor_ring_ok && locality_ring_ok
        p = entries.p
        q = entries.q
        q0 = _sl3_local_constant_coefficient(q, var_idx, R)
        p0 = _sl3_local_constant_coefficient(p, var_idx, R)
        lambda = reduction.right_e21_coefficient
        determinant_ok = det(target) == one(R) && det(eliminated_target) == one(R)
        constants_ok =
            reduction.q0 == q0 &&
            reduction.p0 == p0 &&
            reduction.q0 * reduction.q0_inverse == one(R) &&
            lambda == -reduction.q0_inverse * reduction.p0
        elimination_factor_ok =
            reduction.elimination_factor == elementary_matrix(3, 2, 1, lambda, R) &&
            reduction.inverse_elimination_factor == elementary_matrix(3, 2, 1, -lambda, R)
        eliminated_target_ok =
            target * reduction.elimination_factor == eliminated_target &&
            eliminated_entries.p == p + lambda * q &&
            eliminated_entries.q == entries.q &&
            eliminated_entries.r == entries.r + lambda * entries.s &&
            eliminated_entries.s == entries.s
        divisibility_ok = eliminated_entries.p == reduction.selected_variable * reduction.p_prime
        degree_ok =
            reduction.degree_p == degree(p, var_idx) &&
            reduction.degree_p_prime == degree(reduction.p_prime, var_idx) &&
            reduction.degree_p_prime < reduction.degree_p
        locality_witness_ok =
            locality_witness.eliminated_p == eliminated_entries.p &&
            locality_witness.eliminated_r == eliminated_entries.r &&
            locality_witness.eliminated_s == eliminated_entries.s
        split_certificate = reduction.split_certificate
        split_certificate_ok =
            verify_sl3_local_realization(split_certificate) &&
            split_certificate.target == eliminated_target &&
            split_certificate.selected_variable == reduction.selected_variable &&
            split_certificate.branch == :murthy_split_lemma &&
            split_certificate.witness.split.witness.a == reduction.selected_variable &&
            split_certificate.witness.split.witness.a_prime == reduction.p_prime &&
            split_certificate.witness.split.witness.b == eliminated_entries.q &&
            split_certificate.witness.split.first_child_target[1, 1] == reduction.selected_variable &&
            split_certificate.witness.split.second_child_target[1, 1] == reduction.p_prime
        final_factors_ok =
            split_certificate_ok &&
            _sl3_local_factor_product(
                vcat(split_certificate.factors, [reduction.inverse_elimination_factor]),
                R,
            ) == target
    end

    overall_ok = size_ok && eliminated_ring_ok && ordinary_univariate_ok && shape_ok &&
        variable_ok && scalar_ring_ok && factor_size_ok && factor_ring_ok &&
        locality_keys_ok && locality_ring_ok && determinant_ok && constants_ok &&
        elimination_factor_ok && eliminated_target_ok && divisibility_ok &&
        degree_ok && locality_witness_ok && split_certificate_ok && final_factors_ok
    return (;
        overall_ok,
        size_ok,
        eliminated_ring_ok,
        ordinary_univariate_ok,
        shape_ok,
        variable_ok,
        scalar_ring_ok,
        factor_size_ok,
        factor_ring_ok,
        locality_keys_ok,
        locality_ring_ok,
        determinant_ok,
        constants_ok,
        elimination_factor_ok,
        eliminated_target_ok,
        divisibility_ok,
        degree_ok,
        locality_witness_ok,
        split_certificate_ok,
        final_factors_ok,
    )
end

function _sl3_local_murthy_q0_nonunit_reduction_verification(reduction)
    target = reduction.target
    bezout_target = reduction.bezout_target
    child_link_target = reduction.child_link_target
    target_size_ok = nrows(target) == 3 && ncols(target) == 3
    bezout_size_ok = nrows(bezout_target) == 3 && ncols(bezout_target) == 3
    child_link_size_ok = nrows(child_link_target) == 3 && ncols(child_link_target) == 3
    size_ok = target_size_ok && bezout_size_ok && child_link_size_ok
    R = size_ok ? base_ring(target) : nothing
    bezout_ring_ok = size_ok && _same_base_ring(base_ring(bezout_target), R)
    child_link_ring_ok = size_ok && _same_base_ring(base_ring(child_link_target), R)
    selected_var_idx = size_ok ? findfirst(isequal(reduction.selected_variable), collect(gens(R))) : nothing
    ordinary_univariate_ok = child_link_ring_ok &&
        _sl3_local_supports_murthy_q0_unit_branch(R, selected_var_idx === nothing ? 0 : selected_var_idx)
    entries = ordinary_univariate_ok ? _sl3_local_target_entries(target) : nothing
    bezout_entries = ordinary_univariate_ok ? _sl3_local_target_entries(bezout_target) : nothing
    child_entries = ordinary_univariate_ok ? _sl3_local_target_entries(child_link_target) : nothing
    shape_ok = entries !== nothing && bezout_entries !== nothing && child_entries !== nothing
    variable_ok = ordinary_univariate_ok &&
        parent(reduction.selected_variable) === R &&
        reduction.selected_variable in collect(gens(R))
    var_idx = variable_ok ? findfirst(isequal(reduction.selected_variable), collect(gens(R))) : nothing
    scalar_ring_ok = ordinary_univariate_ok &&
        parent(reduction.p0) === R &&
        parent(reduction.q0) === R &&
        parent(reduction.p_prime) === R &&
        parent(reduction.q_prime) === R &&
        parent(reduction.resultant) === R &&
        parent(reduction.branch_unit) === R &&
        parent(reduction.branch_unit_inverse) === R
    factor_size_ok = ordinary_univariate_ok &&
        nrows(reduction.left_factor) == 3 &&
        ncols(reduction.left_factor) == 3 &&
        nrows(reduction.first_elementary_factor) == 3 &&
        ncols(reduction.first_elementary_factor) == 3
    factor_ring_ok = factor_size_ok &&
        _same_base_ring(base_ring(reduction.left_factor), R) &&
        _same_base_ring(base_ring(reduction.first_elementary_factor), R)
    witness_source_ok =
        reduction.witness_source in (:supplied_bezout_witness, :extracted_bezout_witness)

    determinant_ok = false
    constants_ok = false
    bezout_ok = false
    degree_ok = false
    branch_unit_ok = false
    bezout_target_ok = false
    child_link_target_ok = false
    factors_ok = false
    replay_ok = false
    child_certificate_ok = false
    final_factors_ok = false
    if shape_ok && variable_ok && scalar_ring_ok && factor_ring_ok && witness_source_ok
        p = entries.p
        q = entries.q
        r = entries.r
        s = entries.s
        determinant_ok =
            det(target) == one(R) &&
            det(bezout_target) == one(R) &&
            det(child_link_target) == one(R)
        constants_ok =
            reduction.p0 == _sl3_local_constant_coefficient(p, var_idx, R) &&
            reduction.q0 == _sl3_local_constant_coefficient(q, var_idx, R)
        bezout_ok =
            reduction.resultant == reduction.p_prime * p - reduction.q_prime * q &&
            reduction.resultant == one(R)
        degree_ok =
            reduction.degree_p == degree(p, var_idx) &&
            reduction.degree_q == degree(q, var_idx) &&
            reduction.degree_p_prime == degree(reduction.p_prime, var_idx) &&
            reduction.degree_q_prime == degree(reduction.q_prime, var_idx) &&
            reduction.degree_p_prime < reduction.degree_q &&
            reduction.degree_q_prime < reduction.degree_p
        branch_unit_ok =
            reduction.branch_unit == _sl3_local_constant_coefficient(child_entries.q, var_idx, R) &&
            reduction.branch_unit ==
                reduction.q0 + _sl3_local_constant_coefficient(reduction.p_prime, var_idx, R) &&
            reduction.branch_unit * reduction.branch_unit_inverse == one(R)
        bezout_target_ok =
            bezout_entries.p == p &&
            bezout_entries.q == q &&
            bezout_entries.r == reduction.q_prime &&
            bezout_entries.s == reduction.p_prime
        child_link_target_ok =
            child_entries.p == p + reduction.q_prime &&
            child_entries.q == q + reduction.p_prime &&
            child_entries.r == reduction.q_prime &&
            child_entries.s == reduction.p_prime
        factors_ok =
            reduction.left_factor ==
                elementary_matrix(3, 2, 1, r * reduction.p_prime - s * reduction.q_prime, R) &&
            reduction.first_elementary_factor == elementary_matrix(3, 1, 2, -one(R), R)
        replay_ok =
            target == reduction.left_factor * bezout_target &&
            bezout_target == reduction.first_elementary_factor * child_link_target
        child_certificate = reduction.child_certificate
        child_certificate_ok =
            verify_sl3_local_realization(child_certificate) &&
            child_certificate.target == child_link_target &&
            child_certificate.selected_variable == reduction.selected_variable &&
            child_certificate.branch == :murthy_q0_unit
        final_factors_ok =
            child_certificate_ok &&
            _sl3_local_factor_product(
                vcat(
                    [reduction.left_factor, reduction.first_elementary_factor],
                    child_certificate.factors,
                ),
                R,
            ) == target
    end

    overall_ok = size_ok && bezout_ring_ok && child_link_ring_ok && ordinary_univariate_ok &&
        shape_ok && variable_ok && scalar_ring_ok && factor_size_ok && factor_ring_ok &&
        witness_source_ok && determinant_ok && constants_ok && bezout_ok &&
        degree_ok && branch_unit_ok && bezout_target_ok && child_link_target_ok &&
        factors_ok && replay_ok && child_certificate_ok && final_factors_ok
    return (;
        overall_ok,
        size_ok,
        bezout_ring_ok,
        child_link_ring_ok,
        ordinary_univariate_ok,
        shape_ok,
        variable_ok,
        scalar_ring_ok,
        factor_size_ok,
        factor_ring_ok,
        witness_source_ok,
        determinant_ok,
        constants_ok,
        bezout_ok,
        degree_ok,
        branch_unit_ok,
        bezout_target_ok,
        child_link_target_ok,
        factors_ok,
        replay_ok,
        child_certificate_ok,
        final_factors_ok,
    )
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
    elseif certificate.branch == :q_unit
        return [
            elementary_matrix(3, 2, 1, (s - one(R)) * witness.pivot_inverse, R),
            elementary_matrix(3, 1, 2, q, R),
            elementary_matrix(3, 2, 1, (p - one(R)) * witness.pivot_inverse, R),
        ]
    elseif certificate.branch == :murthy_split_lemma
        return vcat(
            witness.split.prefix_factors,
            witness.first_child_certificate.factors,
            witness.split.middle_factors,
            witness.second_child_certificate.factors,
            witness.split.suffix_factors,
        )
    elseif certificate.branch == :murthy_q_degree_normalization
        return [
            witness.normalization.normalized_target,
            witness.normalization.elementary_correction,
        ]
    elseif certificate.branch == :murthy_q0_unit
        if witness.reduction === nothing
            return vcat(
                witness.normalized_certificate.factors,
                [witness.normalization.elementary_correction],
            )
        end

        return vcat(
            witness.reduction.split_certificate.factors,
            [witness.reduction.inverse_elimination_factor],
        )
    elseif certificate.branch == :murthy_q0_nonunit_bezout_resultant
        return vcat(
            [witness.reduction.left_factor, witness.reduction.first_elementary_factor],
            witness.reduction.child_certificate.factors,
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
    elseif certificate.branch == :q_unit
        return witness.pivot == q && witness.pivot * witness.pivot_inverse == one(R)
    elseif certificate.branch == :murthy_split_lemma
        return _sl3_local_split_certificate_witness_ok(certificate)
    elseif certificate.branch == :murthy_q_degree_normalization
        return _sl3_local_q_degree_certificate_witness_ok(certificate)
    elseif certificate.branch == :murthy_q0_unit
        return _sl3_local_murthy_q0_unit_certificate_witness_ok(certificate)
    elseif certificate.branch == :murthy_q0_nonunit_bezout_resultant
        return _sl3_local_murthy_q0_nonunit_certificate_witness_ok(certificate)
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

function _sl3_local_q_degree_certificate_witness_ok(certificate)
    normalization = certificate.witness.normalization
    verify_sl3_local_q_degree_normalization(normalization) || return false
    normalization.target == certificate.target || return false
    normalization.selected_variable == certificate.selected_variable || return false
    return true
end

function _sl3_local_murthy_q0_unit_certificate_witness_ok(certificate)
    witness = certificate.witness
    if witness.reduction === nothing
        witness.normalization === nothing && return false
        witness.normalized_certificate === nothing && return false
        normalization = witness.normalization
        normalized_certificate = witness.normalized_certificate
        verify_sl3_local_q_degree_normalization(normalization) || return false
        normalization.target == certificate.target || return false
        normalization.selected_variable == certificate.selected_variable || return false
        normalized_certificate.target == normalization.normalized_target || return false
        normalized_certificate.selected_variable == certificate.selected_variable || return false
        verify_sl3_local_realization(normalized_certificate) || return false
        return true
    end

    witness.normalization === nothing || return false
    witness.normalized_certificate === nothing || return false
    reduction = witness.reduction
    verify_sl3_local_murthy_q_unit_reduction(reduction) || return false
    reduction.target == certificate.target || return false
    reduction.selected_variable == certificate.selected_variable || return false
    return true
end

function _sl3_local_murthy_q0_nonunit_certificate_witness_ok(certificate)
    reduction = certificate.witness.reduction
    verify_sl3_local_murthy_q0_nonunit_reduction(reduction) || return false
    reduction.target == certificate.target || return false
    reduction.selected_variable == certificate.selected_variable || return false
    return true
end

function _sl3_local_witness_keys_ok(branch, witness)
    witness isa NamedTuple || return false
    if branch in (:open_s_one, :open_p_one)
        return keys(witness) == (:q, :r)
    elseif branch in (:s_unit, :p_unit, :q_unit)
        return keys(witness) == (:pivot, :pivot_inverse)
    elseif branch == :murthy_split_lemma
        return keys(witness) == (:split, :first_child_certificate, :second_child_certificate)
    elseif branch == :murthy_q_degree_normalization
        return keys(witness) == (:normalization,)
    elseif branch == :murthy_q0_unit
        return keys(witness) == (:normalization, :normalized_certificate, :reduction)
    elseif branch == :murthy_q0_nonunit_bezout_resultant
        return keys(witness) == (:reduction,)
    end

    return false
end

function _throw_staged_sl3_local_failure(reason::AbstractString)
    throw(ArgumentError("staged local SL_3 solver failure: $(reason)"))
end

function _sl3_local_divrem_monic_in_variable(q, p, var_idx::Int, R)
    degree_p = degree(p, var_idx)
    degree_p >= 0 || throw(ArgumentError("p must have nonnegative degree in X"))
    _is_monic_in_variable(p, var_idx, R) || throw(ArgumentError("p must be monic in X"))

    quotient = zero(R)
    remainder = q
    X = collect(gens(R))[var_idx]
    while !iszero(remainder) && degree(remainder, var_idx) >= degree_p
        degree_remainder = degree(remainder, var_idx)
        degree_gap = degree_remainder - degree_p
        leading = _sl3_local_coefficient_in_variable_degree(remainder, var_idx, degree_remainder, R)
        term = leading * X^degree_gap
        quotient += term
        remainder -= term * p
    end

    return quotient, remainder
end

function _sl3_local_coefficient_in_variable_degree(value, var_idx::Int, target_degree::Int, R)
    ring_gens = collect(gens(R))
    total = zero(R)
    for (coeff, exponents) in zip(AbstractAlgebra.coefficients(value), AbstractAlgebra.exponent_vectors(value))
        exponents[var_idx] == target_degree || continue
        term = R(coeff)
        for idx in eachindex(ring_gens)
            idx == var_idx && continue
            exponent = exponents[idx]
            exponent == 0 && continue
            term *= ring_gens[idx]^exponent
        end
        total += term
    end

    return total
end

function _sl3_local_constant_coefficient(value, var_idx::Int, R)
    return _sl3_local_coefficient_in_variable_degree(value, var_idx, 0, R)
end

function _is_monic_in_variable(p, var_idx::Int, R)
    iszero(p) && return false

    target_degree = degree(p, var_idx)
    target_degree < 0 && return false

    return _sl3_local_coefficient_in_variable_degree(p, var_idx, target_degree, R) == one(R)
end

function _sl3_local_supports_murthy_q0_unit_branch(R, var_idx::Int)
    var_idx < 1 && return false
    _is_laurent_polynomial_ring(R) && return false
    return length(collect(gens(R))) == 1
end
