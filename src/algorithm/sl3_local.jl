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
    local_factor_replay
    witness_source::Symbol
end

struct SL3LocalMurthyInputContext
    R
    X
    var_idx::Int
    entries::NamedTuple
    target
    determinant
    degree_p::Int
    degree_q::Int
    p0
    q0
    p_monic::Bool
    global_units::NamedTuple
    local_units::NamedTuple
    local_unit_witnesses::NamedTuple
    split_witness
    bezout_witness
end

struct SL3LocalElementaryFactor
    R
    n::Int
    row::Int
    col::Int
    numerator
    denominator
    selected_variable
    local_unit_witness
end

struct SL3LocalElementaryFactorReplay
    target
    factors::Vector{SL3LocalElementaryFactor}
    selected_variable
    mode::Symbol
    denominator_product
    cleared_product
    materialized_factors
end

struct SL3LocalMurthyQUnitLocalReduction
    target
    context::SL3LocalMurthyInputContext
    q0
    q0_inverse
    p0
    right_e21_coefficient
    elimination_factor::SL3LocalElementaryFactor
    inverse_elimination_factor::SL3LocalElementaryFactor
    source_certificate
    split_certificate
    local_factor_replay::SL3LocalElementaryFactorReplay
    selected_variable
    degree_p::Int
    degree_q::Int
end

function Base.:(==)(left::SL3LocalQDegreeNormalization, right::SL3LocalQDegreeNormalization)
    return left.target == right.target &&
        left.quotient == right.quotient &&
        left.remainder == right.remainder &&
        left.normalized_target == right.normalized_target &&
        left.elementary_correction == right.elementary_correction &&
        left.selected_variable == right.selected_variable
end

function Base.:(==)(left::SL3LocalElementaryFactor, right::SL3LocalElementaryFactor)
    return _same_base_ring(left.R, right.R) &&
        left.n == right.n &&
        left.row == right.row &&
        left.col == right.col &&
        left.numerator == right.numerator &&
        left.denominator == right.denominator &&
        left.selected_variable == right.selected_variable &&
        left.local_unit_witness == right.local_unit_witness
end

function sl3_local_murthy_input_context(A, X; witness=nothing, local_unit_witnesses=(;),
        split_witness=nothing, bezout_witness=nothing)
    nrows(A) == 3 || throw(ArgumentError("Murthy local input context requires a 3x3 special-form matrix"))
    ncols(A) == 3 || throw(ArgumentError("Murthy local input context requires a 3x3 special-form matrix"))
    entries = _sl3_local_target_entries(A)
    entries === nothing && throw(ArgumentError("Murthy local input context requires a special-form SL_3 target"))
    return sl3_local_murthy_input_context(
        entries.p,
        entries.q,
        entries.r,
        entries.s,
        X;
        witness,
        local_unit_witnesses,
        split_witness,
        bezout_witness,
    )
end

function sl3_local_murthy_input_context(p, q, r, s, X; witness=nothing,
        local_unit_witnesses=(;), split_witness=nothing, bezout_witness=nothing)
    return _sl3_local_murthy_input_context(
        p,
        q,
        r,
        s,
        X;
        witness,
        local_unit_witnesses,
        split_witness,
        bezout_witness,
    )
end

function verify_sl3_local_murthy_input_context(context)::Bool
    try
        return _sl3_local_murthy_input_context_verification(context).overall_ok
    catch err
        err isa InterruptException && rethrow()
        return false
    end
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

function realize_sl3_local_certificate(context::SL3LocalMurthyInputContext)
    verify_sl3_local_murthy_input_context(context) ||
        throw(ArgumentError("local Murthy realization requires a verified input context"))
    if context.degree_q >= context.degree_p
        return sl3_local_q_degree_normalization_certificate(context)
    end
    if context.global_units.q0 || context.local_units.q0
        if context.global_units.q0 &&
                _sl3_local_supports_murthy_q0_unit_branch(context.R, context.var_idx)
            return realize_sl3_local_certificate(context.target, context.X)
        end
        return _realize_sl3_local_murthy_q0_unit_local_certificate(context)
    end
    if _sl3_local_supports_murthy_q0_unit_branch(context.R, context.var_idx)
        return realize_sl3_local_certificate(
            context.target,
            context.X;
            murthy_q0_nonunit_witness = context.bezout_witness,
        )
    end
    context.bezout_witness === nothing &&
        _throw_staged_sl3_local_failure("Murthy q(0)-nonunit local Bezout/resultant extraction is unsupported")
    return _realize_sl3_local_murthy_q0_nonunit_local_certificate(context)
end

function sl3_local_q_degree_normalization(A, X; check_monic::Bool=true)
    nrows(A) == 3 || throw(ArgumentError("q-degree normalization requires a local SL_3 special-form matrix"))
    ncols(A) == 3 || throw(ArgumentError("q-degree normalization requires a local SL_3 special-form matrix"))
    entries = _sl3_local_target_entries(A)
    entries === nothing && throw(ArgumentError("q-degree normalization requires a local SL_3 special-form matrix"))
    return sl3_local_q_degree_normalization(entries.p, entries.q, entries.r, entries.s, X; check_monic)
end

function sl3_local_q_degree_normalization(context::SL3LocalMurthyInputContext)
    verify_sl3_local_murthy_input_context(context) ||
        throw(ArgumentError("Murthy q-degree normalization requires a verified local input context"))
    context.degree_q >= context.degree_p ||
        throw(ArgumentError("Murthy q-degree normalization context requires deg(q) >= deg(p)"))
    return sl3_local_q_degree_normalization(
        context.entries.p,
        context.entries.q,
        context.entries.r,
        context.entries.s,
        context.X,
    )
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

function sl3_local_q_degree_normalization_certificate(context::SL3LocalMurthyInputContext)
    return sl3_local_q_degree_normalization_certificate(
        sl3_local_q_degree_normalization(context),
    )
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
        _is_monic_in_variable(p, var_idx, R) || _throw_staged_sl3_local_failure("p must be monic in X")
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

function _realize_sl3_local_murthy_q0_unit_local_certificate(context::SL3LocalMurthyInputContext)
    reduction = _sl3_local_murthy_q0_unit_local_reduction(context)
    certificate = SL3LocalRealizationCertificate(
        context.target,
        :murthy_q0_unit,
        reduction.local_factor_replay.factors,
        context.X,
        (; normalization = nothing, normalized_certificate = nothing, reduction),
    )
    verify_sl3_local_realization(certificate) ||
        error("internal local Murthy q(0)-unit certificate verification failed")
    return certificate
end

function _sl3_local_murthy_q0_unit_local_reduction(context::SL3LocalMurthyInputContext)
    model = _sl3_local_fraction_model(context)
    fraction_target = _sl3_local_to_fraction_matrix(context.target, model)
    source_entries = _sl3_local_target_entries(fraction_target)
    source_entries === nothing &&
        error("internal local q(0)-unit replay source left special form")
    source_form = (;
        family = :murthy_q0_unit,
        R = model.S,
        p = source_entries.p,
        q = source_entries.q,
        r = source_entries.r,
        s = source_entries.s,
        X = model.Y,
        target = fraction_target,
        var_idx = 1,
        murthy_q0_nonunit_witness = nothing,
    )
    source_certificate = _realize_sl3_local_murthy_q0_unit_certificate(source_form)
    source_reduction = _sl3_local_source_q0_reduction(source_certificate)

    records = _sl3_local_local_factors_from_fraction_matrices(
        source_certificate.factors,
        context,
    )
    replay = sl3_local_elementary_factor_replay(context.target, records, context.X)
    elimination_factor = first(_sl3_local_local_factors_from_fraction_matrices(
        [source_reduction.elimination_factor],
        context,
    ))
    inverse_elimination_factor = first(_sl3_local_local_factors_from_fraction_matrices(
        [source_reduction.inverse_elimination_factor],
        context,
    ))

    reduction = SL3LocalMurthyQUnitLocalReduction(
        context.target,
        context,
        context.q0,
        source_reduction.q0_inverse,
        context.p0,
        source_reduction.right_e21_coefficient,
        elimination_factor,
        inverse_elimination_factor,
        source_certificate,
        source_reduction.split_certificate,
        replay,
        context.X,
        context.degree_p,
        context.degree_q,
    )
    verify_sl3_local_murthy_q_unit_reduction(reduction) ||
        error("internal local Murthy q(0)-unit reduction verification failed")
    return reduction
end

function _sl3_local_fraction_model(context::SL3LocalMurthyInputContext)
    R = context.R
    ring_gens = collect(gens(R))
    coefficient_indices = [idx for idx in eachindex(ring_gens) if idx != context.var_idx]
    isempty(coefficient_indices) &&
        throw(ArgumentError("unsupported local-unit denominator witness"))
    coefficient_names = [string(ring_gens[idx]) for idx in coefficient_indices]
    C, coefficient_gens = polynomial_ring(base_ring(R), coefficient_names)
    K = fraction_field(C)
    S, (Y,) = polynomial_ring(K, [string(context.X)])
    return (;
        R,
        S,
        C,
        K,
        Y,
        coefficient_indices,
        coefficient_gens = collect(coefficient_gens),
        coefficient_names,
        var_idx = context.var_idx,
    )
end

function _sl3_local_to_fraction_matrix(target, model)
    return matrix(model.S, [
        _sl3_local_to_fraction_polynomial(target[i, j], model)
        for i in 1:nrows(target), j in 1:ncols(target)
    ])
end

function _sl3_local_to_fraction_polynomial(value, model)
    parent(value) === model.R ||
        throw(ArgumentError("fraction model translation requires a polynomial in the context ring"))
    result = zero(model.S)
    for (coefficient, exponents) in
            zip(AbstractAlgebra.coefficients(value), AbstractAlgebra.exponent_vectors(value))
        coefficient_value = model.C(coefficient)
        for (coefficient_position, original_idx) in enumerate(model.coefficient_indices)
            exponent = exponents[original_idx]
            exponent == 0 && continue
            coefficient_value *= model.coefficient_gens[coefficient_position]^exponent
        end
        result += model.S(model.K(coefficient_value)) * model.Y^exponents[model.var_idx]
    end
    return result
end

function _sl3_local_fraction_polynomial_to_ratio(value, context::SL3LocalMurthyInputContext)
    S = parent(value)
    length(collect(gens(S))) == 1 ||
        throw(ArgumentError("local fraction polynomial translation requires a univariate source"))
    coefficients = collect(AbstractAlgebra.coefficients(value))
    exponents = collect(AbstractAlgebra.exponent_vectors(value))
    isempty(coefficients) && return (; numerator = zero(context.R), denominator = one(context.R))

    coefficient_denominators = [denominator(coefficient) for coefficient in coefficients]
    common_denominator = one(parent(first(coefficient_denominators)))
    for coefficient_denominator in coefficient_denominators
        common_denominator *= coefficient_denominator
    end

    numerator_value = zero(context.R)
    for (coefficient, exponent_vector) in zip(coefficients, exponents)
        coefficient_numerator = numerator(coefficient)
        coefficient_denominator = denominator(coefficient)
        multiplier = divexact(common_denominator, coefficient_denominator)
        numerator_value +=
            _sl3_local_coefficient_model_to_ring(
                coefficient_numerator * multiplier,
                context,
                nothing,
            ) * context.X^first(exponent_vector)
    end
    denominator_value =
        _sl3_local_coefficient_model_to_ring(common_denominator, context, nothing)
    return (; numerator = numerator_value, denominator = denominator_value)
end

function _sl3_local_coefficient_model_to_ring(
        value,
        context::SL3LocalMurthyInputContext,
        coefficient_names,
)
    R = context.R
    coefficient_indices = [idx for idx in eachindex(collect(gens(R))) if idx != context.var_idx]
    source_gens = collect(gens(parent(value)))
    if coefficient_names !== nothing && length(coefficient_names) != length(source_gens)
        throw(ArgumentError("coefficient model variable-name mismatch"))
    end
    length(source_gens) == length(coefficient_indices) ||
        throw(ArgumentError("coefficient model variable count mismatch"))

    ring_gens = collect(gens(R))
    result = zero(R)
    for (coefficient, exponents) in
            zip(AbstractAlgebra.coefficients(value), AbstractAlgebra.exponent_vectors(value))
        term = R(coefficient)
        for (coefficient_position, original_idx) in enumerate(coefficient_indices)
            exponent = exponents[coefficient_position]
            exponent == 0 && continue
            term *= ring_gens[original_idx]^exponent
        end
        result += term
    end
    return result
end

function _sl3_local_derive_local_unit_witness(context::SL3LocalMurthyInputContext, unit)
    R = context.R
    unit = _coerce_into_ring(R, unit, "local elementary factor denominator")
    unit == one(R) && return nothing
    hasproperty(context.local_unit_witnesses, :q0) ||
        throw(ArgumentError("unsupported local-unit denominator witness"))
    template = context.local_unit_witnesses.q0
    generators = tuple(template.maximal_ideal_generators...)
    length(generators) == 1 ||
        throw(ArgumentError("unsupported local-unit denominator witness"))
    generator = first(generators)
    generator_idx = findfirst(isequal(generator), collect(gens(R)))
    generator_idx === nothing &&
        throw(ArgumentError("unsupported local-unit denominator witness"))

    residue_unit = _sl3_local_constant_coefficient(unit, generator_idx, R)
    residue_inverse = _unit_inverse_or_nothing(residue_unit)
    residue_inverse === nothing &&
        throw(ArgumentError("unsupported local-unit denominator witness"))
    difference = unit - residue_unit
    residue_difference_coefficient = try
        divexact(difference, generator)
    catch err
        err isa InterruptException && rethrow() # COV_EXCL_LINE
        throw(ArgumentError("unsupported local-unit denominator witness")) # COV_EXCL_LINE
    end
    generator * residue_difference_coefficient == difference ||
        throw(ArgumentError("unsupported local-unit denominator witness"))

    return (;
        context = template.context,
        unit,
        residue_unit,
        residue_inverse,
        maximal_ideal_generators = generators,
        residue_difference_coefficients = (residue_difference_coefficient,),
        global_unit = is_unit(unit),
    )
end

function _sl3_local_local_factor_from_fraction_matrix(
        factor,
        context::SL3LocalMurthyInputContext,
)
    nrows(factor) == 3 || throw(ArgumentError("local fraction elementary factor must be 3x3"))
    ncols(factor) == 3 || throw(ArgumentError("local fraction elementary factor must be 3x3"))
    S = base_ring(factor)
    row = 0
    col = 0
    coefficient = zero(S)
    for i in 1:3, j in 1:3
        if i == j
            factor[i, j] == one(S) ||
                throw(ArgumentError("local fraction factor diagonal is not identity"))
        elseif factor[i, j] != zero(S)
            row == 0 ||
                throw(ArgumentError("local fraction factor has more than one off-diagonal entry"))
            row = i
            col = j
            coefficient = factor[i, j]
        end
    end
    row != 0 || throw(ArgumentError("local fraction identity factor has no elementary data"))

    ratio = _sl3_local_fraction_polynomial_to_ratio(coefficient, context)
    local_unit_witness = _sl3_local_derive_local_unit_witness(context, ratio.denominator)
    return sl3_local_elementary_factor(
        row,
        col,
        ratio.numerator,
        ratio.denominator,
        context.X;
        local_unit_witness,
    )
end

function _sl3_local_local_factors_from_fraction_matrices(
        factors,
        context::SL3LocalMurthyInputContext,
)
    return SL3LocalElementaryFactor[
        _sl3_local_local_factor_from_fraction_matrix(factor, context)
        for factor in factors
        if factor != identity_matrix(base_ring(factor), nrows(factor))
    ]
end

function _sl3_local_source_q0_reduction(source_certificate)
    source_certificate.branch == :murthy_q0_unit ||
        throw(ArgumentError("local q(0)-unit replay requires a q(0)-unit source certificate"))
    witness = source_certificate.witness
    if witness.reduction !== nothing
        return witness.reduction
    end
    if witness.normalized_certificate !== nothing
        return _sl3_local_source_q0_reduction(witness.normalized_certificate)
    end
    throw(ArgumentError("local q(0)-unit replay source certificate has no q(0)-unit reduction"))
end

function _realize_sl3_local_murthy_q0_nonunit_local_certificate(
        context::SL3LocalMurthyInputContext,
)
    reduction = _sl3_local_murthy_q0_nonunit_local_reduction(context)
    certificate = SL3LocalRealizationCertificate(
        context.target,
        :murthy_q0_nonunit_bezout_resultant,
        reduction.local_factor_replay.factors,
        context.X,
        (; reduction),
    )
    verify_sl3_local_realization(certificate) ||
        error("internal local Murthy q(0)-nonunit Bezout/resultant certificate verification failed")
    return certificate
end

function _sl3_local_murthy_q0_nonunit_local_reduction(
        context::SL3LocalMurthyInputContext,
)
    verify_sl3_local_murthy_input_context(context) ||
        throw(ArgumentError("local Murthy q(0)-nonunit reduction requires a verified input context"))
    bezout_data = _sl3_local_murthy_bezout_data(
        context.R,
        context.X,
        context.entries.p,
        context.entries.q,
        context.p0,
        context.q0,
        context.degree_p,
        context.degree_q,
        context.bezout_witness,
    )
    if bezout_data === nothing
        context.bezout_witness === nothing && # COV_EXCL_LINE
            _throw_staged_sl3_local_failure("Murthy q(0)-nonunit local Bezout/resultant extraction is unsupported")
        error("internal local Murthy q(0)-nonunit Bezout/resultant witness did not verify") # COV_EXCL_LINE
    end

    p_prime = bezout_data.p_prime
    q_prime = bezout_data.q_prime
    resultant = bezout_data.resultant
    branch_unit = bezout_data.branch_unit
    branch_unit_inverse = _unit_inverse_or_nothing(branch_unit)

    bezout_target = _sl3_local_special_form_target(
        context.R,
        context.entries.p,
        context.entries.q,
        q_prime,
        p_prime,
    )
    child_link_target = _sl3_local_special_form_target(
        context.R,
        context.entries.p + q_prime,
        context.entries.q + p_prime,
        q_prime,
        p_prime,
    )
    left_factor = elementary_matrix(
        3,
        2,
        1,
        context.entries.r * p_prime - context.entries.s * q_prime,
        context.R,
    )
    first_elementary_factor = elementary_matrix(3, 1, 2, -one(context.R), context.R)
    child_context = _sl3_local_murthy_q0_nonunit_child_context(context, child_link_target)
    child_certificate = realize_sl3_local_certificate(child_context)

    context.target == left_factor * bezout_target ||
        error("internal local Murthy q(0)-nonunit Bezout reduction first equality failed")
    bezout_target == first_elementary_factor * child_link_target ||
        error("internal local Murthy q(0)-nonunit Bezout reduction q0-unit equality failed")
    verify_sl3_local_realization(child_certificate) ||
        error("internal local Murthy q(0)-nonunit child certificate verification failed")

    reduction_without_replay = SL3LocalMurthyQ0NonunitReduction(
        context.target,
        context.p0,
        context.q0,
        p_prime,
        q_prime,
        resultant,
        bezout_target,
        child_link_target,
        left_factor,
        first_elementary_factor,
        child_certificate,
        context.X,
        context.degree_p,
        context.degree_q,
        bezout_data.degree_p_prime,
        bezout_data.degree_q_prime,
        branch_unit,
        branch_unit_inverse,
        nothing,
        bezout_data.source,
    )
    local_factor_replay = sl3_local_elementary_factor_replay(
        context.target,
        _sl3_local_q0_nonunit_local_factor_records(reduction_without_replay),
        context.X,
    )
    reduction = SL3LocalMurthyQ0NonunitReduction(
        context.target,
        context.p0,
        context.q0,
        p_prime,
        q_prime,
        resultant,
        bezout_target,
        child_link_target,
        left_factor,
        first_elementary_factor,
        child_certificate,
        context.X,
        context.degree_p,
        context.degree_q,
        bezout_data.degree_p_prime,
        bezout_data.degree_q_prime,
        branch_unit,
        branch_unit_inverse,
        local_factor_replay,
        bezout_data.source,
    )
    verify_sl3_local_murthy_q0_nonunit_reduction(reduction) ||
        error("internal local Murthy q(0)-nonunit Bezout/resultant reduction verification failed")
    return reduction
end

function _sl3_local_murthy_q0_nonunit_child_context(
        context::SL3LocalMurthyInputContext,
        child_link_target,
)
    child_local_unit_witnesses =
        hasproperty(context.local_unit_witnesses, :branch_unit) ?
        (; q0 = context.local_unit_witnesses.branch_unit) :
        (;)
    return sl3_local_murthy_input_context(
        child_link_target,
        context.X;
        local_unit_witnesses = child_local_unit_witnesses,
    )
end

function _sl3_local_q0_nonunit_local_factor_records(reduction)
    prefix_records = sl3_local_denominator_one_records_from_matrices(
        [reduction.left_factor, reduction.first_elementary_factor],
        reduction.selected_variable,
    )
    child_factors = reduction.child_certificate.factors
    all(factor -> factor isa SL3LocalElementaryFactor, child_factors) ||
        throw(ArgumentError("local q(0)-nonunit replay requires child local elementary factors"))
    child_records = SL3LocalElementaryFactor[factor for factor in child_factors]
    return vcat(prefix_records, child_records)
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
        nothing,
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

function _sl3_local_murthy_input_context(p, q, r, s, X; witness, local_unit_witnesses,
        split_witness, bezout_witness)
    R = parent(p)
    parent(q) === R || throw(ArgumentError("Murthy local input context entries and X must lie in the same polynomial ring"))
    parent(r) === R || throw(ArgumentError("Murthy local input context entries and X must lie in the same polynomial ring"))
    parent(s) === R || throw(ArgumentError("Murthy local input context entries and X must lie in the same polynomial ring"))
    parent(X) === R || throw(ArgumentError("Murthy local input context entries and X must lie in the same polynomial ring"))

    ring_gens = collect(gens(R))
    var_idx = findfirst(isequal(X), ring_gens)
    var_idx === nothing && throw(ArgumentError("Murthy local input context variable must be a ring generator"))
    _is_laurent_polynomial_ring(R) &&
        throw(ArgumentError("Murthy local input context is only supported for ordinary polynomial rings"))

    target = _sl3_local_special_form_target(R, p, q, r, s)
    determinant = det(target)
    determinant == one(R) || throw(ArgumentError("Murthy local input context target must have determinant 1"))
    p_monic = _is_monic_in_variable(p, var_idx, R)
    p_monic || throw(ArgumentError("Murthy local input context p must be monic in X"))

    normalized_witnesses = _sl3_local_murthy_normalize_witness_data(
        witness,
        local_unit_witnesses,
        split_witness,
        bezout_witness,
    )
    p0 = _sl3_local_constant_coefficient(p, var_idx, R)
    q0 = _sl3_local_constant_coefficient(q, var_idx, R)
    degree_p = degree(p, var_idx)
    degree_q = degree(q, var_idx)
    bezout_data = _sl3_local_murthy_bezout_data(
        R,
        X,
        p,
        q,
        p0,
        q0,
        degree_p,
        degree_q,
        normalized_witnesses.bezout_witness,
    )
    global_units = _sl3_local_murthy_global_units(
        p,
        q,
        r,
        s,
        p0,
        q0;
        resultant = _sl3_local_murthy_optional_bezout_field(bezout_data, :resultant),
        branch_unit = _sl3_local_murthy_optional_bezout_field(bezout_data, :branch_unit),
    )
    local_units = _sl3_local_murthy_local_units(
        R,
        X,
        p,
        q,
        p0,
        q0,
        degree_p,
        degree_q,
        normalized_witnesses.local_unit_witnesses,
        normalized_witnesses.bezout_witness,
    )
    _sl3_local_murthy_validate_required_local_evidence(
        R,
        var_idx,
        degree_p,
        degree_q,
        q0,
        global_units,
        local_units,
        normalized_witnesses.bezout_witness,
    )
    _sl3_local_murthy_verify_split_witness(R, normalized_witnesses.split_witness, target)

    context = SL3LocalMurthyInputContext(
        R,
        X,
        var_idx,
        (; p, q, r, s),
        target,
        determinant,
        degree_p,
        degree_q,
        p0,
        q0,
        p_monic,
        global_units,
        local_units,
        normalized_witnesses.local_unit_witnesses,
        normalized_witnesses.split_witness,
        normalized_witnesses.bezout_witness,
    )
    verify_sl3_local_murthy_input_context(context) ||
        error("internal Murthy local input context verification failed")
    return context
end

function _sl3_local_murthy_normalize_witness_data(
        witness,
        local_unit_witnesses,
        split_witness,
        bezout_witness,
)
    local_unit_witnesses isa NamedTuple ||
        throw(ArgumentError("Murthy local input context local-unit witnesses must be a named tuple"))

    extracted_local_unit_witnesses = (;)
    extracted_split_witness = nothing
    extracted_bezout_witness = nothing

    if witness !== nothing
        witness isa NamedTuple ||
            throw(ArgumentError("Murthy local input context witness must be a named tuple"))
        if hasproperty(witness, :local_unit_witness)
            extracted_local_unit_witnesses =
                merge(extracted_local_unit_witnesses, (; q0 = witness.local_unit_witness))
        end
        if hasproperty(witness, :branch_unit_witness)
            extracted_local_unit_witnesses =
                merge(extracted_local_unit_witnesses, (; branch_unit = witness.branch_unit_witness))
        end
        if hasproperty(witness, :resultant_unit_witness)
            extracted_local_unit_witnesses =
                merge(extracted_local_unit_witnesses, (; resultant = witness.resultant_unit_witness))
        end
        if hasproperty(witness, :split)
            extracted_split_witness = witness.split
        elseif _sl3_local_murthy_split_witness_fields_ok(witness)
            extracted_split_witness = witness
        end
        if hasproperty(witness, :p_prime) || hasproperty(witness, :q_prime)
            extracted_bezout_witness = witness
        end
    end

    return (;
        local_unit_witnesses = merge(extracted_local_unit_witnesses, local_unit_witnesses),
        split_witness = split_witness === nothing ? extracted_split_witness : split_witness,
        bezout_witness = bezout_witness === nothing ? extracted_bezout_witness : bezout_witness,
    )
end

function _sl3_local_murthy_global_units(p, q, r, s, p0, q0; resultant=nothing, branch_unit=nothing)
    return (;
        p = is_unit(p),
        q = is_unit(q),
        r = is_unit(r),
        s = is_unit(s),
        p0 = is_unit(p0),
        q0 = is_unit(q0),
        resultant = resultant !== nothing && is_unit(resultant),
        branch_unit = branch_unit !== nothing && is_unit(branch_unit),
    )
end

function _sl3_local_murthy_local_units(
        R,
        X,
        p,
        q,
        p0,
        q0,
        degree_p,
        degree_q,
        local_unit_witnesses,
        bezout_witness,
)
    local_unit_witnesses isa NamedTuple ||
        throw(ArgumentError("Murthy local input context local-unit witnesses must be a named tuple"))
    allowed_keys = (:p0, :q0, :resultant, :branch_unit)
    for key in keys(local_unit_witnesses)
        key in allowed_keys ||
            throw(ArgumentError("Murthy local input context has an unsupported local-unit witness $(key)"))
    end

    bezout_data = _sl3_local_murthy_bezout_data(
        R,
        X,
        p,
        q,
        p0,
        q0,
        degree_p,
        degree_q,
        bezout_witness,
    )
    resultant = _sl3_local_murthy_optional_bezout_field(bezout_data, :resultant)
    branch_unit = _sl3_local_murthy_optional_bezout_field(bezout_data, :branch_unit)

    return (;
        p0 = _sl3_local_murthy_local_unit_status(
            R,
            X,
            local_unit_witnesses,
            :p0,
            p0,
            "p0 local-unit witness",
        ),
        q0 = _sl3_local_murthy_local_unit_status(
            R,
            X,
            local_unit_witnesses,
            :q0,
            q0,
            "q0 local-unit witness",
        ),
        resultant = _sl3_local_murthy_local_unit_status(
            R,
            X,
            local_unit_witnesses,
            :resultant,
            resultant,
            "resultant local-unit witness",
        ),
        branch_unit = _sl3_local_murthy_local_unit_status(
            R,
            X,
            local_unit_witnesses,
            :branch_unit,
            branch_unit,
            "branch-unit local-unit witness",
        ),
    )
end

function _sl3_local_murthy_local_unit_status(
        R,
        X,
        local_unit_witnesses,
        key::Symbol,
        expected_unit,
        label::AbstractString,
)
    expected_unit === nothing && !hasproperty(local_unit_witnesses, key) && return false
    expected_unit === nothing &&
        throw(ArgumentError("Murthy local input context cannot verify $(label) without its expected unit"))
    if hasproperty(local_unit_witnesses, key)
        _sl3_local_murthy_verify_local_unit_witness(
            R,
            X,
            getproperty(local_unit_witnesses, key),
            expected_unit;
            label,
        )
        return true
    end
    return is_unit(expected_unit)
end

function _sl3_local_murthy_validate_required_local_evidence(
        R,
        var_idx,
        degree_p,
        degree_q,
        q0,
        global_units,
        local_units,
        bezout_witness,
)
    degree_q >= degree_p && return true
    (global_units.q0 || local_units.q0) && return true

    if bezout_witness === nothing &&
            !_sl3_local_supports_murthy_q0_unit_branch(R, var_idx) &&
            !is_unit(_sl3_local_constant_coefficient_outside_variable(q0, var_idx, R))
        throw(ArgumentError(
            "Murthy local input context has unsupported local Bezout/resultant extraction without a supplied Bezout witness",
        ))
    end

    if !(global_units.resultant || local_units.resultant)
        throw(ArgumentError("Murthy local input context requires a local-unit witness for q0 or Bezout resultant evidence"))
    end
    if global_units.branch_unit || local_units.branch_unit
        return true
    end

    if bezout_witness === nothing
        throw(ArgumentError("Murthy local input context requires a local-unit witness for q0 or a Bezout witness with branch-unit evidence"))
    end
    throw(ArgumentError("Murthy local input context requires a local-unit witness for q0 or branch_unit"))
end

function _sl3_local_murthy_verify_local_unit_witness(
        R,
        X,
        witness,
        expected_unit;
        label::AbstractString = "local-unit witness",
)
    witness isa NamedTuple ||
        throw(ArgumentError("Murthy local input context $(label) must be a named tuple"))
    required_fields = (
        :context,
        :unit,
        :residue_unit,
        :residue_inverse,
        :maximal_ideal_generators,
        :residue_difference_coefficients,
    )
    for field in required_fields
        hasproperty(witness, field) ||
            throw(ArgumentError("Murthy local input context $(label) is missing $(field)"))
    end

    context = witness.context
    context isa NamedTuple ||
        throw(ArgumentError("Murthy local input context $(label) context must be a named tuple"))
    hasproperty(context, :kind) && context.kind == :localization_at_maximal_ideal ||
        throw(ArgumentError("Murthy local input context $(label) has unsupported local context kind"))
    hasproperty(context, :selected_variable) && context.selected_variable == X ||
        throw(ArgumentError("Murthy local input context $(label) local context variable mismatch"))

    unit = _coerce_into_ring(R, witness.unit, "Murthy local input context $(label) unit")
    residue_unit = _coerce_into_ring(R, witness.residue_unit, "Murthy local input context $(label) residue_unit")
    residue_inverse = _coerce_into_ring(R, witness.residue_inverse, "Murthy local input context $(label) residue_inverse")
    expected = _coerce_into_ring(R, expected_unit, "Murthy local input context $(label) expected unit")
    unit == expected ||
        throw(ArgumentError("Murthy local input context $(label) unit does not match expected value"))

    generators = tuple(witness.maximal_ideal_generators...)
    coefficients = tuple(witness.residue_difference_coefficients...)
    hasproperty(context, :maximal_ideal_generators) &&
            tuple(context.maximal_ideal_generators...) == generators ||
        throw(ArgumentError("Murthy local input context $(label) local context maximal_ideal_generators mismatch"))
    length(generators) == length(coefficients) ||
        throw(ArgumentError("Murthy local input context $(label) generator/coefficient length mismatch"))
    for generator in generators
        parent(generator) === R ||
            throw(ArgumentError("Murthy local input context $(label) generator ring mismatch"))
    end
    for coefficient in coefficients
        parent(coefficient) === R ||
            throw(ArgumentError("Murthy local input context $(label) coefficient ring mismatch"))
    end

    difference = zero(R)
    for (coefficient, generator) in zip(coefficients, generators)
        difference += coefficient * generator
    end
    unit - residue_unit == difference ||
        throw(ArgumentError("Murthy local input context $(label) local residue equation failed"))
    residue_unit * residue_inverse == one(R) ||
        throw(ArgumentError("Murthy local input context $(label) local residue inverse equation failed"))
    if hasproperty(witness, :global_unit)
        is_unit(unit) == witness.global_unit ||
            throw(ArgumentError("Murthy local input context $(label) global unit flag is incorrect"))
    end
    return true
end

function _sl3_local_murthy_bezout_data(
        R,
        X,
        p,
        q,
        p0,
        q0,
        degree_p,
        degree_q,
        bezout_witness,
)
    parent(X) === R ||
        throw(ArgumentError("Murthy local input context Bezout variable must lie in the target ring"))
    var_idx = findfirst(isequal(X), collect(gens(R)))
    var_idx === nothing &&
        throw(ArgumentError("Murthy local input context Bezout variable must be a ring generator"))

    supplied = bezout_witness !== nothing
    if supplied
        bezout_witness isa NamedTuple ||
            throw(ArgumentError("Murthy local input context Bezout witness must be a named tuple"))
        hasproperty(bezout_witness, :p_prime) ||
            throw(ArgumentError("Murthy local input context Bezout witness must provide p_prime"))
        hasproperty(bezout_witness, :q_prime) ||
            throw(ArgumentError("Murthy local input context Bezout witness must provide q_prime"))
        p_prime = _coerce_into_ring(R, bezout_witness.p_prime, "Murthy local input context p_prime witness")
        q_prime = _coerce_into_ring(R, bezout_witness.q_prime, "Murthy local input context q_prime witness")
        resultant = p_prime * p - q_prime * q
        if hasproperty(bezout_witness, :resultant)
            supplied_resultant =
                _coerce_into_ring(R, bezout_witness.resultant, "Murthy local input context resultant witness")
            supplied_resultant == resultant ||
                throw(ArgumentError("Murthy local input context Bezout equality p_prime*p - q_prime*q failed"))
        end
        if hasproperty(bezout_witness, :p0)
            _coerce_into_ring(R, bezout_witness.p0, "Murthy local input context p0 witness") == p0 ||
                throw(ArgumentError("Murthy local input context Bezout p0 witness is incorrect"))
        end
        if hasproperty(bezout_witness, :q0)
            _coerce_into_ring(R, bezout_witness.q0, "Murthy local input context q0 witness") == q0 ||
                throw(ArgumentError("Murthy local input context Bezout q0 witness is incorrect"))
        end
    else
        (degree_q >= degree_p || is_unit(q0)) && return nothing
        _sl3_local_supports_murthy_q0_unit_branch(R, var_idx) || return nothing
        g, a, b = gcdx(p, q)
        g_inverse = _unit_inverse_or_nothing(g)
        g_inverse === nothing && return nothing
        p_prime = g_inverse * a
        q_prime = -g_inverse * b
        resultant = p_prime * p - q_prime * q
        resultant == one(R) || return nothing
    end

    degree_p_prime = degree(p_prime, var_idx)
    degree_q_prime = degree(q_prime, var_idx)
    degree_guards_ok = degree_p_prime < degree_q && degree_q_prime < degree_p
    if !degree_guards_ok
        supplied &&
            throw(ArgumentError("Murthy local input context Bezout degree guards failed"))
        return nothing
    end
    if supplied && hasproperty(bezout_witness, :p_prime_degree)
        bezout_witness.p_prime_degree == degree_p_prime ||
            throw(ArgumentError("Murthy local input context p_prime degree witness is incorrect"))
    end
    if supplied && hasproperty(bezout_witness, :q_prime_degree)
        bezout_witness.q_prime_degree == degree_q_prime ||
            throw(ArgumentError("Murthy local input context q_prime degree witness is incorrect"))
    end

    branch_unit = q0 + _sl3_local_constant_coefficient(p_prime, var_idx, R)
    if supplied && hasproperty(bezout_witness, :branch_unit)
        _coerce_into_ring(R, bezout_witness.branch_unit, "Murthy local input context branch unit witness") ==
                branch_unit ||
            throw(ArgumentError("Murthy local input context branch unit witness is incorrect"))
    end
    if supplied && hasproperty(bezout_witness, :case1_entries)
        case1_entries = bezout_witness.case1_entries
        case1_entries isa NamedTuple ||
            throw(ArgumentError("Murthy local input context case1_entries witness must be a named tuple"))
        for field in (:p, :q, :r, :s)
            hasproperty(case1_entries, field) ||
                throw(ArgumentError("Murthy local input context case1_entries witness is missing $(field)"))
        end
        case1_entries.p == p + q_prime ||
            throw(ArgumentError("Murthy local input context Case 2 reduction p entry is incorrect"))
        case1_entries.q == q + p_prime ||
            throw(ArgumentError("Murthy local input context Case 2 reduction q entry is incorrect"))
        case1_entries.r == q_prime ||
            throw(ArgumentError("Murthy local input context Case 2 reduction r entry is incorrect"))
        case1_entries.s == p_prime ||
            throw(ArgumentError("Murthy local input context Case 2 reduction s entry is incorrect"))
        _sl3_local_constant_coefficient(case1_entries.q, var_idx, R) == branch_unit ||
            throw(ArgumentError("Murthy local input context Case 2 branch constant is incorrect"))
        det(_sl3_local_special_form_target(
            R,
            case1_entries.p,
            case1_entries.q,
            case1_entries.r,
            case1_entries.s,
        )) == one(R) ||
            throw(ArgumentError("Murthy local input context Case 2 target determinant is not one"))
    end

    return (;
        p_prime,
        q_prime,
        resultant,
        degree_p_prime,
        degree_q_prime,
        branch_unit,
        source = supplied ? :supplied_bezout_witness : :extracted_bezout_witness,
    )
end

function _sl3_local_murthy_verify_split_witness(R, split_witness, expected_target=nothing)
    split_witness === nothing && return true
    _sl3_local_murthy_split_witness_fields_ok(split_witness) ||
        throw(ArgumentError("Murthy local input context split witness has invalid fields"))
    _sl3_local_murthy_split_witness_ring_ok(split_witness, R) ||
        throw(ArgumentError("Murthy local input context split witness ring mismatch"))

    replay = sl3_local_split_lemma_replay(
        split_witness.a,
        split_witness.a_prime,
        split_witness.b,
        split_witness.c,
        split_witness.d,
        split_witness.c1,
        split_witness.d1,
        split_witness.c2,
        split_witness.d2;
        split_id = :murthy_input_context_split,
    )
    if expected_target !== nothing
        replay.original_target == expected_target ||
            throw(ArgumentError("Murthy local input context split witness does not reconstruct the target"))
    end
    return true
end

function _sl3_local_murthy_split_witness_fields_ok(split_witness)
    split_witness isa NamedTuple || return false
    for field in (:a, :a_prime, :b, :c, :d, :c1, :d1, :c2, :d2)
        hasproperty(split_witness, field) || return false
    end
    return true
end

function _sl3_local_murthy_split_witness_ring_ok(split_witness, R)
    for field in (:a, :a_prime, :b, :c, :d, :c1, :d1, :c2, :d2)
        parent(getproperty(split_witness, field)) === R || return false
    end
    return true
end

function _sl3_local_murthy_optional_bezout_field(bezout_data, field::Symbol)
    bezout_data === nothing && return nothing
    return getproperty(bezout_data, field)
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

function sl3_local_elementary_factor(row, col, numerator, denominator, X;
        local_unit_witness=nothing, n::Int=3)
    R = parent(X)
    record = SL3LocalElementaryFactor(
        R,
        n,
        Int(row),
        Int(col),
        _coerce_into_ring(R, numerator, "local elementary factor numerator"),
        _coerce_into_ring(R, denominator, "local elementary factor denominator"),
        X,
        local_unit_witness,
    )
    _require_sl3_local_elementary_factor(record)
    return record
end

function verify_sl3_local_elementary_factor(record)::Bool
    try
        _require_sl3_local_elementary_factor(record)
        return true
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end

function sl3_local_materialize_elementary_factor(record)
    _require_sl3_local_elementary_factor(record)
    record.denominator == one(record.R) ||
        throw(ArgumentError("local elementary factor cannot materialize over the ordinary base ring"))
    return elementary_matrix(record.n, record.row, record.col, record.numerator, record.R)
end

function sl3_local_denominator_one_records_from_matrices(factors, X)
    return [
        _sl3_local_denominator_one_record_from_matrix(factor, X)
        for factor in factors
    ]
end

function sl3_local_elementary_factor_replay(target, records, X)
    collected = SL3LocalElementaryFactor[records...]
    replay = _sl3_local_elementary_factor_replay(target, collected, X)
    verify_sl3_local_elementary_factor_replay(replay) ||
        error("internal local elementary factor replay verification failed")
    return replay
end

function verify_sl3_local_elementary_factor_replay(replay)::Bool
    try
        return _sl3_local_elementary_factor_replay_verification(replay).overall_ok
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end

function _require_sl3_local_elementary_factor(record)
    record.n == 3 || throw(ArgumentError("local elementary factor size must be 3"))
    1 <= record.row <= record.n || throw(ArgumentError("local elementary factor row is out of bounds"))
    1 <= record.col <= record.n || throw(ArgumentError("local elementary factor column is out of bounds"))
    record.row != record.col || throw(ArgumentError("local elementary factor row and column must differ"))
    parent(record.selected_variable) === record.R ||
        throw(ArgumentError("local elementary factor variable must lie in the factor ring"))
    record.selected_variable in collect(gens(record.R)) ||
        throw(ArgumentError("local elementary factor variable must be a ring generator"))
    parent(record.numerator) === record.R ||
        throw(ArgumentError("local elementary factor numerator ring mismatch"))
    parent(record.denominator) === record.R ||
        throw(ArgumentError("local elementary factor denominator ring mismatch"))
    iszero(record.denominator) &&
        throw(ArgumentError("local elementary factor denominator must be nonzero"))
    if record.denominator == one(record.R)
        record.local_unit_witness === nothing ||
            _sl3_local_murthy_verify_local_unit_witness(
                record.R,
                record.selected_variable,
                record.local_unit_witness,
                record.denominator;
                label = "local elementary factor denominator witness",
            )
    else
        record.local_unit_witness === nothing &&
            throw(ArgumentError("local elementary factor denominator requires a local-unit witness"))
        _sl3_local_murthy_verify_local_unit_witness(
            record.R,
            record.selected_variable,
            record.local_unit_witness,
            record.denominator;
            label = "local elementary factor denominator witness",
        )
    end
    return record
end

function _sl3_local_cleared_elementary_factor(record)
    _require_sl3_local_elementary_factor(record)
    cleared = record.denominator * identity_matrix(record.R, record.n)
    cleared[record.row, record.col] += record.numerator
    return cleared
end

function _sl3_local_denominator_one_record_from_matrix(factor, X)
    nrows(factor) == 3 || throw(ArgumentError("ordinary elementary factor must be 3x3"))
    ncols(factor) == 3 || throw(ArgumentError("ordinary elementary factor must be 3x3"))
    R = base_ring(factor)
    parent(X) === R || throw(ArgumentError("ordinary elementary factor variable ring mismatch"))
    row = 0
    col = 0
    coefficient = zero(R)
    for i in 1:3, j in 1:3
        if i == j
            factor[i, j] == one(R) || throw(ArgumentError("ordinary factor diagonal is not identity"))
        elseif factor[i, j] != zero(R)
            row == 0 || throw(ArgumentError("ordinary factor has more than one off-diagonal entry"))
            row = i
            col = j
            coefficient = factor[i, j]
        end
    end
    row != 0 || throw(ArgumentError("ordinary identity factor has no elementary row/column data"))
    return sl3_local_elementary_factor(row, col, coefficient, one(R), X)
end

function _sl3_local_elementary_factor_replay(target, records::Vector{SL3LocalElementaryFactor}, X)
    nrows(target) == 3 || throw(ArgumentError("local elementary factor replay target must be 3x3"))
    ncols(target) == 3 || throw(ArgumentError("local elementary factor replay target must be 3x3"))
    R = base_ring(target)
    parent(X) === R || throw(ArgumentError("local elementary factor replay variable ring mismatch"))
    denominator_product = one(R)
    cleared_product = identity_matrix(R, 3)
    all_materializable = true
    materialized = Any[]
    for record in records
        _require_sl3_local_elementary_factor(record)
        record.R === R || throw(ArgumentError("local elementary factor replay ring mismatch"))
        record.selected_variable == X ||
            throw(ArgumentError("local elementary factor replay variable mismatch"))
        denominator_product *= record.denominator
        cleared_product *= _sl3_local_cleared_elementary_factor(record)
        if record.denominator == one(R)
            push!(materialized, sl3_local_materialize_elementary_factor(record))
        else
            all_materializable = false
        end
    end
    mode = all_materializable ? :ordinary : :denominator_cleared
    return SL3LocalElementaryFactorReplay(
        target,
        records,
        X,
        mode,
        denominator_product,
        cleared_product,
        all_materializable ? collect(materialized) : nothing,
    )
end

function _sl3_local_elementary_factor_replay_verification(replay)
    target_size_ok = nrows(replay.target) == 3 && ncols(replay.target) == 3
    cleared_size_ok = nrows(replay.cleared_product) == 3 && ncols(replay.cleared_product) == 3
    size_ok = target_size_ok && cleared_size_ok
    R = size_ok ? base_ring(replay.target) : nothing
    cleared_ring_ok = size_ok && _same_base_ring(base_ring(replay.cleared_product), R)
    variable_ok = size_ok &&
        parent(replay.selected_variable) === R &&
        replay.selected_variable in collect(gens(R))
    factors_ok = replay.factors isa Vector{SL3LocalElementaryFactor} &&
        all(verify_sl3_local_elementary_factor, replay.factors)

    recomputed_ok = false
    denominator_product_ok = false
    cleared_product_ok = false
    mode_ok = false
    materialized_factors_ok = false
    denominator_cleared_ok = false
    ordinary_ok = false
    if cleared_ring_ok && variable_ok && factors_ok
        expected = _sl3_local_elementary_factor_replay(replay.target, replay.factors, replay.selected_variable)
        expected_materialized = expected.materialized_factors
        recomputed_ok = true
        denominator_product_ok = replay.denominator_product == expected.denominator_product
        cleared_product_ok = replay.cleared_product == expected.cleared_product
        mode_ok = replay.mode == expected.mode
        materialized_factors_ok = replay.materialized_factors == expected_materialized

        denominator_cleared_ok =
            expected.cleared_product == replay.denominator_product * replay.target
        ordinary_ok =
            replay.mode == :ordinary &&
            replay.materialized_factors !== nothing &&
            _sl3_local_factor_product(replay.materialized_factors, R) == replay.target
    end

    overall_ok = size_ok && cleared_ring_ok && variable_ok && factors_ok &&
        recomputed_ok && denominator_product_ok && cleared_product_ok && mode_ok &&
        materialized_factors_ok && (
            ordinary_ok ||
            denominator_cleared_ok
        )
    return (;
        overall_ok,
        size_ok,
        cleared_ring_ok,
        variable_ok,
        factors_ok,
        recomputed_ok,
        denominator_product_ok,
        cleared_product_ok,
        mode_ok,
        materialized_factors_ok,
        ordinary_ok,
        denominator_cleared_ok,
    )
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

function verify_sl3_local_murthy_q_unit_reduction(
        reduction::SL3LocalMurthyQUnitLocalReduction,
)::Bool
    try
        return _sl3_local_murthy_q0_unit_local_reduction_verification(reduction).overall_ok
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

function _sl3_local_murthy_input_context_verification(context)
    target = context.target
    size_ok = nrows(target) == 3 && ncols(target) == 3
    R = size_ok ? base_ring(target) : nothing
    ring_ok = size_ok && context.R === R
    ordinary_polynomial_ok = ring_ok && !_is_laurent_polynomial_ring(R)
    entries = ordinary_polynomial_ok ? _sl3_local_target_entries(target) : nothing
    shape_ok = entries !== nothing
    entries_ok = shape_ok && context.entries == entries

    variable_parent_ok = ordinary_polynomial_ok && parent(context.X) === R
    expected_var_idx =
        variable_parent_ok ? findfirst(isequal(context.X), collect(gens(R))) : nothing
    variable_ok = expected_var_idx !== nothing && context.var_idx == expected_var_idx

    determinant_ok = false
    degree_ok = false
    constants_ok = false
    monic_ok = false
    global_units_ok = false
    local_unit_witnesses_ok = context.local_unit_witnesses isa NamedTuple
    local_units_ok = false
    split_witness_ok = false
    bezout_witness_ok = false
    required_evidence_ok = false
    target_replay_ok = false
    if shape_ok && entries_ok && variable_ok && local_unit_witnesses_ok
        p = entries.p
        q = entries.q
        r = entries.r
        s = entries.s
        determinant = det(target)
        degree_p = degree(p, expected_var_idx)
        degree_q = degree(q, expected_var_idx)
        p0 = _sl3_local_constant_coefficient(p, expected_var_idx, R)
        q0 = _sl3_local_constant_coefficient(q, expected_var_idx, R)
        p_monic = _is_monic_in_variable(p, expected_var_idx, R)
        bezout_data = _sl3_local_murthy_bezout_data(
            R,
            context.X,
            p,
            q,
            p0,
            q0,
            degree_p,
            degree_q,
            context.bezout_witness,
        )
        expected_global_units = _sl3_local_murthy_global_units(
            p,
            q,
            r,
            s,
            p0,
            q0;
            resultant = _sl3_local_murthy_optional_bezout_field(bezout_data, :resultant),
            branch_unit = _sl3_local_murthy_optional_bezout_field(bezout_data, :branch_unit),
        )
        expected_local_units = _sl3_local_murthy_local_units(
            R,
            context.X,
            p,
            q,
            p0,
            q0,
            degree_p,
            degree_q,
            context.local_unit_witnesses,
            context.bezout_witness,
        )

        determinant_ok = determinant == one(R) && context.determinant == determinant
        degree_ok = context.degree_p == degree_p && context.degree_q == degree_q
        constants_ok = context.p0 == p0 && context.q0 == q0
        monic_ok = context.p_monic == p_monic && p_monic
        global_units_ok = context.global_units == expected_global_units
        local_units_ok = context.local_units == expected_local_units
        replayed_target = _sl3_local_special_form_target(R, p, q, r, s)
        split_witness_ok = _sl3_local_murthy_verify_split_witness(
            R,
            context.split_witness,
            replayed_target,
        )
        bezout_witness_ok = context.bezout_witness === nothing || bezout_data !== nothing
        required_evidence_ok = _sl3_local_murthy_validate_required_local_evidence(
            R,
            expected_var_idx,
            degree_p,
            degree_q,
            q0,
            expected_global_units,
            expected_local_units,
            context.bezout_witness,
        )
        target_replay_ok = replayed_target == target
    end

    overall_ok = size_ok && ring_ok && ordinary_polynomial_ok && shape_ok &&
        entries_ok && variable_parent_ok && variable_ok && determinant_ok &&
        degree_ok && constants_ok && monic_ok && global_units_ok &&
        local_unit_witnesses_ok && local_units_ok && split_witness_ok &&
        bezout_witness_ok && required_evidence_ok && target_replay_ok
    return (;
        overall_ok,
        size_ok,
        ring_ok,
        ordinary_polynomial_ok,
        shape_ok,
        entries_ok,
        variable_parent_ok,
        variable_ok,
        determinant_ok,
        degree_ok,
        constants_ok,
        monic_ok,
        global_units_ok,
        local_unit_witnesses_ok,
        local_units_ok,
        split_witness_ok,
        bezout_witness_ok,
        required_evidence_ok,
        target_replay_ok,
    )
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

function _sl3_local_murthy_q0_unit_local_reduction_verification(
        reduction::SL3LocalMurthyQUnitLocalReduction,
)
    context = reduction.context
    context_ok = verify_sl3_local_murthy_input_context(context)
    target_ok = context_ok && reduction.target == context.target
    variable_ok = context_ok && reduction.selected_variable == context.X
    scalar_ok = context_ok &&
        reduction.q0 == context.q0 &&
        reduction.p0 == context.p0 &&
        reduction.degree_p == context.degree_p &&
        reduction.degree_q == context.degree_q &&
        context.degree_q < context.degree_p &&
        (context.global_units.q0 || context.local_units.q0)

    source_certificate = reduction.source_certificate
    source_certificate_ok = false
    source_target_ok = false
    source_reduction = nothing
    source_reduction_ok = false
    split_certificate_ok = false
    source_scalars_ok = false
    local_factor_replay_ok = reduction.local_factor_replay === nothing
    local_factors_ok = false
    elimination_factor_ok = false
    inverse_elimination_factor_ok = false
    if context_ok
        source_certificate_ok =
            verify_sl3_local_realization(source_certificate) &&
            source_certificate.branch == :murthy_q0_unit
        source_target_ok =
            source_certificate_ok &&
            _sl3_local_fraction_matrix_matches_target(
                source_certificate.target,
                context.target,
                context,
            )
        if source_target_ok
            source_reduction = _sl3_local_source_q0_reduction(source_certificate)
            source_reduction_ok = true
            split_certificate_ok =
                reduction.split_certificate == source_reduction.split_certificate &&
                verify_sl3_local_realization(reduction.split_certificate)
            source_scalars_ok =
                reduction.q0_inverse == source_reduction.q0_inverse &&
                reduction.right_e21_coefficient == source_reduction.right_e21_coefficient
            expected_records = _sl3_local_local_factors_from_fraction_matrices(
                source_certificate.factors,
                context,
            )
            local_factor_replay_ok =
                reduction.local_factor_replay.target == context.target &&
                reduction.local_factor_replay.selected_variable == context.X &&
                reduction.local_factor_replay.factors == expected_records &&
                verify_sl3_local_elementary_factor_replay(reduction.local_factor_replay)
            local_factors_ok =
                local_factor_replay_ok &&
                reduction.local_factor_replay.factors isa Vector{SL3LocalElementaryFactor}
            expected_elimination_factor = first(_sl3_local_local_factors_from_fraction_matrices(
                [source_reduction.elimination_factor],
                context,
            ))
            expected_inverse_elimination_factor =
                first(_sl3_local_local_factors_from_fraction_matrices(
                    [source_reduction.inverse_elimination_factor],
                    context,
                ))
            elimination_factor_ok = reduction.elimination_factor == expected_elimination_factor
            inverse_elimination_factor_ok =
                reduction.inverse_elimination_factor == expected_inverse_elimination_factor
        end
    end

    overall_ok = context_ok && target_ok && variable_ok && scalar_ok &&
        source_certificate_ok && source_target_ok && source_reduction_ok &&
        split_certificate_ok && source_scalars_ok && local_factor_replay_ok &&
        local_factors_ok && elimination_factor_ok && inverse_elimination_factor_ok
    return (;
        overall_ok,
        context_ok,
        target_ok,
        variable_ok,
        scalar_ok,
        source_certificate_ok,
        source_target_ok,
        source_reduction_ok,
        split_certificate_ok,
        source_scalars_ok,
        local_factor_replay_ok,
        local_factors_ok,
        elimination_factor_ok,
        inverse_elimination_factor_ok,
    )
end

function _sl3_local_fraction_matrix_matches_target(source_target, target, context)
    nrows(source_target) == nrows(target) || return false
    ncols(source_target) == ncols(target) || return false
    for i in 1:nrows(target), j in 1:ncols(target)
        ratio = _sl3_local_fraction_polynomial_to_ratio(source_target[i, j], context)
        ratio.denominator == one(context.R) || return false
        ratio.numerator == target[i, j] || return false
    end
    return true
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
    ordinary_polynomial_ok = child_link_ring_ok && !_is_laurent_polynomial_ring(R)
    entries = ordinary_polynomial_ok ? _sl3_local_target_entries(target) : nothing
    bezout_entries = ordinary_polynomial_ok ? _sl3_local_target_entries(bezout_target) : nothing
    child_entries = ordinary_polynomial_ok ? _sl3_local_target_entries(child_link_target) : nothing
    shape_ok = entries !== nothing && bezout_entries !== nothing && child_entries !== nothing
    variable_ok = ordinary_polynomial_ok &&
        parent(reduction.selected_variable) === R &&
        reduction.selected_variable in collect(gens(R))
    var_idx = variable_ok ? findfirst(isequal(reduction.selected_variable), collect(gens(R))) : nothing
    scalar_ring_ok = ordinary_polynomial_ok &&
        parent(reduction.p0) === R &&
        parent(reduction.q0) === R &&
        parent(reduction.p_prime) === R &&
        parent(reduction.q_prime) === R &&
        parent(reduction.resultant) === R &&
        parent(reduction.branch_unit) === R &&
        (reduction.branch_unit_inverse === nothing || parent(reduction.branch_unit_inverse) === R)
    factor_size_ok = ordinary_polynomial_ok &&
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
    local_factor_replay_ok = reduction.local_factor_replay === nothing
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
        child_local_q0_evidence_ok =
            child_certificate_ok &&
            child_certificate.witness.reduction isa SL3LocalMurthyQUnitLocalReduction &&
            child_certificate.witness.reduction.context.target == child_link_target &&
            child_certificate.witness.reduction.context.q0 == reduction.branch_unit &&
            child_certificate.witness.reduction.context.local_units.q0
        branch_unit_ok =
            reduction.branch_unit == _sl3_local_constant_coefficient(child_entries.q, var_idx, R) &&
            reduction.branch_unit ==
                reduction.q0 + _sl3_local_constant_coefficient(reduction.p_prime, var_idx, R) &&
            (
                (
                    reduction.branch_unit_inverse !== nothing &&
                    reduction.branch_unit * reduction.branch_unit_inverse == one(R)
                ) ||
                (
                    reduction.branch_unit_inverse === nothing &&
                    child_local_q0_evidence_ok
                )
            )
        if reduction.local_factor_replay === nothing
            final_factors_ok =
                child_certificate_ok &&
                _sl3_local_factor_product(
                    vcat(
                        [reduction.left_factor, reduction.first_elementary_factor],
                        child_certificate.factors,
                    ),
                    R,
                ) == target
        else
            expected_records = _sl3_local_q0_nonunit_local_factor_records(reduction)
            local_factor_replay_ok =
                reduction.local_factor_replay.target == target &&
                reduction.local_factor_replay.selected_variable == reduction.selected_variable &&
                reduction.local_factor_replay.factors == expected_records &&
                verify_sl3_local_elementary_factor_replay(reduction.local_factor_replay)
            final_factors_ok =
                child_certificate_ok &&
                local_factor_replay_ok &&
                reduction.local_factor_replay.factors isa Vector{SL3LocalElementaryFactor}
        end
    end

    overall_ok = size_ok && bezout_ring_ok && child_link_ring_ok && ordinary_polynomial_ok &&
        shape_ok && variable_ok && scalar_ring_ok && factor_size_ok && factor_ring_ok &&
        witness_source_ok && determinant_ok && constants_ok && bezout_ok &&
        degree_ok && branch_unit_ok && bezout_target_ok && child_link_target_ok &&
        factors_ok && replay_ok && child_certificate_ok && local_factor_replay_ok &&
        final_factors_ok
    return (;
        overall_ok,
        size_ok,
        bezout_ring_ok,
        child_link_ring_ok,
        ordinary_polynomial_ok,
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
        local_factor_replay_ok,
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
    local_replay =
        witness_ok &&
        certificate.branch == :murthy_q0_unit &&
        certificate.witness.reduction isa SL3LocalMurthyQUnitLocalReduction ?
        certificate.witness.reduction.local_factor_replay :
        witness_ok &&
                certificate.branch == :murthy_q0_nonunit_bezout_resultant &&
                hasproperty(certificate.witness.reduction, :local_factor_replay) ?
            certificate.witness.reduction.local_factor_replay :
            nothing
    factors_ok =
        local_replay === nothing ?
        size_ok && verify_factorization(target, certificate.factors) :
        verify_sl3_local_elementary_factor_replay(local_replay) &&
            local_replay.target == target &&
            local_replay.factors == certificate.factors
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
        if witness.reduction isa SL3LocalMurthyQUnitLocalReduction
            return witness.reduction.local_factor_replay.factors
        end

        return vcat(
            witness.reduction.split_certificate.factors,
            [witness.reduction.inverse_elimination_factor],
        )
    elseif certificate.branch == :murthy_q0_nonunit_bezout_resultant
        if witness.reduction.local_factor_replay !== nothing
            return witness.reduction.local_factor_replay.factors
        end
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

function _sl3_local_constant_coefficient_outside_variable(value, var_idx::Int, R)
    coefficient = value
    for idx in eachindex(collect(gens(R)))
        idx == var_idx && continue
        coefficient = _sl3_local_constant_coefficient(coefficient, idx, R)
    end
    return coefficient
end

function _sl3_local_monicity_witness(p, var_idx::Int, R)
    variable_count = length(collect(gens(R)))
    if var_idx < 1 || var_idx > variable_count || iszero(p)
        return (;
            variable = var_idx >= 1 && var_idx <= variable_count ? collect(gens(R))[var_idx] : nothing,
            variable_index = var_idx,
            degree = -1,
            leading_coefficient = zero(R),
            is_monic = false,
        )
    end

    target_degree = degree(p, var_idx)
    leading = target_degree < 0 ?
        zero(R) :
        _sl3_local_coefficient_in_variable_degree(p, var_idx, target_degree, R)
    return (;
        variable = collect(gens(R))[var_idx],
        variable_index = var_idx,
        degree = target_degree,
        leading_coefficient = leading,
        is_monic = leading == one(R),
    )
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
