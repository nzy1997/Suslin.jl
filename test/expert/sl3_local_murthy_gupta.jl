using Test
using Suslin
using Oscar

const SL3_MURTHY_GUPTA_ACCEPTANCE_FIXTURE_PATH =
    joinpath(@__DIR__, "..", "fixtures", "sl3_murthy_gupta_cases.jl")

function _mg_acceptance_product(factors, R)
    product = identity_matrix(R, 3)
    for factor in factors
        product *= factor
    end
    return product
end

function _mg_acceptance_degree(value, X)
    var_idx = findfirst(isequal(X), collect(gens(parent(value))))
    return degree(value, var_idx)
end

function _mg_acceptance_special_form_target(R, p, q, r, s)
    return matrix(R, [
        p q zero(R);
        r s zero(R);
        zero(R) zero(R) one(R)
    ])
end

function _mg_acceptance_captured_error(f)
    try
        f()
        return nothing
    catch err
        return err
    end
end

function _mg_acceptance_as_namedtuple(value)
    names = propertynames(value)
    return NamedTuple{names}(Tuple(getproperty(value, name) for name in names))
end

function _mg_acceptance_local_unit_witness(unit, residue_unit, residue_inverse, generator, coefficient, X)
    return (;
        context = (;
            kind = :localization_at_maximal_ideal,
            selected_variable = X,
            maximal_ideal_generators = (generator,),
        ),
        unit,
        residue_unit,
        residue_inverse,
        maximal_ideal_generators = (generator,),
        residue_difference_coefficients = (coefficient,),
        global_unit = is_unit(unit),
    )
end

function _mg_acceptance_assert_local_replay(certificate)
    @test Suslin.verify_sl3_local_realization(certificate)
    reduction = certificate.witness.reduction
    replay = reduction.local_factor_replay
    @test replay.target == certificate.target
    @test replay.factors == certificate.factors
    @test replay.mode == :denominator_cleared
    @test replay.materialized_factors === nothing
    @test all(factor -> factor isa Suslin.SL3LocalElementaryFactor, replay.factors)
    @test Suslin.verify_sl3_local_elementary_factor_replay(replay)

    denominator_index = findfirst(factor -> factor.denominator != one(factor.R), replay.factors)
    @test denominator_index !== nothing
    if denominator_index !== nothing
        @test_throws ArgumentError Suslin.sl3_local_materialize_elementary_factor(
            replay.factors[denominator_index],
        )
    end
    return replay
end

function _mg_acceptance_assert_elementary_sequence(target, factors)
    R = base_ring(target)
    @test !isempty(factors)
    for factor in factors
        @test nrows(factor) == 3
        @test ncols(factor) == 3
        @test Suslin._same_base_ring(base_ring(factor), R)
        nonzero_offdiagonal = 0
        for i in 1:3, j in 1:3
            if i == j
                @test factor[i, j] == one(R)
            elseif factor[i, j] != zero(R)
                nonzero_offdiagonal += 1
            end
        end
        @test nonzero_offdiagonal <= 1
    end
    @test _mg_acceptance_product(factors, R) == target
    @test Suslin.verify_factorization(target, factors)
end

function _mg_acceptance_pre_murthy_open_or_unit_pivot(entry)
    R = entry.ring.object
    p, q, r, s = entry.entries.p, entry.entries.q, entry.entries.r, entry.entries.s
    if s == one(R) && p == one(R) + q * r
        return :open_s_one
    elseif p == one(R) && s == one(R) + q * r
        return :open_p_one
    elseif is_unit(s)
        return :s_unit
    elseif is_unit(p)
        return :p_unit
    end
    throw(ArgumentError("pre-#61 open-slice/unit-pivot solver does not support this target"))
end

function _mg_acceptance_assert_supported_case(entry; kwargs...)
    p, q, r, s = entry.entries.p, entry.entries.q, entry.entries.r, entry.entries.s
    X = entry.variable
    target = entry.target
    R = entry.ring.object

    @test target == _mg_acceptance_special_form_target(R, p, q, r, s)
    @test det(target) == one(R)
    @test Suslin._is_monic_in_variable(p, findfirst(isequal(X), collect(gens(R))), R)
    @test !is_unit(p)
    @test !is_unit(s)

    certificate_from_matrix = Suslin.realize_sl3_local_certificate(target, X; kwargs...)
    certificate_from_entries = Suslin.realize_sl3_local_certificate(p, q, r, s, X; kwargs...)
    @test certificate_from_matrix.target == target
    @test certificate_from_entries.target == target
    @test certificate_from_matrix.branch == certificate_from_entries.branch
    @test Suslin.verify_sl3_local_realization(certificate_from_matrix)
    @test Suslin.verify_sl3_local_realization(certificate_from_entries)
    _mg_acceptance_assert_elementary_sequence(target, certificate_from_matrix.factors)

    factors_from_matrix = Suslin.realize_sl3_local(target, X; kwargs...)
    factors_from_entries = Suslin.realize_sl3_local(p, q, r, s, X; kwargs...)
    @test factors_from_matrix == certificate_from_matrix.factors
    @test factors_from_entries == certificate_from_entries.factors
    _mg_acceptance_assert_elementary_sequence(target, factors_from_matrix)
    return certificate_from_matrix
end

function _mg_acceptance_assert_q0_unit_certificate(certificate)
    @test certificate.branch == :murthy_q0_unit
    @test Suslin.verify_sl3_local_realization(certificate)

    if certificate.witness.normalization !== nothing
        normalization = certificate.witness.normalization
        @test Suslin.verify_sl3_local_q_degree_normalization(normalization)
        @test normalization.target == certificate.target
        @test normalization.selected_variable == certificate.selected_variable
        @test certificate.witness.normalized_certificate !== nothing
        @test certificate.witness.reduction === nothing
        _mg_acceptance_assert_q0_unit_certificate(certificate.witness.normalized_certificate)
        return nothing
    end

    reduction = certificate.witness.reduction
    @test reduction !== nothing
    @test Suslin.verify_sl3_local_murthy_q_unit_reduction(reduction)
    @test reduction.target == certificate.target
    @test reduction.selected_variable == certificate.selected_variable

    if reduction isa Suslin.SL3LocalMurthyQUnitLocalReduction
        @test !isempty(certificate.factors)
        @test all(factor -> factor isa Suslin.SL3LocalElementaryFactor, certificate.factors)
        @test Suslin.verify_sl3_local_realization(reduction.source_certificate)
        @test reduction.split_certificate.branch == :murthy_split_lemma
        @test Suslin.verify_sl3_local_realization(reduction.split_certificate)
        @test Suslin.verify_sl3_local_realization(
            reduction.split_certificate.witness.first_child_certificate,
        )
        @test Suslin.verify_sl3_local_realization(
            reduction.split_certificate.witness.second_child_certificate,
        )
        _mg_acceptance_assert_local_replay(certificate)
        return nothing
    end

    @test reduction.q0 * reduction.q0_inverse == one(base_ring(certificate.target))
    @test _mg_acceptance_degree(reduction.p_prime, reduction.selected_variable) < reduction.degree_p
    @test reduction.split_certificate.branch == :murthy_split_lemma
    @test Suslin.verify_sl3_local_realization(reduction.split_certificate)
    @test Suslin.verify_sl3_local_split_lemma_replay(reduction.split_certificate.witness.split)
    @test reduction.split_certificate.witness.split.split_id == :murthy_q0_unit_split
    return nothing
end

function _mg_acceptance_assert_resultant_certificate(certificate; expected_source::Symbol)
    @test certificate.branch == :murthy_q0_nonunit_bezout_resultant
    @test Suslin.verify_sl3_local_realization(certificate)
    reduction = certificate.witness.reduction
    @test reduction.witness_source == expected_source
    @test Suslin.verify_sl3_local_murthy_q0_nonunit_reduction(reduction)

    local_factor_replay = try
        getproperty(reduction, :local_factor_replay)
    catch err
        err isa FieldError || rethrow()
        nothing
    end

    @test reduction.resultant == one(base_ring(certificate.target))
    @test reduction.p_prime * certificate.target[1, 1] - reduction.q_prime * certificate.target[1, 2] ==
        one(base_ring(certificate.target))
    @test reduction.target == reduction.left_factor * reduction.bezout_target
    @test reduction.bezout_target == reduction.first_elementary_factor * reduction.child_link_target
    @test _mg_acceptance_degree(reduction.p_prime, reduction.selected_variable) < reduction.degree_q
    @test _mg_acceptance_degree(reduction.q_prime, reduction.selected_variable) < reduction.degree_p
    _mg_acceptance_assert_q0_unit_certificate(reduction.child_certificate)

    if local_factor_replay === nothing
        @test reduction.branch_unit * reduction.branch_unit_inverse == one(base_ring(certificate.target))
        return nothing
    end

    child_reduction = reduction.child_certificate.witness.reduction
    @test child_reduction isa Suslin.SL3LocalMurthyQUnitLocalReduction
    @test child_reduction.context.q0 == reduction.branch_unit
    @test child_reduction.context.local_units.q0
    _mg_acceptance_assert_local_replay(certificate)
    return nothing
end

@testset "Issue 61 Murthy-Gupta local SL3 acceptance" begin
    include(SL3_MURTHY_GUPTA_ACCEPTANCE_FIXTURE_PATH)
    catalog = SL3MurthyGuptaFixtureCatalog.catalog()
    by_id = Dict(entry.id => entry for entry in catalog.cases)

    normalization_entry = by_id["mg-q-degree-normalization"]
    q0_unit_entry = by_id["mg-q0-unit-recursion"]
    supplied_entry = by_id["mg-q0-nonunit-normalized-bezout-resultant"]
    extracted_entry = by_id["mg-q0-nonunit-extracted-bezout-resultant"]

    acceptance_entries = (normalization_entry, q0_unit_entry, supplied_entry, extracted_entry)
    @test length(acceptance_entries) >= 3
    @test_throws ArgumentError _mg_acceptance_pre_murthy_open_or_unit_pivot(normalization_entry)
    @test_throws ArgumentError _mg_acceptance_pre_murthy_open_or_unit_pivot(supplied_entry)

    normalization_certificate = _mg_acceptance_assert_supported_case(normalization_entry)
    @test normalization_certificate.branch == :murthy_q0_unit
    @test normalization_certificate.witness.normalization !== nothing
    _mg_acceptance_assert_q0_unit_certificate(normalization_certificate)

    q0_unit_certificate = _mg_acceptance_assert_supported_case(q0_unit_entry)
    @test q0_unit_certificate.witness.normalization === nothing
    _mg_acceptance_assert_q0_unit_certificate(q0_unit_certificate)

    supplied_certificate = _mg_acceptance_assert_supported_case(
        supplied_entry;
        murthy_q0_nonunit_witness = first(supplied_entry.witnesses),
    )
    _mg_acceptance_assert_resultant_certificate(
        supplied_certificate;
        expected_source = :supplied_bezout_witness,
    )

    extracted_certificate = _mg_acceptance_assert_supported_case(extracted_entry)
    _mg_acceptance_assert_resultant_certificate(
        extracted_certificate;
        expected_source = :extracted_bezout_witness,
    )
end

@testset "Issue 182 Murthy local SL3 closeout acceptance" begin
    R, (X,) = Oscar.polynomial_ring(QQ, ["X"])

    p_qdegree = X^2 + X + 1
    q_qdegree = X * p_qdegree + one(R)
    qdegree_target = _mg_acceptance_special_form_target(
        R,
        p_qdegree,
        q_qdegree,
        -one(R),
        -X,
    )
    qdegree_cert = Suslin.realize_sl3_local_certificate(qdegree_target, X)
    @test qdegree_cert.branch == :murthy_q0_unit
    @test qdegree_cert.witness.normalization !== nothing
    _mg_acceptance_assert_q0_unit_certificate(qdegree_cert)
    _mg_acceptance_assert_elementary_sequence(
        qdegree_target,
        Suslin.realize_sl3_local(qdegree_target, X),
    )

    p_q0_unit = X^2 + X + 1
    q0_unit_target = _mg_acceptance_special_form_target(
        R,
        p_q0_unit,
        one(R),
        p_q0_unit^2 - one(R),
        p_q0_unit,
    )
    q0_unit_cert = Suslin.realize_sl3_local_certificate(q0_unit_target, X)
    @test q0_unit_cert.branch == :murthy_q0_unit
    @test q0_unit_cert.witness.normalization === nothing
    _mg_acceptance_assert_q0_unit_certificate(q0_unit_cert)
    _mg_acceptance_assert_elementary_sequence(
        q0_unit_target,
        Suslin.realize_sl3_local(q0_unit_target, X),
    )

    p_nonunit = X^2 + X + 1
    q_nonunit = X
    nonunit_witness = (;
        p0 = one(R),
        q0 = zero(R),
        p_prime = one(R),
        q_prime = X + 1,
        resultant = one(R),
        p_prime_degree = 0,
        q_prime_degree = 1,
        branch_unit = one(R),
        case1_entries = (;
            p = p_nonunit + X + 1,
            q = q_nonunit + one(R),
            r = X + 1,
            s = one(R),
        ),
    )
    nonunit_target = _mg_acceptance_special_form_target(
        R,
        p_nonunit,
        q_nonunit,
        X + 1 + p_nonunit,
        X + 1,
    )
    supplied_nonunit_cert = Suslin.realize_sl3_local_certificate(
        nonunit_target,
        X;
        murthy_q0_nonunit_witness = nonunit_witness,
    )
    _mg_acceptance_assert_resultant_certificate(
        supplied_nonunit_cert;
        expected_source = :supplied_bezout_witness,
    )
    _mg_acceptance_assert_elementary_sequence(
        nonunit_target,
        Suslin.realize_sl3_local(
            nonunit_target,
            X;
            murthy_q0_nonunit_witness = nonunit_witness,
        ),
    )

    RU, (u, Y) = Oscar.polynomial_ring(QQ, ["u", "X"])

    local_q0_unit = Y + u + 2
    local_q0_unit_p = Y * local_q0_unit + one(RU)
    local_q0_unit_target = _mg_acceptance_special_form_target(
        RU,
        local_q0_unit_p,
        local_q0_unit,
        Y + local_q0_unit_p * Y,
        local_q0_unit_p,
    )
    local_q0_unit_witness = (;
        p0 = one(RU),
        q0 = u + 2,
        local_unit_witness = _mg_acceptance_local_unit_witness(
            u + 2,
            RU(2),
            RU(QQ(1) // QQ(2)),
            u,
            one(RU),
            Y,
        ),
        formal_right_e21_coefficient = "-1/(u + 2)",
    )
    local_q0_unit_context = Suslin.sl3_local_murthy_input_context(
        local_q0_unit_target,
        Y;
        witness = local_q0_unit_witness,
    )
    local_q0_unit_cert = Suslin.realize_sl3_local_certificate(local_q0_unit_context)
    @test local_q0_unit_cert.branch == :murthy_q0_unit
    _mg_acceptance_assert_local_replay(local_q0_unit_cert)

    local_nonunit_q = Y + 2 * u
    local_nonunit_p = Y * local_nonunit_q + one(RU)
    local_nonunit_target = _mg_acceptance_special_form_target(
        RU,
        local_nonunit_p,
        local_nonunit_q,
        Y + local_nonunit_p * Y,
        local_nonunit_p,
    )
    local_nonunit_witness = (;
        p0 = one(RU),
        q0 = 2 * u,
        p_prime = one(RU),
        q_prime = Y,
        resultant = one(RU),
        p_prime_degree = 0,
        q_prime_degree = 1,
        branch_unit = one(RU) + 2 * u,
        branch_unit_witness = _mg_acceptance_local_unit_witness(
            one(RU) + 2 * u,
            one(RU),
            one(RU),
            u,
            RU(2),
            Y,
        ),
        case1_entries = (;
            p = local_nonunit_p + Y,
            q = local_nonunit_q + one(RU),
            r = Y,
            s = one(RU),
        ),
    )
    local_nonunit_context = Suslin.sl3_local_murthy_input_context(
        local_nonunit_target,
        Y;
        witness = local_nonunit_witness,
    )
    local_nonunit_cert = Suslin.realize_sl3_local_certificate(local_nonunit_context)
    @test local_nonunit_cert.branch == :murthy_q0_nonunit_bezout_resultant
    _mg_acceptance_assert_resultant_certificate(
        local_nonunit_cert;
        expected_source = :supplied_bezout_witness,
    )
    _mg_acceptance_assert_local_replay(local_nonunit_cert)
end

@testset "Issue 182 Murthy local SL3 closeout negative controls" begin
    R, (X,) = Oscar.polynomial_ring(QQ, ["X"])
    determinant_bad = _mg_acceptance_special_form_target(
        R,
        X + one(R),
        zero(R),
        zero(R),
        one(R),
    )
    determinant_err = _mg_acceptance_captured_error(
        () -> Suslin.realize_sl3_local_certificate(determinant_bad, X),
    )
    @test determinant_err isa ArgumentError
    @test occursin("determinant", sprint(showerror, determinant_err))

    nonmonic = _mg_acceptance_special_form_target(R, 2 * X + one(R), X, R(2), one(R))
    nonmonic_err = _mg_acceptance_captured_error(
        () -> Suslin.realize_sl3_local_certificate(nonmonic, X),
    )
    @test nonmonic_err isa ArgumentError
    @test occursin("p must be monic", sprint(showerror, nonmonic_err))

    RU, (u, Y) = Oscar.polynomial_ring(QQ, ["u", "X"])
    q_local = Y + u + 2
    p_local = Y * q_local + one(RU)
    missing_witness_target = _mg_acceptance_special_form_target(
        RU,
        p_local,
        q_local,
        Y + p_local * Y,
        p_local,
    )
    missing_witness_err = _mg_acceptance_captured_error(
        () -> Suslin.sl3_local_murthy_input_context(missing_witness_target, Y),
    )
    @test missing_witness_err isa ArgumentError
    @test occursin("local-unit witness", sprint(showerror, missing_witness_err))

    q_nonunit = Y + 2 * u
    p_nonunit = Y * q_nonunit + one(RU)
    unsupported_extraction_target = _mg_acceptance_special_form_target(
        RU,
        p_nonunit,
        q_nonunit,
        Y + p_nonunit * Y,
        p_nonunit,
    )
    unsupported_extraction_err = _mg_acceptance_captured_error(
        () -> Suslin.sl3_local_murthy_input_context(unsupported_extraction_target, Y),
    )
    @test unsupported_extraction_err isa ArgumentError
    @test occursin("unsupported local Bezout/resultant extraction", sprint(showerror, unsupported_extraction_err)) ||
        occursin("staged local SL_3 solver failure", sprint(showerror, unsupported_extraction_err))

    p_supported = X^2 + X + 1
    q_supported = one(R)
    supported_target = _mg_acceptance_special_form_target(
        R,
        p_supported,
        q_supported,
        p_supported^2 - one(R),
        p_supported,
    )
    supported_cert = Suslin.realize_sl3_local_certificate(supported_target, X)
    corrupted_factors = copy(supported_cert.factors)
    corrupted_factors[1] =
        corrupted_factors[1] * elementary_matrix(3, 1, 3, one(R), R)
    corrupted_cert = Suslin.SL3LocalRealizationCertificate(
        supported_cert.target,
        supported_cert.branch,
        corrupted_factors,
        supported_cert.selected_variable,
        supported_cert.witness,
    )
    @test !Suslin.verify_sl3_local_realization(corrupted_cert)
end

@testset "Issue 61 staged local SL3 unsupported boundary" begin
    R, (X,) = Oscar.polynomial_ring(QQ, ["X"])
    p = 2 * X + one(R)
    q = X
    r = R(2)
    s = one(R)
    target = _mg_acceptance_special_form_target(R, p, q, r, s)
    @test det(target) == one(R)
    @test !Suslin._is_monic_in_variable(p, 1, R)

    err = try
        Suslin.realize_sl3_local(target, X)
        nothing
    catch caught
        caught
    end
    @test err isa ArgumentError
    @test occursin("staged local SL_3 solver failure", sprint(showerror, err))
    @test occursin("p must be monic in X", sprint(showerror, err))
end
