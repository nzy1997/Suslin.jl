using Test
using Suslin
using Oscar

const SL3_RESULTANT_FIXTURE_PATH =
    joinpath(@__DIR__, "..", "fixtures", "sl3_murthy_gupta_cases.jl")

function _resultant_product(factors, R)
    product = identity_matrix(R, 3)
    for factor in factors
        product *= factor
    end
    return product
end

function _degree_in_variable(value, X)
    var_idx = findfirst(isequal(X), collect(gens(parent(value))))
    return degree(value, var_idx)
end

function _as_namedtuple(value)
    names = propertynames(value)
    return NamedTuple{names}(Tuple(getproperty(value, name) for name in names))
end

function _local_q0_reduction_copy(reduction;
        target = reduction.target,
        context = reduction.context,
        q0 = reduction.q0,
        q0_inverse = reduction.q0_inverse,
        p0 = reduction.p0,
        right_e21_coefficient = reduction.right_e21_coefficient,
        elimination_factor = reduction.elimination_factor,
        inverse_elimination_factor = reduction.inverse_elimination_factor,
        source_certificate = reduction.source_certificate,
        split_certificate = reduction.split_certificate,
        local_factor_replay = reduction.local_factor_replay,
        selected_variable = reduction.selected_variable,
        degree_p = reduction.degree_p,
        degree_q = reduction.degree_q)
    return Suslin.SL3LocalMurthyQUnitLocalReduction(
        target,
        context,
        q0,
        q0_inverse,
        p0,
        right_e21_coefficient,
        elimination_factor,
        inverse_elimination_factor,
        source_certificate,
        split_certificate,
        local_factor_replay,
        selected_variable,
        degree_p,
        degree_q,
    )
end

function _assert_resultant_certificate(cert; expected_source::Symbol)
    @test cert.branch == :murthy_q0_nonunit_bezout_resultant
    @test Suslin.verify_sl3_local_realization(cert)
    reduction = cert.witness.reduction
    @test reduction.witness_source == expected_source
    @test Suslin.verify_sl3_local_murthy_q0_nonunit_reduction(reduction)
    local_factor_replay = try
        getproperty(reduction, :local_factor_replay)
    catch err
        err isa FieldError || rethrow()
        nothing
    end
    if local_factor_replay === nothing
        @test Suslin.verify_factorization(cert.target, cert.factors)
        @test _resultant_product(cert.factors, base_ring(cert.target)) == cert.target
        @test reduction.child_certificate.branch == :murthy_q0_unit
        @test Suslin.verify_sl3_local_realization(reduction.child_certificate)
        @test reduction.bezout_target == reduction.first_elementary_factor * reduction.child_link_target
        @test reduction.target == reduction.left_factor * reduction.bezout_target
        @test reduction.branch_unit * reduction.branch_unit_inverse == one(base_ring(cert.target))
        @test _degree_in_variable(reduction.p_prime, reduction.selected_variable) < reduction.degree_q
        @test _degree_in_variable(reduction.q_prime, reduction.selected_variable) < reduction.degree_p
    else
        @test all(factor -> factor isa Suslin.SL3LocalElementaryFactor, cert.factors)
        @test reduction.local_factor_replay.target == cert.target
        @test reduction.local_factor_replay.factors == cert.factors
        @test Suslin.verify_sl3_local_elementary_factor_replay(reduction.local_factor_replay)
        child_reduction = reduction.child_certificate.witness.reduction
        @test child_reduction isa Suslin.SL3LocalMurthyQUnitLocalReduction
        @test child_reduction.context.q0 == reduction.branch_unit
        @test child_reduction.context.local_units.q0
    end
end

@testset "Murthy q(0)-nonunit Bezout/resultant branch for local SL3" begin
    include(SL3_RESULTANT_FIXTURE_PATH)
    catalog = SL3MurthyGuptaFixtureCatalog.catalog()
    by_id = Dict(entry.id => entry for entry in catalog.cases)

    supplied_fixture = by_id["mg-q0-nonunit-normalized-bezout-resultant"]
    supplied_cert = Suslin.realize_sl3_local_certificate(
        supplied_fixture.entries.p,
        supplied_fixture.entries.q,
        supplied_fixture.entries.r,
        supplied_fixture.entries.s,
        supplied_fixture.variable;
        murthy_q0_nonunit_witness = first(supplied_fixture.witnesses),
    )
    _assert_resultant_certificate(supplied_cert; expected_source = :supplied_bezout_witness)
    @test Suslin.realize_sl3_local(
        supplied_fixture.entries.p,
        supplied_fixture.entries.q,
        supplied_fixture.entries.r,
        supplied_fixture.entries.s,
        supplied_fixture.variable;
        murthy_q0_nonunit_witness = first(supplied_fixture.witnesses),
    ) == supplied_cert.factors

    extracted_fixture = by_id["mg-q0-nonunit-extracted-bezout-resultant"]
    extracted_cert = Suslin.realize_sl3_local_certificate(
        extracted_fixture.entries.p,
        extracted_fixture.entries.q,
        extracted_fixture.entries.r,
        extracted_fixture.entries.s,
        extracted_fixture.variable,
    )
    _assert_resultant_certificate(extracted_cert; expected_source = :extracted_bezout_witness)
    extracted_context = Suslin.sl3_local_murthy_input_context(
        extracted_fixture.target,
        extracted_fixture.variable,
    )
    @test extracted_context.bezout_witness === nothing
    @test Suslin.verify_sl3_local_murthy_input_context(extracted_context)
    extracted_context_cert = Suslin.realize_sl3_local_certificate(extracted_context)
    _assert_resultant_certificate(extracted_context_cert; expected_source = :extracted_bezout_witness)
    @test extracted_context_cert.factors == extracted_cert.factors
    @test Suslin.realize_sl3_local(
        extracted_fixture.entries.p,
        extracted_fixture.entries.q,
        extracted_fixture.entries.r,
        extracted_fixture.entries.s,
        extracted_fixture.variable,
    ) == extracted_cert.factors

    local_fixture = by_id["mg-local-q0-nonunit-bezout-at-u"]
    local_context = Suslin.sl3_local_murthy_input_context(
        local_fixture.target,
        local_fixture.variable;
        witness = first(local_fixture.witnesses),
    )
    local_cert = Suslin.realize_sl3_local_certificate(local_context)
    _assert_resultant_certificate(local_cert; expected_source = :supplied_bezout_witness)
    @test local_cert.target == local_fixture.target

    local_reduction = _as_namedtuple(local_cert.witness.reduction)
    corrupt_local_p_prime = merge(
        local_reduction,
        (; p_prime = local_reduction.p_prime + one(base_ring(local_cert.target))),
    )
    @test !Suslin.verify_sl3_local_murthy_q0_nonunit_reduction(corrupt_local_p_prime)

    corrupt_local_q_prime = merge(
        local_reduction,
        (; q_prime = local_reduction.q_prime + one(base_ring(local_cert.target))),
    )
    @test !Suslin.verify_sl3_local_murthy_q0_nonunit_reduction(corrupt_local_q_prime)

    corrupt_local_resultant = merge(
        local_reduction,
        (; resultant = local_reduction.resultant + one(base_ring(local_cert.target))),
    )
    @test !Suslin.verify_sl3_local_murthy_q0_nonunit_reduction(corrupt_local_resultant)

    corrupt_local_degree_p_prime = merge(
        local_reduction,
        (; degree_p_prime = local_reduction.degree_p_prime + 1),
    )
    @test !Suslin.verify_sl3_local_murthy_q0_nonunit_reduction(corrupt_local_degree_p_prime)

    child_reduction = local_reduction.child_certificate.witness.reduction
    bad_child_replay = Suslin.SL3LocalElementaryFactorReplay(
        child_reduction.local_factor_replay.target,
        copy(child_reduction.local_factor_replay.factors),
        child_reduction.local_factor_replay.selected_variable,
        child_reduction.local_factor_replay.mode,
        child_reduction.local_factor_replay.denominator_product,
        child_reduction.local_factor_replay.cleared_product,
        child_reduction.local_factor_replay.materialized_factors,
    )
    bad_child_replay_factors = copy(bad_child_replay.factors)
    bad_child_replay_factors[1] = Suslin.SL3LocalElementaryFactor(
        bad_child_replay_factors[1].R,
        bad_child_replay_factors[1].n,
        bad_child_replay_factors[1].row,
        bad_child_replay_factors[1].col,
        bad_child_replay_factors[1].numerator + one(bad_child_replay_factors[1].R),
        bad_child_replay_factors[1].denominator,
        bad_child_replay_factors[1].selected_variable,
        bad_child_replay_factors[1].local_unit_witness,
    )
    bad_child_reduction = _local_q0_reduction_copy(
        child_reduction;
        local_factor_replay = Suslin.SL3LocalElementaryFactorReplay(
            bad_child_replay.target,
            bad_child_replay_factors,
            bad_child_replay.selected_variable,
            bad_child_replay.mode,
            bad_child_replay.denominator_product,
            bad_child_replay.cleared_product,
            bad_child_replay.materialized_factors,
        ),
    )
    bad_child_certificate = Suslin.SL3LocalRealizationCertificate(
        local_reduction.child_certificate.target,
        local_reduction.child_certificate.branch,
        local_reduction.child_certificate.factors,
        local_reduction.child_certificate.selected_variable,
        merge(
            _as_namedtuple(local_reduction.child_certificate.witness),
            (; reduction = bad_child_reduction),
        ),
    )
    corrupt_child_witness = merge(
        local_reduction,
        (; child_certificate = bad_child_certificate),
    )
    @test !Suslin.verify_sl3_local_murthy_q0_nonunit_reduction(corrupt_child_witness)

    wrong_child_target = Suslin.SL3LocalRealizationCertificate(
        local_reduction.child_certificate.target + identity_matrix(
            base_ring(local_cert.target),
            3,
        ),
        local_reduction.child_certificate.branch,
        local_reduction.child_certificate.factors,
        local_reduction.child_certificate.selected_variable,
        local_reduction.child_certificate.witness,
    )
    corrupt_wrong_target = merge(
        local_reduction,
        (; child_certificate = wrong_child_target),
    )
    @test !Suslin.verify_sl3_local_murthy_q0_nonunit_reduction(corrupt_wrong_target)

    supplied_reduction = _as_namedtuple(supplied_cert.witness.reduction)
    corrupt_supplied_reduction = merge(
        supplied_reduction,
        (;
            p_prime = supplied_reduction.p_prime + one(base_ring(supplied_cert.target)),
        ),
    )
    @test !Suslin.verify_sl3_local_murthy_q0_nonunit_reduction(corrupt_supplied_reduction)

    shifted_supplied_witness = merge(
        first(supplied_fixture.witnesses),
        (;
            p_prime = first(supplied_fixture.witnesses).p_prime + supplied_fixture.entries.q,
            q_prime = first(supplied_fixture.witnesses).q_prime + supplied_fixture.entries.p,
        ),
    )
    @test shifted_supplied_witness.p_prime * supplied_fixture.entries.p -
        shifted_supplied_witness.q_prime * supplied_fixture.entries.q == one(supplied_fixture.ring.object)
    @test _degree_in_variable(shifted_supplied_witness.p_prime, supplied_fixture.variable) >=
        _degree_in_variable(supplied_fixture.entries.q, supplied_fixture.variable)
    @test _degree_in_variable(shifted_supplied_witness.q_prime, supplied_fixture.variable) >=
        _degree_in_variable(supplied_fixture.entries.p, supplied_fixture.variable)
    degree_bound_err = try
        Suslin.realize_sl3_local_certificate(
            supplied_fixture.entries.p,
            supplied_fixture.entries.q,
            supplied_fixture.entries.r,
            supplied_fixture.entries.s,
            supplied_fixture.variable;
            murthy_q0_nonunit_witness = shifted_supplied_witness,
        )
        nothing
    catch err
        err
    end
    @test degree_bound_err isa ArgumentError
    @test occursin("degree guard", sprint(showerror, degree_bound_err))

    extracted_reduction = _as_namedtuple(extracted_cert.witness.reduction)
    tampered_bezout_reduction = merge(
        extracted_reduction,
        (;
            bezout_target = extracted_reduction.bezout_target + identity_matrix(base_ring(extracted_cert.target), 3),
        ),
    )
    @test !Suslin.verify_sl3_local_murthy_q0_nonunit_reduction(tampered_bezout_reduction)

    tampered_child_reduction = merge(
        extracted_reduction,
        (;
            child_certificate = merge(
                _as_namedtuple(extracted_reduction.child_certificate),
                (; target = extracted_reduction.child_certificate.target + identity_matrix(base_ring(extracted_cert.target), 3)),
            ),
        ),
    )
    @test !Suslin.verify_sl3_local_murthy_q0_nonunit_reduction(tampered_child_reduction)
    @test !Suslin.verify_sl3_local_murthy_q0_nonunit_reduction((;))
end
