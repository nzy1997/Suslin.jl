using Test
using Suslin
using Oscar

const SL3_Q0_UNIT_FIXTURE_PATH =
    joinpath(@__DIR__, "..", "fixtures", "sl3_murthy_gupta_cases.jl")

function _q0_target(R, p, q, r, s)
    return matrix(R, [
        p q zero(R);
        r s zero(R);
        zero(R) zero(R) one(R)
    ])
end

function _q0_product(factors, R)
    product = identity_matrix(R, 3)
    for factor in factors
        product *= factor
    end
    return product
end

function _q0_degree(value, X)
    var_idx = findfirst(isequal(X), collect(gens(parent(value))))
    return degree(value, var_idx)
end

function _assert_elementary_sequence(target, factors)
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
    @test _q0_product(factors, R) == target
    @test Suslin.verify_factorization(target, factors)
end

function _assert_reduction_replay(reduction)
    @test Suslin.verify_sl3_local_murthy_q_unit_reduction(reduction)
    split_cert = reduction.split_certificate
    split = split_cert.witness.split
    @test split.original_target == reduction.eliminated_target
    @test split.witness.a == reduction.selected_variable
    @test split.witness.a_prime == reduction.p_prime
    @test split.witness.b == reduction.eliminated_target[1, 2]
    @test _q0_degree(reduction.p_prime, reduction.selected_variable) < reduction.degree_p
    @test Suslin.verify_sl3_local_realization(split_cert)
    @test Suslin.verify_sl3_local_split_lemma_replay(split)
    @test Suslin.verify_sl3_local_realization(split_cert.witness.first_child_certificate)
    @test Suslin.verify_sl3_local_realization(split_cert.witness.second_child_certificate)
end

function _assert_q0_certificate(cert; normalization_expected::Bool)
    @test cert.branch == :murthy_q0_unit
    @test Suslin.verify_sl3_local_realization(cert)
    _assert_elementary_sequence(cert.target, cert.factors)
    if normalization_expected
        @test cert.witness.normalization !== nothing
        @test Suslin.verify_sl3_local_q_degree_normalization(cert.witness.normalization)
        @test cert.witness.normalized_certificate !== nothing
        @test cert.witness.reduction === nothing
        @test cert.witness.normalized_certificate.branch == :murthy_q0_unit
        _assert_q0_certificate(cert.witness.normalized_certificate; normalization_expected = false)
    else
        @test cert.witness.normalization === nothing
        @test cert.witness.normalized_certificate === nothing
        @test cert.witness.reduction !== nothing
        _assert_reduction_replay(cert.witness.reduction)
    end
end

@testset "Murthy q(0)-unit recursive branch for local SL3" begin
    include(SL3_Q0_UNIT_FIXTURE_PATH)
    catalog = SL3MurthyGuptaFixtureCatalog.catalog()
    by_id = Dict(entry.id => entry for entry in catalog.cases)

    fixture = by_id["mg-q0-unit-recursion"]
    fixture_cert = Suslin.realize_sl3_local_certificate(
        fixture.entries.p,
        fixture.entries.q,
        fixture.entries.r,
        fixture.entries.s,
        fixture.variable,
    )
    _assert_q0_certificate(fixture_cert; normalization_expected = false)
    @test fixture_cert.target == fixture.target
    @test Suslin.realize_sl3_local(
        fixture.entries.p,
        fixture.entries.q,
        fixture.entries.r,
        fixture.entries.s,
        fixture.variable,
    ) == fixture_cert.factors
    fixture_reduction = fixture_cert.witness.reduction
    replacement_split_certificate = Suslin.realize_sl3_local_certificate(
        fixture_reduction.eliminated_target,
        fixture.variable,
    )
    @test replacement_split_certificate.branch == :q_unit
    tampered_reduction = Suslin.SL3LocalMurthyQUnitReduction(
        fixture_reduction.target,
        fixture_reduction.q0,
        fixture_reduction.q0_inverse,
        fixture_reduction.p0,
        fixture_reduction.right_e21_coefficient,
        fixture_reduction.eliminated_target,
        fixture_reduction.elimination_factor,
        fixture_reduction.inverse_elimination_factor,
        fixture_reduction.p_prime,
        replacement_split_certificate,
        fixture_reduction.selected_variable,
        fixture_reduction.degree_p,
        fixture_reduction.degree_p_prime,
        fixture_reduction.locality_witness,
    )
    @test !Suslin.verify_sl3_local_murthy_q_unit_reduction(tampered_reduction)

    R, (X,) = Oscar.polynomial_ring(QQ, ["X"])

    p = X^3 + X + 1
    q = X^2 + X + 1
    _, s, minus_r = gcdx(p, q)
    r = -minus_r
    normalized_target = _q0_target(R, p, q, r, s)
    @test det(normalized_target) == one(R)
    @test _q0_degree(q, X) < _q0_degree(p, X)
    @test !is_unit(p)
    @test !is_unit(s)
    normalized_cert = Suslin.realize_sl3_local_certificate(normalized_target, X)
    _assert_q0_certificate(normalized_cert; normalization_expected = false)

    normalizing_fixture = by_id["mg-q-degree-normalization"]
    normalizing_cert = Suslin.realize_sl3_local_certificate(
        normalizing_fixture.entries.p,
        normalizing_fixture.entries.q,
        normalizing_fixture.entries.r,
        normalizing_fixture.entries.s,
        normalizing_fixture.variable,
    )
    _assert_q0_certificate(normalizing_cert; normalization_expected = true)
    @test normalizing_cert.target == normalizing_fixture.target

    nonunit_q0_p = X^2 + one(R)
    nonunit_q0_q = X
    nonunit_q0_s = X + one(R)
    nonunit_q0_r = div(nonunit_q0_p * nonunit_q0_s - one(R), X)
    nonunit_target = _q0_target(R, nonunit_q0_p, nonunit_q0_q, nonunit_q0_r, nonunit_q0_s)
    @test det(nonunit_target) == one(R)
    @test _q0_degree(nonunit_q0_q, X) < _q0_degree(nonunit_q0_p, X)
    @test_throws ArgumentError Suslin.realize_sl3_local_certificate(nonunit_target, X)
end
