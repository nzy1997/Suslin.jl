using Test
using Suslin
using Oscar

const SL3_SPLIT_FIXTURE_PATH =
    joinpath(@__DIR__, "..", "fixtures", "sl3_murthy_gupta_cases.jl")

function _split_target(R, p, q, r, s)
    return matrix(R, [
        p q zero(R);
        r s zero(R);
        zero(R) zero(R) one(R)
    ])
end

function _split_product(factors, R)
    product = identity_matrix(R, 3)
    for factor in factors
        product *= factor
    end
    return product
end

function _split_wrapped_product(replay)
    R = base_ring(replay.original_target)
    return _split_product(replay.prefix_factors, R) *
        replay.first_child_target *
        _split_product(replay.middle_factors, R) *
        replay.second_child_target *
        _split_product(replay.suffix_factors, R)
end

function _assert_elementary_3x3_factor(factor, R)
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

function _assert_split_replay(replay)
    R = base_ring(replay.original_target)
    witness = replay.witness

    @test witness.a * witness.a_prime * witness.d - witness.b * witness.c == one(R)
    @test witness.a * witness.d1 - witness.b * witness.c1 == one(R)
    @test witness.a_prime * witness.d2 - witness.b * witness.c2 == one(R)
    @test replay.original_target ==
        _split_target(R, witness.a * witness.a_prime, witness.b, witness.c, witness.d)
    @test replay.first_child_target ==
        _split_target(R, witness.a, witness.b, witness.c1, witness.d1)
    @test replay.second_child_target ==
        _split_target(R, witness.a_prime, witness.b, witness.c2, witness.d2)
    @test det(replay.original_target) == one(R)
    @test det(replay.first_child_target) == one(R)
    @test det(replay.second_child_target) == one(R)
    @test replay.wrapper_factors == vcat(
        replay.prefix_factors,
        replay.middle_factors,
        replay.suffix_factors,
    )
    for factor in replay.wrapper_factors
        _assert_elementary_3x3_factor(factor, R)
    end
    @test replay.reassembled_product == _split_wrapped_product(replay)
    @test replay.reassembled_product == replay.original_target
    @test Suslin.verify_sl3_local_split_lemma_replay(replay)
end

function _assert_split_certificate(replay, X)
    first_cert = Suslin.realize_sl3_local_certificate(replay.first_child_target, X)
    second_cert = Suslin.realize_sl3_local_certificate(replay.second_child_target, X)
    cert = Suslin.sl3_local_split_lemma_certificate(replay, first_cert, second_cert, X)
    R = base_ring(replay.original_target)

    @test cert.target == replay.original_target
    @test cert.branch == :murthy_split_lemma
    @test cert.witness.split == replay
    @test cert.witness.first_child_certificate == first_cert
    @test cert.witness.second_child_certificate == second_cert
    @test _split_product(cert.factors, R) == replay.original_target
    @test Suslin.verify_factorization(cert.target, cert.factors)
    @test Suslin.verify_sl3_local_realization(cert)

    @test_throws ArgumentError Suslin.sl3_local_split_lemma_certificate(
        replay,
        second_cert,
        first_cert,
        X,
    )
end

function _split_replay_from_witness(witness; split_id = :murthy_split_lemma)
    return Suslin.sl3_local_split_lemma_replay(
        witness.a,
        witness.a_prime,
        witness.b,
        witness.c,
        witness.d,
        witness.c1,
        witness.d1,
        witness.c2,
        witness.d2;
        split_id,
    )
end

@testset "Murthy split lemma replay for local SL3" begin
    include(SL3_SPLIT_FIXTURE_PATH)
    catalog = SL3MurthyGuptaFixtureCatalog.catalog()
    by_id = Dict(entry.id => entry for entry in catalog.cases)

    fixture_entry = by_id["mg-split-lemma-x-square"]
    fixture_witness = first(fixture_entry.witnesses)
    fixture_replay = _split_replay_from_witness(
        fixture_witness;
        split_id = Symbol(fixture_entry.id),
    )
    @test fixture_replay.split_id == Symbol("mg-split-lemma-x-square")
    @test fixture_replay.original_target == fixture_entry.target
    _assert_split_replay(fixture_replay)

    R, (X,) = Oscar.polynomial_ring(QQ, ["X"])

    open_witness = (;
        a = X + 1,
        a_prime = X + 2,
        b = one(R),
        c = (X + 1) * (X + 2) - one(R),
        d = one(R),
        c1 = X,
        d1 = one(R),
        c2 = X + 1,
        d2 = one(R),
    )
    open_replay = _split_replay_from_witness(
        open_witness;
        split_id = :mg_split_lemma_open_children,
    )
    _assert_split_replay(open_replay)
    _assert_split_certificate(open_replay, X)

    unit_witness = (;
        a = X,
        a_prime = X + 1,
        b = one(R),
        c = X * (X + 1) - one(R),
        d = one(R),
        c1 = 2 * X - one(R),
        d1 = R(2),
        c2 = 3 * (X + 1) - one(R),
        d2 = R(3),
    )
    unit_replay = _split_replay_from_witness(
        unit_witness;
        split_id = :mg_split_lemma_unit_children,
    )
    _assert_split_replay(unit_replay)
    _assert_split_certificate(unit_replay, X)

    @test_throws ArgumentError Suslin.sl3_local_split_lemma_replay(
        unit_witness.a,
        unit_witness.a_prime,
        unit_witness.b,
        unit_witness.c,
        unit_witness.d,
        unit_witness.c1 + one(R),
        unit_witness.d1,
        unit_witness.c2,
        unit_witness.d2,
    )

    tampered_replay = Suslin.SL3LocalSplitLemmaReplay(
        unit_replay.split_id,
        unit_replay.original_target,
        unit_replay.first_child_target,
        unit_replay.second_child_target,
        vcat(unit_replay.prefix_factors[1:(end - 1)], [identity_matrix(R, 3)]),
        unit_replay.middle_factors,
        unit_replay.suffix_factors,
        unit_replay.wrapper_factors,
        unit_replay.reassembled_product,
        unit_replay.witness,
    )
    @test !Suslin.verify_sl3_local_split_lemma_replay(tampered_replay)
end
