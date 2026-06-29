using Test
using Suslin
using Oscar

const SL3_Q_DEGREE_FIXTURE_PATH =
    joinpath(@__DIR__, "..", "fixtures", "sl3_murthy_gupta_cases.jl")

function _qdegree_target(R, p, q, r, s)
    return matrix(R, [
        p q zero(R);
        r s zero(R);
        zero(R) zero(R) one(R)
    ])
end

function _qdegree_product(factors, R)
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

function _assert_qdegree_record(record, p, q, r, s, X, expected_f, expected_g)
    R = parent(p)
    @test record.target == _qdegree_target(R, p, q, r, s)
    @test record.quotient == expected_f
    @test record.remainder == expected_g
    @test q == record.quotient * p + record.remainder
    @test _degree_in_variable(record.remainder, X) < _degree_in_variable(p, X)
    @test record.normalized_target ==
        _qdegree_target(R, p, record.remainder, r, s - record.quotient * r)
    @test record.elementary_correction ==
        Suslin.elementary_matrix(3, 1, 2, record.quotient, R)
    @test record.normalized_target * record.elementary_correction == record.target
    @test Suslin.verify_sl3_local_q_degree_normalization(record)
end

@testset "Murthy q-degree normalization replay for local SL3" begin
    include(SL3_Q_DEGREE_FIXTURE_PATH)
    catalog = SL3MurthyGuptaFixtureCatalog.catalog()
    by_id = Dict(entry.id => entry for entry in catalog.cases)

    fixture = by_id["mg-q-degree-normalization"]
    witness = first(fixture.witnesses)
    record = Suslin.sl3_local_q_degree_normalization(
        fixture.entries.p,
        fixture.entries.q,
        fixture.entries.r,
        fixture.entries.s,
        fixture.variable,
    )
    _assert_qdegree_record(
        record,
        fixture.entries.p,
        fixture.entries.q,
        fixture.entries.r,
        fixture.entries.s,
        fixture.variable,
        witness.quotient,
        witness.remainder,
    )
    @test record.normalized_target[2, 2] == witness.normalized_s

    cert = Suslin.sl3_local_q_degree_normalization_certificate(record)
    @test cert.branch == :murthy_q_degree_normalization
    @test cert.witness.normalization == record
    @test _qdegree_product(cert.factors, fixture.ring.object) == record.target
    @test Suslin.verify_sl3_local_realization(cert)

    local_fixture = by_id["mg-local-q-degree-qq-u-x"]
    local_context = Suslin.sl3_local_murthy_input_context(
        local_fixture.target,
        local_fixture.variable;
        witness = first(local_fixture.witnesses),
    )
    local_record = Suslin.sl3_local_q_degree_normalization(local_context)
    local_witness = first(local_fixture.witnesses)
    _assert_qdegree_record(
        local_record,
        local_fixture.entries.p,
        local_fixture.entries.q,
        local_fixture.entries.r,
        local_fixture.entries.s,
        local_fixture.variable,
        local_witness.quotient,
        local_witness.remainder,
    )
    local_cert = Suslin.sl3_local_q_degree_normalization_certificate(local_context)
    @test local_cert.branch == :murthy_q_degree_normalization
    @test local_cert.target == local_fixture.target
    @test Suslin.verify_sl3_local_realization(local_cert)

    matrix_record = Suslin.sl3_local_q_degree_normalization(fixture.target, fixture.variable)
    @test matrix_record == record

    R, (X,) = Oscar.polynomial_ring(QQ, ["X"])
    p = X^2 + X + 1
    r = one(R)
    s = X^2
    q = p * s - one(R)
    second = Suslin.sl3_local_q_degree_normalization(p, q, r, s, X)
    _assert_qdegree_record(second, p, q, r, s, X, X^2, -one(R))
    second_cert = Suslin.sl3_local_q_degree_normalization_certificate(p, q, r, s, X)
    @test Suslin.verify_sl3_local_realization(second_cert)

    local_q0_unit = by_id["mg-local-q0-unit-at-u"]
    local_q0_context = Suslin.sl3_local_murthy_input_context(
        local_q0_unit.target,
        local_q0_unit.variable;
        witness = first(local_q0_unit.witnesses),
    )
    @test_throws ArgumentError Suslin.sl3_local_q_degree_normalization(local_q0_context)

    L, (T,) = Oscar.laurent_polynomial_ring(QQ, ["T"])
    laurent_p = T
    laurent_q = one(L)
    laurent_r = zero(L)
    laurent_s = inv(T)
    @test det(_qdegree_target(L, laurent_p, laurent_q, laurent_r, laurent_s)) == one(L)
    @test_throws ArgumentError Suslin.sl3_local_q_degree_normalization(
        laurent_p,
        laurent_q,
        laurent_r,
        laurent_s,
        T,
    )

    nonmonic_p = 2 * X + one(R)
    nonmonic_q = 2 * X
    @test det(_qdegree_target(R, nonmonic_p, nonmonic_q, one(R), one(R))) == one(R)
    @test_throws ArgumentError Suslin.sl3_local_q_degree_normalization(
        nonmonic_p,
        nonmonic_q,
        one(R),
        one(R),
        X,
    )

    bad_record = Suslin.SL3LocalQDegreeNormalization(
        second.target,
        second.quotient + one(R),
        second.remainder,
        second.normalized_target,
        second.elementary_correction,
        second.selected_variable,
    )
    @test det(bad_record.target) == one(R)
    @test !Suslin.verify_sl3_local_q_degree_normalization(bad_record)
    bad_cert = Suslin.SL3LocalRealizationCertificate(
        bad_record.target,
        :murthy_q_degree_normalization,
        [bad_record.normalized_target, bad_record.elementary_correction],
        bad_record.selected_variable,
        (; normalization = bad_record),
    )
    @test !Suslin.verify_sl3_local_realization(bad_cert)

    bad_remainder = Suslin.SL3LocalQDegreeNormalization(
        second.target,
        second.quotient,
        second.remainder + one(R),
        second.normalized_target,
        second.elementary_correction,
        second.selected_variable,
    )
    @test !Suslin.verify_sl3_local_q_degree_normalization(bad_remainder)

    bad_normalized_target = copy(second.normalized_target)
    bad_normalized_target[1, 2] += one(R)
    bad_normalized_record = Suslin.SL3LocalQDegreeNormalization(
        second.target,
        second.quotient,
        second.remainder,
        bad_normalized_target,
        second.elementary_correction,
        second.selected_variable,
    )
    @test !Suslin.verify_sl3_local_q_degree_normalization(bad_normalized_record)
    bad_normalized_cert = Suslin.SL3LocalRealizationCertificate(
        bad_normalized_record.target,
        :murthy_q_degree_normalization,
        [bad_normalized_record.normalized_target, bad_normalized_record.elementary_correction],
        bad_normalized_record.selected_variable,
        (; normalization = bad_normalized_record),
    )
    @test !Suslin.verify_sl3_local_realization(bad_normalized_cert)

    bad_correction = Suslin.SL3LocalQDegreeNormalization(
        second.target,
        second.quotient,
        second.remainder,
        second.normalized_target,
        identity_matrix(R, 3),
        second.selected_variable,
    )
    @test !Suslin.verify_sl3_local_q_degree_normalization(bad_correction)

    @test !Suslin.verify_sl3_local_q_degree_normalization(Suslin.SL3LocalQDegreeNormalization(
        second.target,
        second.quotient,
        second.remainder,
        second.normalized_target,
        second.elementary_correction,
        X + one(R),
    ))

    S, (Y,) = Oscar.polynomial_ring(QQ, ["Y"])
    @test !Suslin.verify_sl3_local_q_degree_normalization(Suslin.SL3LocalQDegreeNormalization(
        second.target,
        second.quotient,
        second.remainder,
        second.normalized_target,
        second.elementary_correction,
        Y,
    ))
    @test !Suslin.verify_sl3_local_q_degree_normalization((; target = second.target))
end
