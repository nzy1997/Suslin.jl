using Test
using Suslin
using Oscar

const SL3_Q_DEGREE_FIXTURE_PATH =
    joinpath(@__DIR__, "..", "fixtures", "sl3_murthy_gupta_cases.jl")

if !isdefined(@__MODULE__, :SL3MurthyGuptaFixtureCatalog)
    include(SL3_Q_DEGREE_FIXTURE_PATH)
end

function _q_degree_target(R, p, q, r, s)
    return matrix(R, [
        p q zero(R);
        r s zero(R);
        zero(R) zero(R) one(R)
    ])
end

function _q_degree_case_data_from_fixture()
    catalog = SL3MurthyGuptaFixtureCatalog.catalog()
    entry = only(filter(case -> case.id == "mg-q-degree-normalization", catalog.cases))
    witness = first(entry.witnesses)
    return (;
        id = entry.id,
        R = entry.ring.object,
        X = entry.variable,
        p = entry.entries.p,
        q = entry.entries.q,
        r = entry.entries.r,
        s = entry.entries.s,
        quotient = witness.quotient,
        remainder = witness.remainder,
        normalized_s = witness.normalized_s,
        target = entry.target,
    )
end

function _q_degree_hand_checked_case()
    R, (X,) = Oscar.polynomial_ring(QQ, ["X"])
    p = X^2 + 1
    f = X^2
    g = X
    r = X
    normalized_s = one(R)
    s = normalized_s + f * r
    q = f * p + g
    target = _q_degree_target(R, p, q, r, s)
    @assert det(target) == one(R)
    return (;
        id = "hand-checked-q-degree",
        R,
        X,
        p,
        q,
        r,
        s,
        quotient = f,
        remainder = g,
        normalized_s,
        target,
    )
end

function _assert_q_degree_certificate(case)
    cert = Suslin.normalize_sl3_local_q_degree(
        case.p,
        case.q,
        case.r,
        case.s,
        case.X,
        check_monic = true,
    )
    correction = Suslin.elementary_matrix(3, 1, 2, case.quotient, case.R)
    normalized_target = _q_degree_target(
        case.R,
        case.p,
        case.remainder,
        case.r,
        case.normalized_s,
    )

    @test cert isa Suslin.SL3LocalRealizationCertificate
    @test cert.branch == :q_degree_normalization
    @test cert.target == case.target
    @test cert.selected_variable == case.X
    @test cert.factors == [correction]
    @test cert.witness.p == case.p
    @test cert.witness.q == case.q
    @test cert.witness.r == case.r
    @test cert.witness.s == case.s
    @test cert.witness.quotient == case.quotient
    @test cert.witness.remainder == case.remainder
    @test cert.witness.normalized_target == normalized_target
    @test cert.witness.correction == correction
    @test cert.witness.p_degree == degree(case.p, 1)
    @test cert.witness.remainder_degree == degree(case.remainder, 1)
    @test case.q == case.quotient * case.p + case.remainder
    @test degree(case.remainder, 1) < degree(case.p, 1)
    @test normalized_target * correction == case.target
    @test Suslin.verify_sl3_local_realization(cert)

    matrix_cert = Suslin.normalize_sl3_local_q_degree(case.target, case.X; check_monic = true)
    @test matrix_cert.target == cert.target
    @test matrix_cert.witness == cert.witness
    @test Suslin.verify_sl3_local_realization(matrix_cert)
end

@testset "local SL3 q-degree normalization certificate" begin
    cases = (_q_degree_case_data_from_fixture(), _q_degree_hand_checked_case())

    for case in cases
        @test degree(case.q, 1) >= degree(case.p, 1)
        _assert_q_degree_certificate(case)
    end

    case = first(cases)
    non_monic = 2 * case.X^2 + one(case.R)
    non_monic_q = one(case.R)
    non_monic_r = -one(case.R)
    non_monic_s = zero(case.R)
    non_monic_target = _q_degree_target(case.R, non_monic, non_monic_q, non_monic_r, non_monic_s)
    @test det(non_monic_target) == one(case.R)
    @test_throws ArgumentError Suslin.normalize_sl3_local_q_degree(
        non_monic,
        non_monic_q,
        non_monic_r,
        non_monic_s,
        case.X,
        check_monic = true,
    )

    cert = Suslin.normalize_sl3_local_q_degree(case.p, case.q, case.r, case.s, case.X)
    corrupted = Suslin.SL3LocalRealizationCertificate(
        cert.target,
        cert.branch,
        cert.factors,
        cert.selected_variable,
        merge(cert.witness, (; remainder = cert.witness.remainder + one(case.R))),
    )
    @test det(corrupted.target) == one(case.R)
    @test !Suslin.verify_sl3_local_realization(corrupted)

    extra_field = Suslin.SL3LocalRealizationCertificate(
        cert.target,
        cert.branch,
        cert.factors,
        cert.selected_variable,
        merge(cert.witness, (; extra = one(case.R))),
    )
    @test !Suslin.verify_sl3_local_realization(extra_field)
end
