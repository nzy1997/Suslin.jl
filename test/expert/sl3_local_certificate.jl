using Test
using Suslin
using Oscar

function _sl3_certificate_target(p, q, r, s, R)
    return matrix(R, [
        p q zero(R);
        r s zero(R);
        zero(R) zero(R) one(R)
    ])
end

function _sl3_certificate_product(factors, R)
    product = identity_matrix(R, 3)
    for factor in factors
        product *= factor
    end
    return product
end

function _assert_sl3_certificate_replays(cert)
    R = base_ring(cert.target)
    @test _sl3_certificate_product(cert.factors, R) == cert.target
    @test Suslin.verify_factorization(cert.target, cert.factors)
    @test Suslin.verify_sl3_local_realization(cert)
end

function _tamper_first_factor(cert)
    R = base_ring(cert.target)
    factors = copy(cert.factors)
    factors[1] = identity_matrix(R, 3)
    return Suslin.SL3LocalRealizationCertificate(
        cert.target,
        cert.branch,
        factors,
        cert.selected_variable,
        cert.witness,
    )
end

function _tamper_witness_q(cert)
    R = base_ring(cert.target)
    return Suslin.SL3LocalRealizationCertificate(
        cert.target,
        cert.branch,
        cert.factors,
        cert.selected_variable,
        merge(cert.witness, (; q = cert.witness.q + one(R))),
    )
end

function _tamper_extra_witness_field(cert)
    R = base_ring(cert.target)
    return Suslin.SL3LocalRealizationCertificate(
        cert.target,
        cert.branch,
        cert.factors,
        cert.selected_variable,
        merge(cert.witness, (; extra = one(R))),
    )
end

@testset "local SL3 realization certificates" begin
    R, (X,) = Oscar.polynomial_ring(QQ, ["X"])

    q_open = one(R)
    r_open = X
    p_open = one(R) + q_open * r_open
    s_open = one(R)
    open_target = _sl3_certificate_target(p_open, q_open, r_open, s_open, R)
    open_cert = Suslin.realize_sl3_local_certificate(p_open, q_open, r_open, s_open, X)
    @test open_cert.target == open_target
    @test open_cert.branch == :open_s_one
    @test open_cert.selected_variable == X
    @test open_cert.witness.q == q_open
    @test open_cert.witness.r == r_open
    _assert_sl3_certificate_replays(open_cert)

    legacy_factors = Suslin.realize_sl3_local(p_open, q_open, r_open, s_open, X)
    @test legacy_factors isa Vector
    @test _sl3_certificate_product(legacy_factors, R) == open_target
    @test Suslin.verify_factorization(open_target, legacy_factors)

    q_dual = X + one(R)
    r_dual = X
    p_dual = one(R)
    s_dual = one(R) + q_dual * r_dual
    dual_cert = Suslin.realize_sl3_local_certificate(p_dual, q_dual, r_dual, s_dual, X)
    @test dual_cert.branch == :open_p_one
    _assert_sl3_certificate_replays(dual_cert)

    s_unit = R(2)
    q_s_unit = one(R)
    r_s_unit = 2 * X
    p_s_unit = X + R(1 // 2)
    s_unit_target = _sl3_certificate_target(p_s_unit, q_s_unit, r_s_unit, s_unit, R)
    s_unit_cert = Suslin.realize_sl3_local_certificate(s_unit_target, X)
    @test s_unit_cert.branch == :s_unit
    @test s_unit_cert.witness.pivot == s_unit
    @test s_unit_cert.witness.pivot_inverse == inv(s_unit)
    _assert_sl3_certificate_replays(s_unit_cert)

    p_unit = R(2)
    q_p_unit = X
    r_p_unit = one(R)
    s_p_unit = (X + one(R)) * R(1 // 2)
    p_unit_target = _sl3_certificate_target(p_unit, q_p_unit, r_p_unit, s_p_unit, R)
    p_unit_cert = Suslin.realize_sl3_local_certificate(
        p_unit,
        q_p_unit,
        r_p_unit,
        s_p_unit,
        X;
        check_monic = false,
    )
    @test p_unit_cert.target == p_unit_target
    @test p_unit_cert.branch == :p_unit
    @test p_unit_cert.witness.pivot == p_unit
    @test p_unit_cert.witness.pivot_inverse == inv(p_unit)
    _assert_sl3_certificate_replays(p_unit_cert)

    @test !Suslin.verify_sl3_local_realization(_tamper_first_factor(open_cert))
    @test !Suslin.verify_sl3_local_realization(_tamper_witness_q(open_cert))
    @test !Suslin.verify_sl3_local_realization(_tamper_extra_witness_field(open_cert))

    S, (Y,) = Oscar.polynomial_ring(QQ, ["Y"])
    foreign_variable_cert = Suslin.SL3LocalRealizationCertificate(
        open_cert.target,
        open_cert.branch,
        open_cert.factors,
        Y,
        open_cert.witness,
    )
    @test parent(Y) === S
    @test base_ring(foreign_variable_cert.target) === R
    @test !Suslin.verify_sl3_local_realization(foreign_variable_cert)
end
