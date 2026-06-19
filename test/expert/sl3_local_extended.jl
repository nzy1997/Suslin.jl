using Test
using Suslin
using Oscar

function _sl3_extended_target(p, q, r, s, R)
    return matrix(R, [
        p q zero(R);
        r s zero(R);
        zero(R) zero(R) one(R)
    ])
end

function _sl3_extended_product(factors, R)
    product = identity_matrix(R, 3)
    for factor in factors
        product *= factor
    end
    return product
end

function _sl3_extended_captured_error(f)
    try
        f()
        return nothing
    catch err
        return err
    end
end

function _test_exact_sl3_factorization(target, factors)
    R = base_ring(target)
    @test target == _sl3_extended_product(factors, R)
    @test verify_factorization(target, factors)
end

@testset "extended local SL3 special-form realization" begin
    R, (X,) = Oscar.polynomial_ring(QQ, ["X"])

    q_open = one(R)
    r_open = X
    p_open = one(R) + q_open * r_open
    s_open = one(R)
    open_target = _sl3_extended_target(p_open, q_open, r_open, s_open, R)
    open_factors = Suslin.realize_sl3_local(p_open, q_open, r_open, s_open, X)
    @test length(open_factors) == 2
    _test_exact_sl3_factorization(open_target, open_factors)

    p_dual = one(R)
    q_dual = X + one(R)
    r_dual = X
    s_dual = one(R) + q_dual * r_dual
    dual_target = _sl3_extended_target(p_dual, q_dual, r_dual, s_dual, R)
    dual_factors = Suslin.realize_sl3_local(p_dual, q_dual, r_dual, s_dual, X)
    @test length(dual_factors) == 2
    _test_exact_sl3_factorization(dual_target, dual_factors)

    s_unit = R(2)
    q_unit = one(R)
    r_unit = 2 * X
    p_unit = X + R(1 // 2)
    unit_s_target = _sl3_extended_target(p_unit, q_unit, r_unit, s_unit, R)
    unit_s_factors = Suslin.realize_sl3_local(unit_s_target, X)
    @test length(unit_s_factors) > 2
    _test_exact_sl3_factorization(unit_s_target, unit_s_factors)

    L, (x, y) = suslin_laurent_polynomial_ring(QQ, ["x", "y"])
    p_laurent = x
    q_laurent = x * y
    r_laurent = one(L)
    s_laurent = x^-1 + y
    laurent_target = _sl3_extended_target(p_laurent, q_laurent, r_laurent, s_laurent, L)
    laurent_factors = Suslin.realize_sl3_local(
        p_laurent,
        q_laurent,
        r_laurent,
        s_laurent,
        x;
        check_monic=false,
    )
    @test length(laurent_factors) > 2
    _test_exact_sl3_factorization(laurent_target, laurent_factors)

    unsupported = _sl3_extended_target(X, -one(R), one(R), zero(R), R)
    unsupported_err = _sl3_extended_captured_error(() -> Suslin.realize_sl3_local(unsupported, X))
    @test unsupported_err isa ArgumentError
    @test occursin("staged local SL_3 solver failure", sprint(showerror, unsupported_err))
end
