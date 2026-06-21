using Test
using Suslin
using Oscar

function product_of_factors(factors, R, n)
    product = identity_matrix(R, n)
    for factor in factors
        product *= factor
    end
    return product
end

@testset "local SL3 special-form realization" begin
    F = QQ
    R, (X,) = Oscar.polynomial_ring(F, ["X"])

    q = one(R)
    r = X
    p = one(R) + q * r
    s = one(R)
    target = matrix(R, [
        p q 0;
        r s 0;
        0 0 1
    ])

    factors = Suslin.realize_sl3_local(p, q, r, s, X)

    @test target == product_of_factors(factors, R, 3)

    p_dual = one(R)
    q_dual = X + one(R)
    r_dual = X
    s_dual = one(R) + r_dual * q_dual
    dual_target = matrix(R, [
        p_dual q_dual 0;
        r_dual s_dual 0;
        0      0      1
    ])
    dual_factors = Suslin.realize_sl3_local(p_dual, q_dual, r_dual, s_dual, X)

    @test dual_target == product_of_factors(dual_factors, R, 3)

    q_nonmonic = 2 * X
    r_nonmonic = one(R)
    p_nonmonic = one(R) + q_nonmonic * r_nonmonic
    nonmonic_target = matrix(R, [
        p_nonmonic q_nonmonic 0;
        r_nonmonic one(R)     0;
        0         0           1
    ])
    nonmonic_factors = Suslin.realize_sl3_local(
        p_nonmonic,
        q_nonmonic,
        r_nonmonic,
        one(R),
        X;
        check_monic=false,
    )

    @test nonmonic_target == product_of_factors(nonmonic_factors, R, 3)

    normalized_target = matrix(R, [
        one(R) + X X          0;
        -X         one(R) - X 0;
        0          0          1
    ])
    normalized_factors = Suslin.realize_sl3_local(
        one(R) + X,
        X,
        -X,
        one(R) - X,
        X,
    )
    @test normalized_target == product_of_factors(normalized_factors, R, 3)

    monic_err = try
        Suslin.realize_sl3_local(
            one(R) + 2 * X^2,
            2 * X,
            X,
            one(R),
            X,
        )
        nothing
    catch err
        err
    end
    @test monic_err isa ArgumentError
    @test occursin("p must be monic in X", sprint(showerror, monic_err))

    S, (Y,) = Oscar.polynomial_ring(QQ, ["Y"])
    @test_throws ArgumentError Suslin.realize_sl3_local(
        p,
        q,
        r,
        s,
        Y,
    )
end
