using Test
using Suslin
using Oscar

@testset "elementary factorization small examples" begin
    R, (X,) = Oscar.polynomial_ring(QQ, ["X"])

    A = matrix(R, [
        one(R)      one(R) + X      zero(R);
        X           one(R) + X + X^2 zero(R);
        zero(R)     zero(R) one(R)
    ])

    factors = elementary_factorization(A)

    @test verify_factorization(A, factors)
    @test !verify_factorization(A, eltype(factors)[])

    bad = matrix(R, [
        one(R) + X  one(R)  zero(R);
        X           one(R) + X  zero(R);
        zero(R)     zero(R) one(R)
    ])
    @test_throws ArgumentError elementary_factorization(bad)
end
