using Test
using Suslin
using Oscar

@testset "documentation smoke" begin
    R, (X,) = Oscar.polynomial_ring(QQ, ["X"])
    A = matrix(R, [
        one(R)      one(R) + X      zero(R);
        X           one(R) + X + X^2 zero(R);
        zero(R)     zero(R)         one(R)
    ])

    factors = elementary_factorization(A)
    @test verify_factorization(A, factors)
end
