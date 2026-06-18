using Suslin
using Test
using Oscar

@testset "api surface" begin
    @test isdefined(Suslin, :suslin_polynomial_ring)
    @test isdefined(Suslin, :elementary_matrix)
    @test isdefined(Suslin, :elementary_factorization)
    @test isdefined(Suslin, :realize_cohn_type)
    @test isdefined(Suslin, :realize_conjugate_elementary)
    @test isdefined(Suslin, :verify_factorization)
    @test Suslin.elementary_matrix === elementary_matrix
    @test Suslin.elementary_factorization === elementary_factorization
    @test Suslin.realize_cohn_type === realize_cohn_type
    @test Suslin.realize_conjugate_elementary === realize_conjugate_elementary
    @test Suslin.verify_factorization === verify_factorization

    R, (X,) = Oscar.polynomial_ring(QQ, ["X"])
    A = matrix(R, [
        one(R) + X  one(R)  zero(R);
        X           one(R)  zero(R);
        zero(R)     zero(R) one(R)
    ])

    factors = elementary_factorization(A)
    @test verify_factorization(A, factors)
end
