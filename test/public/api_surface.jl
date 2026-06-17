using SuslinStability
using Test
using Oscar

@testset "api surface" begin
    @test isdefined(SuslinStability, :suslin_polynomial_ring)
    @test isdefined(SuslinStability, :elementary_matrix)
    @test isdefined(SuslinStability, :elementary_factorization)
    @test isdefined(SuslinStability, :realize_cohn_type)
    @test isdefined(SuslinStability, :realize_conjugate_elementary)
    @test isdefined(SuslinStability, :verify_factorization)
    @test SuslinStability.elementary_matrix === elementary_matrix
    @test SuslinStability.elementary_factorization === elementary_factorization
    @test SuslinStability.realize_cohn_type === realize_cohn_type
    @test SuslinStability.realize_conjugate_elementary === realize_conjugate_elementary
    @test SuslinStability.verify_factorization === verify_factorization

    R, (X,) = Oscar.polynomial_ring(QQ, ["X"])
    A = matrix(R, [
        one(R) + X  one(R)  zero(R);
        X           one(R)  zero(R);
        zero(R)     zero(R) one(R)
    ])

    factors = elementary_factorization(A)
    @test verify_factorization(A, factors)
end
