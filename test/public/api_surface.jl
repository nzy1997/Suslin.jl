using Suslin
using Test
using Oscar

@testset "api surface" begin
    @test isdefined(Suslin, :suslin_polynomial_ring)
    @test isdefined(Suslin, :suslin_laurent_polynomial_ring)
    @test isdefined(Suslin, :elementary_matrix)
    @test isdefined(Suslin, :block_embedding)
    @test isdefined(Suslin, :embed_factor_sequence)
    @test isdefined(Suslin, :compose_factor_sequences)
    @test isdefined(Suslin, :elementary_factorization)
    @test isdefined(Suslin, :realize_cohn_type)
    @test isdefined(Suslin, :realize_conjugate_elementary)
    @test isdefined(Suslin, :verify_factorization)
    @test Suslin.elementary_matrix === elementary_matrix
    @test Suslin.block_embedding === block_embedding
    @test Suslin.embed_factor_sequence === embed_factor_sequence
    @test Suslin.compose_factor_sequences === compose_factor_sequences
    @test Suslin.suslin_laurent_polynomial_ring === suslin_laurent_polynomial_ring
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
