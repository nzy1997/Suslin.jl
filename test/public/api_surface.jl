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
    @test isdefined(Suslin, :classify_laurent_determinant)
    @test isdefined(Suslin, :normalize_laurent_gl_matrix)
    @test isdefined(Suslin, :normalize_laurent_object)
    @test isdefined(Suslin, :solve_laurent_linear)
    @test isdefined(Suslin, :realize_cohn_type)
    @test isdefined(Suslin, :realize_conjugate_elementary)
    @test isdefined(Suslin, :lift_laurent_normalization)
    @test isdefined(Suslin, :verify_laurent_gl_normalization)
    @test isdefined(Suslin, :verify_laurent_normalization)
    @test isdefined(Suslin, :verify_factorization)
    @test isdefined(Suslin, :SL3LocalObligation)
    @test isdefined(Suslin, :SLNToSL3Reduction)
    @test isdefined(Suslin, :reduce_sln_to_sl3)
    @test isdefined(Suslin, :verify_sln_to_sl3_reduction)
    @test isdefined(Suslin, :LocalCertificate)
    @test isdefined(Suslin, :common_denominator_factor)
    @test isdefined(Suslin, :patched_substitution)
    @test isdefined(Suslin, :QuillenDenominatorData)
    @test isdefined(Suslin, :QuillenElementaryCorrection)
    @test isdefined(Suslin, :QuillenLocalContribution)
    @test isdefined(Suslin, :QuillenPatchVerification)
    @test isdefined(Suslin, :QuillenPatch)
    @test isdefined(Suslin, :construct_quillen_patch)
    @test isdefined(Suslin, :verify_quillen_patch)
    @test Suslin.classify_laurent_determinant === classify_laurent_determinant
    @test Suslin.elementary_matrix === elementary_matrix
    @test Suslin.block_embedding === block_embedding
    @test Suslin.embed_factor_sequence === embed_factor_sequence
    @test Suslin.compose_factor_sequences === compose_factor_sequences
    @test Suslin.lift_laurent_normalization === lift_laurent_normalization
    @test Suslin.normalize_laurent_gl_matrix === normalize_laurent_gl_matrix
    @test Suslin.normalize_laurent_object === normalize_laurent_object
    @test Suslin.suslin_laurent_polynomial_ring === suslin_laurent_polynomial_ring
    @test Suslin.solve_laurent_linear === solve_laurent_linear
    @test Suslin.elementary_factorization === elementary_factorization
    @test Suslin.realize_cohn_type === realize_cohn_type
    @test Suslin.realize_conjugate_elementary === realize_conjugate_elementary
    @test Suslin.verify_laurent_gl_normalization === verify_laurent_gl_normalization
    @test Suslin.verify_laurent_normalization === verify_laurent_normalization
    @test Suslin.verify_factorization === verify_factorization
    @test Suslin.SL3LocalObligation === SL3LocalObligation
    @test Suslin.SLNToSL3Reduction === SLNToSL3Reduction
    @test Suslin.reduce_sln_to_sl3 === reduce_sln_to_sl3
    @test Suslin.verify_sln_to_sl3_reduction === verify_sln_to_sl3_reduction
    @test Suslin.LocalCertificate === LocalCertificate
    @test Suslin.common_denominator_factor === common_denominator_factor
    @test Suslin.patched_substitution === patched_substitution
    @test Suslin.QuillenDenominatorData === QuillenDenominatorData
    @test Suslin.QuillenElementaryCorrection === QuillenElementaryCorrection
    @test Suslin.QuillenLocalContribution === QuillenLocalContribution
    @test Suslin.QuillenPatchVerification === QuillenPatchVerification
    @test Suslin.QuillenPatch === QuillenPatch
    @test Suslin.construct_quillen_patch === construct_quillen_patch
    @test Suslin.verify_quillen_patch === verify_quillen_patch

    R, (X,) = Oscar.polynomial_ring(QQ, ["X"])
    A = matrix(R, [
        one(R) + X  one(R)  zero(R);
        X           one(R)  zero(R);
        zero(R)     zero(R) one(R)
    ])

    factors = elementary_factorization(A)
    @test verify_factorization(A, factors)
end
