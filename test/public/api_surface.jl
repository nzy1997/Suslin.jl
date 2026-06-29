using Suslin
using Test
using Oscar

@testset "api surface" begin
    @test isdefined(Suslin, :suslin_polynomial_ring)
    @test isdefined(Suslin, :suslin_laurent_polynomial_ring)
    @test isdefined(Suslin, :elementary_matrix)
    @test isdefined(Suslin, :max_elementary_factor_monomial_degree)
    @test isdefined(Suslin, :total_elementary_factor_offdiagonal_monomials)
    @test isdefined(Suslin, :elementary_preconditioning_step)
    @test isdefined(Suslin, :replay_elementary_preconditioning)
    @test isdefined(Suslin, :verify_elementary_preconditioning)
    @test isdefined(Suslin, :block_embedding)
    @test isdefined(Suslin, :embed_factor_sequence)
    @test isdefined(Suslin, :compose_factor_sequences)
    @test isdefined(Suslin, :elementary_factorization)
    @test isdefined(Suslin, :classify_laurent_determinant)
    @test isdefined(Suslin, :LaurentGLFactorizationCertificate)
    @test isdefined(Suslin, :LaurentLazyGLHoistCertificate)
    @test isdefined(Suslin, :laurent_gl_factorization_certificate)
    @test isdefined(Suslin, :normalize_laurent_gl_matrix)
    @test isdefined(Suslin, :normalize_laurent_object)
    @test isdefined(Suslin, :solve_laurent_linear)
    @test isdefined(Suslin, :realize_cohn_type)
    @test isdefined(Suslin, :realize_conjugate_elementary)
    @test isdefined(Suslin, :lift_laurent_normalization)
    @test isdefined(Suslin, :verify_laurent_gl_factorization_certificate)
    @test isdefined(Suslin, :verify_laurent_gl_normalization)
    @test isdefined(Suslin, :verify_laurent_normalization)
    @test isdefined(Suslin, :verify_factorization)
    @test isdefined(Suslin, :ConjugatedElementaryNormalityCertificate)
    @test isdefined(Suslin, :realize_conjugate_elementary_certificate)
    @test isdefined(Suslin, :verify_conjugate_elementary_certificate)
    @test isdefined(Suslin, :SL3LocalObligation)
    @test isdefined(Suslin, :SLNToSL3Reduction)
    @test isdefined(Suslin, :SL3LocalReductionDiagnostic)
    @test isdefined(Suslin, :SLNToSL3ReductionDiagnostic)
    @test isdefined(Suslin, :reduce_sln_to_sl3)
    @test isdefined(Suslin, :diagnose_sln_to_sl3_reduction)
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
    @test Suslin.LaurentGLFactorizationCertificate === LaurentGLFactorizationCertificate
    @test Suslin.LaurentLazyGLHoistCertificate === LaurentLazyGLHoistCertificate
    @test Suslin.laurent_gl_factorization_certificate === laurent_gl_factorization_certificate
    @test Suslin.elementary_matrix === elementary_matrix
    @test Suslin.max_elementary_factor_monomial_degree === max_elementary_factor_monomial_degree
    @test Suslin.total_elementary_factor_offdiagonal_monomials === total_elementary_factor_offdiagonal_monomials
    @test Suslin.elementary_preconditioning_step === elementary_preconditioning_step
    @test Suslin.replay_elementary_preconditioning === replay_elementary_preconditioning
    @test Suslin.verify_elementary_preconditioning === verify_elementary_preconditioning
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
    @test Suslin.verify_laurent_gl_factorization_certificate === verify_laurent_gl_factorization_certificate
    @test Suslin.verify_laurent_gl_normalization === verify_laurent_gl_normalization
    @test Suslin.verify_laurent_normalization === verify_laurent_normalization
    @test Suslin.verify_factorization === verify_factorization
    @test Suslin.ConjugatedElementaryNormalityCertificate === ConjugatedElementaryNormalityCertificate
    @test Suslin.realize_conjugate_elementary_certificate === realize_conjugate_elementary_certificate
    @test Suslin.verify_conjugate_elementary_certificate === verify_conjugate_elementary_certificate
    @test Suslin.SL3LocalObligation === SL3LocalObligation
    @test Suslin.SLNToSL3Reduction === SLNToSL3Reduction
    @test Suslin.SL3LocalReductionDiagnostic === SL3LocalReductionDiagnostic
    @test Suslin.SLNToSL3ReductionDiagnostic === SLNToSL3ReductionDiagnostic
    @test Suslin.reduce_sln_to_sl3 === reduce_sln_to_sl3
    @test Suslin.diagnose_sln_to_sl3_reduction === diagnose_sln_to_sl3_reduction
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

    LR, (u,) = suslin_laurent_polynomial_ring(QQ, ["u"])
    L = matrix(LR, [
        one(LR) u       zero(LR);
        zero(LR) one(LR) zero(LR);
        zero(LR) zero(LR) one(LR)
    ])

    one_arg_certificate = laurent_gl_factorization_certificate(L)
    @test one_arg_certificate isa LaurentGLFactorizationCertificate
    @test verify_laurent_gl_factorization_certificate(one_arg_certificate)

    lazy_keyword_certificate = laurent_gl_factorization_certificate(
        L;
        determinant_strategy = :lazy,
        correction_side = :row,
    )
    @test lazy_keyword_certificate isa LaurentLazyGLHoistCertificate
    @test lazy_keyword_certificate.determinant_source == :deferred_submatrix
    @test lazy_keyword_certificate.correction_side == :row
    @test verify_laurent_gl_factorization_certificate(lazy_keyword_certificate)

    rejected_strategy = try
        laurent_gl_factorization_certificate(L; determinant_strategy = :unsupported)
        nothing
    catch err
        err
    end
    @test rejected_strategy isa ArgumentError
    @test occursin(":eager", sprint(showerror, rejected_strategy))
    @test occursin(":lazy", sprint(showerror, rejected_strategy))

    rejected_side_without_lazy = try
        laurent_gl_factorization_certificate(L; correction_side = :column)
        nothing
    catch err
        err
    end
    @test rejected_side_without_lazy isa ArgumentError
    @test occursin("determinant_strategy = :lazy", sprint(showerror, rejected_side_without_lazy))
end
