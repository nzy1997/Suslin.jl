module Suslin

using Oscar

export suslin_polynomial_ring
export suslin_laurent_polynomial_ring
export elementary_matrix
export elementary_preconditioning_step
export replay_elementary_preconditioning
export verify_elementary_preconditioning
export block_embedding
export embed_factor_sequence
export compose_factor_sequences
export elementary_factorization
export SL3LocalObligation
export SLNToSL3Reduction
export SL3LocalReductionDiagnostic
export SLNToSL3ReductionDiagnostic
export reduce_sln_to_sl3
export diagnose_sln_to_sl3_reduction
export verify_sln_to_sl3_reduction
export classify_laurent_determinant
export LaurentGLFactorizationCertificate
export laurent_gl_factorization_certificate
export normalize_laurent_object
export normalize_laurent_gl_matrix
export solve_laurent_linear
export realize_cohn_type
export realize_conjugate_elementary
export lift_laurent_normalization
export verify_laurent_gl_factorization_certificate
export verify_laurent_gl_normalization
export verify_laurent_normalization
export verify_factorization
export LocalCertificate
export common_denominator_factor
export patched_substitution
export QuillenDenominatorData
export QuillenElementaryCorrection
export QuillenLocalContribution
export QuillenPatchVerification
export QuillenPatch
export construct_quillen_patch
export verify_quillen_patch

include("core/rings.jl")
include("core/polynomials.jl")
include("core/groebner_tools.jl")
include("core/laurent_linear_solve.jl")
include("core/elementary_matrices.jl")
include("core/unimodular.jl")
include("core/gl_laurent_normalization.jl")
include("algorithm/cohn_type.jl")
include("algorithm/column_reduction.jl")
include("algorithm/normality.jl")
include("algorithm/sl3_local.jl")
include("algorithm/laurent_column_peel.jl")
include("algorithm/laurent_gl_certificate.jl")
include("algorithm/factorization.jl")
include("algorithm/sln_to_sl3_reduction.jl")
include("algorithm/quillen_induction.jl")

function _coerce_into_ring(R, value, label::AbstractString)
    try
        return R(value)
    catch err
        if err isa ArgumentError || err isa MethodError
            throw(ArgumentError("$label must be coercible into the target ring"))
        end
        if err isa ErrorException && (
            occursin("Coercion not supported", err.msg) ||
            occursin("Unable to coerce polynomial", err.msg)
        )
            throw(ArgumentError("$label must be coercible into the target ring"))
        end
        rethrow()
    end
end

end
