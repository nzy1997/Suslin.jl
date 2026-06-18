module Suslin

using Oscar

export suslin_polynomial_ring
export elementary_matrix
export elementary_factorization
export realize_cohn_type
export realize_conjugate_elementary
export verify_factorization

include("core/rings.jl")
include("core/polynomials.jl")
include("core/groebner_tools.jl")
include("core/elementary_matrices.jl")
include("core/unimodular.jl")
include("algorithm/cohn_type.jl")
include("algorithm/column_reduction.jl")
include("algorithm/normality.jl")
include("algorithm/sl3_local.jl")
include("algorithm/factorization.jl")
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
