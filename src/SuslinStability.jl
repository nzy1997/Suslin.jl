module SuslinStability

using Oscar

export suslin_polynomial_ring
export elementary_matrix
export realize_cohn_type

include("core/rings.jl")
include("core/polynomials.jl")
include("core/elementary_matrices.jl")
include("algorithm/cohn_type.jl")

function _coerce_into_ring(R, value, label::AbstractString)
    try
        return R(value)
    catch err
        if err isa ArgumentError || err isa MethodError || err isa ErrorException
            throw(ArgumentError("$label must be coercible into the target ring"))
        end
        rethrow()
    end
end

end
