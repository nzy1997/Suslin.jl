module SuslinStability

using Oscar

export suslin_polynomial_ring
export elementary_matrix

include("core/rings.jl")
include("core/polynomials.jl")
include("core/elementary_matrices.jl")

end
