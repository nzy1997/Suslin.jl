using SuslinStability
using Test

@testset "public" begin
    include("public/api_surface.jl")
end

@testset "expert" begin
    include("expert/elementary_matrices.jl")
    include("expert/cohn_type.jl")
    include("expert/normality.jl")
    include("expert/quillen_induction.jl")
    include("expert/unimodular_columns.jl")
end

@testset "internal" begin
    include("internal/rings.jl")
end
