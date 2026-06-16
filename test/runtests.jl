using SuslinStability
using Test

@testset "public" begin
    include("public/api_surface.jl")
end

@testset "expert" begin
    include("expert/elementary_matrices.jl")
end

@testset "internal" begin
    include("internal/rings.jl")
end
