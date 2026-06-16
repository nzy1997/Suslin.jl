using SuslinStability
using Test

@testset "public" begin
    include("public/api_surface.jl")
end

@testset "expert" begin
end

@testset "internal" begin
    include("internal/rings.jl")
end
