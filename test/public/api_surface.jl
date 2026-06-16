using SuslinStability
using Test

@testset "api surface" begin
    @test isdefined(SuslinStability, :suslin_polynomial_ring)
    @test isdefined(SuslinStability, :elementary_matrix)
    @test SuslinStability.elementary_matrix === elementary_matrix
end
