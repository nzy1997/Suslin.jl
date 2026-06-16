using SuslinStability
using Test

@testset "api surface" begin
    @test isdefined(SuslinStability, :suslin_polynomial_ring)
end
