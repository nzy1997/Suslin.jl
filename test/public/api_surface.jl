using SuslinStability
using Test

@testset "api surface" begin
    @test isdefined(SuslinStability, :suslin_polynomial_ring)
    @test isdefined(SuslinStability, :elementary_matrix)
    @test isdefined(SuslinStability, :realize_cohn_type)
    @test isdefined(SuslinStability, :realize_conjugate_elementary)
    @test SuslinStability.elementary_matrix === elementary_matrix
    @test SuslinStability.realize_cohn_type === realize_cohn_type
    @test SuslinStability.realize_conjugate_elementary === realize_conjugate_elementary
end
