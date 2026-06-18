using Test
using Suslin
using Oscar

@testset "suslin polynomial ring" begin
    R, vars = suslin_polynomial_ring(GF(2), ["x", "y"])
    @test ngens(R) == 2
    @test parent(vars[1]) === R
    @test length(vars) == 2
    @test string(vars[1]) == "x"
    @test string(vars[2]) == "y"
end
