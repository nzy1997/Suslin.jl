using Test
using SuslinStability
using Oscar

@testset "elementary matrix" begin
    F = GF(2)
    R, (x,) = Oscar.polynomial_ring(F, ["x"])
    E = elementary_matrix(3, 1, 2, x + 1, R)
    expected = identity_matrix(R, 3)
    expected[1, 2] = x + 1

    S, (y,) = Oscar.polynomial_ring(F, ["y"])

    @test size(E) == (3, 3)
    @test E == expected
    @test E[1, 2] == x + 1
    @test det(E) == one(R)
    @test_throws ArgumentError elementary_matrix(3, 1, 1, x + 1, R)
    @test_throws ArgumentError elementary_matrix(3, 1, 2, y + 1, R)
end
