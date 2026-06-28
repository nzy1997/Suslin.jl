using Test
using Suslin
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

@testset "elementary factor sequence analysis" begin
    L, (x, y) = suslin_laurent_polynomial_ring(QQ, ["x", "y"])
    laurent_factor = elementary_matrix(3, 1, 2, x^5 * y^-7 + x^-2 + one(L), L)
    diagonal_noise = identity_matrix(L, 3)
    diagonal_noise[1, 1] = x^99 * y^-99

    @test max_elementary_factor_monomial_degree([laurent_factor, diagonal_noise]) == 12
    @test total_elementary_factor_offdiagonal_monomials([laurent_factor, diagonal_noise]) == 3

    R, (u, v) = suslin_polynomial_ring(QQ, ["u", "v"])
    polynomial_factor = elementary_matrix(3, 2, 3, u^2 * v + u + 3, R)
    constant_factor = elementary_matrix(3, 3, 1, R(7), R)
    polynomial_diagonal_noise = identity_matrix(R, 3)
    polynomial_diagonal_noise[2, 2] = u^100

    @test max_elementary_factor_monomial_degree([
        polynomial_factor,
        constant_factor,
        polynomial_diagonal_noise,
    ]) == 3
    @test total_elementary_factor_offdiagonal_monomials([
        polynomial_factor,
        constant_factor,
        polynomial_diagonal_noise,
    ]) == 4
    @test max_elementary_factor_monomial_degree(Matrix{Any}[]) == 0
    @test total_elementary_factor_offdiagonal_monomials(Matrix{Any}[]) == 0
end
