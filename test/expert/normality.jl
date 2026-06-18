using Test
using Suslin
using Oscar

function product_of_factors(factors)
    R = base_ring(first(factors))
    product = identity_matrix(R, size(first(factors), 1))
    for factor in factors
        product *= factor
    end
    return product
end

@testset "constructive normality" begin
    F = QQ
    R, (x,) = Oscar.polynomial_ring(F, ["x"])

    B = matrix(R, [
        1 0 0;
        x 1 0;
        0 1 1
    ])
    E = elementary_matrix(3, 1, 3, x + 1, R)

    factors = Suslin.realize_conjugate_elementary(B, 1, 3, x + 1)

    @test B * E * inv(B) == product_of_factors(factors)

    B4 = matrix(R, [
        1 0 0 0;
        x 1 0 0;
        0 1 1 0;
        0 0 x 1
    ])
    E4 = elementary_matrix(4, 1, 4, x + 1, R)
    factors4 = Suslin.realize_conjugate_elementary(B4, 1, 4, x + 1)

    @test B4 * E4 * inv(B4) == product_of_factors(factors4)

    zero_factors = Suslin.realize_conjugate_elementary(B, 1, 3, zero(R))
    @test isempty(zero_factors)

    singular = matrix(R, [
        1 0 0;
        0 0 0;
        0 0 1
    ])

    @test_throws ArgumentError Suslin.realize_conjugate_elementary(B, 1, 1, x + 1)
    @test_throws ArgumentError Suslin.realize_conjugate_elementary(matrix(R, [1 0; x 1]), 1, 2, x + 1)
    @test_throws ErrorException Suslin.realize_conjugate_elementary(singular, 1, 3, x + 1)
end
