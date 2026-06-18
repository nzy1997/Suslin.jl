using Test
using Suslin
using Oscar

@testset "suslin Laurent polynomial ring" begin
    R, vars = suslin_laurent_polynomial_ring(GF(2), ["x", "y"])
    x, y = vars

    @test length(vars) == 2
    @test ngens(R) == 2
    @test parent(x) === R
    @test parent(y) === R
    @test string(x) == "x"
    @test string(y) == "y"

    x_inv = x^-1
    mixed = x^-1 * y
    @test parent(x_inv) === R
    @test parent(mixed) === R
    @test string(x_inv) == "x^-1"
    @test string(mixed) == "x^-1*y"
end

@testset "Laurent validators" begin
    R, (x, y) = suslin_laurent_polynomial_ring(GF(2), ["x", "y"])
    S, (u, v) = suslin_laurent_polynomial_ring(GF(2), ["u", "v"])
    P, (px, py) = suslin_polynomial_ring(GF(2), ["x", "y"])

    @test Suslin._is_laurent_polynomial_ring(R)
    @test !Suslin._is_laurent_polynomial_ring(P)
    @test Suslin._require_laurent_polynomial_ring(R) === R
    @test_throws ArgumentError Suslin._require_laurent_polynomial_ring(P)

    value = x^-1 * y
    @test Suslin._require_laurent_element(value) === value
    @test Suslin._require_laurent_element(value, R) === value
    @test_throws ArgumentError Suslin._require_laurent_element(u, R)
    @test_throws ArgumentError Suslin._require_laurent_element(px)
    @test_throws ArgumentError Suslin._require_laurent_element(px, R)

    @test Suslin._require_same_laurent_parent([x, value, y^-1]) === R
    @test_throws ArgumentError Suslin._require_same_laurent_parent([x, u])
    @test_throws ArgumentError Suslin._require_same_laurent_parent([x, px])
    @test_throws ArgumentError Suslin._require_same_laurent_parent([])
end
