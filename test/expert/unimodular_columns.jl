using Test
using SuslinStability
using Oscar

function unimodular_product_of_factors(factors, R, n)
    product = identity_matrix(R, n)
    for factor in factors
        product *= factor
    end
    return product
end

function apply_factors_to_column(factors, v, R)
    n = length(v)
    column = matrix(R, n, 1, collect(v))
    return unimodular_product_of_factors(factors, R, n) * column
end

@testset "unimodular column reduction" begin
    F2 = GF(2)
    R2, (x, y) = Oscar.polynomial_ring(F2, ["x", "y"])

    v = [x, y, one(R2)]
    factors = SuslinStability.reduce_unimodular_column(v, R2)
    reduced = apply_factors_to_column(factors, v, R2)

    @test SuslinStability.is_unimodular_column(v, R2)
    @test reduced == matrix(R2, 3, 1, [zero(R2), zero(R2), one(R2)])

    vperm = [x, one(R2), y]
    factors_perm = SuslinStability.reduce_unimodular_column(vperm, R2)
    reduced_perm = apply_factors_to_column(factors_perm, vperm, R2)

    @test SuslinStability.is_unimodular_column(vperm, R2)
    @test reduced_perm == matrix(R2, 3, 1, [zero(R2), zero(R2), one(R2)])

    F3 = GF(3)
    R3, (t,) = Oscar.polynomial_ring(F3, ["t"])

    v3 = [t + one(R3), t, R3(2)]
    factors3 = SuslinStability.reduce_unimodular_column(v3, R3)
    reduced3 = apply_factors_to_column(factors3, v3, R3)

    @test SuslinStability.is_unimodular_column(v3, R3)
    @test reduced3 == matrix(R3, 3, 1, [zero(R3), zero(R3), one(R3)])

    v3perm = [t + one(R3), R3(2), t]
    factors3perm = SuslinStability.reduce_unimodular_column(v3perm, R3)
    reduced3perm = apply_factors_to_column(factors3perm, v3perm, R3)

    @test SuslinStability.is_unimodular_column(v3perm, R3)
    @test reduced3perm == matrix(R3, 3, 1, [zero(R3), zero(R3), one(R3)])

    monic_gap = [x, y, x + one(R2)]
    monic_factors = SuslinStability.reduce_unimodular_column(monic_gap, R2)
    monic_reduced = apply_factors_to_column(monic_factors, monic_gap, R2)

    @test SuslinStability.is_unimodular_column(monic_gap, R2)
    @test monic_reduced == matrix(R2, 3, 1, [zero(R2), zero(R2), one(R2)])

    variable_change_gap = [x + y^2, x * y + x + one(R2), x^2 + x * y + y + one(R2)]
    variable_change_factors = SuslinStability.reduce_unimodular_column(variable_change_gap, R2)
    variable_change_reduced = apply_factors_to_column(variable_change_factors, variable_change_gap, R2)

    @test SuslinStability.is_unimodular_column(variable_change_gap, R2)
    @test variable_change_reduced == matrix(R2, 3, 1, [zero(R2), zero(R2), one(R2)])

    heuristic_gap = [zero(R2), x^2, x * y + one(R2)]
    @test SuslinStability.is_unimodular_column(heuristic_gap, R2)
    @test_throws ArgumentError SuslinStability.reduce_unimodular_column(heuristic_gap, R2)

    bad = [x, y, x * y]

    @test !SuslinStability.is_unimodular_column(bad, R2)
    @test_throws ArgumentError SuslinStability.reduce_unimodular_column(bad, R2)
end
