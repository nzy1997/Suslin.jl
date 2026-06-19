using Test
using Suslin
using Oscar

const LAURENT_FIXTURE_CATALOG_PATH = joinpath(@__DIR__, "..", "fixtures", "laurent_cases.jl")
include(LAURENT_FIXTURE_CATALOG_PATH)

function laurent_product_of_factors(factors)
    isempty(factors) && throw(ArgumentError("factor list must be nonempty"))

    R = base_ring(first(factors))
    product = identity_matrix(R, size(first(factors), 1))
    for factor in factors
        product *= factor
    end
    return product
end

function laurent_cohn_type_target(n::Int, i::Int, j::Int, a, v::AbstractVector, R)
    target = identity_matrix(R, n)
    for row in 1:n
        target[row, i] += a * v[row] * v[j]
        target[row, j] -= a * v[row] * v[i]
    end
    return target
end

@testset "Laurent elementary core" begin
    catalog = LaurentFixtureCatalog.catalog()
    fixture_ring = catalog.ring.object
    fx, fy = catalog.ring.generators

    fixture_entry = only(filter(entry -> entry.id == "laurent-negative-exponent-normalization", catalog.cases))
    elementary_entry = fixture_entry.inputs.vector[1, 1]
    E = elementary_matrix(3, 1, 2, elementary_entry, fixture_ring)

    @test base_ring(E) === fixture_ring
    @test E[1, 2] == elementary_entry
    @test E[1, 1] == one(fixture_ring)
    @test E[2, 1] == zero(fixture_ring)
    @test parent(fx^-1 * fy) === fixture_ring

    R, (x, y) = suslin_laurent_polynomial_ring(QQ, ["x", "y"])
    a = x * y^-1 - y

    v3 = [one(R), x^-1, y]
    target3 = laurent_cohn_type_target(3, 1, 2, a, v3, R)
    factors3 = Suslin.realize_cohn_type(3, 1, 2, a, v3, R)

    @test target3 == laurent_product_of_factors(factors3)
    @test verify_factorization(target3, factors3)
    @test all(factor -> base_ring(factor) === R, factors3)

    v4 = [one(R), x^-1, y, x * y^-1]
    target4 = laurent_cohn_type_target(4, 2, 4, a, v4, R)
    factors4 = Suslin.realize_cohn_type(4, 2, 4, a, v4, R)

    @test target4 == laurent_product_of_factors(factors4)
    @test verify_factorization(target4, factors4)
    @test all(factor -> base_ring(factor) === R, factors4)

    B = matrix(R, [
        1      0  0;
        x^-1   1  0;
        y      x  1
    ])
    Binv = matrix(R, [
        1          0   0;
        -x^-1      1   0;
        1 - y     -x   1
    ])
    normality_entry = x + y^-1
    normality_elementary = elementary_matrix(3, 1, 3, normality_entry, R)
    normality_target = B * normality_elementary * Binv
    normality_factors = Suslin.realize_conjugate_elementary(B, 1, 3, normality_entry)

    @test B * Binv == identity_matrix(R, 3)
    @test Binv * B == identity_matrix(R, 3)
    @test normality_target == laurent_product_of_factors(normality_factors)
    @test verify_factorization(normality_target, normality_factors)
    @test all(factor -> base_ring(factor) === R, normality_factors)

    S, (u, _) = suslin_laurent_polynomial_ring(QQ, ["u", "v"])
    wrong_parent_factor = elementary_matrix(3, 1, 2, u, S)

    @test_throws ArgumentError verify_factorization(identity_matrix(R, 3), [wrong_parent_factor])
end
