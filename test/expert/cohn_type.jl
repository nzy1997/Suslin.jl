using Test
using SuslinStability
using Oscar

struct ZeroBasedPair{T} <: AbstractVector{T}
    data::NTuple{2, T}
end

Base.size(::ZeroBasedPair) = (2,)
Base.axes(::ZeroBasedPair) = (0:1,)
Base.IndexStyle(::Type{<:ZeroBasedPair}) = IndexLinear()
Base.getindex(v::ZeroBasedPair, i::Int) = v.data[i + 1]

function product_of_factors(factors)
    R = base_ring(first(factors))
    product = identity_matrix(R, size(first(factors), 1))
    for factor in factors
        product *= factor
    end
    return product
end

function cohn_type_target(n::Int, i::Int, j::Int, a, v::AbstractVector, R)
    target = identity_matrix(R, n)
    for row in 1:n
        target[row, i] += a * v[row] * v[j]
        target[row, j] -= a * v[row] * v[i]
    end
    return target
end

@testset "cohn-type realization" begin
    F = QQ
    R, (a, v1, v2, v3, v4) = Oscar.polynomial_ring(F, ["a", "v1", "v2", "v3", "v4"])

    v = [v1, v2, v3]
    target = cohn_type_target(3, 1, 2, a, v, R)
    factors = SuslinStability.realize_cohn_type(3, 1, 2, a, v, R)

    @test target == product_of_factors(factors)

    vfull = [v1, v2, v3, v4]
    target4 = cohn_type_target(4, 2, 4, a, vfull, R)
    factors4 = SuslinStability.realize_cohn_type(4, 2, 4, a, vfull, R)

    @test target4 == product_of_factors(factors4)

    @test_throws ArgumentError SuslinStability.realize_cohn_type(2, 1, 2, a, v, R)
    @test_throws ArgumentError SuslinStability.realize_cohn_type(3, 1, 1, a, v, R)
    @test_throws ArgumentError SuslinStability.realize_cohn_type(3, 1, 2, a, [v1], R)
    @test_throws ArgumentError SuslinStability.realize_cohn_type(3, 1, 2, a, ZeroBasedPair((v1, v2)), R)
end
