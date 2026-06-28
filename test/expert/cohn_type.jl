using Test
using Suslin
using Oscar

const PARK_WOODBURN_POLYNOMIAL_NORMALITY_CATALOG_PATH =
    joinpath(@__DIR__, "..", "fixtures", "polynomial_normality_cases.jl")

if !isdefined(Main, :ParkWoodburnPolynomialNormalityFixtureCatalog)
    include(PARK_WOODBURN_POLYNOMIAL_NORMALITY_CATALOG_PATH)
end

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

function replace_cohn_certificate_field(cert, field::Symbol, value)
    names = fieldnames(typeof(cert))
    idx = findfirst(==(field), names)
    idx === nothing && throw(ArgumentError("certificate field $(field) not found"))
    values = ntuple(k -> k == idx ? value : getfield(cert, names[k]), length(names))
    return typeof(cert)(values...)
end

@testset "cohn-type realization" begin
    F = QQ
    R, (a, v1, v2, v3, v4) = Oscar.polynomial_ring(F, ["a", "v1", "v2", "v3", "v4"])

    v = [v1, v2, v3]
    target = cohn_type_target(3, 1, 2, a, v, R)
    factors = Suslin.realize_cohn_type(3, 1, 2, a, v, R)

    @test target == product_of_factors(factors)

    vfull = [v1, v2, v3, v4]
    target4 = cohn_type_target(4, 2, 4, a, vfull, R)
    factors4 = Suslin.realize_cohn_type(4, 2, 4, a, vfull, R)

    @test target4 == product_of_factors(factors4)

    @test_throws ArgumentError Suslin.realize_cohn_type(2, 1, 2, a, v, R)
    @test_throws ArgumentError Suslin.realize_cohn_type(3, 1, 1, a, v, R)
    @test_throws ArgumentError Suslin.realize_cohn_type(3, 1, 2, a, [v1], R)
    @test_throws ArgumentError Suslin.realize_cohn_type(3, 1, 2, a, ZeroBasedPair((v1, v2)), R)
end

@testset "cohn-type replay certificate" begin
    cases = Main.ParkWoodburnPolynomialNormalityFixtureCatalog.cases_by_id()
    fixture = cases["pw-section2-cohn-type-qq"]
    R = fixture.ring.object
    inputs = fixture.inputs
    n = nrows(fixture.target_matrix)

    cert = Suslin.realize_cohn_type_certificate(n, inputs.i, inputs.j, inputs.a, inputs.v, R)

    @test cert isa Suslin.CohnTypeRealizationCertificate
    @test cert.n == n
    @test cert.i == inputs.i
    @test cert.j == inputs.j
    @test cert.a == R(inputs.a)
    @test cert.v == [R(value) for value in inputs.v]
    @test cert.target == fixture.target_matrix
    @test cert.product == cert.target
    @test cert.verification.target_matches_product_ok
    @test Suslin.verify_cohn_type_certificate(cert)
    @test Suslin.realize_cohn_type(n, inputs.i, inputs.j, inputs.a, inputs.v, R) == cert.factors

    tampered_factor_cert = Suslin.realize_cohn_type_certificate(n, inputs.i, inputs.j, inputs.a, inputs.v, R)
    tampered_factor_cert.factors[1][1, 1] += one(R)
    @test !Suslin.verify_cohn_type_certificate(tampered_factor_cert)

    tampered_target_cert = Suslin.realize_cohn_type_certificate(n, inputs.i, inputs.j, inputs.a, inputs.v, R)
    tampered_target_cert.target[1, 1] += one(R)
    @test !Suslin.verify_cohn_type_certificate(tampered_target_cert)

    changed_a_cert = replace_cohn_certificate_field(cert, :a, cert.a + one(R))
    @test !Suslin.verify_cohn_type_certificate(changed_a_cert)

    L, (lx, ly) = Oscar.laurent_polynomial_ring(QQ, ["lx", "ly"])
    @test_throws ArgumentError Suslin.realize_cohn_type_certificate(
        3,
        1,
        2,
        lx,
        [one(L), lx, ly],
        L,
    )
end
