using Test
using Suslin
using Oscar

const PARK_WOODBURN_POLYNOMIAL_NORMALITY_CATALOG_PATH =
    joinpath(@__DIR__, "..", "fixtures", "polynomial_normality_cases.jl")

if !isdefined(Main, :ParkWoodburnPolynomialNormalityFixtureCatalog)
    include(PARK_WOODBURN_POLYNOMIAL_NORMALITY_CATALOG_PATH)
end

function rank_one_product_of_factors(factors, R, n::Int)
    product = identity_matrix(R, n)
    for factor in factors
        product *= factor
    end
    return product
end

function rank_one_target_from_vectors(v, w, R, n::Int)
    target = identity_matrix(R, n)
    for row in 1:n, col in 1:n
        target[row, col] += v[row] * w[col]
    end
    return target
end

function expected_rank_one_cohn_coefficients(v, w, g)
    n = length(v)
    return [
        (; i, j, a = w[i] * g[j] - w[j] * g[i])
        for i in 1:(n - 1) for j in (i + 1):n
    ]
end

function replace_rank_one_certificate_field(cert, field::Symbol, value)
    names = fieldnames(typeof(cert))
    idx = findfirst(==(field), names)
    idx === nothing && throw(ArgumentError("certificate field $(field) not found"))
    values = ntuple(k -> k == idx ? value : getfield(cert, names[k]), length(names))
    return typeof(cert)(values...)
end

@testset "orthogonal rank-one normality replay certificate" begin
    cases = Main.ParkWoodburnPolynomialNormalityFixtureCatalog.cases_by_id()
    fixture = cases["pw-section2-orthogonal-rank-one-qq"]
    R = fixture.ring.object
    inputs = fixture.inputs
    n = nrows(fixture.target_matrix)

    cert = Suslin.realize_rank_one_normality_certificate(inputs.v, inputs.w, inputs.g, R)

    @test cert isa Suslin.RankOneNormalityCertificate
    @test cert.n == n
    @test cert.v == [R(value) for value in inputs.v]
    @test cert.w == [R(value) for value in inputs.w]
    @test cert.g == [R(value) for value in inputs.g]
    @test cert.orthogonality == zero(R)
    @test cert.bezout == one(R)
    @test cert.cohn_coefficients == expected_rank_one_cohn_coefficients(cert.v, cert.w, cert.g)
    @test length(cert.child_certificates) == count(entry -> entry.a != zero(R), cert.cohn_coefficients)
    @test all(Suslin.verify_cohn_type_certificate, cert.child_certificates)
    @test cert.factors == reduce(vcat, [child.factors for child in cert.child_certificates]; init = Any[])
    @test cert.product == rank_one_product_of_factors(cert.factors, R, n)
    @test cert.target == fixture.target_matrix
    @test cert.target == rank_one_target_from_vectors(cert.v, cert.w, R, n)
    @test cert.product == cert.target
    @test cert.verification.target_matches_product_ok
    @test Suslin.verify_rank_one_normality_certificate(cert)

    bad_w = copy(inputs.w)
    bad_w[end] += one(R)
    @test_throws ArgumentError Suslin.realize_rank_one_normality_certificate(inputs.v, bad_w, inputs.g, R)

    bad_g = copy(inputs.g)
    bad_g[1] += one(R)
    @test_throws ArgumentError Suslin.realize_rank_one_normality_certificate(inputs.v, inputs.w, bad_g, R)

    tampered_child_cert = Suslin.realize_rank_one_normality_certificate(inputs.v, inputs.w, inputs.g, R)
    tampered_child_cert.child_certificates[1].factors[1][1, 1] += one(R)
    @test !Suslin.verify_rank_one_normality_certificate(tampered_child_cert)

    tampered_factor_cert = Suslin.realize_rank_one_normality_certificate(inputs.v, inputs.w, inputs.g, R)
    tampered_factor_cert.factors[1][1, 1] += one(R)
    @test !Suslin.verify_rank_one_normality_certificate(tampered_factor_cert)

    tampered_verification = replace_rank_one_certificate_field(
        cert.verification,
        :target_matches_product_ok,
        false,
    )
    @test !Suslin.verify_rank_one_normality_certificate(
        replace_rank_one_certificate_field(cert, :verification, tampered_verification),
    )

    changed_table = copy(cert.cohn_coefficients)
    changed_table[1] = (; changed_table[1].i, changed_table[1].j, a = changed_table[1].a + one(R))
    @test !Suslin.verify_rank_one_normality_certificate(
        replace_rank_one_certificate_field(cert, :cohn_coefficients, changed_table),
    )
end
