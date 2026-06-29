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

function product_of_factors(factors, R, n::Int)
    product = identity_matrix(R, n)
    for factor in factors
        product *= factor
    end
    return product
end

function replace_conjugate_certificate_field(cert, field::Symbol, value)
    names = fieldnames(typeof(cert))
    idx = findfirst(==(field), names)
    idx === nothing && throw(ArgumentError("certificate field $(field) not found"))
    values = ntuple(k -> k == idx ? value : getfield(cert, names[k]), length(names))
    if cert isa NamedTuple
        return typeof(cert)(values)
    end
    return typeof(cert)(values...)
end

function assert_conjugate_certificate(cert, A, i::Int, j::Int, a)
    R = base_ring(A)
    n = nrows(A)
    inverse_A = inv(A)
    elementary = elementary_matrix(n, i, j, a, R)
    target = A * elementary * inverse_A
    coerced_a = R(a)
    expected_v = [A[row, i] for row in 1:n]
    expected_w = [coerced_a * inverse_A[j, col] for col in 1:n]
    expected_g = [inverse_A[i, col] for col in 1:n]

    @test cert isa Suslin.ConjugatedElementaryNormalityCertificate
    @test cert.n == n
    @test cert.A == A
    @test cert.i == i
    @test cert.j == j
    @test cert.a == coerced_a
    @test cert.ring == R
    @test cert.determinant == one(R)
    @test cert.inverse_A == inverse_A
    @test cert.elementary_matrix == elementary
    @test cert.conjugation_convention == :A_E_invA
    @test cert.conjugation_target == target
    @test cert.v == expected_v
    @test cert.w == expected_w
    @test cert.g == expected_g
    @test Suslin.verify_rank_one_normality_certificate(cert.rank_one_certificate)
    @test cert.rank_one_certificate.v == expected_v
    @test cert.rank_one_certificate.w == expected_w
    @test cert.rank_one_certificate.g == expected_g
    @test cert.rank_one_certificate.target == target
    @test cert.factors == cert.rank_one_certificate.factors
    @test cert.product == product_of_factors(cert.factors, R, n)
    @test cert.product == cert.conjugation_target
    @test cert.verification.target_matches_product_ok
    @test Suslin.verify_conjugate_elementary_certificate(cert)
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

    cert = Suslin.realize_conjugate_elementary_certificate(B, 1, 3, x + 1)
    assert_conjugate_certificate(cert, B, 1, 3, x + 1)
    @test Suslin.realize_conjugate_elementary(B, 1, 3, x + 1) == cert.factors

    B4 = matrix(R, [
        1 0 0 0;
        x 1 0 0;
        0 1 1 0;
        0 0 x 1
    ])
    E4 = elementary_matrix(4, 1, 4, x + 1, R)
    factors4 = Suslin.realize_conjugate_elementary(B4, 1, 4, x + 1)

    @test B4 * E4 * inv(B4) == product_of_factors(factors4)

    cert4 = Suslin.realize_conjugate_elementary_certificate(B4, 1, 4, x + 1)
    assert_conjugate_certificate(cert4, B4, 1, 4, x + 1)
    @test Suslin.realize_conjugate_elementary(B4, 1, 4, x + 1) == cert4.factors

    zero_factors = Suslin.realize_conjugate_elementary(B, 1, 3, zero(R))
    @test isempty(zero_factors)

    fixture_path = joinpath(@__DIR__, "..", "fixtures", "polynomial_normality_cases.jl")
    if !isdefined(Main, :ParkWoodburnPolynomialNormalityFixtureCatalog)
        include(fixture_path)
    end
    fixture = Main.ParkWoodburnPolynomialNormalityFixtureCatalog.cases_by_id()["pw-section2-conjugated-elementary-qq"]
    fixture_cert = Suslin.realize_conjugate_elementary_certificate(
        fixture.inputs.B,
        fixture.inputs.i,
        fixture.inputs.j,
        fixture.inputs.a,
    )
    assert_conjugate_certificate(
        fixture_cert,
        fixture.inputs.B,
        fixture.inputs.i,
        fixture.inputs.j,
        fixture.inputs.a,
    )
    @test fixture_cert.conjugation_target == fixture.target_matrix

    singular = matrix(R, [
        1 0 0;
        0 0 0;
        0 0 1
    ])

    @test_throws ArgumentError Suslin.realize_conjugate_elementary(B, 1, 1, x + 1)
    @test_throws ArgumentError Suslin.realize_conjugate_elementary(matrix(R, [1 0; x 1]), 1, 2, x + 1)
    @test_throws ErrorException Suslin.realize_conjugate_elementary(singular, 1, 3, x + 1)

    non_sl = matrix(R, [
        2 0 0;
        0 1 0;
        0 0 1
    ])
    @test_throws ArgumentError Suslin.realize_conjugate_elementary_certificate(non_sl, 1, 2, one(R))
    @test_throws ArgumentError Suslin.realize_conjugate_elementary_certificate(B, 1, 1, x + 1)

    L, (u,) = Oscar.laurent_polynomial_ring(QQ, ["u"])
    @test_throws ArgumentError Suslin.realize_conjugate_elementary_certificate(B, 1, 3, u)

    tampered_factor_cert = Suslin.realize_conjugate_elementary_certificate(B, 1, 3, x + 1)
    tampered_factor_cert.factors[1][1, 1] += one(R)
    @test !Suslin.verify_conjugate_elementary_certificate(tampered_factor_cert)

    tampered_convention_cert = replace_conjugate_certificate_field(cert, :conjugation_convention, :invA_E_A)
    @test !Suslin.verify_conjugate_elementary_certificate(tampered_convention_cert)

    tampered_verification = replace_conjugate_certificate_field(cert.verification, :target_matches_product_ok, false)
    tampered_verification_cert = replace_conjugate_certificate_field(cert, :verification, tampered_verification)
    @test !Suslin.verify_conjugate_elementary_certificate(tampered_verification_cert)
end
