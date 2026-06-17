using Test
using SuslinStability
using Oscar

function product_of_factors(factors, R, n)
    product = identity_matrix(R, n)
    for factor in factors
        product *= factor
    end
    return product
end

@testset "quillen induction scaffolding" begin
    F = QQ
    R, (X, r, g) = Oscar.polynomial_ring(F, ["X", "r", "g"])
    K = fraction_field(R)

    local_entries = [
        K((X + 1) // (r + 1)),
        K((g + 2) // (X^2 + 1)),
    ]
    certificate = SuslinStability.LocalCertificate(
        [1, 3],
        [denominator(entry) for entry in local_entries],
    )

    @test certificate.indices == [1, 3]
    @test certificate.denominators == [r + 1, X^2 + 1]
    @test SuslinStability.common_denominator_factor(local_entries) == (r + 1) * (X^2 + 1)
    @test SuslinStability.common_denominator_factor([local_entries[1], X + g]) == (r + 1)

    @test_throws ArgumentError SuslinStability.LocalCertificate([1], [r + 1, X^2 + 1])
    @test_throws ArgumentError SuslinStability.common_denominator_factor(typeof(local_entries)())

    A = matrix(R, [
        X^2 + r  g + 1  0;
        r * X    X + g  1;
        0        r      1
    ])
    shift = X + r^2 * g
    patched = SuslinStability.patched_substitution(A, X, r, 2, g)
    expected = matrix(R, [
        shift^2 + r  g + 1      0;
        r * shift    shift + g  1;
        0            r          1
    ])

    @test patched == expected

    patch_factor = SuslinStability.patched_substitution(
        elementary_matrix(3, 1, 2, X, R),
        X,
        r,
        2,
        g,
    )
    expected_factor = elementary_matrix(3, 1, 2, shift, R)
    realized = SuslinStability.realize_cohn_type(
        3,
        2,
        1,
        shift,
        [one(R), zero(R), zero(R)],
        R,
    )

    @test patch_factor == expected_factor
    @test product_of_factors(realized, R, 3) == expected_factor

    @test_throws ArgumentError SuslinStability.patched_substitution(A, X, r, -1, g)
    @test_throws ArgumentError SuslinStability.patched_substitution(A, one(R), r, 2, g)

    S, (y,) = Oscar.polynomial_ring(QQ, ["y"])
    @test_throws ArgumentError SuslinStability.patched_substitution(A, X, y, 2, g)
end
