using Test
using Suslin
using Oscar

function quillen_patch_product(factors, R, n)
    product = identity_matrix(R, n)
    for factor in factors
        product *= factor
    end
    return product
end

@testset "constructive quillen patching exact" begin
    R, (X, r, g) = Oscar.polynomial_ring(QQ, ["X", "r", "g"])
    n = 3
    target_entry = X + g + 1
    target = elementary_matrix(n, 1, 2, target_entry, R)
    contributions = [
        QuillenLocalContribution(
            LocalCertificate([1, 2], [r, r]),
            r,
            one(R),
            QuillenElementaryCorrection(1, 2, target_entry),
        ),
        QuillenLocalContribution(
            LocalCertificate([1, 2], [one(R) - r, one(R) - r]),
            one(R) - r,
            one(R),
            QuillenElementaryCorrection(1, 2, target_entry),
        ),
    ]

    patch = construct_quillen_patch(n, X, contributions; target)

    @test verify_quillen_patch(patch)
    @test patch.ring == R
    @test patch.size == n
    @test patch.substitution_variable == X
    @test patch.product == target
    @test quillen_patch_product(patch.factors, R, n) == target
    @test patch.verification.product == target
    @test patch.verification.target == target
    @test patch.verification.coverage_sum == one(R)
    @test patch.verification.coverage_ok
    @test patch.verification.product_ok
    @test patch.verification.denominator_data_ok

    base_matrix = matrix(R, [
        one(R)  X       g;
        zero(R) one(R)  r;
        zero(R) zero(R) one(R)
    ])
    @test base_matrix * patch.product == base_matrix * target

    @test length(patch.denominator_data) == 2
    @test [data.denominator for data in patch.denominator_data] == [r, one(R) - r]
    @test [data.coverage_multiplier for data in patch.denominator_data] == [one(R), one(R)]
    @test [contribution.certificate.indices for contribution in patch.local_contributions] == [[1, 2], [1, 2]]
    @test [contribution.certificate.denominators for contribution in patch.local_contributions] == [
        [r, r],
        [one(R) - r, one(R) - r],
    ]

    uncovered_contributions = [
        contributions[1],
        QuillenLocalContribution(
            LocalCertificate([1, 2], [one(R) - r, one(R) - r]),
            one(R) - r,
            r,
            QuillenElementaryCorrection(1, 2, target_entry),
        ),
    ]
    @test_throws ArgumentError construct_quillen_patch(n, X, uncovered_contributions; target)

    tampered_denominator_data = copy(patch.denominator_data)
    tampered_denominator_data[2] = QuillenDenominatorData(one(R) - r + X, one(R))
    tampered_patch = QuillenPatch(
        patch.ring,
        patch.size,
        patch.substitution_variable,
        tampered_denominator_data,
        patch.local_contributions,
        patch.factors,
        patch.product,
        patch.target,
        patch.verification,
    )
    @test !verify_quillen_patch(tampered_patch)

    certificate_mismatch = [
        QuillenLocalContribution(
            LocalCertificate([1, 2], [r + 1, r + 1]),
            r,
            one(R),
            QuillenElementaryCorrection(1, 2, target_entry),
        ),
        contributions[2],
    ]
    @test_throws ArgumentError construct_quillen_patch(n, X, certificate_mismatch; target)
    @test_throws ArgumentError construct_quillen_patch(n, one(R), contributions; target)

    L, (x, u) = suslin_laurent_polynomial_ring(QQ, ["x", "u"])
    laurent_target_entry = x^-1 + u
    laurent_target = elementary_matrix(n, 2, 3, laurent_target_entry, L)
    laurent_contributions = [
        QuillenLocalContribution(
            LocalCertificate([2, 3], [u, u]),
            u,
            one(L),
            QuillenElementaryCorrection(2, 3, laurent_target_entry),
        ),
        QuillenLocalContribution(
            LocalCertificate([2, 3], [one(L) - u, one(L) - u]),
            one(L) - u,
            one(L),
            QuillenElementaryCorrection(2, 3, laurent_target_entry),
        ),
    ]
    laurent_patch = construct_quillen_patch(n, x, laurent_contributions; target = laurent_target)
    @test verify_quillen_patch(laurent_patch)
    @test laurent_patch.product == laurent_target
    @test quillen_patch_product(laurent_patch.factors, L, n) == laurent_target
    @test laurent_patch.verification.coverage_sum == one(L)
end
