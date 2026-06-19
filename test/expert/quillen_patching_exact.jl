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

    certificate_pair_mismatch = [
        QuillenLocalContribution(
            LocalCertificate([1, 2, 3], [r + 1, r, r]),
            r,
            one(R),
            QuillenElementaryCorrection(1, 3, X + r),
        ),
        QuillenLocalContribution(
            LocalCertificate([1, 3], [-r, -r]),
            -r,
            one(R),
            QuillenElementaryCorrection(1, 3, X + r),
        ),
    ]
    pair_target = elementary_matrix(n, 1, 3, X + r, R)
    @test_throws ArgumentError construct_quillen_patch(n, X, certificate_pair_mismatch; target = pair_target)

    tampered_local_contributions = copy(patch.local_contributions)
    tampered_local_contributions[1] = QuillenLocalContribution(
        LocalCertificate([1, 2, 3], [r + 1, r, r]),
        tampered_local_contributions[1].denominator,
        tampered_local_contributions[1].coverage_multiplier,
        QuillenElementaryCorrection(1, 3, tampered_local_contributions[1].correction.entry),
    )
    tampered_certificate_patch = QuillenPatch(
        patch.ring,
        patch.size,
        patch.substitution_variable,
        patch.denominator_data,
        tampered_local_contributions,
        patch.factors,
        patch.product,
        patch.target,
        patch.verification,
    )
    @test !verify_quillen_patch(tampered_certificate_patch)

    tampered_verification = QuillenPatchVerification(
        patch.verification.denominator_data_ok,
        patch.verification.coverage_sum,
        patch.verification.coverage_ok,
        patch.verification.product + elementary_matrix(n, 1, 3, one(R), R),
        patch.verification.target,
        patch.verification.product_ok,
    )
    tampered_verification_patch = QuillenPatch(
        patch.ring,
        patch.size,
        patch.substitution_variable,
        patch.denominator_data,
        patch.local_contributions,
        patch.factors,
        patch.product,
        patch.target,
        tampered_verification,
    )
    @test !verify_quillen_patch(tampered_verification_patch)

    RR, (Y, s) = Oscar.polynomial_ring(RealField(), ["Y", "s"])
    inexact_target_entry = Y + 1
    inexact_target = elementary_matrix(n, 1, 2, inexact_target_entry, RR)
    inexact_contributions = [
        QuillenLocalContribution(
            LocalCertificate([1, 2], [s, s]),
            s,
            one(RR),
            QuillenElementaryCorrection(1, 2, inexact_target_entry),
        ),
        QuillenLocalContribution(
            LocalCertificate([1, 2], [one(RR) - s, one(RR) - s]),
            one(RR) - s,
            one(RR),
            QuillenElementaryCorrection(1, 2, inexact_target_entry),
        ),
    ]
    @test_throws ArgumentError construct_quillen_patch(n, Y, inexact_contributions; target = inexact_target)

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

    S, t = Oscar.polynomial_ring(QQ, "t")
    univariate_target_entry = t + 1
    univariate_target = elementary_matrix(n, 1, 2, univariate_target_entry, S)
    univariate_contributions = [
        QuillenLocalContribution(
            LocalCertificate([1, 2], [t, t]),
            t,
            one(S),
            QuillenElementaryCorrection(1, 2, univariate_target_entry),
        ),
        QuillenLocalContribution(
            LocalCertificate([1, 2], [one(S) - t, one(S) - t]),
            one(S) - t,
            one(S),
            QuillenElementaryCorrection(1, 2, univariate_target_entry),
        ),
    ]
    univariate_patch = construct_quillen_patch(n, t, univariate_contributions; target = univariate_target)
    @test verify_quillen_patch(univariate_patch)
    @test univariate_patch.product == univariate_target
    @test quillen_patch_product(univariate_patch.factors, S, n) == univariate_target
end
