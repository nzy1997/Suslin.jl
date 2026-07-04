using Test
using Suslin
using Oscar

const STEINBERG_OPTIMIZATION_FIXTURE_PATH =
    joinpath(@__DIR__, "..", "fixtures", "steinberg_optimization_cases.jl")

if !isdefined(Main, :SteinbergOptimizationFixtureCatalog)
    include(STEINBERG_OPTIMIZATION_FIXTURE_PATH)
end

@testset "Steinberg canonical elementary factor records" begin
    R, (x, y) = Oscar.polynomial_ring(QQ, ["x", "y"])
    coefficient = x + y + one(R)

    elementary_factor = elementary_matrix(3, 1, 2, coefficient, R)
    elementary_record = Suslin._canonical_elementary_factor_record(elementary_factor)

    @test elementary_record.kind == :elementary
    @test elementary_record.n == 3
    @test Suslin._same_base_ring(elementary_record.ring, R)
    @test elementary_record.row == 1
    @test elementary_record.col == 2
    @test elementary_record.coefficient == coefficient
    @test Suslin._elementary_factor_record_matrix(elementary_record) == elementary_factor

    zero_elementary_factor = elementary_matrix(3, 1, 2, zero(R), R)
    identity_record = Suslin._canonical_elementary_factor_record(zero_elementary_factor)

    @test identity_record.kind == :identity
    @test identity_record.n == 3
    @test Suslin._same_base_ring(identity_record.ring, R)
    @test Suslin._elementary_factor_record_matrix(identity_record) == zero_elementary_factor

    nonsquare_factor = zero_matrix(R, 2, 3)
    bad_diagonal = identity_matrix(R, 3)
    bad_diagonal[2, 2] = x
    two_offdiagonal = identity_matrix(R, 3)
    two_offdiagonal[1, 2] = x
    two_offdiagonal[2, 3] = y

    @test_throws DimensionMismatch Suslin._canonical_elementary_factor_record(nonsquare_factor)
    @test_throws ArgumentError Suslin._canonical_elementary_factor_record(bad_diagonal)
    @test_throws ArgumentError Suslin._canonical_elementary_factor_record(two_offdiagonal)
    @test_throws ArgumentError Suslin._elementary_factor_record_matrix((; kind = :unsupported, n = 3, ring = R))
end

function _assert_valid_commutator_certificate(
    certificate,
    original_factors,
    expected_optimized_factors,
    expected_rule_names,
)
    @test certificate isa Suslin.SteinbergOptimizationCertificate
    @test Suslin._verify_steinberg_optimization_certificate(certificate)
    @test certificate.original_factors == original_factors
    @test certificate.optimized_factors == expected_optimized_factors
    @test certificate.original_product == certificate.optimized_product
    @test certificate.comparison_summary.products_equal
    @test certificate.comparison_summary.original_factor_count == length(original_factors)
    @test certificate.comparison_summary.optimized_factor_count == length(expected_optimized_factors)
    @test certificate.comparison_summary.factor_count_delta ==
          length(expected_optimized_factors) - length(original_factors)
    @test [rewrite.rule_name for rewrite in certificate.applied_rewrites] == expected_rule_names
    @test all(rewrite -> rewrite.metadata.local_products_equal, certificate.applied_rewrites)
    return certificate
end

@testset "Steinberg conservative commutator optimizer positives" begin
    R, (x, y) = Oscar.polynomial_ring(QQ, ["x", "y"])
    a = x + one(R)
    b = y

    forward_factors = [
        elementary_matrix(3, 1, 2, a, R),
        elementary_matrix(3, 2, 3, b, R),
        elementary_matrix(3, 1, 2, -a, R),
        elementary_matrix(3, 2, 3, -b, R),
    ]
    forward_expected = [
        elementary_matrix(3, 1, 3, a * b, R),
    ]
    forward_certificate =
        Suslin._steinberg_commutator_rewrite_optimization_certificate(forward_factors)
    _assert_valid_commutator_certificate(
        forward_certificate,
        forward_factors,
        forward_expected,
        [:commutator_forward],
    )
    @test forward_certificate.applied_rewrites[1].original_span == (start = 1, stop = 4)
    @test forward_certificate.applied_rewrites[1].optimized_span == (start = 1, stop = 1)
    @test forward_certificate.applied_rewrites[1].metadata.indices == (i = 1, j = 2, l = 3)

    reverse_factors = [
        elementary_matrix(3, 2, 3, x, R),
        elementary_matrix(3, 1, 2, y + one(R), R),
        elementary_matrix(3, 2, 3, -x, R),
        elementary_matrix(3, 1, 2, -(y + one(R)), R),
    ]
    reverse_expected = [
        elementary_matrix(3, 1, 3, -(x * (y + one(R))), R),
    ]
    reverse_certificate =
        Suslin._steinberg_commutator_rewrite_optimization_certificate(reverse_factors)
    _assert_valid_commutator_certificate(
        reverse_certificate,
        reverse_factors,
        reverse_expected,
        [:commutator_reverse],
    )
    @test reverse_certificate.applied_rewrites[1].metadata.indices == (l = 1, i = 2, j = 3)

    disjoint_factors = [
        elementary_matrix(4, 1, 2, x, R),
        elementary_matrix(4, 3, 4, y + one(R), R),
        elementary_matrix(4, 1, 2, -x, R),
        elementary_matrix(4, 3, 4, -(y + one(R)), R),
    ]
    disjoint_certificate =
        Suslin._steinberg_commutator_rewrite_optimization_certificate(disjoint_factors)
    _assert_valid_commutator_certificate(
        disjoint_certificate,
        disjoint_factors,
        typeof(first(disjoint_factors))[],
        [:disjoint_commutator_identity],
    )
    @test disjoint_certificate.applied_rewrites[1].optimized_span == (start = 1, stop = 0)
    @test disjoint_certificate.applied_rewrites[1].metadata.indices == (i = 1, j = 2, l = 3, p = 4)
end

@testset "Steinberg commutator optimizer fixture catalog positives" begin
    entries = SteinbergOptimizationFixtureCatalog.cases_by_id()
    for id in (
        "steinberg-commutator-forward-qq",
        "steinberg-commutator-reverse-qq",
        "steinberg-disjoint-commutator-identity-qq",
    )
        entry = entries[id]
        original_factors = collect(entry.factors)
        expected_factors = collect(entry.expected_rewrite_factors)
        certificate =
            Suslin._steinberg_commutator_rewrite_optimization_certificate(original_factors)

        _assert_valid_commutator_certificate(
            certificate,
            original_factors,
            expected_factors,
            [entry.rule_name],
        )
        @test certificate.optimized_product == entry.rewritten_product
        @test certificate.original_product == entry.original_product
    end
end

@testset "Steinberg commutator optimizer negative controls" begin
    R, (x, y) = Oscar.polynomial_ring(QQ, ["x", "y"])
    a = x + one(R)
    b = y + one(R)

    reordered_factors = [
        elementary_matrix(3, 1, 2, a, R),
        elementary_matrix(3, 1, 2, -a, R),
        elementary_matrix(3, 2, 3, b, R),
        elementary_matrix(3, 2, 3, -b, R),
    ]
    reordered_certificate =
        Suslin._steinberg_commutator_rewrite_optimization_certificate(reordered_factors)
    _assert_valid_commutator_certificate(reordered_certificate, reordered_factors, reordered_factors, Symbol[])

    wrong_inverse_factors = [
        elementary_matrix(3, 1, 2, a, R),
        elementary_matrix(3, 2, 3, b, R),
        elementary_matrix(3, 1, 2, -(a + one(R)), R),
        elementary_matrix(3, 2, 3, -b, R),
    ]
    wrong_inverse_certificate =
        Suslin._steinberg_commutator_rewrite_optimization_certificate(wrong_inverse_factors)
    _assert_valid_commutator_certificate(
        wrong_inverse_certificate,
        wrong_inverse_factors,
        wrong_inverse_factors,
        Symbol[],
    )

    invalid_forward_indices = [
        elementary_matrix(3, 1, 2, a, R),
        elementary_matrix(3, 2, 1, b, R),
        elementary_matrix(3, 1, 2, -a, R),
        elementary_matrix(3, 2, 1, -b, R),
    ]
    invalid_forward_certificate =
        Suslin._steinberg_commutator_rewrite_optimization_certificate(invalid_forward_indices)
    _assert_valid_commutator_certificate(
        invalid_forward_certificate,
        invalid_forward_indices,
        invalid_forward_indices,
        Symbol[],
    )

    invalid_disjoint_indices = [
        elementary_matrix(3, 1, 2, a, R),
        elementary_matrix(3, 3, 1, b, R),
        elementary_matrix(3, 1, 2, -a, R),
        elementary_matrix(3, 3, 1, -b, R),
    ]
    invalid_disjoint_certificate =
        Suslin._steinberg_commutator_rewrite_optimization_certificate(invalid_disjoint_indices)
    _assert_valid_commutator_certificate(
        invalid_disjoint_certificate,
        invalid_disjoint_indices,
        invalid_disjoint_indices,
        Symbol[],
    )
end

@testset "Steinberg optimization certificate replay" begin
    entries = SteinbergOptimizationFixtureCatalog.cases_by_id()
    entry = entries["steinberg-same-position-merge-qq"]
    original_factors = collect(entry.factors)

    certificate = Suslin._steinberg_optimization_certificate(original_factors, copy(original_factors), ())
    summary = certificate.comparison_summary

    @test certificate isa Suslin.SteinbergOptimizationCertificate
    @test Suslin._verify_steinberg_optimization_certificate(certificate)
    @test certificate.original_product == certificate.optimized_product
    @test summary.original_product == certificate.original_product
    @test summary.optimized_product == certificate.optimized_product
    @test summary.products_equal
    @test summary.verification_status
    @test isempty(summary.applied_rewrites)
    @test summary.original_factor_count == length(original_factors)
    @test summary.optimized_factor_count == length(original_factors)
    @test summary.factor_count_delta == 0
    @test summary.original_metrics.max_elementary_factor_monomial_degree ==
          max_elementary_factor_monomial_degree(original_factors)
    @test summary.optimized_metrics.max_elementary_factor_monomial_degree ==
          max_elementary_factor_monomial_degree(original_factors)
    @test summary.original_metrics.total_elementary_factor_offdiagonal_monomials ==
          total_elementary_factor_offdiagonal_monomials(original_factors)
    @test summary.optimized_metrics.total_elementary_factor_offdiagonal_monomials ==
          total_elementary_factor_offdiagonal_monomials(original_factors)
    @test certificate.verification.products_equal
    @test certificate.verification.overall_ok

    R = base_ring(first(original_factors))
    n = nrows(first(original_factors))
    tampered_optimized = copy(original_factors)
    tampered_optimized[1] = elementary_matrix(n, 1, 3, one(R), R)
    tampered_certificate =
        Suslin._steinberg_optimization_certificate(original_factors, tampered_optimized, ())
    @test !Suslin._verify_steinberg_optimization_certificate(tampered_certificate)
    @test !tampered_certificate.comparison_summary.products_equal
    @test !tampered_certificate.comparison_summary.verification_status

    stale_rule_log = [(
        rule_name = :same_position_merge,
        original_factor_count = length(original_factors),
        optimized_factor_count = length(original_factors) + 1,
    )]
    stale_log_certificate = Suslin.SteinbergOptimizationCertificate(
        certificate.original_factors,
        certificate.optimized_factors,
        stale_rule_log,
        certificate.comparison_summary,
        certificate.original_product,
        certificate.optimized_product,
        certificate.verification,
    )
    @test !Suslin._verify_steinberg_optimization_certificate(stale_log_certificate)

    impossible_same_delta_log = [(
        rule_name = :same_position_merge,
        original_factor_count = length(original_factors) + 100,
        optimized_factor_count = length(original_factors) + 100,
    )]
    impossible_same_delta_certificate = Suslin._steinberg_optimization_certificate(
        original_factors,
        copy(original_factors),
        impossible_same_delta_log,
    )
    @test !Suslin._verify_steinberg_optimization_certificate(impossible_same_delta_certificate)

    span_rule_log = [(
        rule_name = :span_replay_probe,
        original_factor_count = 1,
        optimized_factor_count = 1,
        original_span = (start = 1, stop = 1),
        optimized_span = (start = 1, stop = 1),
        metadata = (source = :test,),
    )]
    span_certificate = Suslin._steinberg_optimization_certificate(
        original_factors,
        copy(original_factors),
        span_rule_log,
    )
    @test Suslin._verify_steinberg_optimization_certificate(span_certificate)
    @test span_certificate.applied_rewrites[1].metadata.source == :test

    empty_optimized_rule_log = [(
        rule_name = :empty_optimized_probe,
        original_factor_count = length(original_factors),
        optimized_factor_count = 0,
        original_span = (start = 1, stop = length(original_factors)),
        optimized_span = (start = 1, stop = 0),
    )]
    empty_optimized_certificate = Suslin._steinberg_optimization_certificate(
        original_factors,
        original_factors[1:0],
        empty_optimized_rule_log,
    )
    @test !Suslin._verify_steinberg_optimization_certificate(empty_optimized_certificate)
    @test empty_optimized_certificate.comparison_summary.optimized_factor_count == 0

    malformed_certificate = Suslin.SteinbergOptimizationCertificate(
        Any[],
        Any[],
        Any[],
        nothing,
        nothing,
        nothing,
        nothing,
    )
    @test !Suslin._verify_steinberg_optimization_certificate(malformed_certificate)

    @test_throws ArgumentError Suslin._steinberg_optimization_certificate(
        original_factors,
        copy(original_factors),
        [(original_factor_count = 1, optimized_factor_count = 1)],
    )
    @test_throws ArgumentError Suslin._steinberg_optimization_certificate(
        original_factors,
        copy(original_factors),
        [(
            rule_name = "span_replay_probe",
            original_factor_count = 1,
            optimized_factor_count = 1,
        )],
    )
    @test_throws ArgumentError Suslin._steinberg_optimization_certificate(
        original_factors,
        copy(original_factors),
        [(
            rule_name = :span_replay_probe,
            original_factor_count = 1.0,
            optimized_factor_count = 1,
        )],
    )
    @test_throws ArgumentError Suslin._steinberg_optimization_certificate(
        original_factors,
        copy(original_factors),
        [(
            rule_name = :span_replay_probe,
            original_factor_count = -1,
            optimized_factor_count = 1,
        )],
    )
    @test_throws ArgumentError Suslin._steinberg_optimization_certificate(
        original_factors,
        copy(original_factors),
        [(
            rule_name = :span_replay_probe,
            original_factor_count = 1,
            optimized_factor_count = 1,
            original_span = (start = 1,),
        )],
    )
    @test_throws ArgumentError Suslin._steinberg_optimization_certificate(
        original_factors,
        copy(original_factors),
        [(
            rule_name = :span_replay_probe,
            original_factor_count = 1,
            optimized_factor_count = 1,
            original_span = (start = 1.0, stop = 1),
        )],
    )
    @test_throws ArgumentError Suslin._steinberg_optimization_certificate(
        original_factors,
        copy(original_factors),
        [(
            rule_name = :span_replay_probe,
            original_factor_count = 1,
            optimized_factor_count = 1,
            original_span = (start = 0, stop = 0),
        )],
    )
    @test_throws ArgumentError Suslin._steinberg_optimization_certificate(
        original_factors,
        copy(original_factors),
        [(
            rule_name = :span_replay_probe,
            original_factor_count = 1,
            optimized_factor_count = 1,
            original_span = (start = 1, stop = 2),
        )],
    )
    @test_throws ArgumentError Suslin._steinberg_optimization_certificate(
        original_factors,
        copy(original_factors),
        [(
            rule_name = :span_replay_probe,
            original_factor_count = 1,
            optimized_factor_count = 1,
            optimized_span = (start = 1, stop = 2),
        )],
    )

    R_alt, (u, v) = Oscar.polynomial_ring(QQ, ["u", "v"])
    mixed_ring_factors = copy(original_factors)
    mixed_ring_factors[2] = elementary_matrix(n, 1, 2, u + one(R_alt), R_alt)
    @test_throws ArgumentError Suslin._steinberg_optimization_certificate(
        original_factors,
        mixed_ring_factors,
        (),
    )
end

@testset "Steinberg optimization certificate accepts univariate ordinary polynomial rings" begin
    R, x = Oscar.polynomial_ring(QQ, "x")
    original_factors = [
        elementary_matrix(3, 1, 2, x + one(R), R),
        elementary_matrix(3, 2, 3, x^2 + one(R), R),
    ]

    certificate = Suslin._steinberg_optimization_certificate(original_factors, copy(original_factors), ())

    @test certificate isa Suslin.SteinbergOptimizationCertificate
    @test Suslin._verify_steinberg_optimization_certificate(certificate)
    @test certificate.original_product == certificate.optimized_product
    @test certificate.comparison_summary.original_factor_count == length(original_factors)
    @test certificate.comparison_summary.optimized_factor_count == length(original_factors)
    @test certificate.comparison_summary.factor_count_delta == 0
    @test isempty(certificate.applied_rewrites)
end

@testset "Steinberg adjacent identity merge cancellation optimizer" begin
    R, (x, y) = Oscar.polynomial_ring(QQ, ["x", "y"])
    inverse_coefficient = x * y + one(R)
    original_factors = [
        elementary_matrix(3, 1, 2, zero(R), R),
        elementary_matrix(3, 2, 3, inverse_coefficient, R),
        elementary_matrix(3, 2, 3, -inverse_coefficient, R),
        elementary_matrix(3, 1, 3, x, R),
        elementary_matrix(3, 1, 3, y + one(R), R),
    ]

    certificate = Suslin._steinberg_adjacent_rewrite_optimization_certificate(original_factors)

    @test certificate isa Suslin.SteinbergOptimizationCertificate
    @test Suslin._verify_steinberg_optimization_certificate(certificate)
    @test length(certificate.optimized_factors) < length(original_factors)
    @test certificate.comparison_summary.original_factor_count == length(original_factors)
    @test certificate.comparison_summary.optimized_factor_count == 1
    @test certificate.comparison_summary.factor_count_delta == -4
    @test certificate.original_product == certificate.optimized_product
    @test certificate.comparison_summary.products_equal
    @test certificate.optimized_factors == [
        elementary_matrix(3, 1, 3, x + y + one(R), R),
    ]
    @test [rewrite.rule_name for rewrite in certificate.applied_rewrites] == [
        :identity_removal,
        :inverse_cancellation,
        :same_position_merge,
    ]
    @test certificate.applied_rewrites[1].original_span == (start = 1, stop = 1)
    @test certificate.applied_rewrites[2].original_span == (start = 2, stop = 3)
    @test certificate.applied_rewrites[3].original_span == (start = 4, stop = 5)

    different_position_factors = [
        elementary_matrix(3, 1, 2, x, R),
        elementary_matrix(3, 1, 3, -x, R),
    ]
    different_position_certificate =
        Suslin._steinberg_adjacent_rewrite_optimization_certificate(different_position_factors)

    @test Suslin._verify_steinberg_optimization_certificate(different_position_certificate)
    @test isempty(different_position_certificate.applied_rewrites)
    @test different_position_certificate.optimized_factors == different_position_factors
    @test different_position_certificate.comparison_summary.factor_count_delta == 0

    tampered_optimized = copy(certificate.optimized_factors)
    tampered_optimized[1] = elementary_matrix(3, 1, 3, x + y + 2 * one(R), R)
    tampered_certificate = Suslin.SteinbergOptimizationCertificate(
        certificate.original_factors,
        tampered_optimized,
        certificate.applied_rewrites,
        certificate.comparison_summary,
        certificate.original_product,
        certificate.optimized_product,
        certificate.verification,
    )
    @test !Suslin._verify_steinberg_optimization_certificate(tampered_certificate)
end
