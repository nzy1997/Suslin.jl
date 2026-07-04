using Test
using Suslin
using Oscar

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

const STEINBERG_OPTIMIZATION_FIXTURE_PATH =
    joinpath(@__DIR__, "..", "fixtures", "steinberg_optimization_cases.jl")

if !isdefined(Main, :SteinbergOptimizationFixtureCatalog)
    include(STEINBERG_OPTIMIZATION_FIXTURE_PATH)
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
