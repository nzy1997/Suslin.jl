using Test

const TORICBUILDER_CACHE_STATUS_REPORT_SCRIPT =
    joinpath(@__DIR__, "..", "..", "scripts", "report_toricbuilder_cache_q_blocks.jl")
const TORICBUILDER_CACHE_STATUS_REPORT_PATH =
    joinpath(@__DIR__, "..", "..", "docs", "audits", "2026-06-24-toricbuilder-cache-q-block-status.md")

@testset "ToricBuilder cache Q-block status report" begin
    @test isfile(TORICBUILDER_CACHE_STATUS_REPORT_SCRIPT)

    include(TORICBUILDER_CACHE_STATUS_REPORT_SCRIPT)
    report = ToricBuilderCacheQBlockStatusReport.build_report()

    @test length(report.rows) == 12
    @test report.source_fixture == "test/fixtures/toricbuilder_cache_q_blocks.jl"

    by_id = Dict(row.case_id => row for row in report.rows)
    @test Set(keys(by_id)) == Set(["case_$(lpad(string(idx), 3, "0"))" for idx in 1:12])

    @test by_id["case_001"].route_status == :gl_certificate_pass
    @test by_id["case_001"].public_elementary_status == :staged_boundary
    @test by_id["case_001"].determinant_class == :laurent_monomial_unit
    @test by_id["case_001"].verified == true
    @test by_id["case_001"].factor_count > 0
    @test by_id["case_001"].runtime_seconds > 0

    @test by_id["case_002"].route_status == :gl_certificate_pass
    @test by_id["case_002"].public_elementary_status == :staged_boundary
    @test by_id["case_002"].verified == true
    @test by_id["case_002"].decomposed_base_matrix_count == 186
    @test by_id["case_002"].runtime_seconds > 0

    @test by_id["case_003"].route_status == :gl_certificate_pass
    @test by_id["case_003"].public_elementary_status == :staged_boundary
    @test by_id["case_003"].verified == true
    @test by_id["case_003"].factor_count > 0
    @test by_id["case_003"].runtime_seconds > 0

    @test by_id["case_004"].route_status == :gl_certificate_pass
    @test by_id["case_004"].decomposed_base_matrix_count == 189
    @test by_id["case_005"].route_status == :gl_certificate_pass
    @test by_id["case_005"].decomposed_base_matrix_count == 168
    @test by_id["case_006"].route_status == :gl_certificate_pass
    @test by_id["case_006"].decomposed_base_matrix_count == 212

    @test by_id["case_010"].route_status == :gl_certificate_pass
    @test by_id["case_010"].public_elementary_status == :staged_boundary
    @test by_id["case_010"].determinant_class == :laurent_monomial_unit
    @test by_id["case_010"].verified == true
    @test by_id["case_010"].decomposed_base_matrix_count > 0
    @test by_id["case_010"].runtime_seconds > 0
    @test by_id["case_010"].error_details == "none"

    @test by_id["case_011"].matrix_size == (288, 288)
    @test by_id["case_011"].route_status == :not_exercised_in_default_report
    @test by_id["case_011"].public_elementary_status == :not_run
    @test by_id["case_011"].decomposed_base_matrix_count == :not_run
    @test by_id["case_011"].runtime_seconds == :not_run

    markdown = ToricBuilderCacheQBlockStatusReport.render_markdown(report)
    @test occursin("# ToricBuilder Cache Q-Block Status Report", markdown)
    @test occursin("Matrix size", markdown)
    @test occursin("Decomposed base matrices", markdown)
    @test occursin("Runtime seconds", markdown)
    @test occursin(r"\| case_001 \| 6x6 \| 30 \| default_contract \| gl_certificate_pass \| staged_boundary \| laurent_monomial_unit \| 50 \| [0-9]+\.[0-9]{3} \|", markdown)
    @test occursin(r"\| case_002 \| 14x14 \| 79 \| default_contract \| gl_certificate_pass \| staged_boundary \| laurent_monomial_unit \| 186 \| [0-9]+\.[0-9]{3} \|", markdown)
    @test occursin(r"\| case_003 \| 6x6 \| 27 \| default_contract \| gl_certificate_pass \| staged_boundary \| laurent_monomial_unit \| 49 \| [0-9]+\.[0-9]{3} \|", markdown)
    @test occursin(r"\| case_004 \| 18x18 \| 73 \| default_contract \| gl_certificate_pass \| staged_boundary \| laurent_monomial_unit \| 189 \| [0-9]+\.[0-9]{3} \|", markdown)
    @test occursin(r"\| case_005 \| 14x14 \| 90 \| default_contract \| gl_certificate_pass \| staged_boundary \| laurent_monomial_unit \| 168 \| [0-9]+\.[0-9]{3} \|", markdown)
    @test occursin(r"\| case_006 \| 18x18 \| 99 \| default_contract \| gl_certificate_pass \| staged_boundary \| laurent_monomial_unit \| 212 \| [0-9]+\.[0-9]{3} \|", markdown)
    @test occursin(r"\| case_010 \| 6x6 \| 34 \| default_contract \| gl_certificate_pass \| staged_boundary \| laurent_monomial_unit \| [1-9][0-9]* \| [0-9]+\.[0-9]{3} \|", markdown)
    @test occursin("| case_011 | 288x288 | 14713 | optional_slow | not_exercised_in_default_report | not_run | not_run | not_run | not_run |", markdown)
    @test hasproperty(by_id["case_001"], :stage_timings)
    @test by_id["case_001"].stage_timings.determinant_classification.status == :pass
    @test by_id["case_001"].stage_timings.normalization.status == :pass
    @test by_id["case_001"].stage_timings.certificate_construction.status == :pass
    @test by_id["case_001"].stage_timings.verification.status == :pass
    @test by_id["case_001"].stage_timings.determinant_classification.elapsed_seconds >= 0

    @test hasproperty(by_id["case_007"], :stage_timings)
    @test by_id["case_007"].stage_timings.determinant_classification.status == :not_run
    @test by_id["case_007"].stage_timings.normalization.status == :not_run
    @test by_id["case_007"].stage_timings.certificate_construction.status == :not_run
    @test by_id["case_007"].stage_timings.verification.status == :not_run

    @test occursin("## Stage Timing Details", markdown)
    @test occursin("Determinant classification", markdown)
    @test occursin("Certificate construction", markdown)
    @test !occursin("## Route Error Details", markdown)
    @test !occursin("unsupported exact unimodular column reduction", markdown)
    @test occursin("julia --project=. scripts/report_toricbuilder_cache_q_blocks.jl", markdown)

    parsed_timeout = ToricBuilderCacheQBlockStatusReport._parse_args([
        "--exercise=case_007",
        "--timeout-seconds=1.5",
        "--output=/tmp/qblock-timeout.md",
    ])
    @test parsed_timeout.exercised == ["case_007"]
    @test parsed_timeout.timeout_seconds == 1.5

    @test_throws ArgumentError ToricBuilderCacheQBlockStatusReport._parse_args([
        "--timeout-seconds=0",
    ])
    @test_throws ArgumentError ToricBuilderCacheQBlockStatusReport._parse_args([
        "--timeout-seconds=not-a-number",
    ])

    @test_throws ArgumentError ToricBuilderCacheQBlockStatusReport.build_report(;
        exercised_case_ids = ("case_999",),
    )

    unknown_output = tempname()
    @test_throws ArgumentError ToricBuilderCacheQBlockStatusReport.main([
        "--exercise=case_999",
        "--output=$(unknown_output)",
    ])
    @testset "_record_stage! failure paths" begin
        entry = (;
            id = "case_failure",
            dimensions = (; matrix = (3, 3)),
            sparse_entry_count = 9,
            expected_test_level = :default_contract,
        )

        start_ns = time_ns()
        boundary_timings = Dict{Symbol, Any}()
        boundary_result = ToricBuilderCacheQBlockStatusReport._record_stage!(
            boundary_timings,
            :determinant_classification,
            () -> throw(ArgumentError("unsupported Laurent GL_n determinant test boundary")),
        )
        @test boundary_result.status == :certified_algorithm_boundary
        boundary_row = ToricBuilderCacheQBlockStatusReport._stage_failure_row(
            entry,
            boundary_timings,
            boundary_result,
            start_ns,
        )
        @test boundary_row.route_status == :certified_algorithm_boundary
        @test boundary_row.stage_timings.determinant_classification.status ==
              :certified_algorithm_boundary

        start_ns = time_ns()
        route_timings = Dict{Symbol, Any}()
        route_result = ToricBuilderCacheQBlockStatusReport._record_stage!(
            route_timings,
            :determinant_classification,
            () -> throw(ErrorException("unexpected qblock test failure")),
        )
        @test route_result.status == :route_error
        route_row = ToricBuilderCacheQBlockStatusReport._stage_failure_row(
            entry,
            route_timings,
            route_result,
            start_ns,
        )
        @test route_row.route_status == :route_error
        @test route_row.stage_timings.determinant_classification.status == :route_error
    end
    @test !isfile(unknown_output)
    @test isfile(TORICBUILDER_CACHE_STATUS_REPORT_PATH)
end
