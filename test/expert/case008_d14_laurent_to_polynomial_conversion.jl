using Test
using Suslin
using Oscar

if !(@isdefined case008_d14_laurent_to_polynomial_conversion_report)
    include(joinpath(
        @__DIR__,
        "..",
        "..",
        "scripts",
        "report_case008_d14_laurent_to_polynomial_conversion.jl",
    ))
end

@testset "case_008 d=14 Laurent-to-polynomial conversion report" begin
    fixture = ToricBuilderCase008D14ColumnBoundary.boundary_fixture()
    @test ToricBuilderCase008D14ColumnBoundary.validate_boundary_fixture(fixture) == :ok

    report = case008_d14_laurent_to_polynomial_conversion_report(
        fixture;
        timeout_seconds = 300,
    )
    @test _case008_d14_ltp_report_has_required_fields(report)
    @test report.case_id == "case_008"
    @test report.source_dimension == 14
    @test report.ring_generators == ("u", "v")
    @test report.source_boundary == :ToricBuilderCase008D14ColumnBoundary
    @test report.selected_entry == 9
    @test report.selected_generator == "v"
    @test report.conversion_status == :ok
    @test report.conversion_replay == :pass
    @test report.conversion_factors == 17
    @test report.source_terms == 50311
    @test report.polynomial_terms == 2774551
    @test report.polynomialization_power == 93
    @test report.polynomial_unimodular
    @test report.polynomial_witness_replay == :pass
    @test report.polynomial_witness_terms isa Integer
    @test report.polynomial_witness_terms >= 0
    @test report.polynomial_witness_provenance isa Symbol
    @test isfinite(report.elapsed_seconds)
    @test report.elapsed_seconds <= report.timeout_seconds
    @test report.final_status == :ok
    @test report.status == CASE008_D14_LTP_REPORT_STATUS
    @test validate_case008_d14_laurent_to_polynomial_conversion_report(
        report,
        fixture,
    ) == :ok
end

@testset "case_008 d=14 Laurent-to-polynomial conversion negative controls" begin
    report = case008_d14_laurent_to_polynomial_conversion_report(;
        timeout_seconds = 300,
    )
    controls = report.negative_controls
    @test report.negative_controls_pass
    @test controls.nonunimodular_source == :rejected
    @test controls.modified_preprocessing_factor != :ok
    @test controls.stale_intermediate_column != :ok
    @test controls.modified_polynomialization_power != :ok
    @test controls.noninvertible_inverse_map != :ok
    @test controls.copied_polynomial_column != :ok
    @test controls.modified_witness_coefficient != :ok
    @test controls.copied_witness_from_different_column != :ok
end

@testset "case_008 d=14 Laurent-to-polynomial conversion markdown" begin
    report = case008_d14_laurent_to_polynomial_conversion_report(;
        timeout_seconds = 300,
    )
    markdown = render_case008_d14_laurent_to_polynomial_conversion_markdown(report)
    @test occursin("case_id = case_008", markdown)
    @test occursin("source_dimension = 14", markdown)
    @test occursin("selected_entry = 9", markdown)
    @test occursin("selected_generator = v", markdown)
    @test occursin("conversion_status = ok", markdown)
    @test occursin("conversion_replay = pass", markdown)
    @test occursin("conversion_factors = 17", markdown)
    @test occursin("source_terms = 50311", markdown)
    @test occursin("polynomial_terms = 2774551", markdown)
    @test occursin("polynomialization_power = 93", markdown)
    @test occursin("polynomial_unimodular = true", markdown)
    @test occursin("polynomial_witness_replay = pass", markdown)
end
