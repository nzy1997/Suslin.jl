using Test
using Suslin
using Oscar

include(joinpath(@__DIR__, "..", "fixtures", "toricbuilder_case008_d14_column_boundary.jl"))

const CASE008_D14_LAURENT_TO_POLYNOMIAL_REQUIRED_FIELDS = (
    :case_id,
    :source_dimension,
    :ring_generators,
    :source_boundary,
    :conversion_certificate_status,
    :ordinary_ecp_child_status,
    :lifted_factor_count,
    :source_term_count,
    :converted_term_count,
    :degree_growth_before,
    :degree_growth_after,
    :elapsed_seconds,
    :final_target_column,
    :verification_status,
)

@testset "case_008 d=14 Laurent-to-polynomial bridge contract" begin
    fixture = ToricBuilderCase008D14ColumnBoundary.boundary_fixture()
    @test ToricBuilderCase008D14ColumnBoundary.validate_boundary_fixture(fixture) == :ok

    manifest_path = joinpath(@__DIR__, "..", "ci", "shards.toml")
    manifest = read(manifest_path, String)
    @test occursin(
        "path = \"expert/case008_d14_laurent_to_polynomial_bridge.jl\"",
        manifest,
    )
    manifest_expectations = read(joinpath(@__DIR__, "..", "ci", "runtests.jl"), String)
    @test occursin(
        "(\"expert/case008_d14_laurent_to_polynomial_bridge.jl\", \"expert\", \"expert-core\")",
        manifest_expectations,
    )

    report = case008_d14_laurent_to_polynomial_bridge_report(; timeout_seconds = 600)
    @test all(field -> hasproperty(report, field), CASE008_D14_LAURENT_TO_POLYNOMIAL_REQUIRED_FIELDS)
    @test report.case_id == "case_008"
    @test report.source_dimension == 14
    @test report.ring_generators == ("u", "v")
    @test report.source_boundary == :ToricBuilderCase008D14ColumnBoundary
    @test report.conversion_certificate_status == :ok
    @test report.ordinary_ecp_child_status == :ok
    @test report.lifted_factor_count isa Integer
    @test report.source_term_count isa Integer
    @test report.converted_term_count isa Integer
    @test report.degree_growth_before isa Integer
    @test report.degree_growth_after isa Integer
    @test report.elapsed_seconds isa Real
    @test report.verification_status == :laurent_to_polynomial_bridge_pass
    @test validate_case008_d14_laurent_to_polynomial_bridge_report(report) ==
          :laurent_to_polynomial_bridge_pass
end
