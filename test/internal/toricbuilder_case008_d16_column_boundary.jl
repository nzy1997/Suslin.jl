using Test
using Suslin
using Oscar

const TORICBUILDER_CASE008_D16_COLUMN_BOUNDARY_PATH =
    joinpath(@__DIR__, "..", "fixtures", "toricbuilder_case008_d16_column_boundary.jl")

function _case008_d16_stage_detail(diagnostic, stage::Symbol)
    idx = findfirst(detail -> detail.stage == stage, diagnostic.stage_details)
    return idx === nothing ? nothing : diagnostic.stage_details[idx]
end

function _test_case008_d16_stage_details_shape(diagnostic)
    hasproperty(diagnostic, :stage_details) || return
    stage_details = diagnostic.stage_details
    @test stage_details isa Tuple
    @test all(detail -> detail isa NamedTuple, stage_details)
end

@testset "ToricBuilder case_008 d=16 Laurent column boundary" begin
    @test isfile(TORICBUILDER_CASE008_D16_COLUMN_BOUNDARY_PATH)

    include(TORICBUILDER_CASE008_D16_COLUMN_BOUNDARY_PATH)
    fixture = ToricBuilderCase008D16ColumnBoundary.boundary_fixture()

    @test fixture.case_id == "case_008"
    @test fixture.source_matrix_dimensions == (30, 30)
    @test fixture.source_column_transformation_dimensions == (60, 60)
    @test fixture.first_failing_peel_dimension == 16
    @test fixture.passed_peel_dimensions == (30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17)
    @test length(fixture.failing_column) == 16
    @test fixture.ring_description == "GF(2)[u^+/-1, v^+/-1]"
    @test Tuple(string.(gens(fixture.ring))) == ("u", "v")
    @test Suslin.is_unimodular_column(fixture.failing_column, fixture.ring)
    @test count(is_unit, fixture.failing_column) == 0

    diagnostic = Suslin.diagnose_unimodular_column_reduction(
        fixture.failing_column,
        fixture.ring,
    )
    @test diagnostic.status == :supported
    @test diagnostic.failure_code === nothing
    @test diagnostic.column_length == 16
    @test diagnostic.ring_profile.kind == :laurent_polynomial
    @test diagnostic.ring_profile.generators == ("u", "v")
    @test hasproperty(diagnostic, :stage_details)
    _test_case008_d16_stage_details_shape(diagnostic)
    @test length(diagnostic.stage_details) == length(diagnostic.attempted_stages)

    witness_detail = _case008_d16_stage_detail(diagnostic, :laurent_witness_unit)
    @test witness_detail !== nothing
    @test witness_detail.outcome == :witness_without_unit
    @test witness_detail.witness_unit_index === nothing

    preconditioning_detail =
        _case008_d16_stage_detail(diagnostic, :laurent_elementary_row_preconditioning)
    @test preconditioning_detail !== nothing
    @test preconditioning_detail.outcome == :supported
    @test preconditioning_detail.target_index == 1
    @test preconditioning_detail.source_index == 10
    @test preconditioning_detail.coefficient == one(fixture.ring)
    @test preconditioning_detail.transformed_stage == :witness_unit

    @test ToricBuilderCase008D16ColumnBoundary.validate_boundary_fixture(fixture) == :ok

    corrupted = ToricBuilderCase008D16ColumnBoundary.non_unimodular_negative_control(fixture)
    @test !Suslin.is_unimodular_column(corrupted.failing_column, corrupted.ring)
    @test ToricBuilderCase008D16ColumnBoundary.validate_boundary_fixture(corrupted) == :not_unimodular

    wrong_dimension = merge(fixture, (; first_failing_peel_dimension = 15))
    @test ToricBuilderCase008D16ColumnBoundary.validate_boundary_fixture(wrong_dimension) == :wrong_peel_dimension

    wrong_passed = merge(fixture, (; passed_peel_dimensions = (30, 29)))
    @test ToricBuilderCase008D16ColumnBoundary.validate_boundary_fixture(wrong_passed) == :wrong_passed_peel_dimensions

    wrong_ring, (wrong_u, wrong_v) = Suslin.suslin_laurent_polynomial_ring(GF(3), ["u", "v"])
    wrong_ring_fixture = merge(
        fixture,
        (;
            ring = wrong_ring,
            failing_column = [idx == 1 ? wrong_u : idx == 2 ? wrong_v : zero(wrong_ring) for idx in 1:16],
        ),
    )
    @test ToricBuilderCase008D16ColumnBoundary.validate_boundary_fixture(wrong_ring_fixture) == :wrong_ring
end
