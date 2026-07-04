using Test
using Suslin
using Oscar

const TORICBUILDER_CASE008_D15_COLUMN_BOUNDARY_PATH =
    joinpath(@__DIR__, "..", "fixtures", "toricbuilder_case008_d15_column_boundary.jl")

@testset "ToricBuilder case_008 d=15 Laurent column boundary" begin
    @test isfile(TORICBUILDER_CASE008_D15_COLUMN_BOUNDARY_PATH)

    include(TORICBUILDER_CASE008_D15_COLUMN_BOUNDARY_PATH)
    fixture = ToricBuilderCase008D15ColumnBoundary.boundary_fixture()

    @test fixture.case_id == "case_008"
    @test fixture.source_matrix_dimensions == (30, 30)
    @test fixture.source_column_transformation_dimensions == (60, 60)
    @test fixture.first_failing_peel_dimension == 15
    @test fixture.passed_peel_dimensions == (30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16)
    @test length(fixture.failing_column) == 15
    @test fixture.ring_description == "GF(2)[u^+/-1, v^+/-1]"
    @test Tuple(string.(gens(fixture.ring))) == ("u", "v")
    @test Suslin.is_unimodular_column(fixture.failing_column, fixture.ring)
    @test count(is_unit, fixture.failing_column) == 0

    diagnostic = Suslin.diagnose_unimodular_column_reduction(
        fixture.failing_column,
        fixture.ring,
    )
    @test diagnostic.status == :unsupported
    @test diagnostic.failure_code == :unsupported_laurent_column_family
    @test diagnostic.column_length == 15
    @test diagnostic.ring_profile.kind == :laurent_polynomial
    @test diagnostic.ring_profile.generators == ("u", "v")
    @test ToricBuilderCase008D15ColumnBoundary.validate_boundary_fixture(fixture) == :ok
    @test_throws ArgumentError Suslin.reduce_unimodular_column(
        fixture.failing_column,
        fixture.ring,
    )

    corrupted = ToricBuilderCase008D15ColumnBoundary.non_unimodular_negative_control(fixture)
    @test !Suslin.is_unimodular_column(corrupted.failing_column, corrupted.ring)
    @test ToricBuilderCase008D15ColumnBoundary.validate_boundary_fixture(corrupted) == :not_unimodular
    @test_throws ArgumentError Suslin.reduce_unimodular_column(
        corrupted.failing_column,
        corrupted.ring,
    )

    wrong_dimension = merge(fixture, (; first_failing_peel_dimension = 16))
    @test ToricBuilderCase008D15ColumnBoundary.validate_boundary_fixture(wrong_dimension) == :wrong_peel_dimension

    wrong_passed = merge(fixture, (; passed_peel_dimensions = (30, 29)))
    @test ToricBuilderCase008D15ColumnBoundary.validate_boundary_fixture(wrong_passed) == :wrong_passed_peel_dimensions

    wrong_ring, (wrong_u, wrong_v) = Suslin.suslin_laurent_polynomial_ring(GF(3), ["u", "v"])
    wrong_ring_fixture = merge(
        fixture,
        (;
            ring = wrong_ring,
            failing_column = [idx == 1 ? wrong_u : idx == 2 ? wrong_v : zero(wrong_ring) for idx in 1:15],
        ),
    )
    @test ToricBuilderCase008D15ColumnBoundary.validate_boundary_fixture(wrong_ring_fixture) == :wrong_ring
end
