using Test
using Suslin
using Oscar

const TORICBUILDER_CASE008_D14_COLUMN_BOUNDARY_PATH =
    joinpath(@__DIR__, "..", "fixtures", "toricbuilder_case008_d14_column_boundary.jl")

@testset "ToricBuilder case_008 d=14 Laurent column boundary" begin
    @test isfile(TORICBUILDER_CASE008_D14_COLUMN_BOUNDARY_PATH)

    include(TORICBUILDER_CASE008_D14_COLUMN_BOUNDARY_PATH)
    fixture = ToricBuilderCase008D14ColumnBoundary.boundary_fixture()

    @test fixture.case_id == "case_008"
    @test fixture.source_block == :column_transformation_upper_left_q_block
    @test fixture.source_matrix_dimensions == (30, 30)
    @test fixture.source_column_transformation_dimensions == (60, 60)
    @test fixture.first_failing_peel_dimension == 14
    @test fixture.passed_peel_dimensions == (30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15)
    @test length(fixture.failing_column) == 14
    @test fixture.ring_description == "GF(2)[u^+/-1, v^+/-1]"
    @test Tuple(string.(gens(fixture.ring))) == ("u", "v")
    @test coefficient_ring(fixture.ring) == GF(2)
    @test Suslin.is_unimodular_column(fixture.failing_column, fixture.ring)
    @test count(!iszero, fixture.failing_column) == 14
    @test count(is_unit, fixture.failing_column) == 0
    @test ToricBuilderCase008D14ColumnBoundary.column_statistics(fixture.failing_column) ==
          (; last_column_nonzero_count = 14, max_entry_term_count = 3734)

    provenance = fixture.boundary_provenance
    @test provenance.source == :explicit_bounded_case008_report
    @test provenance.stage == :certificate_construction
    @test provenance.current_peel_dimension == 14
    @test provenance.last_completed_peel_dimension == 15
    @test provenance.failure_code == :unsupported_laurent_column_family
    @test provenance.old_d15_boundary_cleared == true
    @test ToricBuilderCase008D14ColumnBoundary.validate_boundary_fixture(fixture) == :ok

    corrupted = ToricBuilderCase008D14ColumnBoundary.corrupted_column_negative_control(fixture)
    @test ToricBuilderCase008D14ColumnBoundary.validate_boundary_fixture(corrupted) == :wrong_column

    non_unimodular = ToricBuilderCase008D14ColumnBoundary.non_unimodular_negative_control(fixture)
    @test !Suslin.is_unimodular_column(non_unimodular.failing_column, non_unimodular.ring)
    @test ToricBuilderCase008D14ColumnBoundary.validate_boundary_fixture(non_unimodular) == :not_unimodular

    old_d15 = merge(
        fixture,
        (;
            first_failing_peel_dimension = 15,
            boundary_provenance = merge(
                fixture.boundary_provenance,
                (;
                    current_peel_dimension = 15,
                    last_completed_peel_dimension = 16,
                    old_d15_boundary_cleared = false,
                ),
            ),
        ),
    )
    @test ToricBuilderCase008D14ColumnBoundary.validate_boundary_fixture(old_d15) == :old_d15_boundary

    wrong_passed = merge(fixture, (; passed_peel_dimensions = (30, 29)))
    @test ToricBuilderCase008D14ColumnBoundary.validate_boundary_fixture(wrong_passed) == :wrong_passed_peel_dimensions

    wrong_ring, (wrong_u, wrong_v) = Suslin.suslin_laurent_polynomial_ring(GF(3), ["u", "v"])
    wrong_ring_fixture = merge(
        fixture,
        (;
            ring = wrong_ring,
            failing_column = [idx == 1 ? wrong_u : idx == 2 ? wrong_v : zero(wrong_ring) for idx in 1:14],
        ),
    )
    @test ToricBuilderCase008D14ColumnBoundary.validate_boundary_fixture(wrong_ring_fixture) == :wrong_ring
end
