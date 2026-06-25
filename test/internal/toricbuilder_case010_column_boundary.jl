using Test
using Suslin
using Oscar

const TORICBUILDER_CASE010_COLUMN_BOUNDARY_PATH =
    joinpath(@__DIR__, "..", "fixtures", "toricbuilder_case010_column_boundary.jl")

@testset "ToricBuilder case_010 Laurent column boundary" begin
    @test isfile(TORICBUILDER_CASE010_COLUMN_BOUNDARY_PATH)

    include(TORICBUILDER_CASE010_COLUMN_BOUNDARY_PATH)
    fixture = ToricBuilderCase010ColumnBoundary.boundary_fixture()

    @test fixture.case_id == "case_010"
    @test (nrows(fixture.original_matrix), ncols(fixture.original_matrix)) == (6, 6)
    @test Suslin.verify_laurent_gl_normalization(fixture.original_matrix, fixture.normalization)
    @test det(fixture.normalized_matrix) == one(fixture.ring)
    @test fixture.first_failing_peel_dimension == 5
    @test length(fixture.failing_column) == 5
    @test fixture.ring_description == "GF(2)[u^+/-1, v^+/-1]"
    @test Tuple(string.(gens(fixture.ring))) == ("u", "v")
    @test Suslin.is_unimodular_column(fixture.failing_column, fixture.ring)
    @test ToricBuilderCase010ColumnBoundary.validate_boundary_fixture(fixture) == :ok

    diagnostic = Suslin.diagnose_unimodular_column_reduction(
        fixture.failing_column,
        fixture.ring,
    )
    @test diagnostic.status == :supported
    @test :laurent_unit_creation in diagnostic.attempted_stages

    negative = ToricBuilderCase010ColumnBoundary.non_unimodular_negative_control(fixture)
    @test !Suslin.is_unimodular_column(negative.failing_column, negative.ring)
    @test ToricBuilderCase010ColumnBoundary.validate_boundary_fixture(negative) == :not_unimodular

    wrong_ring, (wrong_u, wrong_v) = Suslin.suslin_laurent_polynomial_ring(GF(3), ["u", "v"])
    wrong_ring_fixture = merge(
        fixture,
        (;
            ring = wrong_ring,
            failing_column = [wrong_u, wrong_v, one(wrong_ring), zero(wrong_ring), zero(wrong_ring)],
        ),
    )
    @test ToricBuilderCase010ColumnBoundary.validate_boundary_fixture(wrong_ring_fixture) == :wrong_ring

    for perturbed in ToricBuilderCase010ColumnBoundary.single_entry_zero_perturbations(fixture)
        @test Suslin.is_unimodular_column(perturbed, fixture.ring)
    end
end
