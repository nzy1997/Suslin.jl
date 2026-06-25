using Test
using Suslin
using Oscar

include(joinpath(@__DIR__, "..", "fixtures", "toricbuilder_case010_column_boundary.jl"))

@testset "Laurent column reduction diagnostics" begin
    fixture = ToricBuilderCase010ColumnBoundary.boundary_fixture()

    case010 = Suslin.diagnose_unimodular_column_reduction(
        fixture.failing_column,
        fixture.ring,
    )
    @test case010.status == :supported
    @test case010.failure_code === nothing
    @test case010.column_length == 5
    @test case010.ring_profile.kind == :laurent_polynomial
    @test case010.ring_profile.generators == ("u", "v")
    @test :laurent_unit_creation in case010.attempted_stages

    unsupported_ring, (x, y) = Suslin.suslin_laurent_polynomial_ring(GF(2), ["x", "y"])
    unsupported_column = [x + y, x * y + x + one(unsupported_ring), x^2 + x * y + y^2 + one(unsupported_ring)]
    unsupported = Suslin.diagnose_unimodular_column_reduction(
        unsupported_column,
        unsupported_ring,
    )
    @test unsupported.status == :unsupported
    @test unsupported.failure_code == :unsupported_laurent_column_family
    @test unsupported.column_length == 3
    @test unsupported.ring_profile.kind == :laurent_polynomial
    @test unsupported.ring_profile.generators == ("x", "y")
    @test occursin("unsupported exact unimodular column reduction", unsupported.message)
    for stage in (:unit_entry, :laurent_unit_creation, :laurent_normalization, :witness_unit, :monicity_normalization)
        @test stage in unsupported.attempted_stages
    end

    d = nrows(fixture.normalized_matrix)
    supported_column = [fixture.normalized_matrix[row, d] for row in 1:d]
    supported = Suslin.diagnose_unimodular_column_reduction(supported_column, fixture.ring)
    @test supported.status == :supported
    @test supported.failure_code === nothing
    @test supported.column_length == d
    @test !isempty(supported.attempted_stages)
    Suslin.reduce_unimodular_column(supported_column, fixture.ring)

    negative = ToricBuilderCase010ColumnBoundary.non_unimodular_negative_control(fixture)
    precondition = Suslin.diagnose_unimodular_column_reduction(
        negative.failing_column,
        negative.ring,
    )
    @test precondition.status == :precondition_failed
    @test precondition.failure_code == :not_unimodular
    @test precondition.failure_code != :unsupported_laurent_column_family
    @test isempty(precondition.attempted_stages)
end
