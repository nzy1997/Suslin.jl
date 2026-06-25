using Test
using Suslin
using Oscar

include(joinpath(@__DIR__, "..", "fixtures", "toricbuilder_case010_column_boundary.jl"))

@testset "Laurent column reduction diagnostics" begin
    fixture = ToricBuilderCase010ColumnBoundary.boundary_fixture()

    unsupported = Suslin.diagnose_unimodular_column_reduction(
        fixture.failing_column,
        fixture.ring,
    )
    @test unsupported.status == :unsupported
    @test unsupported.failure_code == :unsupported_laurent_column_family
    @test unsupported.column_length == 5
    @test unsupported.ring_profile.kind == :laurent_polynomial
    @test unsupported.ring_profile.generators == ("u", "v")
    @test occursin("unsupported exact unimodular column reduction", unsupported.message)
    @test occursin(
        "no supported unit, witness-unit, monicity-normalized, or 3-entry block reduction stage applies",
        unsupported.message,
    )
    for stage in (:unit_entry, :witness_unit, :monicity_normalization, :three_entry_block)
        @test stage in unsupported.attempted_stages
    end
    @test :laurent_normalization in unsupported.attempted_stages

    thrown = try
        Suslin.reduce_unimodular_column(fixture.failing_column, fixture.ring)
        nothing
    catch err
        err
    end
    @test thrown isa ArgumentError
    @test occursin(
        "no supported unit, witness-unit, monicity-normalized, or 3-entry block reduction stage applies",
        sprint(showerror, thrown),
    )

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
