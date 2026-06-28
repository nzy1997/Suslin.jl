using Test
using Suslin
using Oscar

include(joinpath(@__DIR__, "..", "fixtures", "toricbuilder_case010_column_boundary.jl"))
include(joinpath(@__DIR__, "..", "fixtures", "toricbuilder_case008_d21_column_boundary.jl"))

function _diagnostic_stage_detail(diagnostic, stage::Symbol)
    idx = findfirst(detail -> detail.stage == stage, diagnostic.stage_details)
    return idx === nothing ? nothing : diagnostic.stage_details[idx]
end

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
    @test hasproperty(case010, :stage_details)
    @test length(case010.stage_details) == length(case010.attempted_stages)
    case010_unit_creation = _diagnostic_stage_detail(case010, :laurent_unit_creation)
    @test case010_unit_creation !== nothing
    @test case010_unit_creation.outcome == :supported
    @test case010_unit_creation.pivot_index isa Integer

    case008 = ToricBuilderCase008D21ColumnBoundary.boundary_fixture()
    case008_d21 = Suslin.diagnose_unimodular_column_reduction(
        case008.failing_column,
        case008.ring,
    )
    @test case008_d21.status == :supported
    @test case008_d21.failure_code === nothing
    @test case008_d21.column_length == 21
    @test case008_d21.ring_profile.kind == :laurent_polynomial
    @test case008_d21.ring_profile.generators == ("u", "v")
    @test :laurent_witness_unit in case008_d21.attempted_stages
    @test hasproperty(case008_d21, :stage_details)
    @test length(case008_d21.stage_details) == length(case008_d21.attempted_stages)
    case008_d21_witness = _diagnostic_stage_detail(case008_d21, :laurent_witness_unit)
    @test case008_d21_witness !== nothing
    @test case008_d21_witness.outcome == :supported
    @test case008_d21_witness.witness_unit_index isa Integer

    unsupported_ring, (x, y) = Suslin.suslin_laurent_polynomial_ring(GF(2), ["x", "y"])
    unsupported_column = [x * y + x, x^2 + x + one(unsupported_ring), x * y + y^2 + one(unsupported_ring)]
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
    for stage in (:unit_entry, :laurent_unit_creation, :laurent_witness_unit, :laurent_normalization, :witness_unit, :monicity_normalization)
        @test stage in unsupported.attempted_stages
    end
    @test hasproperty(unsupported, :stage_details)
    @test length(unsupported.stage_details) == length(unsupported.attempted_stages)
    @test !any(detail -> detail.outcome == :supported, unsupported.stage_details)

    d = nrows(fixture.normalized_matrix)
    supported_column = [fixture.normalized_matrix[row, d] for row in 1:d]
    supported = Suslin.diagnose_unimodular_column_reduction(supported_column, fixture.ring)
    @test supported.status == :supported
    @test supported.failure_code === nothing
    @test supported.column_length == d
    @test !isempty(supported.attempted_stages)
    @test hasproperty(supported, :stage_details)
    @test length(supported.stage_details) == length(supported.attempted_stages)
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
    @test hasproperty(precondition, :stage_details)
    @test isempty(precondition.stage_details)
end
