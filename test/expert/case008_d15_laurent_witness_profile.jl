using Test
using Suslin
using Oscar

include(joinpath(@__DIR__, "..", "fixtures", "toricbuilder_case008_d15_column_boundary.jl"))

const CASE008_D15_EXPECTED_ATTEMPTED_STAGES = (
    :unit_entry,
    :laurent_unit_creation,
    :laurent_witness_unit,
    :laurent_normalization,
    :laurent_elementary_row_preconditioning,
)

function _diagnostic_stage_detail(diagnostic, stage::Symbol)
    idx = findfirst(detail -> detail.stage == stage, diagnostic.stage_details)
    return idx === nothing ? nothing : diagnostic.stage_details[idx]
end

function _case008_d15_supported_profile(diagnostic)::Bool
    diagnostic.status == :supported || return false
    diagnostic.failure_code === nothing || return false
    diagnostic.column_length == 15 || return false
    diagnostic.attempted_stages == CASE008_D15_EXPECTED_ATTEMPTED_STAGES || return false
    hasproperty(diagnostic, :stage_details) || return false
    length(diagnostic.stage_details) == length(CASE008_D15_EXPECTED_ATTEMPTED_STAGES) || return false
    any(detail -> detail.outcome == :supported, diagnostic.stage_details) || return false

    unit_entry = _diagnostic_stage_detail(diagnostic, :unit_entry)
    unit_entry !== nothing || return false
    unit_entry.outcome == :no_unit_entry || return false
    unit_entry.pivot_index === nothing || return false

    unit_creation = _diagnostic_stage_detail(diagnostic, :laurent_unit_creation)
    unit_creation !== nothing || return false
    unit_creation.outcome == :no_unit_creation_candidate || return false

    witness = _diagnostic_stage_detail(diagnostic, :laurent_witness_unit)
    witness !== nothing || return false
    witness.outcome == :witness_without_unit || return false
    witness.witness_unit_index === nothing || return false

    normalization = _diagnostic_stage_detail(diagnostic, :laurent_normalization)
    normalization !== nothing || return false
    normalization.outcome == :normalized_not_unimodular || return false
    normalization.normalized_column_length == 15 || return false
    normalization.normalized_status == :precondition_failed || return false
    normalization.normalized_failure_code == :not_unimodular || return false

    row_preconditioning =
        _diagnostic_stage_detail(diagnostic, :laurent_elementary_row_preconditioning)
    row_preconditioning !== nothing || return false
    row_preconditioning.outcome == :supported || return false
    row_preconditioning.target_index == 1 || return false
    row_preconditioning.source_indices == Tuple(2:15) || return false
    row_preconditioning.coefficient_strategy == :target_unit_laurent_linear_synthesis ||
        return false
    row_preconditioning.coefficient_count == 14 || return false
    row_preconditioning.transformed_stage == :unit_entry || return false

    return true
end

@testset "case_008 d=15 Laurent witness profile" begin
    fixture = ToricBuilderCase008D15ColumnBoundary.boundary_fixture()
    @test ToricBuilderCase008D15ColumnBoundary.validate_boundary_fixture(fixture) == :ok

    diagnostic = Suslin.diagnose_unimodular_column_reduction(
        fixture.failing_column,
        fixture.ring,
    )

    @test diagnostic.status == :supported
    @test diagnostic.failure_code === nothing
    @test diagnostic.column_length == 15
    @test diagnostic.ring_profile.kind == :laurent_polynomial
    @test diagnostic.ring_profile.generators == ("u", "v")
    @test diagnostic.attempted_stages == CASE008_D15_EXPECTED_ATTEMPTED_STAGES
    @test hasproperty(diagnostic, :stage_details)
    @test diagnostic.stage_details isa Tuple
    @test all(detail -> detail isa NamedTuple, diagnostic.stage_details)
    @test length(diagnostic.stage_details) == length(CASE008_D15_EXPECTED_ATTEMPTED_STAGES)
    @test any(detail -> detail.outcome == :supported, diagnostic.stage_details)

    unit_entry = _diagnostic_stage_detail(diagnostic, :unit_entry)
    @test unit_entry !== nothing
    @test unit_entry.outcome == :no_unit_entry
    @test unit_entry.pivot_index === nothing

    unit_creation = _diagnostic_stage_detail(diagnostic, :laurent_unit_creation)
    @test unit_creation !== nothing
    @test unit_creation.outcome == :no_unit_creation_candidate

    witness = _diagnostic_stage_detail(diagnostic, :laurent_witness_unit)
    @test witness !== nothing
    @test witness.outcome == :witness_without_unit
    @test witness.witness_unit_index === nothing

    normalization = _diagnostic_stage_detail(diagnostic, :laurent_normalization)
    @test normalization !== nothing
    @test normalization.outcome == :normalized_not_unimodular
    @test normalization.normalized_column_length == 15
    @test normalization.normalized_ring_kind == :polynomial
    @test normalization.normalized_status == :precondition_failed
    @test normalization.normalized_failure_code == :not_unimodular

    row_preconditioning =
        _diagnostic_stage_detail(diagnostic, :laurent_elementary_row_preconditioning)
    @test row_preconditioning !== nothing
    @test row_preconditioning.outcome == :supported
    @test row_preconditioning.target_index == 1
    @test row_preconditioning.source_indices == Tuple(2:15)
    @test row_preconditioning.coefficient_strategy == :target_unit_laurent_linear_synthesis
    @test row_preconditioning.coefficient_count == 14
    @test row_preconditioning.transformed_stage == :unit_entry

    @test _case008_d15_supported_profile(diagnostic)
end

@testset "synthetic direct-unit Laurent profile control" begin
    R, (u, v) = Suslin.suslin_laurent_polynomial_ring(GF(2), ["u", "v"])
    direct_unit_column = [one(R), u + v, zero(R)]

    diagnostic = Suslin.diagnose_unimodular_column_reduction(direct_unit_column, R)

    @test diagnostic.status == :supported
    @test diagnostic.failure_code === nothing
    @test diagnostic.column_length == 3
    @test diagnostic.attempted_stages == (:unit_entry,)
    @test !_case008_d15_supported_profile(diagnostic)

    unit_entry = _diagnostic_stage_detail(diagnostic, :unit_entry)
    @test unit_entry !== nothing
    @test unit_entry.outcome == :supported
    @test unit_entry.pivot_index == 1
end

@testset "synthetic unsupported Laurent profile control" begin
    R, (u, v) = Suslin.suslin_laurent_polynomial_ring(GF(2), ["u", "v"])
    unsupported_column = [u * v + u, u^2 + u + one(R), u * v + v^2 + one(R)]

    diagnostic = Suslin.diagnose_unimodular_column_reduction(unsupported_column, R)

    @test diagnostic.status == :unsupported
    @test diagnostic.failure_code == :unsupported_laurent_column_family
    @test !_case008_d15_supported_profile(diagnostic)
end
