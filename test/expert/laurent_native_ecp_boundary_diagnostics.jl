using Test
using Suslin
using Oscar

include(joinpath(@__DIR__, "..", "fixtures", "toricbuilder_case008_d14_column_boundary.jl"))
include(joinpath(@__DIR__, "..", "fixtures", "toricbuilder_case008_d15_column_boundary.jl"))

struct _EntryWithoutCoefficientIterator end

Base.iszero(::_EntryWithoutCoefficientIterator) = false

function _diagnostic_stage_detail(diagnostic, stage::Symbol)
    idx = findfirst(detail -> detail.stage == stage, diagnostic.stage_details)
    return idx === nothing ? nothing : diagnostic.stage_details[idx]
end

function _assert_laurent_native_ecp_boundary_detail(detail)
    @test detail !== nothing
    @test detail.outcome == :staged_boundary
    @test detail.boundary == :laurent_native_ecp
    @test detail.requires_descent_measure == true
    @test detail.requires_link_witness == true
    @test detail.requires_endpoint_reduction == true
    @test detail.requires_laurent_normality_replay == true
    @test detail.requires_recursive_peel_integration == true
    @test detail.fallback_policy == :diagnostic_only
end

@testset "Laurent native ECP boundary diagnostics" begin
    runtests = read(joinpath(@__DIR__, "..", "runtests.jl"), String)
    @test occursin(
        "\"expert/laurent_native_ecp_boundary_diagnostics.jl\"",
        runtests,
    )

    coefficient_failure_entry = _EntryWithoutCoefficientIterator()
    @test Suslin._column_reduction_entry_term_count(coefficient_failure_entry) === nothing
    @test Suslin._column_reduction_max_entry_term_count([coefficient_failure_entry]) === nothing
    @test Suslin._laurent_diagnostic_large_support_decline([coefficient_failure_entry]) === nothing

    d14_fixture = ToricBuilderCase008D14ColumnBoundary.boundary_fixture()
    @test ToricBuilderCase008D14ColumnBoundary.validate_boundary_fixture(d14_fixture) == :ok

    d14 = Suslin.diagnose_unimodular_column_reduction(
        d14_fixture.failing_column,
        d14_fixture.ring,
        assume_unimodular = true,
        laurent_large_support_diagnostic_decline = true,
    )
    @test d14.status == :unsupported
    @test d14.failure_code == :unsupported_laurent_column_family
    @test d14.column_length == 14
    @test :laurent_elementary_row_preconditioning in d14.attempted_stages
    @test :laurent_native_ecp_boundary in d14.attempted_stages
    preconditioning_idx = findfirst(==(:laurent_elementary_row_preconditioning), d14.attempted_stages)
    boundary_idx = findfirst(==(:laurent_native_ecp_boundary), d14.attempted_stages)
    @test preconditioning_idx !== nothing
    @test boundary_idx !== nothing
    @test boundary_idx > preconditioning_idx
    @test length(d14.stage_details) == length(d14.attempted_stages)
    d14_witness = _diagnostic_stage_detail(d14, :laurent_witness_unit)
    @test d14_witness !== nothing
    @test d14_witness.outcome == :witness_support_too_large
    @test d14_witness.max_entry_term_count == 3734
    d14_normalization = _diagnostic_stage_detail(d14, :laurent_normalization)
    @test d14_normalization !== nothing
    @test d14_normalization.outcome == :delegation_declined_large_support
    @test d14_normalization.normalized_status == :declined
    @test d14_normalization.normalized_failure_code == :support_too_large
    @test d14_normalization.max_entry_term_count == 3734
    _assert_laurent_native_ecp_boundary_detail(
        _diagnostic_stage_detail(d14, :laurent_native_ecp_boundary),
    )

    d15_fixture = ToricBuilderCase008D15ColumnBoundary.boundary_fixture()
    d15 = Suslin.diagnose_unimodular_column_reduction(
        d15_fixture.failing_column,
        d15_fixture.ring,
    )
    @test d15.status == :supported
    @test d15.failure_code === nothing
    @test :laurent_elementary_row_preconditioning in d15.attempted_stages
    @test !(:laurent_native_ecp_boundary in d15.attempted_stages)
    @test _diagnostic_stage_detail(d15, :laurent_native_ecp_boundary) === nothing
    d15_preconditioned =
        _diagnostic_stage_detail(d15, :laurent_elementary_row_preconditioning)
    @test d15_preconditioned !== nothing
    @test d15_preconditioned.outcome == :supported

    R, (u, _) = Suslin.suslin_laurent_polynomial_ring(GF(2), ["u", "v"])
    large_support_entry = sum(u^k for k in 1:1001)
    large_support_witness_column = [
        large_support_entry,
        large_support_entry + one(R),
        zero(R),
    ]
    large_support_witness = Suslin.diagnose_unimodular_column_reduction(
        large_support_witness_column,
        R,
    )
    @test large_support_witness.status == :supported
    @test large_support_witness.failure_code === nothing
    @test :laurent_witness_unit in large_support_witness.attempted_stages
    @test !(:laurent_native_ecp_boundary in large_support_witness.attempted_stages)
    large_support_witness_detail =
        _diagnostic_stage_detail(large_support_witness, :laurent_witness_unit)
    @test large_support_witness_detail !== nothing
    @test large_support_witness_detail.outcome == :supported
    @test large_support_witness_detail.witness_unit_index !== nothing

    non_unimodular =
        ToricBuilderCase008D14ColumnBoundary.non_unimodular_negative_control(d14_fixture)
    non_unimodular_diagnostic = Suslin.diagnose_unimodular_column_reduction(
        non_unimodular.failing_column,
        non_unimodular.ring,
    )
    @test non_unimodular_diagnostic.status == :precondition_failed
    @test non_unimodular_diagnostic.failure_code == :not_unimodular
    @test isempty(non_unimodular_diagnostic.attempted_stages)
    @test isempty(non_unimodular_diagnostic.stage_details)
    @test !(:laurent_native_ecp_boundary in non_unimodular_diagnostic.attempted_stages)
end
