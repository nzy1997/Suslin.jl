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

function _assert_laurent_native_ecp_boundary_detail(
    detail;
    requires_descent_measure::Bool = true,
    certified_descent_scope = nothing,
    requires_link_witness::Bool = true,
    next_boundary = certified_descent_scope === nothing ? nothing : :laurent_link_witness,
)
    @test detail !== nothing
    @test detail.outcome == :staged_boundary
    @test detail.boundary == :laurent_native_ecp
    @test detail.requires_descent_measure == requires_descent_measure
    @test detail.certified_descent_scope == certified_descent_scope
    @test detail.next_boundary == next_boundary
    @test detail.requires_link_witness == requires_link_witness
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
    descent_idx = findfirst(==(:laurent_descent_step_certificate), d14.attempted_stages)
    link_witness_idx =
        findfirst(==(:laurent_link_witness_certificate), d14.attempted_stages)
    @test preconditioning_idx !== nothing
    @test boundary_idx !== nothing
    @test descent_idx !== nothing
    @test link_witness_idx !== nothing
    @test boundary_idx > preconditioning_idx
    @test boundary_idx > descent_idx
    @test link_witness_idx > descent_idx
    @test boundary_idx > link_witness_idx
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
    descent = _diagnostic_stage_detail(d14, :laurent_descent_step_certificate)
    @test descent.outcome == :certified_descent_step
    @test descent.descent_scope == :single_certified_step
    @test descent.operation_family == :entry_addition
    @test descent.target_index == 1
    @test descent.source_index == 2
    @test descent.coefficient == 1
    @test descent.exponent == (-1, 1)
    @test descent.measure_relation == :strict_decrease
    @test descent.replay_status == :ok
    @test descent.next_boundary == :laurent_link_witness
    @test Suslin._strictly_decreases_laurent_measure(
        descent.before_measure,
        descent.after_measure,
    )
    link_witness =
        _diagnostic_stage_detail(d14, :laurent_link_witness_certificate)
    @test link_witness !== nothing
    @test link_witness.outcome == :certified_link_witness
    @test link_witness.witness_family == :two_entry_laurent_combination
    @test link_witness.pivot_index == 10
    @test link_witness.partner_index == 1
    @test link_witness.coefficient == 1
    @test link_witness.exponent == (1, -1)
    @test link_witness.replay_status == :ok
    @test link_witness.identity_status == :verified
    @test link_witness.certificate_status == :link_witness_certificate
    @test link_witness.context_status == :link_witness_context
    @test link_witness.source_endpoint.status == :link_witness_endpoint_metadata
    @test link_witness.target_endpoint.status == :link_witness_endpoint_metadata
    @test link_witness.source_endpoint.case_id == "case_008"
    @test link_witness.target_endpoint.case_id == "case_008"
    @test link_witness.source_endpoint.entry_index == 10
    @test link_witness.target_endpoint.entry_index == 10
    @test link_witness.next_boundary == :laurent_endpoint_reduction
    _assert_laurent_native_ecp_boundary_detail(
        _diagnostic_stage_detail(d14, :laurent_native_ecp_boundary);
        requires_descent_measure = false,
        certified_descent_scope = :single_certified_step,
        requires_link_witness = false,
        next_boundary = :laurent_endpoint_reduction,
    )
    descent_only_boundary = Suslin._laurent_native_ecp_boundary_stage_detail(
        d14_fixture.ring;
        certified_descent_step = true,
        certified_link_witness = false,
    )
    _assert_laurent_native_ecp_boundary_detail(
        descent_only_boundary;
        requires_descent_measure = false,
        certified_descent_scope = :single_certified_step,
        requires_link_witness = true,
        next_boundary = :laurent_link_witness,
    )

    ordinary_R, _ = Oscar.polynomial_ring(GF(2), ["x", "y"])
    @test Suslin._laurent_link_witness_diagnostic_certificate(
        [one(ordinary_R), zero(ordinary_R)],
        ordinary_R,
    ) === nothing

    d15_fixture = ToricBuilderCase008D15ColumnBoundary.boundary_fixture()
    d15 = Suslin.diagnose_unimodular_column_reduction(
        d15_fixture.failing_column,
        d15_fixture.ring,
    )
    @test d15.status == :supported
    @test d15.failure_code === nothing
    @test :laurent_elementary_row_preconditioning in d15.attempted_stages
    @test !(:laurent_descent_step_certificate in d15.attempted_stages)
    @test !(:laurent_link_witness_certificate in d15.attempted_stages)
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
    @test !(:laurent_descent_step_certificate in non_unimodular_diagnostic.attempted_stages)
    @test !(:laurent_link_witness_certificate in non_unimodular_diagnostic.attempted_stages)
    @test !(:laurent_native_ecp_boundary in non_unimodular_diagnostic.attempted_stages)

    tampered_column = copy(d14_fixture.failing_column)
    tampered_column[1] = tampered_column[1] + one(d14_fixture.ring)
    tampered = Suslin.diagnose_unimodular_column_reduction(
        tampered_column,
        d14_fixture.ring;
        assume_unimodular = true,
        laurent_large_support_diagnostic_decline = true,
    )
    @test tampered.status == :unsupported
    @test !(:laurent_descent_step_certificate in tampered.attempted_stages)
    @test !(:laurent_link_witness_certificate in tampered.attempted_stages)
    tampered_boundary = _diagnostic_stage_detail(tampered, :laurent_native_ecp_boundary)
    _assert_laurent_native_ecp_boundary_detail(tampered_boundary)
    @test tampered_boundary.requires_descent_measure == true
end
