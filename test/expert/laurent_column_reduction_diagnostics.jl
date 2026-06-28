using Test
using Suslin
using Oscar

include(joinpath(@__DIR__, "..", "fixtures", "toricbuilder_case010_column_boundary.jl"))
include(joinpath(@__DIR__, "..", "fixtures", "toricbuilder_case008_d21_column_boundary.jl"))
include(joinpath(@__DIR__, "..", "fixtures", "toricbuilder_case008_d16_column_boundary.jl"))

function _diagnostic_stage_detail(diagnostic, stage::Symbol)
    idx = findfirst(detail -> detail.stage == stage, diagnostic.stage_details)
    return idx === nothing ? nothing : diagnostic.stage_details[idx]
end

function _diagnostic_supported_stage(diagnostic)
    supported = filter(detail -> detail.outcome == :supported, diagnostic.stage_details)
    return isempty(supported) ? nothing : only(supported).stage
end

function _test_diagnostic_stage_details_shape(diagnostic)
    hasproperty(diagnostic, :stage_details) || return
    stage_details = diagnostic.stage_details
    @test stage_details isa Tuple
    @test all(detail -> detail isa NamedTuple, stage_details)
end

struct _DiagnosticNonunitEntry end
struct _DiagnosticRing end
struct _DiagnosticNormalizedEntry end
struct _DiagnosticNormalizedRing end

Suslin.is_unit(::_DiagnosticNonunitEntry) = false

function Suslin.normalize_laurent_object(column::Vector{_DiagnosticNonunitEntry})
    return (;
        normalized_object = [_DiagnosticNormalizedEntry() for _ in eachindex(column)],
        metadata = (; polynomial_ring = _DiagnosticNormalizedRing()),
    )
end

function Suslin.is_unimodular_column(::Vector{_DiagnosticNormalizedEntry}, ::_DiagnosticNormalizedRing)
    error("forced diagnostic normalized-column check failure")
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
    _test_diagnostic_stage_details_shape(case010)
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
    _test_diagnostic_stage_details_shape(case008_d21)
    @test length(case008_d21.stage_details) == length(case008_d21.attempted_stages)
    case008_d21_witness = _diagnostic_stage_detail(case008_d21, :laurent_witness_unit)
    @test case008_d21_witness !== nothing
    @test case008_d21_witness.outcome == :supported
    @test case008_d21_witness.witness_unit_index isa Integer

    case008_d16_fixture = ToricBuilderCase008D16ColumnBoundary.boundary_fixture()
    case008_d16 = Suslin.diagnose_unimodular_column_reduction(
        case008_d16_fixture.failing_column,
        case008_d16_fixture.ring,
    )
    @test case008_d16.status == :supported
    @test case008_d16.failure_code === nothing
    @test case008_d16.column_length == 16
    @test :laurent_elementary_row_preconditioning in case008_d16.attempted_stages
    @test !(:case008_special_case in case008_d16.attempted_stages)
    case008_d16_preconditioned =
        _diagnostic_stage_detail(case008_d16, :laurent_elementary_row_preconditioning)
    @test case008_d16_preconditioned !== nothing
    @test case008_d16_preconditioned.outcome == :supported
    @test case008_d16_preconditioned.target_index isa Integer
    @test case008_d16_preconditioned.source_index isa Integer
    @test case008_d16_preconditioned.coefficient == one(case008_d16_fixture.ring)
    case008_d16_cert = Suslin.ecp_column_reduction_certificate(
        case008_d16_fixture.failing_column,
        case008_d16_fixture.ring,
    )
    @test _diagnostic_supported_stage(case008_d16) == case008_d16_cert.stages[end].kind
    @test _diagnostic_stage_detail(case008_d16, :case008_special_case) === nothing

    R, (norm_x, norm_y) = Suslin.suslin_laurent_polynomial_ring(GF(2), ["x", "y"])
    normalization_column = [
        norm_x^-1 * norm_y^-1 * (norm_x + norm_y^2),
        norm_x^-1 * norm_y^-1 * (norm_x * norm_y + norm_x + one(R)),
        norm_x^-1 * norm_y^-1 * (norm_x^2 + norm_x * norm_y + norm_y + one(R)),
        zero(R),
        zero(R),
        zero(R),
    ]
    normalization_cert = Suslin.ecp_column_reduction_certificate(normalization_column, R)
    normalization_diagnostic =
        Suslin.diagnose_unimodular_column_reduction(normalization_column, R)
    @test normalization_cert.stages[end].kind == :laurent_normalization
    @test normalization_diagnostic.status == :supported
    normalization_detail =
        _diagnostic_stage_detail(normalization_diagnostic, :laurent_normalization)
    @test normalization_detail !== nothing
    @test normalization_detail.outcome == :delegated_to_polynomial
    @test normalization_detail.normalized_status == :supported
    @test :laurent_normalization in normalization_diagnostic.attempted_stages
    @test !(:laurent_elementary_row_preconditioning in normalization_diagnostic.attempted_stages)

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
    _test_diagnostic_stage_details_shape(unsupported)
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
    _test_diagnostic_stage_details_shape(supported)
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
    _test_diagnostic_stage_details_shape(precondition)
    @test isempty(precondition.stage_details)
end

@testset "Ordinary column reduction diagnostic stage details" begin
    R, (x, y) = Oscar.polynomial_ring(GF(2), ["x", "y"])

    unit = Suslin.diagnose_unimodular_column_reduction([x, y, one(R)], R)
    @test unit.status == :supported
    unit_detail = _diagnostic_stage_detail(unit, :unit_entry)
    @test unit_detail !== nothing
    @test unit_detail.outcome == :supported
    @test unit_detail.pivot_index == 3

    witness = Suslin.diagnose_unimodular_column_reduction([x, y, x + one(R)], R)
    @test witness.status == :supported
    witness_detail = _diagnostic_stage_detail(witness, :witness_unit)
    @test witness_detail !== nothing
    @test witness_detail.outcome == :supported
    @test witness_detail.witness_unit_index isa Integer

    monic_column = [
        x + y^2,
        x * y + x + one(R),
        x^2 + x * y + y + one(R),
    ]
    monic = Suslin.diagnose_unimodular_column_reduction(monic_column, R)
    @test monic.status == :supported
    monic_detail = _diagnostic_stage_detail(monic, :monicity_normalization)
    @test monic_detail !== nothing
    @test monic_detail.outcome == :supported
    @test monic_detail.normalized_column_length == 3

    block_column = [
        x + y,
        x * y,
        x^2,
        x^2 + x * y + y + one(R),
    ]
    @test Suslin.is_unimodular_column(block_column, R)
    block = Suslin.diagnose_unimodular_column_reduction(block_column, R)
    @test block.status == :supported
    block_detail = _diagnostic_stage_detail(block, :three_entry_block)
    @test block_detail !== nothing
    @test block_detail.outcome == :supported
    @test block_detail.block_indices == (1, 2, 4)
    @test block_detail.pivot_index == 4

    unsupported_column = [zero(R), x^2, x * y + one(R), zero(R)]
    @test Suslin.is_unimodular_column(unsupported_column, R)
    unsupported = Suslin.diagnose_unimodular_column_reduction(unsupported_column, R)
    @test unsupported.status == :unsupported
    no_block_detail = _diagnostic_stage_detail(unsupported, :three_entry_block)
    @test no_block_detail !== nothing
    @test no_block_detail.outcome == :no_supported_three_block
end

@testset "Defensive diagnostic stage details" begin
    unsupported_attempted = Symbol[]
    unsupported_details = Any[]
    unsupported_result = Suslin._diagnose_supported_unimodular_column_reduction(
        [_DiagnosticNonunitEntry(), _DiagnosticNonunitEntry(), _DiagnosticNonunitEntry()],
        _DiagnosticRing(),
        unsupported_attempted,
        unsupported_details,
    )
    @test !unsupported_result.supported
    @test unsupported_attempted == [:unit_entry, :witness_unit]
    unsupported_diagnostic = (; stage_details = tuple(unsupported_details...))
    unavailable_witness = _diagnostic_stage_detail(unsupported_diagnostic, :witness_unit)
    @test unavailable_witness !== nothing
    @test unavailable_witness.outcome == :witness_unavailable
    @test unavailable_witness.witness_unit_index === nothing

    laurent_attempted = Symbol[]
    laurent_details = Any[]
    laurent_result = Suslin._diagnose_laurent_unimodular_column_reduction(
        [_DiagnosticNonunitEntry(), _DiagnosticNonunitEntry(), _DiagnosticNonunitEntry()],
        _DiagnosticRing(),
        laurent_attempted,
        laurent_details,
    )
    @test !laurent_result.supported
    @test laurent_attempted == [
        :unit_entry,
        :laurent_unit_creation,
        :laurent_witness_unit,
        :laurent_normalization,
    ]
    laurent_diagnostic = (; stage_details = tuple(laurent_details...))
    laurent_witness = _diagnostic_stage_detail(laurent_diagnostic, :laurent_witness_unit)
    @test laurent_witness !== nothing
    @test laurent_witness.outcome == :witness_unavailable
    @test laurent_witness.witness_unit_index === nothing
    row_preconditioning =
        _diagnostic_stage_detail(laurent_diagnostic, :laurent_elementary_row_preconditioning)
    @test row_preconditioning === nothing
    normalization = _diagnostic_stage_detail(laurent_diagnostic, :laurent_normalization)
    @test normalization !== nothing
    @test normalization.outcome == :normalized_unimodularity_check_failed
    @test normalization.normalized_status == :precondition_failed
    @test normalization.normalized_failure_code == :unimodularity_check_failed
    @test occursin("forced diagnostic normalized-column check failure", normalization.normalized_message)
end
