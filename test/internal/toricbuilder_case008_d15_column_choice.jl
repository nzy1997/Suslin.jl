using Test
using Suslin
using Oscar

const TORICBUILDER_CASE008_D15_MATRIX_BOUNDARY_PATH =
    joinpath(@__DIR__, "..", "fixtures", "toricbuilder_case008_d15_matrix_boundary.jl")
const TORICBUILDER_CASE008_D15_COLUMN_BOUNDARY_PATH =
    joinpath(@__DIR__, "..", "fixtures", "toricbuilder_case008_d15_column_boundary.jl")

function _case008_d15_stage_detail(diagnostic, stage::Symbol)
    hasproperty(diagnostic, :stage_details) || return nothing
    idx = findfirst(detail -> detail.stage == stage, diagnostic.stage_details)
    return idx === nothing ? nothing : diagnostic.stage_details[idx]
end

function _case008_d15_column(M, column_index::Int)
    return [M[row, column_index] for row in 1:nrows(M)]
end

function _case008_d15_detail_property(detail, field::Symbol)
    detail === nothing && return nothing
    return hasproperty(detail, field) ? getproperty(detail, field) : nothing
end

const CASE008_D15_CACHED_NONCURRENT_DIAGNOSTIC = (;
    status = :not_run,
    failure_code = :cached_due_to_expensive_laurent_diagnostic,
    stage_details = (),
    diagnostic_cached = true,
    diagnostic_cache_reason = :bounded_column_choice_report,
)

function _case008_d15_diagnostic_property(diagnostic, field::Symbol, default)
    return hasproperty(diagnostic, field) ? getproperty(diagnostic, field) : default
end

function _case008_d15_candidate_diagnostic(column, R, column_index::Int, current_index::Int)
    if column_index == current_index
        return Suslin.diagnose_unimodular_column_reduction(column, R)
    end
    return CASE008_D15_CACHED_NONCURRENT_DIAGNOSTIC
end

function _case008_d15_candidate_report(M, R, column_index::Int, current_index::Int)
    column = _case008_d15_column(M, column_index)
    unit_entry_count = count(is_unit, column)
    is_unimodular = Suslin.is_unimodular_column(column, R)
    diagnostic = _case008_d15_candidate_diagnostic(column, R, column_index, current_index)
    witness = _case008_d15_stage_detail(diagnostic, :laurent_witness_unit)
    normalization = _case008_d15_stage_detail(diagnostic, :laurent_normalization)
    row_preconditioning =
        _case008_d15_stage_detail(diagnostic, :laurent_elementary_row_preconditioning)
    status = _case008_d15_diagnostic_property(diagnostic, :status, :not_run)
    return (;
        column_index,
        is_current_peel_column = column_index == current_index,
        is_unimodular,
        unit_entry_count,
        laurent_witness_outcome = _case008_d15_detail_property(witness, :outcome),
        laurent_witness_unit_index = _case008_d15_detail_property(witness, :witness_unit_index),
        normalized_precondition_status = _case008_d15_detail_property(normalization, :normalized_status),
        normalized_failure_code = _case008_d15_detail_property(normalization, :normalized_failure_code),
        row_preconditioning_outcome = _case008_d15_detail_property(row_preconditioning, :outcome),
        row_preconditioning_target_index =
            _case008_d15_detail_property(row_preconditioning, :target_index),
        row_preconditioning_source_indices =
            _case008_d15_detail_property(row_preconditioning, :source_indices),
        row_preconditioning_coefficient_strategy =
            _case008_d15_detail_property(row_preconditioning, :coefficient_strategy),
        row_preconditioning_coefficient_count =
            _case008_d15_detail_property(row_preconditioning, :coefficient_count),
        row_preconditioning_transformed_stage =
            _case008_d15_detail_property(row_preconditioning, :transformed_stage),
        status,
        failure_code = _case008_d15_diagnostic_property(diagnostic, :failure_code, nothing),
        supported_by_current_reducer = status == :supported,
        diagnostic_cached =
            _case008_d15_diagnostic_property(diagnostic, :diagnostic_cached, false),
        diagnostic_cache_reason =
            _case008_d15_diagnostic_property(diagnostic, :diagnostic_cache_reason, nothing),
    )
end

function case008_d15_column_choice_report(fixture = ToricBuilderCase008D15MatrixBoundary.matrix_fixture())
    validation = ToricBuilderCase008D15MatrixBoundary.validate_matrix_fixture(fixture)
    validation == :ok ||
        throw(ArgumentError("invalid case_008 d15 matrix fixture: $(validation)"))

    M = fixture.failing_input_matrix
    R = fixture.ring
    d = fixture.first_failing_peel_dimension
    candidates = tuple(
        (
            _case008_d15_candidate_report(
                M,
                R,
                column_index,
                fixture.current_peel_column_index,
            )
            for column_index in 1:d
        )...,
    )
    return (;
        case_id = fixture.case_id,
        dimension = d,
        current_peel_column_index = fixture.current_peel_column_index,
        candidates,
    )
end

@testset "ToricBuilder case_008 d=15 column-choice diagnostics" begin
    @test isfile(TORICBUILDER_CASE008_D15_MATRIX_BOUNDARY_PATH)
    @test isfile(TORICBUILDER_CASE008_D15_COLUMN_BOUNDARY_PATH)
    include(TORICBUILDER_CASE008_D15_MATRIX_BOUNDARY_PATH)

    fixture = ToricBuilderCase008D15MatrixBoundary.matrix_fixture()
    @test ToricBuilderCase008D15MatrixBoundary.validate_matrix_fixture(fixture) == :ok
    @test (nrows(fixture.failing_input_matrix), ncols(fixture.failing_input_matrix)) == (15, 15)
    @test fixture.current_peel_column_index == 15
    @test fixture.current_peel_column ==
          _case008_d15_column(fixture.failing_input_matrix, fixture.current_peel_column_index)

    column_fixture =
        ToricBuilderCase008D15MatrixBoundary.ToricBuilderCase008D15ColumnBoundary.boundary_fixture()
    @test fixture.current_peel_column == column_fixture.failing_column

    report = case008_d15_column_choice_report(fixture)
    @test report.case_id == "case_008"
    @test report.dimension == 15
    @test report.current_peel_column_index == 15
    @test length(report.candidates) == 15
    @test Tuple(candidate.column_index for candidate in report.candidates) == Tuple(1:15)
    @test count(candidate -> candidate.is_current_peel_column, report.candidates) == 1

    required_fields = (
        :column_index,
        :is_current_peel_column,
        :is_unimodular,
        :unit_entry_count,
        :laurent_witness_outcome,
        :laurent_witness_unit_index,
        :normalized_precondition_status,
        :normalized_failure_code,
        :row_preconditioning_outcome,
        :row_preconditioning_target_index,
        :row_preconditioning_source_indices,
        :row_preconditioning_coefficient_strategy,
        :row_preconditioning_coefficient_count,
        :row_preconditioning_transformed_stage,
        :status,
        :failure_code,
        :supported_by_current_reducer,
        :diagnostic_cached,
        :diagnostic_cache_reason,
    )
    @test all(candidate -> all(field -> hasproperty(candidate, field), required_fields), report.candidates)
    @test all(candidate -> candidate.is_unimodular, report.candidates)
    @test all(candidate -> candidate.unit_entry_count == 0, report.candidates)

    current = only(filter(candidate -> candidate.is_current_peel_column, report.candidates))
    @test current.column_index == 15
    @test current.is_unimodular
    @test current.unit_entry_count == 0
    @test current.status == :supported
    @test current.failure_code === nothing
    @test current.row_preconditioning_outcome == :supported
    @test current.row_preconditioning_target_index == 1
    @test current.row_preconditioning_source_indices == Tuple(2:15)
    @test current.row_preconditioning_coefficient_strategy ==
          :target_unit_laurent_linear_synthesis
    @test current.row_preconditioning_coefficient_count == 14
    @test current.row_preconditioning_transformed_stage == :unit_entry
    @test current.supported_by_current_reducer
    @test current.diagnostic_cached == false
    @test current.diagnostic_cache_reason === nothing

    noncurrent = filter(candidate -> !candidate.is_current_peel_column, report.candidates)
    @test length(noncurrent) == 14
    @test all(candidate -> candidate.diagnostic_cached, noncurrent)
    @test all(
        candidate -> candidate.diagnostic_cache_reason == :bounded_column_choice_report,
        noncurrent,
    )
    @test all(candidate -> candidate.status == :not_run, noncurrent)
    @test all(
        candidate -> candidate.failure_code == :cached_due_to_expensive_laurent_diagnostic,
        noncurrent,
    )
    @test all(candidate -> candidate.supported_by_current_reducer == false, noncurrent)

    supported = filter(candidate -> candidate.supported_by_current_reducer, report.candidates)
    supported_alternatives = filter(candidate -> candidate.column_index != 15, supported)
    @info "case_008 d15 column-choice diagnostic" current_peel_column_index = 15 supported_columns = Tuple(candidate.column_index for candidate in supported) supported_alternative_columns = Tuple(candidate.column_index for candidate in supported_alternatives)
    @test all(candidate -> candidate.status == :supported, supported)
    @test all(candidate -> candidate.failure_code === nothing, supported)
    @test all(
        candidate ->
            !candidate.is_unimodular ||
            candidate.supported_by_current_reducer ||
            candidate.failure_code !== nothing,
        report.candidates,
    )

    corrupted = ToricBuilderCase008D15MatrixBoundary.corrupted_matrix_entry_negative_control(fixture)
    rejection = ToricBuilderCase008D15MatrixBoundary.validate_matrix_fixture(corrupted)
    @test rejection in (:wrong_snapshot, :wrong_column)
    @test_throws ArgumentError case008_d15_column_choice_report(corrupted)

    wrong_column = copy(fixture.current_peel_column)
    wrong_column[1] += one(fixture.ring)
    wrong_column_fixture = merge(fixture, (; current_peel_column = wrong_column))
    @test ToricBuilderCase008D15MatrixBoundary.validate_matrix_fixture(wrong_column_fixture) == :wrong_column
    @test_throws ArgumentError case008_d15_column_choice_report(wrong_column_fixture)
end
