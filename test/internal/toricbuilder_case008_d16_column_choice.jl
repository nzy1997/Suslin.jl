using Test
using Suslin
using Oscar

const TORICBUILDER_CASE008_D16_MATRIX_BOUNDARY_PATH =
    joinpath(@__DIR__, "..", "fixtures", "toricbuilder_case008_d16_matrix_boundary.jl")

function _case008_d16_stage_detail(diagnostic, stage::Symbol)
    hasproperty(diagnostic, :stage_details) || return nothing
    idx = findfirst(detail -> detail.stage == stage, diagnostic.stage_details)
    return idx === nothing ? nothing : diagnostic.stage_details[idx]
end

function _case008_d16_column(M, column_index::Int)
    return [M[row, column_index] for row in 1:nrows(M)]
end

function _case008_d16_detail_property(detail, field::Symbol)
    detail === nothing && return nothing
    return hasproperty(detail, field) ? getproperty(detail, field) : nothing
end

function _case008_d16_candidate_report(M, R, column_index::Int, current_index::Int)
    column = _case008_d16_column(M, column_index)
    unit_entry_count = count(is_unit, column)
    is_unimodular = Suslin.is_unimodular_column(column, R)
    diagnostic = Suslin.diagnose_unimodular_column_reduction(column, R)
    witness = _case008_d16_stage_detail(diagnostic, :laurent_witness_unit)
    normalization = _case008_d16_stage_detail(diagnostic, :laurent_normalization)
    row_preconditioning =
        _case008_d16_stage_detail(diagnostic, :laurent_elementary_row_preconditioning)
    return (;
        column_index,
        is_current_peel_column = column_index == current_index,
        is_unimodular,
        unit_entry_count,
        laurent_witness_outcome = _case008_d16_detail_property(witness, :outcome),
        laurent_witness_unit_index = _case008_d16_detail_property(witness, :witness_unit_index),
        normalized_precondition_status = _case008_d16_detail_property(normalization, :normalized_status),
        normalized_failure_code = _case008_d16_detail_property(normalization, :normalized_failure_code),
        row_preconditioning_outcome = _case008_d16_detail_property(row_preconditioning, :outcome),
        row_preconditioning_transformed_stage =
            _case008_d16_detail_property(row_preconditioning, :transformed_stage),
        status = diagnostic.status,
        failure_code = diagnostic.failure_code,
        supported_by_current_reducer = diagnostic.status == :supported,
    )
end

function case008_d16_column_choice_report(fixture = ToricBuilderCase008D16MatrixBoundary.matrix_fixture())
    validation = ToricBuilderCase008D16MatrixBoundary.validate_matrix_fixture(fixture)
    validation == :ok ||
        throw(ArgumentError("invalid case_008 d16 matrix fixture: $(validation)"))

    M = fixture.failing_input_matrix
    R = fixture.ring
    d = fixture.first_failing_peel_dimension
    candidates = tuple(
        (
            _case008_d16_candidate_report(
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

@testset "ToricBuilder case_008 d=16 column-choice diagnostics" begin
    @test isfile(TORICBUILDER_CASE008_D16_MATRIX_BOUNDARY_PATH)
    include(TORICBUILDER_CASE008_D16_MATRIX_BOUNDARY_PATH)

    fixture = ToricBuilderCase008D16MatrixBoundary.matrix_fixture()
    @test ToricBuilderCase008D16MatrixBoundary.validate_matrix_fixture(fixture) == :ok
    @test (nrows(fixture.failing_input_matrix), ncols(fixture.failing_input_matrix)) == (16, 16)
    @test fixture.current_peel_column_index == 16
    @test fixture.current_peel_column ==
          _case008_d16_column(fixture.failing_input_matrix, fixture.current_peel_column_index)

    report = case008_d16_column_choice_report(fixture)
    @test report.case_id == "case_008"
    @test report.dimension == 16
    @test report.current_peel_column_index == 16
    @test length(report.candidates) == 16
    @test Tuple(candidate.column_index for candidate in report.candidates) == Tuple(1:16)
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
        :row_preconditioning_transformed_stage,
        :status,
        :failure_code,
        :supported_by_current_reducer,
    )
    @test all(candidate -> all(field -> hasproperty(candidate, field), required_fields), report.candidates)

    current = only(filter(candidate -> candidate.is_current_peel_column, report.candidates))
    @test current.column_index == 16
    @test current.is_unimodular
    @test current.unit_entry_count == 0
    @test current.status == :supported
    @test current.failure_code === nothing
    @test current.laurent_witness_outcome == :witness_without_unit
    @test current.normalized_precondition_status == :precondition_failed
    @test current.normalized_failure_code == :not_unimodular
    @test current.row_preconditioning_outcome == :supported
    @test current.row_preconditioning_transformed_stage == :witness_unit
    @test current.supported_by_current_reducer

    supported = filter(candidate -> candidate.supported_by_current_reducer, report.candidates)
    supported_alternatives = filter(candidate -> candidate.column_index != 16, supported)
    @test Tuple(candidate.column_index for candidate in supported) == (8, 9, 10, 16)
    @test Tuple(candidate.column_index for candidate in supported_alternatives) == (8, 9, 10)
    @info "case_008 d16 column-choice diagnostic" current_peel_column_index = 16 supported_columns = Tuple(candidate.column_index for candidate in supported) supported_alternative_columns = Tuple(candidate.column_index for candidate in supported_alternatives)
    @test all(candidate -> candidate.status == :supported, supported)

    if isempty(supported)
        @test all(
            candidate ->
                !candidate.is_unimodular ||
                (!candidate.supported_by_current_reducer && candidate.failure_code !== nothing),
            report.candidates,
        )
    else
        @test all(candidate -> candidate.failure_code === nothing, supported)
    end

    corrupted = ToricBuilderCase008D16MatrixBoundary.corrupted_matrix_entry_negative_control(fixture)
    rejection = ToricBuilderCase008D16MatrixBoundary.validate_matrix_fixture(corrupted)
    @test rejection in (:wrong_snapshot, :wrong_column)
    @test_throws ArgumentError case008_d16_column_choice_report(corrupted)
end
