using Test
using Suslin
using Oscar

include(joinpath(@__DIR__, "..", "fixtures", "toricbuilder_case008_d15_matrix_boundary.jl"))

function _case008_d15_column(M, column_index::Int)
    return [M[row, column_index] for row in 1:nrows(M)]
end

function _checked_nonnegative_integer(value, label::AbstractString)::Int
    value isa Integer || throw(ArgumentError("$label must be an integer"))
    checked = Int(value)
    checked >= 0 || throw(ArgumentError("$label must be nonnegative"))
    return checked
end

function _checked_index(index, limit::Int, label::AbstractString)::Int
    index isa Integer || throw(ArgumentError("$label must be an integer"))
    checked = Int(index)
    1 <= checked <= limit || throw(ArgumentError("$label must be between 1 and $limit"))
    return checked
end

function _checked_side(side)::Symbol
    side isa Symbol || throw(ArgumentError("side must be a symbol"))
    side in (:right, :left) || throw(ArgumentError("unsupported side $(repr(side)); expected :right or :left"))
    return side
end

function _checked_side_order(side_order)
    side_order === nothing && return (:right, :left)
    raw_values = side_order isa Symbol ? (side_order,) : Tuple(side_order)
    return tuple((_checked_side(side) for side in raw_values)...)
end

function _checked_operation_family(operation_family)::Symbol
    operation_family == :column_addition ||
        throw(ArgumentError("only :column_addition search is supported"))
    return operation_family
end

function _checked_source_column_candidates(candidates, column_index::Int, limit::Int)
    candidates === nothing && return Tuple(idx for idx in 1:limit if idx != column_index)
    checked = Int[]
    for candidate in candidates
        source = _checked_index(candidate, limit, "source column")
        source == column_index &&
            throw(ArgumentError("source column candidates must not include the target column"))
        push!(checked, source)
    end
    return tuple(checked...)
end

function _checked_row_synthesis_pivots(pivots, limit::Int)
    pivots === nothing && return Tuple(1:limit)
    checked = Int[]
    for pivot in pivots
        push!(checked, _checked_index(pivot, limit, "row synthesis pivot"))
    end
    return tuple(checked...)
end

function _checked_coefficient_candidates(R, coefficient_candidates)
    coefficient_candidates === nothing && return (one(R),)
    return tuple((R(coefficient) for coefficient in coefficient_candidates)...)
end

function _case008_d15_preconditioning_bounds(
    M,
    R;
    max_depth,
    side_order,
    operation_family,
    coefficient_candidates,
    column_index,
    source_column_candidates,
    right_full_diagnostic_limit,
    row_synthesis_pivots,
    row_synthesis_max_steps,
)
    target_column = _checked_index(column_index, ncols(M), "column_index")
    max_depth_checked = _checked_nonnegative_integer(max_depth, "max_depth")
    max_depth_checked <= 1 ||
        throw(ArgumentError("max_depth values greater than 1 are not supported"))
    return (;
        max_depth = max_depth_checked,
        side_order = _checked_side_order(side_order),
        operation_family = _checked_operation_family(operation_family),
        column_index = target_column,
        source_column_candidates = _checked_source_column_candidates(
            source_column_candidates,
            target_column,
            ncols(M),
        ),
        coefficient_candidates = _checked_coefficient_candidates(R, coefficient_candidates),
        right_full_diagnostic_limit = _checked_nonnegative_integer(
            right_full_diagnostic_limit,
            "right_full_diagnostic_limit",
        ),
        row_synthesis_pivots = _checked_row_synthesis_pivots(row_synthesis_pivots, nrows(M)),
        row_synthesis_max_steps = _checked_nonnegative_integer(
            row_synthesis_max_steps,
            "row_synthesis_max_steps",
        ),
    )
end

function _not_found_preconditioning_result(
    original_matrix,
    R,
    bounds,
    attempt_count::Int,
    progress_summary,
)
    transformed_column = _case008_d15_column(original_matrix, bounds.column_index)
    diagnostic = Suslin.diagnose_unimodular_column_reduction(transformed_column, R)
    return (;
        status = :not_found,
        bounds,
        attempt_count,
        progress_summary = tuple(progress_summary...),
        steps = (),
        transformed_column,
        reducer_diagnostic = diagnostic,
    )
end

function _row_synthesis_steps(M, column, pivot::Int, max_steps::Int)
    R = base_ring(M)
    sources = [idx for idx in 1:length(column) if idx != pivot]
    A = matrix(R, 1, length(sources), [column[idx] for idx in sources])
    B = matrix(R, 1, 1, [one(R) + column[pivot]])
    solution = Suslin.solve_laurent_linear(A, B)
    coefficients = [solution[idx, 1] for idx in 1:nrows(solution)]
    nonzero_pairs = [(source, coefficients[idx]) for (idx, source) in enumerate(sources) if !iszero(coefficients[idx])]

    length(nonzero_pairs) <= max_steps || return (;
        status = :too_many_nonzero_coefficients,
        step_count = length(nonzero_pairs),
    )

    steps = NamedTuple[]
    current_matrix = M
    for (source, coefficient) in nonzero_pairs
        step = Suslin.elementary_preconditioning_step(
            current_matrix,
            :left,
            pivot,
            source,
            coefficient,
        )
        push!(steps, step)
        current_matrix = step.transformed_matrix
    end

    return (;
        status = :ok,
        step_count = length(steps),
        steps = tuple(steps...),
    )
end

function case008_d15_preconditioning_search(
    fixture = ToricBuilderCase008D15MatrixBoundary.matrix_fixture();
    max_depth = 1,
    side_order = (:right, :left),
    operation_family = :column_addition,
    coefficient_candidates = nothing,
    column_index = nothing,
    source_column_candidates = nothing,
    right_full_diagnostic_limit = 0,
    row_synthesis_pivots = nothing,
    row_synthesis_max_steps = 14,
)::NamedTuple
    validation = ToricBuilderCase008D15MatrixBoundary.validate_matrix_fixture(fixture)
    validation == :ok ||
        throw(ArgumentError("invalid case_008 d15 matrix fixture: $(validation)"))

    original_matrix = fixture.failing_input_matrix
    R = fixture.ring
    target_column = column_index === nothing ? fixture.current_peel_column_index : column_index
    bounds = _case008_d15_preconditioning_bounds(
        original_matrix,
        R;
        max_depth,
        side_order,
        operation_family,
        coefficient_candidates,
        column_index = target_column,
        source_column_candidates,
        right_full_diagnostic_limit,
        row_synthesis_pivots,
        row_synthesis_max_steps,
    )

    attempt_count = 1
    progress = NamedTuple[]
    right_diagnostic_budget = bounds.right_full_diagnostic_limit

    if bounds.max_depth == 0
        push!(progress, (
            side = :right,
            source = nothing,
            coefficient_index = nothing,
            outcome = :max_depth_zero,
            direct_unit_rows = (),
            normalized_unimodular = false,
        ))
        return _not_found_preconditioning_result(original_matrix, R, bounds, attempt_count, progress)
    end

    for side in bounds.side_order
        if side == :right
            if isempty(bounds.source_column_candidates)
                push!(progress, (
                    side = :right,
                    source = nothing,
                    coefficient_index = nothing,
                    outcome = :no_source_columns,
                    direct_unit_rows = (),
                    normalized_unimodular = false,
                ))
                continue
            end

            for source in bounds.source_column_candidates
                for (coefficient_index, coefficient) in enumerate(bounds.coefficient_candidates)
                    step = Suslin.elementary_preconditioning_step(
                        original_matrix,
                        :right,
                        bounds.column_index,
                        source,
                        coefficient,
                    )
                    attempt_count += 1

                    transformed_column = _case008_d15_column(
                        step.transformed_matrix,
                        bounds.column_index,
                    )
                    direct_unit_rows = Tuple(findall(is_unit, transformed_column))
                    normalization = Suslin.normalize_laurent_object(transformed_column)
                    normalized_unimodular = Suslin.is_unimodular_column(
                        normalization.normalized_object,
                        normalization.metadata.polynomial_ring,
                    )

                    outcome = :cheap_checks_failed
                    diagnostic = nothing
                    if !isempty(direct_unit_rows) || normalized_unimodular
                        if right_diagnostic_budget > 0
                            diagnostic = Suslin.diagnose_unimodular_column_reduction(
                                transformed_column,
                                R,
                            )
                            right_diagnostic_budget -= 1
                            outcome = diagnostic.status == :supported ?
                                :supported :
                                :full_diagnostic_not_supported
                        else
                            outcome = :expensive_witness_diagnostic_skipped
                        end
                    end

                    push!(progress, (
                        side = :right,
                        source,
                        coefficient_index,
                        outcome,
                        direct_unit_rows,
                        normalized_unimodular,
                    ))

                    if diagnostic !== nothing && diagnostic.status == :supported
                        final_matrix = Suslin.replay_elementary_preconditioning(
                            original_matrix,
                            (step,),
                        )
                        Suslin.verify_elementary_preconditioning(
                            original_matrix,
                            (step,),
                            final_matrix,
                        ) || error("internal preconditioning replay invariant failed")
                        transformed_column = _case008_d15_column(
                            final_matrix,
                            bounds.column_index,
                        )
                        return (;
                        status = :found,
                        bounds,
                        attempt_count,
                        progress_summary = tuple(progress...),
                        steps = (step,),
                        transformed_column,
                        reducer_diagnostic = diagnostic,
                        )
                    end
                end
            end
        elseif side == :left
            if isempty(bounds.row_synthesis_pivots)
                push!(progress, (
                    side = :left,
                    pivot = nothing,
                    outcome = :no_row_synthesis_pivots,
                    reason = :empty_candidates,
                ))
                continue
            end

            for pivot in bounds.row_synthesis_pivots
                attempt_count += 1
                transformed_column = _case008_d15_column(original_matrix, bounds.column_index)

                synthesis = try
                    _row_synthesis_steps(
                        original_matrix,
                        transformed_column,
                        pivot,
                        bounds.row_synthesis_max_steps,
                    )
                catch err
                    err isa InterruptException && rethrow()
                    push!(progress, (
                        side = :left,
                        pivot,
                        outcome = :solve_failed,
                        reason = sprint(showerror, err),
                    ))
                    continue
                end

                if synthesis.status != :ok
                    push!(progress, (
                        side = :left,
                        pivot,
                        outcome = synthesis.status,
                        step_count = synthesis.step_count,
                    ))
                    continue
                end

                steps = synthesis.steps
                final_matrix = Suslin.replay_elementary_preconditioning(original_matrix, steps)
                Suslin.verify_elementary_preconditioning(original_matrix, steps, final_matrix) ||
                    error("internal preconditioning replay invariant failed")

                transformed_column = _case008_d15_column(final_matrix, bounds.column_index)
                diagnostic = Suslin.diagnose_unimodular_column_reduction(transformed_column, R)
                outcome = diagnostic.status == :supported ? :supported : :unsupported_diagnostic
                push!(progress, (
                    side = :left,
                    pivot,
                    outcome,
                    step_count = synthesis.step_count,
                ))

                if diagnostic.status == :supported
                    return (;
                    status = :found,
                    bounds,
                    attempt_count,
                    progress_summary = tuple(progress...),
                    steps,
                    transformed_column,
                    reducer_diagnostic = diagnostic,
                    )
                end
            end
        else
            error("unsupported side order entry $(repr(side))")
        end
    end

    return _not_found_preconditioning_result(original_matrix, R, bounds, attempt_count, progress)
end

function _case008_d15_preconditioning_result_is_verified(original_matrix, result)::Bool
    try
        result.status == :found || return false
        result.attempt_count > 0 || return false
        !isempty(result.steps) || return false
        hasproperty(result, :bounds) || return false
        hasproperty(result.bounds, :column_index) || return false
        hasproperty(result, :transformed_column) || return false
        hasproperty(result, :reducer_diagnostic) || return false
        result.reducer_diagnostic.status == :supported || return false

        replayed = Suslin.replay_elementary_preconditioning(original_matrix, result.steps)
        Suslin.verify_elementary_preconditioning(original_matrix, result.steps, replayed) ||
            return false

        replayed_column = _case008_d15_column(replayed, result.bounds.column_index)
        replayed_column == result.transformed_column || return false

        actual_diagnostic = Suslin.diagnose_unimodular_column_reduction(
            result.transformed_column,
            base_ring(original_matrix),
        )
        actual_diagnostic.status == :supported || return false
        return true
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end

function _tamper_preconditioning_factor(step)
    tampered_factor = copy(step.factor)
    R = base_ring(tampered_factor)
    row, col = step.side == :left ? (step.target, step.source) : (step.source, step.target)
    tampered_factor[row, col] += one(R)
    return merge(step, (; factor = tampered_factor))
end

@testset "case_008 d=15 bounded preconditioning search" begin
    fixture = ToricBuilderCase008D15MatrixBoundary.matrix_fixture()
    @test ToricBuilderCase008D15MatrixBoundary.validate_matrix_fixture(fixture) == :ok

    result = case008_d15_preconditioning_search(fixture)

    @test result.status == :found
    @test result.bounds.max_depth == 1
    @test result.bounds.side_order == (:right, :left)
    @test result.bounds.operation_family == :column_addition
    @test result.bounds.column_index == fixture.current_peel_column_index
    @test result.bounds.source_column_candidates == Tuple(1:14)
    @test result.bounds.coefficient_candidates == (one(fixture.ring),)
    @test result.bounds.right_full_diagnostic_limit == 0
    @test result.bounds.row_synthesis_pivots == Tuple(1:15)
    @test result.bounds.row_synthesis_max_steps == 14
    @test result.attempt_count > 0
    @test !isempty(result.steps)
    @test all(step -> step.side == :left, result.steps)
    @test all(step -> step.target == 1, result.steps)
    @test result.transformed_column[1] |> is_unit
    @test result.reducer_diagnostic.status == :supported
    @test all(record -> hasproperty(record, :outcome), result.progress_summary)
    @test _case008_d15_preconditioning_result_is_verified(
        fixture.failing_input_matrix,
        result,
    )

    final_matrix = Suslin.replay_elementary_preconditioning(
        fixture.failing_input_matrix,
        result.steps,
    )
    @test Suslin.verify_elementary_preconditioning(
        fixture.failing_input_matrix,
        result.steps,
        final_matrix,
    )
    @test _case008_d15_column(final_matrix, result.bounds.column_index) == result.transformed_column

    tampered_steps = collect(result.steps)
    tampered_steps[1] = _tamper_preconditioning_factor(tampered_steps[1])
    @test !Suslin.verify_elementary_preconditioning(
        fixture.failing_input_matrix,
        tuple(tampered_steps...),
        final_matrix,
    )

    unsupported_found = merge(
        result,
        (;
            reducer_diagnostic = (; status = :unsupported),
        ),
    )
    @test !_case008_d15_preconditioning_result_is_verified(
        fixture.failing_input_matrix,
        unsupported_found,
    )
end

@testset "case_008 d=15 max depth bound" begin
    fixture = ToricBuilderCase008D15MatrixBoundary.matrix_fixture()
    max_depth_zero = case008_d15_preconditioning_search(fixture; max_depth = 0)

    @test max_depth_zero.status == :not_found
    @test isempty(max_depth_zero.steps)
    @test any(
        record -> record.outcome == :max_depth_zero,
        max_depth_zero.progress_summary,
    )
    @test_throws ArgumentError case008_d15_preconditioning_search(fixture; max_depth = 2)
end

@testset "case_008 d=15 empty candidate controls" begin
    fixture = ToricBuilderCase008D15MatrixBoundary.matrix_fixture()

    result = case008_d15_preconditioning_search(
        fixture;
        source_column_candidates = (),
        row_synthesis_pivots = (),
    )

    @test result.status == :not_found
    @test result.bounds.source_column_candidates == ()
    @test result.bounds.row_synthesis_pivots == ()
    @test result.attempt_count > 0
    @test any(record -> record.outcome == :no_source_columns, result.progress_summary)
    @test any(record -> record.outcome == :no_row_synthesis_pivots, result.progress_summary)
end
