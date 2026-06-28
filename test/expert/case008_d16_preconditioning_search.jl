using Test
using Suslin
using Oscar

include(joinpath(@__DIR__, "..", "fixtures", "toricbuilder_case008_d16_matrix_boundary.jl"))

function _case008_d16_column(M, column_index::Int)
    return [M[row, column_index] for row in 1:nrows(M)]
end

function _checked_preconditioning_depth(max_depth)::Int
    max_depth isa Integer || throw(ArgumentError("max_depth must be an integer"))
    depth = Int(max_depth)
    depth >= 0 || throw(ArgumentError("max_depth must be nonnegative"))
    return depth
end

function _checked_preconditioning_side(side)::Symbol
    side == :right || throw(ArgumentError("only :right column-addition search is supported"))
    return side
end

function _checked_preconditioning_operation_family(operation_family)::Symbol
    operation_family == :column_addition ||
        throw(ArgumentError("only :column_addition search is supported"))
    return operation_family
end

function _checked_preconditioning_index(index, limit::Int, label::AbstractString)::Int
    index isa Integer || throw(ArgumentError("$label must be an integer"))
    idx = Int(index)
    1 <= idx <= limit || throw(ArgumentError("$label must be between 1 and $limit"))
    return idx
end

function _checked_source_column_candidates(candidates, column_index::Int, limit::Int)
    checked = Int[]
    for candidate in candidates
        idx = _checked_preconditioning_index(candidate, limit, "source column")
        idx == column_index && throw(ArgumentError("source column candidates must not include the target column"))
        push!(checked, idx)
    end
    return tuple(checked...)
end

function _checked_coefficient_candidates(R, coefficient_candidates)
    coefficient_candidates === nothing && return (one(R),)
    return tuple((R(coefficient) for coefficient in coefficient_candidates)...)
end

function _case008_d16_preconditioning_bounds(
    M,
    R;
    max_depth,
    side,
    operation_family,
    column_index,
    source_column_candidates,
    coefficient_candidates,
)
    target_column = _checked_preconditioning_index(column_index, ncols(M), "column_index")
    sources = source_column_candidates === nothing ?
        Tuple(idx for idx in 1:ncols(M) if idx != target_column) :
        _checked_source_column_candidates(source_column_candidates, target_column, ncols(M))
    return (;
        max_depth = _checked_preconditioning_depth(max_depth),
        side = _checked_preconditioning_side(side),
        operation_family = _checked_preconditioning_operation_family(operation_family),
        column_index = target_column,
        source_column_candidates = sources,
        coefficient_candidates = _checked_coefficient_candidates(R, coefficient_candidates),
    )
end

function _diagnose_preconditioned_column(M, R, column_index::Int)
    column = _case008_d16_column(M, column_index)
    return column, Suslin.diagnose_unimodular_column_reduction(column, R)
end

function _not_found_preconditioning_result(M, R, bounds, attempt_count::Int)
    column, diagnostic = _diagnose_preconditioned_column(M, R, bounds.column_index)
    return (;
        status = :not_found,
        bounds,
        attempt_count,
        steps = (),
        transformed_column = column,
        reducer_diagnostic = diagnostic,
    )
end

function case008_d16_preconditioning_search(
    fixture = ToricBuilderCase008D16MatrixBoundary.matrix_fixture();
    max_depth = 1,
    side = :right,
    operation_family = :column_addition,
    coefficient_candidates = nothing,
    column_index = nothing,
    source_column_candidates = nothing,
)::NamedTuple
    validation = ToricBuilderCase008D16MatrixBoundary.validate_matrix_fixture(fixture)
    validation == :ok ||
        throw(ArgumentError("invalid case_008 d16 matrix fixture: $(validation)"))

    M = fixture.failing_input_matrix
    R = fixture.ring
    target_column = column_index === nothing ? fixture.current_peel_column_index : column_index
    bounds = _case008_d16_preconditioning_bounds(
        M,
        R;
        max_depth,
        side,
        operation_family,
        column_index = target_column,
        source_column_candidates,
        coefficient_candidates,
    )

    attempt_count = 1
    initial_column, initial_diagnostic = _diagnose_preconditioned_column(M, R, bounds.column_index)
    if initial_diagnostic.status == :supported
        return (;
            status = :already_supported,
            bounds,
            attempt_count,
            steps = (),
            transformed_column = initial_column,
            reducer_diagnostic = initial_diagnostic,
        )
    end
    frontier = [(matrix = M, steps = ())]

    for _depth in 1:bounds.max_depth
        next_frontier = NamedTuple[]
        for state in frontier
            for source in bounds.source_column_candidates
                for coefficient in bounds.coefficient_candidates
                    step = Suslin.elementary_preconditioning_step(
                        state.matrix,
                        bounds.side,
                        bounds.column_index,
                        source,
                        coefficient,
                    )
                    steps = tuple(state.steps..., step)
                    transformed_matrix = step.transformed_matrix
                    Suslin.verify_elementary_preconditioning(M, steps, transformed_matrix) ||
                        error("internal preconditioning replay invariant failed")

                    attempt_count += 1
                    transformed_column, diagnostic =
                        _diagnose_preconditioned_column(transformed_matrix, R, bounds.column_index)
                    if diagnostic.status == :supported
                        return (;
                            status = :found,
                            bounds,
                            attempt_count,
                            steps,
                            transformed_column,
                            reducer_diagnostic = diagnostic,
                        )
                    end
                    push!(next_frontier, (matrix = transformed_matrix, steps = steps))
                end
            end
        end
        frontier = next_frontier
    end

    return _not_found_preconditioning_result(M, R, bounds, attempt_count)
end

function _preconditioning_result_is_verified(original_matrix, result)::Bool
    try
        result.status == :found || return false
        result.attempt_count > 0 || return false
        !isempty(result.steps) || return false
        hasproperty(result.bounds, :column_index) || return false

        replayed = Suslin.replay_elementary_preconditioning(original_matrix, result.steps)
        Suslin.verify_elementary_preconditioning(original_matrix, result.steps, replayed) ||
            return false
        _case008_d16_column(replayed, result.bounds.column_index) == result.transformed_column ||
            return false
        actual_diagnostic = Suslin.diagnose_unimodular_column_reduction(
            result.transformed_column,
            base_ring(original_matrix),
        )
        actual_diagnostic.status == :supported || return false
        return result.reducer_diagnostic.status == :supported
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

@testset "case_008 d=16 bounded preconditioning search" begin
    fixture = ToricBuilderCase008D16MatrixBoundary.matrix_fixture()
    @test ToricBuilderCase008D16MatrixBoundary.validate_matrix_fixture(fixture) == :ok

    result = case008_d16_preconditioning_search(fixture; max_depth = 1)

    @test result.status in (:already_supported, :found, :not_found)
    @test result.bounds.max_depth == 1
    @test result.bounds.side == :right
    @test result.bounds.operation_family == :column_addition
    @test result.bounds.column_index == fixture.current_peel_column_index
    @test result.bounds.source_column_candidates == Tuple(1:15)
    @test result.bounds.coefficient_candidates == (one(fixture.ring),)
    @test result.attempt_count > 0
    @test hasproperty(result.reducer_diagnostic, :status)

    if result.status == :already_supported
        @test result.steps == ()
        @test result.transformed_column == fixture.current_peel_column
        @test result.reducer_diagnostic.status == :supported
        @test result.reducer_diagnostic.failure_code === nothing
        @test :laurent_elementary_row_preconditioning in
              result.reducer_diagnostic.attempted_stages
    elseif result.status == :found
        final_matrix = Suslin.replay_elementary_preconditioning(
            fixture.failing_input_matrix,
            result.steps,
        )
        @test _preconditioning_result_is_verified(fixture.failing_input_matrix, result)
        @test Suslin.verify_elementary_preconditioning(
            fixture.failing_input_matrix,
            result.steps,
            final_matrix,
        )
        @test _case008_d16_column(final_matrix, result.bounds.column_index) ==
              result.transformed_column
        @test result.reducer_diagnostic.status == :supported

        tampered_steps = collect(result.steps)
        tampered_steps[1] = _tamper_preconditioning_factor(tampered_steps[1])
        @test !Suslin.verify_elementary_preconditioning(
            fixture.failing_input_matrix,
            tuple(tampered_steps...),
            final_matrix,
        )
    else
        @test result.status == :not_found
        @test result.steps == ()
        @test result.transformed_column == fixture.current_peel_column
        @test result.reducer_diagnostic.status == :supported
        @test result.reducer_diagnostic.failure_code === nothing
        @test :laurent_elementary_row_preconditioning in
              result.reducer_diagnostic.attempted_stages
    end
end

@testset "case_008 d=16 preconditioning search short-circuits supported column" begin
    fixture = ToricBuilderCase008D16MatrixBoundary.matrix_fixture()
    already_supported = case008_d16_preconditioning_search(fixture; max_depth = 0)

    @test already_supported.status == :already_supported
    @test already_supported.bounds.max_depth == 0
    @test already_supported.bounds.side == :right
    @test already_supported.bounds.operation_family == :column_addition
    @test already_supported.bounds.column_index == 16
    @test already_supported.bounds.source_column_candidates == Tuple(1:15)
    @test already_supported.bounds.coefficient_candidates == (one(fixture.ring),)
    @test already_supported.attempt_count == 1
    @test already_supported.steps == ()
    @test already_supported.transformed_column == fixture.current_peel_column
    @test already_supported.reducer_diagnostic.status == :supported
    @test already_supported.reducer_diagnostic.failure_code === nothing
    @test :laurent_elementary_row_preconditioning in
          already_supported.reducer_diagnostic.attempted_stages
end

@testset "preconditioning search negative controls reject tampering" begin
    R, (u, v) = Suslin.suslin_laurent_polynomial_ring(GF(2), ["u", "v"])
    A = matrix(R, [
        one(R)  zero(R)  zero(R);
        zero(R) one(R)   zero(R);
        zero(R) zero(R)  one(R)
    ])

    step = Suslin.elementary_preconditioning_step(A, :right, 3, 1, u)
    final_matrix = step.transformed_matrix
    @test Suslin.verify_elementary_preconditioning(A, (step,), final_matrix)

    tampered_step = _tamper_preconditioning_factor(step)
    @test !Suslin.verify_elementary_preconditioning(A, (tampered_step,), final_matrix)

    supported_column = _case008_d16_column(final_matrix, 3)
    supported_diagnostic = Suslin.diagnose_unimodular_column_reduction(supported_column, R)
    @test supported_diagnostic.status == :supported

    unsupported_found = (;
        status = :found,
        bounds = (; column_index = 3),
        attempt_count = 1,
        steps = (step,),
        transformed_column = supported_column,
        reducer_diagnostic = (; status = :unsupported),
    )
    @test !_preconditioning_result_is_verified(A, unsupported_found)

    valid_found = merge(
        unsupported_found,
        (;
            reducer_diagnostic = supported_diagnostic,
        ),
    )
    @test _preconditioning_result_is_verified(A, valid_found)
end
