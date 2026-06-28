# Issue 168 Case008 D16 Preconditioning Search Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add the expert-only bounded elementary preconditioning search requested by issue #168 for the `case_008` `d=16` fixture.

**Architecture:** Keep the search harness local to `test/expert/case008_d16_preconditioning_search.jl`. The helper builds exact side-aware elementary preconditioning steps, diagnoses each transformed target column with the existing reducer diagnostic, and returns either a replay-verified found record or a bounded `:not_found` record.

**Tech Stack:** Julia, Test stdlib, Suslin.jl exported elementary preconditioning helpers, Suslin.jl reducer diagnostics, Oscar Laurent polynomial matrices.

## Global Constraints

- Repository has no `AGENTS.md`, `CLAUDE.md`, or `GEMINI.md`; follow existing Julia test style in `test/expert` and `test/internal`.
- Do not add a production reducer stage, exported API, or new source module.
- Create `test/expert/case008_d16_preconditioning_search.jl`.
- Reuse `elementary_preconditioning_step`, `replay_elementary_preconditioning`, and `verify_elementary_preconditioning`.
- Return records must include `status`, `bounds`, `attempt_count`, `steps`, `transformed_column`, and `reducer_diagnostic`.
- The focused issue verification command is `julia --project=. -e 'include("test/expert/case008_d16_preconditioning_search.jl")'`.
- Required Agent Desk verification command is `julia --project=. -e 'using Pkg; Pkg.test()'`.
- Keep the test expert-only and do not register it in `test/runtests.jl`.

---

## File Structure

- Create `test/expert/case008_d16_preconditioning_search.jl`: local bounded search helpers, result verifier, d16 fixture test, `:not_found` bounds test, and replay tamper negative controls.
- Keep `src/` unchanged.

---

### Task 1: Add Expert Bounded Preconditioning Search Test

**Files:**
- Create: `test/expert/case008_d16_preconditioning_search.jl`

**Interfaces:**
- Consumes: `ToricBuilderCase008D16MatrixBoundary.matrix_fixture()`, `ToricBuilderCase008D16MatrixBoundary.validate_matrix_fixture(fixture)`, `Suslin.elementary_preconditioning_step`, `Suslin.replay_elementary_preconditioning`, `Suslin.verify_elementary_preconditioning`, and `Suslin.diagnose_unimodular_column_reduction`.
- Produces: local helper `case008_d16_preconditioning_search(fixture; max_depth, side, operation_family, coefficient_candidates, column_index, source_column_candidates)::NamedTuple`.
- Produces: local helper `_preconditioning_result_is_verified(original_matrix, result)::Bool`.

- [ ] **Step 1: Write the failing expert test shell**

Create `test/expert/case008_d16_preconditioning_search.jl` with this initial content:

```julia
using Test
using Suslin
using Oscar

include(joinpath(@__DIR__, "..", "fixtures", "toricbuilder_case008_d16_matrix_boundary.jl"))

@testset "case_008 d=16 bounded preconditioning search" begin
    fixture = ToricBuilderCase008D16MatrixBoundary.matrix_fixture()
    @test ToricBuilderCase008D16MatrixBoundary.validate_matrix_fixture(fixture) == :ok

    result = case008_d16_preconditioning_search(fixture; max_depth = 0)

    @test result.status == :not_found
    @test result.attempt_count > 0
    @test result.bounds.max_depth == 0
    @test result.steps == ()
    @test result.transformed_column == fixture.current_peel_column
end
```

- [ ] **Step 2: Run the expert test and verify RED**

Run:

```bash
julia --project=. -e 'include("test/expert/case008_d16_preconditioning_search.jl")'
```

Expected when dependencies are available: FAIL with `UndefVarError: case008_d16_preconditioning_search not defined`.

If the command fails before loading the test because Oscar is not installed, record that environment blocker and continue implementing from the test-first shell. Do not change production code to work around the missing dependency.

- [ ] **Step 3: Replace the shell with the full expert search test**

Replace `test/expert/case008_d16_preconditioning_search.jl` with this complete implementation:

```julia
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
)
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
    _diagnose_preconditioned_column(M, R, bounds.column_index)
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

    @test result.status in (:found, :not_found)
    @test result.bounds.max_depth == 1
    @test result.bounds.side == :right
    @test result.bounds.operation_family == :column_addition
    @test result.bounds.column_index == fixture.current_peel_column_index
    @test result.bounds.source_column_candidates == Tuple(1:15)
    @test result.bounds.coefficient_candidates == (one(fixture.ring),)
    @test result.attempt_count > 0
    @test hasproperty(result.reducer_diagnostic, :status)

    if result.status == :found
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
        @test result.reducer_diagnostic.status == :unsupported
        @test result.reducer_diagnostic.failure_code == :unsupported_laurent_column_family
    end
end

@testset "case_008 d=16 preconditioning not-found bounds are explicit" begin
    fixture = ToricBuilderCase008D16MatrixBoundary.matrix_fixture()
    not_found = case008_d16_preconditioning_search(fixture; max_depth = 0)

    @test not_found.status == :not_found
    @test not_found.bounds.max_depth == 0
    @test not_found.bounds.side == :right
    @test not_found.bounds.operation_family == :column_addition
    @test not_found.bounds.column_index == 16
    @test not_found.bounds.source_column_candidates == Tuple(1:15)
    @test not_found.bounds.coefficient_candidates == (one(fixture.ring),)
    @test not_found.attempt_count == 1
    @test not_found.steps == ()
    @test not_found.transformed_column == fixture.current_peel_column
    @test not_found.reducer_diagnostic.status == :unsupported
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
```

- [ ] **Step 4: Run the focused issue test**

Run:

```bash
julia --project=. -e 'include("test/expert/case008_d16_preconditioning_search.jl")'
```

Expected when dependencies are available: PASS. If a d16 candidate is found, the found branch verifies replay and tampering. If no candidate is found, the result reports `:not_found`, `attempt_count > 0`, and exact bounds. If dependencies are unavailable, record the exact dependency error.

- [ ] **Step 5: Run required package verification**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected when dependencies are available: PASS. If dependencies are unavailable, record the exact dependency error.

- [ ] **Step 6: Commit implementation**

Run:

```bash
git add test/expert/case008_d16_preconditioning_search.jl docs/superpowers/plans/2026-06-28-issue-168-case008-d16-preconditioning-search.md
git commit -m "test: search case008 d16 preconditioning candidates"
```

Expected: commit succeeds and contains only the expert search test and this plan.

---

## Self-Review

- Spec coverage: Task 1 covers the expert-only bounded search, result fields, deterministic bounds, replay verification, `:not_found` reporting, found-result diagnostic rejection, tamper negative controls, and required commands.
- Placeholder scan: no unresolved markers or incomplete implementation sections are present.
- Type consistency: the local helper returns `NamedTuple` records and all later tests use the same fields: `status`, `bounds`, `attempt_count`, `steps`, `transformed_column`, and `reducer_diagnostic`.
