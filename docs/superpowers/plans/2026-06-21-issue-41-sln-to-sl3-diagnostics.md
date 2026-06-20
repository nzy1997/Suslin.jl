# Issue 41 SL_n To SL_3 Diagnostics Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add structured diagnostics for failed `SL_n -> local SL_3` reductions, including Issue #38 local-shape and `3+3` partition-search evidence.

**Architecture:** Keep the existing reduction and staged `ArgumentError` text intact. Add public diagnostic record types and `diagnose_sln_to_sl3_reduction` beside `reduce_sln_to_sl3`, reusing existing normalization, block-location, local-shape, local-solver, and exact-reconstruction helpers. Bound alternative partition search to `6 x 6` matrices and only run it after a primary local failure.

**Tech Stack:** Julia, Oscar, Test standard library, existing Suslin internal helpers.

## Global Constraints

- Do not make the Issue #38 matrix factorize.
- Keep current staged reduction error text recognizable, including `failed to solve local SL_3 obligation on block [1, 2, 3]`.
- Follow TDD: write the failing diagnostics test first, run it and confirm the expected missing-API failure, then write production code.
- Alternative `3+3` partition search is bounded to `6 x 6` matrices and must not run for supported `40 x 40` fixtures.
- The Issue #38 row and column determinant-one cores must report `:determinant_one`, `:local_shape_failure`, `:not_embedded_2x2_with_trailing_identity`, searched partitions, and zero successful partitions.
- The synthetic `40 x 40` block-local Laurent acceptance fixture must report success or no failure diagnostics.
- Run `julia --project=. -e 'include("test/expert/sln_to_sl3_diagnostics.jl")'`.
- Run `julia --project=. -e 'using Pkg; Pkg.test()'`.

---

## File Structure

- Modify `src/Suslin.jl`: export the new diagnostic types and function.
- Modify `src/algorithm/sln_to_sl3_reduction.jl`: add diagnostic record types and helpers near the reduction path.
- Create `test/expert/sln_to_sl3_diagnostics.jl`: focused Issue #41 tests.
- Modify `test/runtests.jl`: register the new expert test.
- Modify `test/public/api_surface.jl`: verify the new public names are exported.

### Task 1: Structured SL_n To SL_3 Diagnostics

**Files:**
- Modify: `src/Suslin.jl`
- Modify: `src/algorithm/sln_to_sl3_reduction.jl`
- Create: `test/expert/sln_to_sl3_diagnostics.jl`
- Modify: `test/runtests.jl`
- Modify: `test/public/api_surface.jl`

**Interfaces:**
- Consumes: `reduce_sln_to_sl3(A; block_locations=nothing)`, `_construct_sln_to_sl3_reduction(A, block_locations)`, `_normalize_reduction_block_locations`, `_principal_submatrix`, `realize_sl3_local`, `normalize_laurent_gl_matrix`, `test/fixtures/toricbuilder_issue38_cases.jl`, and `test/fixtures/laurent_large_acceptance_cases.jl`.
- Produces:
  - `SL3LocalReductionDiagnostic`
  - `SLNToSL3ReductionDiagnostic`
  - `diagnose_sln_to_sl3_reduction(A; block_locations=nothing, search_partitions=true)`

- [ ] **Step 1: Write the failing diagnostics test**

Create `test/expert/sln_to_sl3_diagnostics.jl` with focused tests:

```julia
using Test
using Suslin
using Oscar

include("../fixtures/toricbuilder_issue38_cases.jl")
include("../fixtures/laurent_large_acceptance_cases.jl")

function _issue41_captured_error(f)
    try
        f()
        return nothing
    catch err
        return err
    end
end

function _issue41_failure_diagnostics(diagnostic)
    return filter(block -> block.status == :failure, diagnostic.block_diagnostics)
end

function _issue41_assert_issue38_diagnostic(core, label::AbstractString)
    diagnostic = diagnose_sln_to_sl3_reduction(core)
    @test diagnostic isa SLNToSL3ReductionDiagnostic
    @test diagnostic.status == :failure
    @test diagnostic.failure_code == :local_shape_failure
    @test diagnostic.determinant_status == :determinant_one
    @test diagnostic.determinant_classification == :one

    failures = _issue41_failure_diagnostics(diagnostic)
    @test length(failures) >= 1
    failure = first(failures)
    @test failure isa SL3LocalReductionDiagnostic
    @test failure.block_location == [1, 2, 3]
    @test failure.failure_code == :local_shape_failure
    @test failure.local_shape_reason == :not_embedded_2x2_with_trailing_identity
    @test failure.solver_status == :not_attempted

    @test diagnostic.partition_search.searched
    @test diagnostic.partition_search.status == :no_success
    @test diagnostic.partition_search.attempted_count == 10
    @test isempty(diagnostic.partition_search.successful_partitions)

    err = _issue41_captured_error(() -> reduce_sln_to_sl3(core))
    @test err isa ArgumentError
    message = sprint(showerror, err)
    @test occursin("staged SL_n to local SL_3 reduction failure", message)
    @test occursin("failed to solve local SL_3 obligation on block [1, 2, 3]", message)

    return diagnostic
end

@testset "Issue 38 Laurent core diagnostics" begin
    entry = only(ToricBuilderIssue38Cases.catalog().cases)
    row_core = entry.normalizations.row.core
    column_core = entry.normalizations.column.core

    row_diagnostic = _issue41_assert_issue38_diagnostic(row_core, "row")
    column_diagnostic = _issue41_assert_issue38_diagnostic(column_core, "column")

    @test row_diagnostic.message !== nothing
    @test column_diagnostic.message !== nothing
end

@testset "Supported Laurent block-local diagnostics" begin
    catalog = LaurentLargeAcceptanceCases.acceptance_catalog()
    case = only(filter(entry -> entry.id == "laurent-block-local-40x40", catalog.cases))

    diagnostic = diagnose_sln_to_sl3_reduction(case.matrix; block_locations = case.block_locations)
    @test diagnostic isa SLNToSL3ReductionDiagnostic
    @test diagnostic.status == :success
    @test diagnostic.failure_code === nothing
    @test diagnostic.determinant_status == :determinant_one
    @test isempty(_issue41_failure_diagnostics(diagnostic))
    @test !diagnostic.partition_search.searched
    @test diagnostic.partition_search.status == :not_applicable
end
```

Register the new file in the `"expert"` group in `test/runtests.jl`, after `expert/sln_to_sl3_reduction.jl`.

- [ ] **Step 2: Run the focused test and confirm RED**

Run:

```bash
julia --project=. -e 'include("test/expert/sln_to_sl3_diagnostics.jl")'
```

Expected before implementation: failure because `diagnose_sln_to_sl3_reduction` or `SLNToSL3ReductionDiagnostic` is not defined.

- [ ] **Step 3: Implement diagnostic records and public export**

In `src/Suslin.jl`, add exports:

```julia
export SL3LocalReductionDiagnostic
export SLNToSL3ReductionDiagnostic
export diagnose_sln_to_sl3_reduction
```

In `src/algorithm/sln_to_sl3_reduction.jl`, add record types near the existing reduction types and implement diagnostics with these behaviors:

```julia
struct SL3LocalReductionDiagnostic
    block_location::Vector{Int}
    status::Symbol
    failure_code::Union{Nothing, Symbol}
    determinant_status::Symbol
    local_shape_reason::Symbol
    solver_status::Symbol
    message::Union{Nothing, String}
end

struct SLNToSL3ReductionDiagnostic
    status::Symbol
    failure_code::Union{Nothing, Symbol}
    ring_profile::Symbol
    determinant_status::Symbol
    determinant_classification::Union{Nothing, Symbol}
    block_diagnostics::Vector{SL3LocalReductionDiagnostic}
    partition_search
    message::Union{Nothing, String}
end

function diagnose_sln_to_sl3_reduction(A; block_locations=nothing, search_partitions::Bool=true)
    return _diagnose_sln_to_sl3_reduction(A, block_locations, search_partitions)
end
```

Use helper functions to keep the implementation small:

- `_sln_determinant_status(normalized_A, ring_profile, normalization)` returns `:determinant_one`, `:determinant_not_one`, `:determinant_requires_correction`, or `:determinant_check_failed`.
- `_sl3_local_determinant_status(local_target, R)` returns `:determinant_one`, `:determinant_not_one`, or `:determinant_check_failed`.
- `_sl3_local_shape_reason(local_target, R)` returns `:embedded_2x2_with_trailing_identity` or `:not_embedded_2x2_with_trailing_identity`.
- `_diagnose_sl3_local_obligation(A, R, indices, X, ring_profile)` returns one `SL3LocalReductionDiagnostic`.
- `_three_by_three_partitions(n)` returns the 10 unordered `3+3` partitions only when `n == 6`.
- `_diagnose_three_by_three_partition_search(original_A, normalized_A, search_partitions, primary_failed)` returns a named tuple with `searched`, `status`, `attempted_count`, and `successful_partitions`.

The main diagnostic helper should:

1. Reuse existing validation, normalization, determinant, generator, and block-location checks.
2. Return a determinant failure diagnostic instead of throwing when the determinant-one precondition is not met after normalization.
3. Diagnose non-identity local blocks in order.
4. Stop at the first local failure for the primary path, matching the staged reduction failure behavior.
5. If all local diagnostics succeed, call `_construct_sln_to_sl3_reduction(A, block_locations)` to verify exact reconstruction. Return success if it works, otherwise return `:reassembly_failure`.
6. Run bounded `3+3` partition search only after a primary failure and only for `6 x 6` inputs.

- [ ] **Step 4: Update API surface tests**

In `test/public/api_surface.jl`, add `isdefined` and export-identity checks for:

```julia
SL3LocalReductionDiagnostic
SLNToSL3ReductionDiagnostic
diagnose_sln_to_sl3_reduction
```

- [ ] **Step 5: Run the focused diagnostic test and confirm GREEN**

Run:

```bash
julia --project=. -e 'include("test/expert/sln_to_sl3_diagnostics.jl")'
```

Expected after implementation: all Issue #41 diagnostic tests pass.

- [ ] **Step 6: Run package tests**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: public and internal tests pass.

- [ ] **Step 7: Commit the implementation**

Stage only the intended files and commit:

```bash
git add src/Suslin.jl src/algorithm/sln_to_sl3_reduction.jl test/expert/sln_to_sl3_diagnostics.jl test/runtests.jl test/public/api_surface.jl docs/superpowers/plans/2026-06-21-issue-41-sln-to-sl3-diagnostics.md
git commit -m "Add SLn to SL3 failure diagnostics"
```

## Self Review

- Spec coverage: the task covers public diagnostics, local shape reason, determinant status, bounded partition search, existing error text, Issue #38 row and column cores, and the `40 x 40` negative control.
- Placeholder scan: no TBD/TODO placeholders remain.
- Type consistency: exported names match the test and production signatures.
