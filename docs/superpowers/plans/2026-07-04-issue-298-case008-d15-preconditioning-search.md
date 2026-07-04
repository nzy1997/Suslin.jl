# Issue 298 Case008 D15 Preconditioning Search Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add the expert-only bounded preconditioning search requested by issue #298 for the `case_008 d=15` full-matrix fixture.

**Architecture:** Keep all search helpers local to `test/expert/case008_d15_preconditioning_search.jl`. The helper records cheap right-column probe progress, then performs bounded left row-synthesis using a Laurent linear solve to create a direct unit entry in target column `15`; every accepted result is replayed and verified with the exported elementary preconditioning helpers.

**Tech Stack:** Julia, Test stdlib, Suslin.jl exported elementary preconditioning and reducer diagnostic helpers, Oscar Laurent polynomial matrices.

## Global Constraints

- Repository has no `AGENTS.md`, `CLAUDE.md`, or `GEMINI.md`; follow existing Julia test style in `test/expert` and `test/internal`.
- Do not add a production reducer stage, exported API, or new source module.
- Create `test/expert/case008_d15_preconditioning_search.jl`.
- Keep the expert file out of `test/runtests.jl`.
- Use the #296 fixture `test/fixtures/toricbuilder_case008_d15_matrix_boundary.jl`.
- Verify every accepted candidate with `Suslin.replay_elementary_preconditioning` and `Suslin.verify_elementary_preconditioning`.
- The returned record must include `status`, `bounds`, `attempt_count`, `steps`, `transformed_column`, `reducer_diagnostic`, and `progress_summary`.
- The focused issue verification command is `julia --project=. -e 'include("test/expert/case008_d15_preconditioning_search.jl")'`.
- Required Agent Desk verification command is `julia --project=. -e 'using Pkg; Pkg.test()'`.

---

## File Structure

- Create `test/expert/case008_d15_preconditioning_search.jl`: local bounds validation, cheap right-probe progress, row-synthesis candidate construction, result validation, found-candidate assertions, and negative controls.
- Keep `src/` unchanged.

---

### Task 1: Add Failing Expert Test Shell

**Files:**
- Create: `test/expert/case008_d15_preconditioning_search.jl`

**Interfaces:**
- Consumes: `ToricBuilderCase008D15MatrixBoundary.matrix_fixture()`.
- Produces: an initial failing test that calls `case008_d15_preconditioning_search`.

- [ ] **Step 1: Write the failing test shell**

Create `test/expert/case008_d15_preconditioning_search.jl` with:

```julia
using Test
using Suslin
using Oscar

include(joinpath(@__DIR__, "..", "fixtures", "toricbuilder_case008_d15_matrix_boundary.jl"))

@testset "case_008 d=15 bounded preconditioning search" begin
    fixture = ToricBuilderCase008D15MatrixBoundary.matrix_fixture()
    @test ToricBuilderCase008D15MatrixBoundary.validate_matrix_fixture(fixture) == :ok

    result = case008_d15_preconditioning_search(fixture)

    @test result.status == :found
    @test result.attempt_count > 0
    @test !isempty(result.steps)
    @test result.reducer_diagnostic.status == :supported
end
```

- [ ] **Step 2: Run the expert test and verify RED**

Run:

```bash
julia --project=. -e 'include("test/expert/case008_d15_preconditioning_search.jl")'
```

Expected: FAIL with `UndefVarError: case008_d15_preconditioning_search not defined`.

- [ ] **Step 3: Commit the failing shell only if the repository convention requires red commits**

Do not commit the red shell in this repository. Continue to Task 2 and commit the final green expert test.

---

### Task 2: Implement Local Search Helper And Negative Controls

**Files:**
- Modify: `test/expert/case008_d15_preconditioning_search.jl`

**Interfaces:**
- Consumes: `Suslin.elementary_preconditioning_step`, `Suslin.replay_elementary_preconditioning`, `Suslin.verify_elementary_preconditioning`, `Suslin.diagnose_unimodular_column_reduction`, `Suslin.normalize_laurent_object`, `Suslin.is_unimodular_column`, and `Suslin.solve_laurent_linear`.
- Produces: `case008_d15_preconditioning_search(fixture; max_depth, side_order, operation_family, coefficient_candidates, column_index, source_column_candidates, right_full_diagnostic_limit, row_synthesis_pivots, row_synthesis_max_steps)::NamedTuple`.
- Produces: `_case008_d15_preconditioning_result_is_verified(original_matrix, result)::Bool`.

- [ ] **Step 1: Add local utility helpers**

Implement helpers with these exact responsibilities:

```julia
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
```

Add analogous tuple-normalization helpers for side order, source columns, row synthesis pivots, and coefficient candidates. Reject unsupported sides except `:right` and `:left`, and reject operation families except `:column_addition`.

- [ ] **Step 2: Add bounds and progress records**

Return bounds as a named tuple containing:

```julia
(;
    max_depth,
    side_order,
    operation_family,
    column_index,
    source_column_candidates,
    coefficient_candidates,
    right_full_diagnostic_limit,
    row_synthesis_pivots,
    row_synthesis_max_steps,
)
```

Represent progress as tuples of named tuples. Right probes should include `side = :right`, `source`, `coefficient_index`, `outcome`, `direct_unit_rows`, and `normalized_unimodular`. Row synthesis records should include `side = :left`, `pivot`, `outcome`, and either `step_count` or `reason`.

- [ ] **Step 3: Implement cheap right-side probes**

For each right-side source and coefficient, build one `Suslin.elementary_preconditioning_step`, increment `attempt_count`, collect the transformed target column, check direct unit rows, and check normalized ordinary-polynomial unimodularity:

```julia
normalization = Suslin.normalize_laurent_object(transformed_column)
normalized_unimodular = Suslin.is_unimodular_column(
    normalization.normalized_object,
    normalization.metadata.polynomial_ring,
)
```

If either cheap check is true and the full diagnostic budget allows it, run `Suslin.diagnose_unimodular_column_reduction`. Otherwise record `:expensive_witness_diagnostic_skipped` and continue to row synthesis.

- [ ] **Step 4: Implement row synthesis**

For each pivot in `row_synthesis_pivots`, solve:

```julia
sources = [idx for idx in 1:length(column) if idx != pivot]
A = matrix(R, 1, length(sources), [column[idx] for idx in sources])
B = matrix(R, 1, 1, [one(R) + column[pivot]])
solution = Suslin.solve_laurent_linear(A, B)
```

Build left-side steps with `Suslin.elementary_preconditioning_step(current_matrix, :left, pivot, source, coefficient)` for nonzero coefficients in ascending source order. The next step must use the previous step's transformed matrix. Reject the candidate if the nonzero coefficient count exceeds `row_synthesis_max_steps`.

- [ ] **Step 5: Verify and return the found row-synthesis candidate**

After building row-synthesis steps:

```julia
final_matrix = Suslin.replay_elementary_preconditioning(original_matrix, steps)
Suslin.verify_elementary_preconditioning(original_matrix, steps, final_matrix) ||
    error("internal preconditioning replay invariant failed")
transformed_column = _case008_d15_column(final_matrix, bounds.column_index)
diagnostic = Suslin.diagnose_unimodular_column_reduction(transformed_column, R)
```

Return `status = :found` only when `diagnostic.status == :supported`. The default search must find pivot `1`, produce nonempty left-side steps, and make `transformed_column[1]` a unit.

- [ ] **Step 6: Add validator and negative controls**

Implement `_case008_d15_preconditioning_result_is_verified(original_matrix, result)` to reject non-`:found` results, empty steps, replay mismatches, unsupported diagnostics, and transformed columns not supplied by the replayed matrix.

Add tests that:

- assert the default result is `:found`;
- assert replay verification succeeds and supplies `result.transformed_column`;
- tamper the first found step's factor and assert verification fails;
- construct a fake `:found` result with unsupported diagnostic and assert the validator rejects it;
- run `case008_d15_preconditioning_search(fixture; source_column_candidates = (), row_synthesis_pivots = ())` and assert `status == :not_found`, exact empty bounds are preserved, `attempt_count > 0`, and progress records include `:no_source_columns` and `:no_row_synthesis_pivots`.

- [ ] **Step 7: Run focused verification**

Run:

```bash
julia --project=. -e 'include("test/expert/case008_d15_preconditioning_search.jl")'
```

Expected: PASS. The command should report a found candidate with supported reducer diagnostic.

- [ ] **Step 8: Commit the expert test**

Run:

```bash
git add test/expert/case008_d15_preconditioning_search.jl
git commit -m "test: add case008 d15 preconditioning search"
```

---

### Task 3: Final Verification And PR Preparation

**Files:**
- No new files beyond Task 2.

**Interfaces:**
- Consumes: committed design, plan, and expert test.
- Produces: verified branch ready for PR.

- [ ] **Step 1: Run focused issue command**

Run:

```bash
julia --project=. -e 'include("test/expert/case008_d15_preconditioning_search.jl")'
```

Expected: exit 0.

- [ ] **Step 2: Run required package tests**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: exit 0.

- [ ] **Step 3: Inspect git status**

Run:

```bash
git status --short --branch
```

Expected: branch has committed changes and no unintended untracked files except ignored Julia manifests.

- [ ] **Step 4: Use finishing workflow**

Use `superpowers:verification-before-completion`, then `superpowers:finishing-a-development-branch`. When prompted, choose `Push and create a Pull Request` under the Standing Answer Policy.

## Self-Review

- Spec coverage: Task 2 covers the local helper, found record fields, replay verification, progress/skip summary, and all issue negative controls. Task 3 covers required verification.
- Red-flag scan: no deferred implementation markers remain.
- Type consistency: helper names, result fields, and bounds fields match the design document.
