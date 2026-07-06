# Issue 322 Case008 D14 Laurent Elementary Move Search Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a deterministic expert-only bounded `case_008 d=14` Laurent elementary entry-addition search report.

**Architecture:** Create one expert test/report file that includes the existing d14 measure contract, evaluates the declared operation bounds through exact support-set arithmetic, and replay-verifies every recorded decreasing candidate from the original polynomial column. Register the file in `test/runtests.jl`; do not touch production reducer code.

**Tech Stack:** Julia, Test stdlib, Suslin/Oscar test fixtures, existing d14 profile and measure helpers.

## Global Constraints

- Keep the search report expert-only and test-only.
- Do not modify production code under `src/`.
- Do not add a production Laurent reducer.
- Do not require a local ToricBuilder checkout.
- Use the checked-in `case_008 d=14` boundary fixture from #315.
- Use the descent measure contract from #321.
- Default operation family must be `:entry_addition`.
- Default ordered pairs must cover `target_index, source_index in 1:14` with `target_index != source_index`.
- Default exponent vectors must cover `(a, b)` with `-1 <= a <= 1` and `-1 <= b <= 1`.
- Default coefficient family must be `(1,)` over `GF(2)`.
- Default `checked_operation_count` must be exactly `14 * 13 * 9 == 1638`.
- Operation replay semantics must be `target_entry <- target_entry + coefficient * u^a * v^b * source_entry`.
- Report `case_id` must equal `"case_008"`.
- Report `dimension` must equal `14`.
- Report `input_measure` must equal the d14 baseline measure.
- Report `operation_families` must name `(:entry_addition,)`.
- Report `search_bounds` must name the exact exponent radius and coefficient set.
- Report `status` must be either `:candidate_found` or `:exhausted`.
- Every recorded candidate must contain before/after measures and the elementary operation needed to replay it.
- No measure decrease may be accepted unless the transformed column is recomputed from the recorded elementary operation.
- If `status == :exhausted`, the report must record `candidate_count == 0` and `replay_verified_count == 0`.
- Required focused verification: `julia --project=. -e 'include("test/expert/case008_d14_laurent_elementary_move_search.jl")'`.
- Required package verification: `julia --project=. -e 'using Pkg; Pkg.test()'`.

---

## File Structure

- Create: `test/expert/case008_d14_laurent_elementary_move_search.jl`
  - Owns operation construction, operation replay, support-set measure evaluation, bounded search report generation, replay candidate verification, d14 assertions, and synthetic controls.
- Modify: `test/runtests.jl`
  - Registers `expert/case008_d14_laurent_elementary_move_search.jl` immediately after `expert/case008_d14_laurent_descent_measure_contract.jl`.

### Task 1: Expert Search Report

**Files:**
- Create: `test/expert/case008_d14_laurent_elementary_move_search.jl`
- Modify: `test/runtests.jl`

**Interfaces:**
- Consumes: `case008_d14_laurent_descent_profile(fixture) -> NamedTuple`.
- Consumes: `case008_d14_laurent_descent_measure(profile; fixture) -> NamedTuple`.
- Consumes: `strictly_decreases_laurent_measure(before, after)::Bool`.
- Produces: `replay_laurent_elementary_entry_addition(column, R, operation) -> Vector`.
- Produces: `verify_laurent_elementary_move_candidate(original_column, R, candidate)::Bool`.
- Produces: `case008_d14_laurent_elementary_move_search_report(fixture = boundary_fixture()) -> NamedTuple`.

- [ ] **Step 1: Write the failing tests and runner registration**

Create `test/expert/case008_d14_laurent_elementary_move_search.jl` with a testset that includes the d14 measure contract, calls `case008_d14_laurent_elementary_move_search_report()`, and asserts:

```julia
@test report.case_id == "case_008"
@test report.dimension == 14
@test report.input_measure == baseline_measure
@test report.operation_families == (:entry_addition,)
@test report.search_bounds.exponent_radius == 1
@test report.search_bounds.coefficient_family == (1,)
@test report.checked_operation_count == 1638
@test report.status in (:candidate_found, :exhausted)
```

Also add synthetic tests that call `verify_laurent_elementary_move_candidate` for a known-good decreasing operation, a known-bad non-decreasing operation, a malformed operation, and a stale after-measure candidate.

Register the file in `test/runtests.jl` immediately after:

```julia
"expert/case008_d14_laurent_descent_measure_contract.jl",
```

Run:

```bash
julia --project=. -e 'include("test/expert/case008_d14_laurent_elementary_move_search.jl")'
```

Expected: fail before implementation with `UndefVarError: case008_d14_laurent_elementary_move_search_report not defined` or equivalent missing helper errors.

- [ ] **Step 2: Implement operation and measure helpers**

Add helper functions in the new expert file:

```julia
const CASE008_D14_ELEMENTARY_OPERATION_FAMILIES = (:entry_addition,)
const CASE008_D14_ELEMENTARY_EXPONENT_RADIUS = 1
const CASE008_D14_ELEMENTARY_COEFFICIENT_FAMILY = (1,)

function replay_laurent_elementary_entry_addition(column, R, operation)
    operation.family == :entry_addition ||
        throw(ArgumentError("unsupported operation family $(repr(operation.family))"))
    n = length(column)
    target = _checked_entry_index(operation.target_index, n, "target_index")
    source = _checked_entry_index(operation.source_index, n, "source_index")
    target != source || throw(ArgumentError("target_index and source_index must differ"))
    exponent = _checked_exponent_pair(operation.exponent)
    coefficient = R(operation.coefficient)
    generators = gens(R)
    monomial = coefficient * generators[1]^exponent[1] * generators[2]^exponent[2]
    transformed = copy(column)
    transformed[target] = transformed[target] + monomial * column[source]
    return transformed
end
```

Implement support helpers that convert `_entry_support(entry)` into `Set{Tuple{Int, Int}}`, update a target support by symmetric difference with the shifted source support, and build a #321-compatible measure from support sets with component order `CASE008_D14_MEASURE_COMPONENTS`.

- [ ] **Step 3: Implement candidate verification and report generation**

Implement:

```julia
function verify_laurent_elementary_move_candidate(original_column, R, candidate)::Bool
    # validate candidate fields, replay operation, recompute before/after measures,
    # compare with recorded measures, and require strictly_decreases_laurent_measure.
end
```

Implement:

```julia
function case008_d14_laurent_elementary_move_search_report(
    fixture = ToricBuilderCase008D14ColumnBoundary.boundary_fixture(),
)::NamedTuple
    # validate fixture, compute baseline profile/measure, scan all 1638 default
    # operations, replay-verify decreasing candidates, and return the stable report.
end
```

The scan must loop in deterministic order: target index ascending, source index
ascending with target skipped, exponent `a` ascending, exponent `b` ascending,
then coefficient family order.

- [ ] **Step 4: Run focused red/green verification**

Run:

```bash
julia --project=. -e 'include("test/expert/case008_d14_laurent_elementary_move_search.jl")'
```

Expected: pass. The output must include the new report test and assert `checked_operation_count == 1638`.

- [ ] **Step 5: Run broader verification and commit**

Run:

```bash
git diff --check
julia --project=. test/runtests.jl expert
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: all commands exit 0. Do not commit `Manifest.toml` if it is generated only as a local verification artifact and remains ignored by git.

Commit:

```bash
git add docs/superpowers/specs/2026-07-06-issue-322-case008-d14-laurent-elementary-move-search-design.md docs/superpowers/plans/2026-07-06-issue-322-case008-d14-laurent-elementary-move-search.md test/expert/case008_d14_laurent_elementary_move_search.jl test/runtests.jl
git commit -m "test: add case008 d14 Laurent elementary move search"
```

## Self-Review

- Spec coverage: task covers report shape, default bounds, exact operation count, replay verification, synthetic controls, runner registration, focused verification, and package verification.
- Placeholder scan: no `TBD` or deferred implementation requirements remain.
- Type consistency: produced helper names match the task interfaces and tests.
