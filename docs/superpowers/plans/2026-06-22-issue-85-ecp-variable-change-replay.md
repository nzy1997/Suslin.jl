# Issue 85 ECP Variable-Change Replay Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add explicit replay-checked variable-change data to ECP monicity-normalization certificates.

**Architecture:** Extend the existing `:monicity_normalization` certificate stage instead of adding a new reducer or changing the public API. The stage records deterministic substitution choices, selected monic entry data, first-coordinate status, transformed reduction proof, and inverse-substituted original reduction proof; replay recomputes every field exactly.

**Tech Stack:** Julia, Oscar polynomial rings and matrices, existing Suslin ECP certificate helpers, and Test stdlib.

## Global Constraints

- Preserve existing public behavior: `reduce_unimodular_column(v, R)` returns only the factor sequence.
- Do not add a general Quillen patching step.
- Do not add a broader variable-change search algorithm.
- Keep deterministic choices: record the selected source variable, target variable, shift polynomial, and variable order.
- Every recorded certificate field must participate in exact replay verification.
- The focused issue verification command is `julia --project=. -e 'include("test/expert/ecp_variable_change_replay.jl")'`.
- The required package verification command is `julia --project=. -e 'using Pkg; Pkg.test()'`.

---

## File Structure

- Modify `src/algorithm/column_reduction.jl`: enrich `:monicity_normalization` stage construction and replay checks; add small private helpers for substitution maps and selected monic entry detection.
- Modify `test/fixtures/ecp_column_cases.jl`: add one second catalog case that reaches the variable-change stage.
- Modify `test/internal/ecp_column_fixtures.jl`: require the new fixture id so the catalog test protects it.
- Create `test/expert/ecp_variable_change_replay.jl`: focused issue tests and tamper controls.
- Modify `test/runtests.jl`: register the new expert test.

## Task 1: Add Red Variable-Change Replay Coverage

**Files:**
- Modify: `test/fixtures/ecp_column_cases.jl`
- Modify: `test/internal/ecp_column_fixtures.jl`
- Create: `test/expert/ecp_variable_change_replay.jl`
- Modify: `test/runtests.jl`

**Interfaces:**
- Consumes: `ECPColumnFixtureCatalog.cases_by_id()`, `Suslin.ecp_column_reduction_certificate(v, R)`, `Suslin.verify_ecp_column_reduction(cert)`, and existing private expert helpers in `Suslin`.
- Produces: RED tests for the enriched monicity-normalization stage fields: `variable_order`, `source_variable_index`, `source_variable`, `target_variable_index`, `target_variable`, `shift_polynomial`, `forward_substitution`, `inverse_substitution`, `transformed_column`, `selected_monic_index`, `selected_monic_entry`, `first_coordinate_strategy`, `first_coordinate_move_factors`, `first_coordinate_column`, and `variable_change_verification`.

- [ ] **Step 1: Add the second catalog case**

In `test/fixtures/ecp_column_cases.jl`, after `variable_change_case`, add:

```julia
    variable_change_permuted_case = _case(
        id = "ecp-variable-change-permuted-gf2",
        kind = :variable_change,
        stage_coverage = :supported,
        ring_constructor = gf2_ring_constructor,
        ring = _ring_metadata("GF(2)[x, y]", R2, ("x", "y"), (x, y)),
        variable_order = (:x, :y),
        entries = (
            a = x + y^2,
            b = x * y + x + one(R2),
            c = x^2 + x * y + y + one(R2),
        ),
        column_order = (:b, :a, :c),
        monicity = (;
            selected_entry = :b,
            variable_name = :y,
            substitution = (; x = x + y^2),
            transformed_entry = x * y + x + y^3 + y^2 + one(R2),
        ),
        witnesses = (),
        expected = (; current_status = :passes),
        source_refs = ("Issue 85 variable-change replay coverage",),
        consumer_issue_ids = ("#62", "#85"),
    )
```

Add `variable_change_permuted_case` to the `cases` array immediately after `variable_change_case`.

- [ ] **Step 2: Require the new fixture id**

In `test/internal/ecp_column_fixtures.jl`, add `"ecp-variable-change-permuted-gf2"` to `REQUIRED_ECP_COLUMN_IDS`.

- [ ] **Step 3: Write the failing expert test**

Create `test/expert/ecp_variable_change_replay.jl` with:

```julia
using Test
using Suslin
using Oscar

include(joinpath(@__DIR__, "..", "fixtures", "ecp_column_cases.jl"))

const VARIABLE_CHANGE_REPLAY_CASE_IDS = (
    "ecp-variable-change-monic-gf2",
    "ecp-variable-change-permuted-gf2",
)

function _vc_column(entry)
    return [getproperty(entry.entries, name) for name in entry.column_order]
end

function _vc_target_column(R, n::Int)
    target = zero_matrix(R, n, 1)
    target[n, 1] = one(R)
    return target
end

function _vc_factor_product(factors, R, n::Int)
    product = identity_matrix(R, n)
    for factor in factors
        product *= factor
    end
    return product
end

function _vc_apply_factors(factors, column, R)
    return _vc_factor_product(factors, R, length(column)) * matrix(R, length(column), 1, collect(column))
end

function _vc_stage(cert)
    stages = [stage for stage in cert.stages if stage.kind == :monicity_normalization]
    @test length(stages) == 1
    return only(stages)
end

function _vc_map_values(substitution_map)
    return tuple((entry.value for entry in substitution_map)...)
end

function _vc_map_variables(substitution_map)
    return tuple((entry.variable for entry in substitution_map)...)
end

function _tamper_inverse_substitution_map(cert)
    stages = collect(cert.stages)
    stage_idx = findfirst(stage -> stage.kind == :monicity_normalization, stages)
    stage_idx === nothing && error("certificate has no variable-change stage")
    stage = stages[stage_idx]
    inverse_substitution = collect(stage.inverse_substitution)
    changed = inverse_substitution[stage.source_variable_index]
    inverse_substitution[stage.source_variable_index] = merge(
        changed,
        (; value = changed.value + one(cert.ring)),
    )
    stages[stage_idx] = merge(stage, (; inverse_substitution = tuple(inverse_substitution...)))
    return Suslin.ECPColumnReductionCertificate(
        cert.original_column,
        cert.ring,
        tuple(stages...),
        cert.factors,
        cert.final_column,
        cert.verification,
    )
end

function _tamper_selected_monic_index(cert)
    stages = collect(cert.stages)
    stage_idx = findfirst(stage -> stage.kind == :monicity_normalization, stages)
    stage_idx === nothing && error("certificate has no variable-change stage")
    stage = stages[stage_idx]
    bad_index = stage.selected_monic_index == length(stage.transformed_column) ? 1 : stage.selected_monic_index + 1
    stages[stage_idx] = merge(stage, (; selected_monic_index = bad_index))
    return Suslin.ECPColumnReductionCertificate(
        cert.original_column,
        cert.ring,
        tuple(stages...),
        cert.factors,
        cert.final_column,
        cert.verification,
    )
end

function _assert_variable_change_stage(entry)
    column = _vc_column(entry)
    R = entry.ring.object
    cert = Suslin.ecp_column_reduction_certificate(column, R)
    stage = _vc_stage(cert)
    ring_gens = tuple(gens(R)...)
    target = _vc_target_column(R, length(column))

    @test !any(is_unit, column)
    @test Suslin._reduce_supported_unimodular_column_certificate(column, R) === nothing
    @test Suslin.verify_ecp_column_reduction(cert)
    @test stage.variable_order == ring_gens
    @test Tuple(Symbol(string(gen)) for gen in stage.variable_order) == entry.variable_order
    @test stage.source_variable_index == stage.variable_index
    @test stage.source_variable == ring_gens[stage.source_variable_index]
    @test stage.target_variable_index == stage.last_variable_index
    @test stage.target_variable == ring_gens[end]
    @test stage.shift_polynomial == stage.shift_sign * stage.target_variable^stage.shift_power
    @test _vc_map_variables(stage.forward_substitution) == ring_gens
    @test _vc_map_variables(stage.inverse_substitution) == ring_gens
    @test _vc_map_values(stage.forward_substitution) == stage.forward_values
    @test _vc_map_values(stage.inverse_substitution) == stage.inverse_values

    recomputed_transformed = tuple((
        R(evaluate(entry, collect(stage.forward_values))) for entry in column
    )...)
    @test stage.transformed_column == recomputed_transformed
    @test 1 <= stage.selected_monic_index <= length(stage.transformed_column)
    @test stage.selected_monic_entry == stage.transformed_column[stage.selected_monic_index]
    @test Suslin._is_monic_in_last_variable(stage.selected_monic_entry, R)

    expected_strategy = stage.selected_monic_index == 1 ? :already_first : :not_moved
    @test stage.first_coordinate_strategy == expected_strategy
    @test isempty(stage.first_coordinate_move_factors)
    @test stage.first_coordinate_column == stage.transformed_column

    @test _vc_apply_factors(stage.transformed_factors, collect(stage.transformed_column), R) == target
    @test all(factor -> base_ring(factor) == R, stage.inverse_substituted_factors)
    @test _vc_apply_factors(stage.inverse_substituted_factors, column, R) == target
    @test stage.variable_change_verification.selected_monic_ok == true
    @test stage.variable_change_verification.transformed_reduction_ok == true
    @test stage.variable_change_verification.original_reduction_ok == true
    @test !Suslin.verify_ecp_column_reduction(_tamper_inverse_substitution_map(cert))
    @test !Suslin.verify_ecp_column_reduction(_tamper_selected_monic_index(cert))
end

@testset "ECP variable-change replay records" begin
    cases = ECPColumnFixtureCatalog.cases_by_id()
    @test Set(VARIABLE_CHANGE_REPLAY_CASE_IDS) ⊆ Set(keys(cases))
    for id in VARIABLE_CHANGE_REPLAY_CASE_IDS
        _assert_variable_change_stage(cases[id])
    end
end
```

- [ ] **Step 4: Register the expert test**

In `test/runtests.jl`, add `"expert/ecp_variable_change_replay.jl"` immediately after `"expert/ecp_column_certificate.jl"`.

- [ ] **Step 5: Run the focused test to verify RED**

Run:

```bash
julia --project=. -e 'include("test/expert/ecp_variable_change_replay.jl")'
```

Expected: FAIL because the current stage does not yet have `variable_order` or another new variable-change field.

- [ ] **Step 6: Commit the red coverage**

```bash
git add test/fixtures/ecp_column_cases.jl test/internal/ecp_column_fixtures.jl test/expert/ecp_variable_change_replay.jl test/runtests.jl docs/superpowers/plans/2026-06-22-issue-85-ecp-variable-change-replay.md
git commit -m "test: cover ECP variable-change replay records"
```

## Task 2: Implement Replay-Checked Variable-Change Stage Fields

**Files:**
- Modify: `src/algorithm/column_reduction.jl`
- Modify: `test/expert/ecp_column_certificate.jl` if a helper must refer to renamed fields

**Interfaces:**
- Consumes: existing `_reduce_after_monicity_normalization_certificate`, `_ecp_replay_stage`, `_substitute_matrix_entries`, `_is_monic_in_last_variable`, `_target_reduced_column`, and `_apply_reduction_factors`.
- Produces: enriched `:monicity_normalization` stages whose new fields exactly match Task 1 tests.

- [ ] **Step 1: Add private helper functions**

Near `_ecp_column_tuple`, add:

```julia
function _ecp_substitution_map_tuple(variables, values)
    return tuple(((; variable = variables[idx], value = values[idx]) for idx in eachindex(variables))...)
end

function _ecp_first_monic_entry_index(column, R)
    return findfirst(entry -> _is_monic_in_last_variable(entry, R), column)
end

function _ecp_first_coordinate_strategy(selected_monic_index::Int)
    return selected_monic_index == 1 ? :already_first : :not_moved
end
```

- [ ] **Step 2: Enrich stage construction**

In `_reduce_after_monicity_normalization_certificate`, replace the `any(...) || continue` check with:

```julia
        selected_monic_index = _ecp_first_monic_entry_index(transformed, R)
        selected_monic_index === nothing && continue
```

After computing `factors`, compute:

```julia
        target = _target_reduced_column(R, length(column))
        transformed_output = _apply_reduction_factors(transformed_result.factors, transformed, R)
        original_output = _apply_reduction_factors(factors, column, R)
        forward_substitution = _ecp_substitution_map_tuple(ring_gens, forward_values)
        inverse_substitution = _ecp_substitution_map_tuple(ring_gens, inverse_values)
        selected_monic_entry = transformed[selected_monic_index]
        first_coordinate_move_factors = typeof(identity_matrix(R, length(column)))[]
        first_coordinate_column = tuple(transformed...)
        variable_change_verification = (;
            selected_monic_ok = _is_monic_in_last_variable(selected_monic_entry, R),
            transformed_reduction_ok = transformed_output == target,
            original_reduction_ok = original_output == target,
        )
```

Add these fields to the stage record in this order after `input_column` and before `transformed_stage`, while preserving existing `variable_index`, `last_variable_index`, `forward_values`, `inverse_values`, `transformed_factors`, `inverse_substituted_factors`, `factors`, and `output_column`:

```julia
            variable_order = tuple(ring_gens...),
            variable_index = var_idx,
            source_variable_index = var_idx,
            source_variable = ring_gens[var_idx],
            last_variable_index = last_var_idx,
            target_variable_index = last_var_idx,
            target_variable = last_var,
            shift_power,
            shift_sign,
            shift_polynomial = shift_sign * last_var^shift_power,
            forward_values = tuple(forward_values...),
            inverse_values = tuple(inverse_values...),
            forward_substitution,
            inverse_substitution,
            transformed_column = tuple(transformed...),
            selected_monic_index,
            selected_monic_entry,
            first_coordinate_strategy = _ecp_first_coordinate_strategy(selected_monic_index),
            first_coordinate_move_factors,
            first_coordinate_column,
```

Keep `output_column = original_output` and add `variable_change_verification` at the end of the stage record.

- [ ] **Step 3: Enrich replay checks**

In `_ecp_replay_stage` for `stage.kind == :monicity_normalization`, recompute:

```julia
        forward_substitution = _ecp_substitution_map_tuple(ring_gens, forward_values)
        inverse_substitution = _ecp_substitution_map_tuple(ring_gens, inverse_values)
        selected_index_ok = stage.selected_monic_index isa Integer &&
            1 <= stage.selected_monic_index <= length(substituted_column)
        selected_monic_entry = selected_index_ok ? substituted_column[stage.selected_monic_index] : zero(R)
        selected_monic_ok = selected_index_ok &&
            stage.selected_monic_entry == selected_monic_entry &&
            _is_monic_in_last_variable(selected_monic_entry, R)
        first_coordinate_strategy = selected_index_ok ?
            _ecp_first_coordinate_strategy(stage.selected_monic_index) :
            :invalid
        transformed_output = _apply_reduction_factors(transformed_replay.factors, substituted_column, R)
        target = _target_reduced_column(R, length(input_column))
        variable_change_verification = (;
            selected_monic_ok,
            transformed_reduction_ok = transformed_output == target,
            original_reduction_ok = expected_output == target,
        )
```

Replace the expected key tuple for `:monicity_normalization` with:

```julia
(:kind, :input_column, :variable_order, :variable_index, :source_variable_index, :source_variable, :last_variable_index, :target_variable_index, :target_variable, :shift_power, :shift_sign, :shift_polynomial, :forward_values, :inverse_values, :forward_substitution, :inverse_substitution, :transformed_column, :selected_monic_index, :selected_monic_entry, :first_coordinate_strategy, :first_coordinate_move_factors, :first_coordinate_column, :transformed_stage, :transformed_factors, :inverse_substituted_factors, :factors, :output_column, :variable_change_verification)
```

Replace `stage.substituted_column == tuple(substituted_column...)` with `stage.transformed_column == tuple(substituted_column...)`, and add replay checks:

```julia
            stage.variable_order == tuple(ring_gens...) &&
            stage.source_variable_index == stage.variable_index &&
            stage.source_variable == ring_gens[stage.source_variable_index] &&
            stage.target_variable_index == stage.last_variable_index &&
            stage.target_variable == last_var &&
            stage.shift_polynomial == stage.shift_sign * last_var^stage.shift_power &&
            stage.forward_substitution == forward_substitution &&
            stage.inverse_substitution == inverse_substitution &&
            selected_monic_ok &&
            stage.first_coordinate_strategy == first_coordinate_strategy &&
            isempty(stage.first_coordinate_move_factors) &&
            stage.first_coordinate_column == tuple(substituted_column...) &&
            stage.variable_change_verification == variable_change_verification &&
            all(values(variable_change_verification)) &&
```

- [ ] **Step 4: Run the focused test to verify GREEN**

Run:

```bash
julia --project=. -e 'include("test/expert/ecp_variable_change_replay.jl")'
```

Expected: PASS.

- [ ] **Step 5: Run the existing certificate regression**

Run:

```bash
julia --project=. -e 'include("test/expert/ecp_column_certificate.jl")'
```

Expected: PASS.

- [ ] **Step 6: Commit implementation**

```bash
git add src/algorithm/column_reduction.jl test/expert/ecp_column_certificate.jl
git commit -m "feat: record ECP variable-change replay stages"
```

## Task 3: Expert Registration and Regression Sweep

**Files:**
- Modify only if needed by failures: `test/expert/ecp_variable_change_replay.jl`, `test/runtests.jl`, `src/algorithm/column_reduction.jl`

**Interfaces:**
- Consumes: Task 1 tests and Task 2 implementation.
- Produces: clean focused, expert, and package verification results.

- [ ] **Step 1: Run focused issue verification**

Run:

```bash
julia --project=. -e 'include("test/expert/ecp_variable_change_replay.jl")'
```

Expected: PASS.

- [ ] **Step 2: Run package verification**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: PASS with public and internal groups.

- [ ] **Step 3: Run expert regression**

Run:

```bash
julia --project=. test/runtests.jl expert
```

Expected: PASS.

- [ ] **Step 4: Run whitespace check**

Run:

```bash
git diff --check
```

Expected: no output and exit 0.

- [ ] **Step 5: Commit any regression fixes**

If Step 1-4 required changes, commit them:

```bash
git add src/algorithm/column_reduction.jl test/fixtures/ecp_column_cases.jl test/internal/ecp_column_fixtures.jl test/expert/ecp_variable_change_replay.jl test/expert/ecp_column_certificate.jl test/runtests.jl
git commit -m "test: register ECP variable-change replay regression"
```

If no files changed, do not create an empty commit.

## Self-Review

- Spec coverage: Task 1 covers the required two catalog cases and negative controls; Task 2 records and replays the requested stage data; Task 3 runs the issue, package, and expert verification commands.
- Placeholder scan: no `TBD`, `TODO`, or unspecified implementation steps remain.
- Type consistency: all tasks use the same field names for the enriched `:monicity_normalization` stage.

## Execution Choice

Plan complete and saved to `docs/superpowers/plans/2026-06-22-issue-85-ecp-variable-change-replay.md`. Two execution options:

1. Subagent-Driven (recommended) - dispatch a fresh subagent per task, review between tasks, fast iteration
2. Inline Execution - execute tasks in this session using executing-plans, batch execution with checkpoints

Use option 1 under the Standing Answer Policy because it is marked recommended.
