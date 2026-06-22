# Issue 86 ECP Monicity Search Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a deterministic, variable-order-driven ECP monicity search helper that returns replayable success records or staged search-exhaustion failures.

**Architecture:** The implementation keeps `reduce_unimodular_column(v, R)` factor-returning and routes its monicity normalization stage through a new internal search helper. Successful searches reuse the #85 `:monicity_normalization` stage shape; exhausted searches return a separate failure struct that records the bounded search space without claiming non-unimodularity.

**Tech Stack:** Julia, Oscar/AbstractAlgebra polynomial rings, existing `src/algorithm/column_reduction.jl` internals, fixture-backed expert tests.

## Global Constraints

- Repository has no `AGENTS.md`.
- Base branch is `main`; worker branch is `agent/issue-86-generalize-the-deterministic-ecp-monicity-search-run-1`.
- Preserve `reduce_unimodular_column(v, R)` public behavior: it returns only factors.
- Do not export the new helper or record types from `src/Suslin.jl`.
- Search must be deterministic over the supplied variable order, bounded shift powers, and supplied shift signs.
- Successful search records must reuse the #85 replayable monicity stage fields.
- Staged failure must name the searched variable order and bound and must not pretend the column is non-unimodular.
- Do not use randomness.
- Do not implement Quillen/local-to-global patching.
- Do not claim support for arbitrary Noether-normalization choices.

---

## File Structure

- Modify `src/algorithm/column_reduction.jl`: add `ECPMonicitySearchResult`, `ECPMonicitySearchFailure`, variable-order normalization helpers, deterministic search helper, and replay adjustments for non-native variable order.
- Create `test/expert/ecp_monicity_search.jl`: fixture-backed expert coverage for success, non-native variable order, staged exhaustion, and negative controls.

### Task 1: Deterministic ECP Monicity Search

**Files:**
- Modify: `src/algorithm/column_reduction.jl`
- Create: `test/expert/ecp_monicity_search.jl`

**Interfaces:**
- Consumes: `_reduce_supported_unimodular_column_certificate(column, R)`, `_is_monic_in_last_variable(entry, R)`, `_substitute_matrix_entries(M, values, R)`, `_ecp_substitution_map_tuple(variables, values)`, `_checked_reduction_factors(factors, column, R, stage)`.
- Produces: `ECPMonicitySearchResult`, `ECPMonicitySearchFailure`, `_deterministic_ecp_monicity_search(column, R; variable_order, max_shift_power, shift_signs)`, and a reducer path where `_reduce_after_monicity_normalization_certificate(column, R)` returns `result.stage` on success and `nothing` on staged failure.

- [ ] **Step 1: Write the failing expert test**

Create `test/expert/ecp_monicity_search.jl` with focused tests that call the wished-for helper before it exists.

```julia
using Test
using Suslin
using Oscar

include(joinpath(@__DIR__, "..", "fixtures", "ecp_column_cases.jl"))

function _search_column(entry)
    return [getproperty(entry.entries, name) for name in entry.column_order]
end

function _search_target_column(R, n::Int)
    target = zero_matrix(R, n, 1)
    target[n, 1] = one(R)
    return target
end

function _search_factor_product(factors, R, n::Int)
    product = identity_matrix(R, n)
    for factor in factors
        product *= factor
    end
    return product
end

function _search_apply_factors(factors, column, R)
    return _search_factor_product(factors, R, length(column)) * matrix(R, length(column), 1, collect(column))
end

function _search_stage(result)
    @test result isa Suslin.ECPMonicitySearchResult
    @test result.stage.kind == :monicity_normalization
    return result.stage
end

function _assert_success_reduces(entry; variable_order = tuple(gens(entry.ring.object)...), max_shift_power = 3)
    column = _search_column(entry)
    R = entry.ring.object
    result = Suslin._deterministic_ecp_monicity_search(
        column,
        R;
        variable_order,
        max_shift_power,
    )
    stage = _search_stage(result)
    @test result.original_column == tuple(column...)
    @test result.ring == R
    @test result.variable_order == tuple(Suslin._ecp_normalize_variable_order(R, variable_order)...)
    @test result.max_shift_power == max_shift_power
    @test result.stage === stage
    @test stage.variable_order == result.variable_order
    @test stage.source_variable == result.source_variable
    @test stage.target_variable == result.target_variable
    @test stage.shift_power == result.shift_power
    @test stage.shift_sign == result.shift_sign
    @test stage.shift_polynomial == result.shift_polynomial
    @test stage.selected_monic_index == result.selected_monic_index
    @test stage.selected_monic_entry == result.selected_monic_entry
    @test Suslin._is_monic_in_last_variable(stage.selected_monic_entry, R)
    @test _search_apply_factors(result.factors, column, R) == _search_target_column(R, length(column))
    @test _search_apply_factors(stage.factors, column, R) == _search_target_column(R, length(column))
    cert = Suslin.ecp_column_reduction_certificate(column, R)
    @test Suslin.verify_ecp_column_reduction(cert)
    @test any(cert_stage -> cert_stage.kind == :monicity_normalization, cert.stages)
    return result
end

function _target_x_fixture()
    R, (x, y) = Oscar.polynomial_ring(GF(2), ["x", "y"])
    return (;
        id = "ecp-variable-change-target-x-gf2",
        ring = (; object = R),
        variable_order = (:y, :x),
        entries = (
            a = y + x^2,
            b = x * y + y + one(R),
            c = y^2 + x * y + x + one(R),
        ),
        column_order = (:a, :b, :c),
    )
end

@testset "deterministic ECP monicity search" begin
    cases = ECPColumnFixtureCatalog.cases_by_id()

    old_bounded = _assert_success_reduces(cases["ecp-variable-change-monic-gf2"])
    @test old_bounded.source_variable == gens(old_bounded.ring)[1]
    @test old_bounded.target_variable == gens(old_bounded.ring)[2]
    @test old_bounded.shift_power == 2

    target_x_fixture = _target_x_fixture()
    target_x = _assert_success_reduces(
        target_x_fixture;
        variable_order = reverse(tuple(gens(target_x_fixture.ring.object)...)),
    )
    @test target_x.variable_order == reverse(tuple(gens(target_x.ring)...))
    @test target_x.source_variable == gens(target_x.ring)[2]
    @test target_x.target_variable == gens(target_x.ring)[1]
    @test target_x.shift_power == 2
    @test target_x.stage.target_variable_index == 1
    @test target_x.stage.source_variable_index == 2

    exhausted_entry = cases["ecp-unsupported-unimodular-gf2"]
    exhausted_column = _search_column(exhausted_entry)
    exhausted = Suslin._deterministic_ecp_monicity_search(
        exhausted_column,
        exhausted_entry.ring.object;
        variable_order = tuple(gens(exhausted_entry.ring.object)...),
        max_shift_power = 3,
    )
    @test exhausted isa Suslin.ECPMonicitySearchFailure
    @test exhausted.kind == :monicity_search_exhausted
    @test exhausted.variable_order == tuple(gens(exhausted_entry.ring.object)...)
    @test exhausted.max_shift_power == 3
    @test exhausted.attempted_candidates == 3
    @test occursin("exhausted deterministic ECP monicity search", exhausted.message)

    low_bound = Suslin._deterministic_ecp_monicity_search(
        _search_column(cases["ecp-variable-change-monic-gf2"]),
        cases["ecp-variable-change-monic-gf2"].ring.object;
        variable_order = tuple(gens(cases["ecp-variable-change-monic-gf2"].ring.object)...),
        max_shift_power = 1,
    )
    @test low_bound isa Suslin.ECPMonicitySearchFailure
    @test low_bound.max_shift_power == 1

    missing_required_source = Suslin._deterministic_ecp_monicity_search(
        _search_column(cases["ecp-variable-change-monic-gf2"]),
        cases["ecp-variable-change-monic-gf2"].ring.object;
        variable_order = (last(gens(cases["ecp-variable-change-monic-gf2"].ring.object)),),
        max_shift_power = 3,
    )
    @test missing_required_source isa Suslin.ECPMonicitySearchFailure
    @test missing_required_source.variable_order == (last(gens(cases["ecp-variable-change-monic-gf2"].ring.object)),)
    @test missing_required_source.attempted_candidates == 0
end
```

- [ ] **Step 2: Run the red test**

Run:

```bash
julia --project=. -e 'include("test/expert/ecp_monicity_search.jl")'
```

Expected: FAIL because `Suslin.ECPMonicitySearchResult` or `Suslin._deterministic_ecp_monicity_search` is not defined.

- [ ] **Step 3: Add the minimal search records and helper**

In `src/algorithm/column_reduction.jl`, add `ECPMonicitySearchResult` and `ECPMonicitySearchFailure` near `ECPColumnReductionCertificate`. Refactor `_reduce_after_monicity_normalization_certificate` so it calls `_deterministic_ecp_monicity_search(column, R)` and returns `result.stage`/`result.factors` only when the result is a success.

The helper must:

```julia
function _deterministic_ecp_monicity_search(
    column::AbstractVector,
    R;
    variable_order = tuple(gens(R)...),
    max_shift_power::Integer = 3,
    shift_signs = (one(R), -one(R)),
)
    normalized_order = _ecp_normalize_variable_order(R, variable_order)
    max_shift_power < 0 && throw(ArgumentError("max_shift_power must be nonnegative"))
    signs = tuple((_coerce_into_ring(R, sign, "shift sign") for sign in shift_signs)...)
    isempty(normalized_order) && return _ecp_monicity_search_failure(column, R, normalized_order, max_shift_power, signs, 0)
    length(normalized_order) < 2 && return _ecp_monicity_search_failure(column, R, normalized_order, max_shift_power, signs, 0)

    ring_gens = collect(gens(R))
    target_variable = normalized_order[end]
    target_variable_index = _ecp_generator_index(R, target_variable)
    attempted = 0
    for source_variable in normalized_order[1:(end - 1)]
        source_variable_index = _ecp_generator_index(R, source_variable)
        for shift_power in 1:max_shift_power, shift_sign in signs
            attempted += 1
            candidate = _ecp_monicity_candidate_stage(
                column,
                R,
                tuple(normalized_order...),
                source_variable_index,
                target_variable_index,
                shift_power,
                shift_sign,
            )
            candidate === nothing && continue
            return candidate
        end
    end

    return _ecp_monicity_search_failure(column, R, normalized_order, max_shift_power, signs, attempted)
end
```

Use helper functions for `_ecp_normalize_variable_order`, `_ecp_generator_index`, `_ecp_monicity_candidate_stage`, and `_ecp_monicity_search_failure` so the replay branch can use the same metadata. The candidate stage should build the same #85 fields as the old loop, except `variable_order`, `source_variable_index`, and `target_variable_index` come from the requested order rather than assuming `gens(R)[end]` is the target.

- [ ] **Step 4: Update monicity replay for requested variable order**

In `_ecp_replay_stage(stage, input_column, R)` for `stage.kind == :monicity_normalization`, replace the hard-coded `last_var = ring_gens[end]` logic with:

```julia
variable_order = collect(stage.variable_order)
target_variable = stage.target_variable
target_variable_index = stage.target_variable_index
source_variable_index = stage.source_variable_index
source_variable = stage.source_variable
```

Then verify:

- `stage.variable_order` is a duplicate-free subset of `tuple(gens(R)...)`.
- `stage.target_variable == variable_order[end]`.
- `stage.source_variable` appears before the target in `stage.variable_order`.
- `stage.variable_index == stage.source_variable_index`.
- `stage.last_variable_index == stage.target_variable_index`.
- `stage.shift_polynomial == stage.shift_sign * stage.target_variable^stage.shift_power`.
- `forward_values` and `inverse_values` are native-generator-order tuples passed to `evaluate`.
- `selected_monic_entry` is monic in `stage.target_variable`, not always in `gens(R)[end]`.

Add a helper `_is_monic_in_variable(p, R, variable_index::Int)` and make `_is_monic_in_last_variable(p, R)` delegate to it with `ngens(R)` so existing tests keep working.

- [ ] **Step 5: Run focused verification**

Run:

```bash
julia --project=. -e 'include("test/expert/ecp_monicity_search.jl")'
julia --project=. -e 'include("test/expert/ecp_variable_change_replay.jl")'
julia --project=. -e 'include("test/expert/ecp_column_certificate.jl")'
```

Expected: all commands exit 0.

- [ ] **Step 6: Run full verification and commit**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
git diff --check
```

Expected: both commands exit 0. Then commit:

```bash
git add src/algorithm/column_reduction.jl test/expert/ecp_monicity_search.jl docs/superpowers/plans/2026-06-22-issue-86-ecp-monicity-search-plan.md
git commit -m "feat: add deterministic ECP monicity search"
```

## Self-Review

- Spec coverage: the task adds success and failure records, deterministic variable-order search, reducer integration, replay support, and fixture-backed tests.
- Placeholder scan: no `TBD`, `TODO`, or "fill in later" placeholders are present.
- Type consistency: the produced helper names and record names match the tests and reducer integration.
