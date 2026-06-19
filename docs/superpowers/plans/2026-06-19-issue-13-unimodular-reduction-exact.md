# Issue 13 Exact Unimodular Column Reduction Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace `reduce_unimodular_column`'s heuristic top-level behavior with a staged exact pipeline that reduces supported ordinary and Laurent-normalized unimodular columns to `e_n`.

**Architecture:** Keep the public function in `src/algorithm/column_reduction.jl`, split validation, ordinary reduction, 3-entry block reduction, Laurent normalization/lift-back, exact verification, and staged unsupported failure into private helpers. Add one focused expert test file and register it in the expert test group.

**Tech Stack:** Julia, Oscar exact polynomial and Laurent rings, existing Suslin elementary matrix helpers, Laurent normalization helpers from issue #8, `Test`.

## Global Constraints

- Work in the existing linked worktree on branch `agent/issue-13-replace-heuristic-unimodular-column-reduction-wi-run-1`; do not create another worktree.
- Follow TDD: write the failing issue test first, run it and confirm the expected failure, then write production code.
- Keep the existing public function name `reduce_unimodular_column(v, R)`.
- Do not add new public exports.
- Non-unimodular columns must fail immediately with `ArgumentError("v must be a unimodular column")`.
- Unsupported but unimodular columns must fail with a staged unsupported message that does not claim the column is non-unimodular.
- Preserve the old successful small examples from `test/expert/unimodular_columns.jl`.
- Laurent support must use the existing `normalize_laurent_object` metadata and lift factors back; do not claim general Laurent factorization.
- The issue-specific verification command is `julia --project=. -e 'include("test/expert/unimodular_reduction_exact.jl")'`.
- The Agent Desk package verification command is `julia --project=. -e 'using Pkg; Pkg.test()'`.
- The documented full suite from #21 is `julia --project=. test/runtests.jl all`.
- Do not commit `Manifest.toml`.

---

## File Structure

- Modify `src/algorithm/column_reduction.jl`: refactor the public reducer and add private exact pipeline helpers.
- Create `test/expert/unimodular_reduction_exact.jl`: focused issue #13 tests.
- Modify `test/runtests.jl`: register the focused expert test.

### Task 1: Add Exact Reduction Expert Tests

**Files:**
- Create: `test/expert/unimodular_reduction_exact.jl`
- Modify: `test/runtests.jl`

**Interfaces:**
- Consumes: existing `Suslin.reduce_unimodular_column(v, R)` and `Suslin.is_unimodular_column(v, R)`.
- Produces: focused expert coverage that fails before the exact block and Laurent-normalized stages are implemented.

- [ ] **Step 1: Write the failing expert test file**

Create `test/expert/unimodular_reduction_exact.jl` with this content:

```julia
using Test
using Suslin
using Oscar

function exact_reduction_product(factors, R, n::Int)
    product = identity_matrix(R, n)
    for factor in factors
        product *= factor
    end
    return product
end

function exact_apply_factors(factors, v, R)
    n = length(v)
    column = matrix(R, n, 1, collect(v))
    return exact_reduction_product(factors, R, n) * column
end

function exact_target_column(R, n::Int)
    target = zero_matrix(R, n, 1)
    target[n, 1] = one(R)
    return target
end

function assert_reduces_to_last_unit(v, R)
    factors = Suslin.reduce_unimodular_column(v, R)
    @test exact_apply_factors(factors, v, R) == exact_target_column(R, length(v))
    return factors
end

function captured_reduction_error(v, R)
    try
        Suslin.reduce_unimodular_column(v, R)
        return nothing
    catch err
        return err
    end
end

@testset "exact unimodular reduction supports longer ordinary columns" begin
    F2 = GF(2)
    R, (x, y) = Oscar.polynomial_ring(F2, ["x", "y"])

    hard_slice = [
        x + y^2,
        x * y + x + one(R),
        x^2 + x * y + y + one(R),
    ]

    length6 = vcat(hard_slice, [x^2, x * y, y^2 + x])
    @test Suslin.is_unimodular_column(length6, R)
    assert_reduces_to_last_unit(length6, R)

    length8 = [
        x^2 + y,
        hard_slice[1],
        x * y,
        y^2 + one(R),
        hard_slice[2],
        x^3 + y,
        x * y + y,
        hard_slice[3],
    ]
    @test Suslin.is_unimodular_column(length8, R)
    assert_reduces_to_last_unit(length8, R)

    witness_slice = [x, y, x + one(R)]
    length12 = [
        x^2,
        y^2,
        x * y,
        witness_slice[1],
        x^2 + y,
        y^3 + x,
        witness_slice[2],
        x * y + y,
        x^3 + y^2,
        y^2 + y,
        x^2 * y + x,
        witness_slice[3],
    ]
    @test Suslin.is_unimodular_column(length12, R)
    assert_reduces_to_last_unit(length12, R)
end

@testset "exact unimodular reduction supports Laurent-normalized columns" begin
    R, (x, y) = Suslin.suslin_laurent_polynomial_ring(GF(2), ["x", "y"])

    v = [
        x^-1,
        x^-2 * y,
        x^-1 + x^-2,
        x^-1 * y,
        x^-2 * y^2 + x^-1,
        x^-2,
    ]

    @test Suslin.is_unimodular_column(v, R)
    assert_reduces_to_last_unit(v, R)
end

@testset "exact unimodular reduction preserves old small cases" begin
    F2 = GF(2)
    R2, (x, y) = Oscar.polynomial_ring(F2, ["x", "y"])

    assert_reduces_to_last_unit([x, y, one(R2)], R2)
    assert_reduces_to_last_unit([x, one(R2), y], R2)
    assert_reduces_to_last_unit([x, y, x + one(R2)], R2)
    assert_reduces_to_last_unit([
        x + y^2,
        x * y + x + one(R2),
        x^2 + x * y + y + one(R2),
    ], R2)

    F3 = GF(3)
    R3, (t,) = Oscar.polynomial_ring(F3, ["t"])
    assert_reduces_to_last_unit([t + one(R3), t, R3(2)], R3)
    assert_reduces_to_last_unit([t + one(R3), R3(2), t], R3)
end

@testset "exact unimodular reduction staged failures" begin
    F2 = GF(2)
    R, (x, y) = Oscar.polynomial_ring(F2, ["x", "y"])

    non_unimodular = [x, y, x * y]
    non_unimodular_err = captured_reduction_error(non_unimodular, R)
    @test non_unimodular_err isa ArgumentError
    @test occursin("v must be a unimodular column", sprint(showerror, non_unimodular_err))

    unsupported = [zero(R), x^2, x * y + one(R)]
    @test Suslin.is_unimodular_column(unsupported, R)
    unsupported_err = captured_reduction_error(unsupported, R)
    @test unsupported_err isa ArgumentError
    @test occursin("unsupported exact unimodular column reduction", sprint(showerror, unsupported_err))
    @test !occursin("not unimodular", sprint(showerror, unsupported_err))
end
```

- [ ] **Step 2: Register the expert test**

In `test/runtests.jl`, add the new expert file after `expert/unimodular_columns.jl`:

```julia
"expert/unimodular_reduction_exact.jl",
```

- [ ] **Step 3: Run the focused test to verify RED**

Run:

```bash
julia --project=. -e 'include("test/expert/unimodular_reduction_exact.jl")'
```

Expected before Task 2: failure in the longer ordinary or Laurent-normalized
test because the current reducer has no exact block or Laurent-normalization
pipeline. If it fails earlier due to a test typo, fix the test and re-run until
the failure demonstrates missing production behavior.

- [ ] **Step 4: Commit the red test**

Commit only the new test and test registration:

```bash
git add test/expert/unimodular_reduction_exact.jl test/runtests.jl
git commit -m "test: cover exact unimodular reduction pipeline"
```

### Task 2: Implement the Staged Exact Reduction Pipeline

**Files:**
- Modify: `src/algorithm/column_reduction.jl`

**Interfaces:**
- Consumes: `elementary_matrix`, `block_embedding`, `normalize_laurent_object`, `lift_laurent_normalization` metadata shape, `_is_laurent_polynomial_ring`, `_coerce_into_ring`, and existing unit/witness/monicity helpers.
- Produces: `reduce_unimodular_column(v, R)` with immediate validation, exact ordinary pipeline, exact 3-entry block pipeline, Laurent normalization/lift-back pipeline, and staged unsupported failures.

- [ ] **Step 1: Refactor public validation and dispatch**

Replace the top of `src/algorithm/column_reduction.jl` through the old
`reduce_unimodular_column` body with helper-based dispatch shaped like this:

```julia
function reduce_unimodular_column(v::AbstractVector, R)
    column = _validated_unimodular_column(v, R)

    factors = _is_laurent_polynomial_ring(R) ?
        _reduce_laurent_unimodular_column(column, R) :
        _reduce_polynomial_unimodular_column_exact(column, R)

    factors !== nothing && return _checked_reduction_factors(factors, column, R, "public reducer")
    _throw_unsupported_unimodular_column_reduction(column, R)
end

function _validated_unimodular_column(v::AbstractVector, R)
    Base.require_one_based_indexing(v)
    n = length(v)
    n >= 3 || throw(ArgumentError("v must have length at least 3"))

    column = [_coerce_into_ring(R, v[idx], "v[$idx]") for idx in 1:n]
    is_unimodular_column(column, R) || throw(ArgumentError("v must be a unimodular column"))
    return column
end
```

- [ ] **Step 2: Add exact product verification helpers**

Add these private helpers near the public reducer:

```julia
function _factor_sequence_product(factors, R, n::Int)
    product = identity_matrix(R, n)
    for factor in factors
        nrows(factor) == n || throw(ArgumentError("reduction factor has wrong row count"))
        ncols(factor) == n || throw(ArgumentError("reduction factor has wrong column count"))
        _same_base_ring(base_ring(factor), R) || throw(ArgumentError("reduction factor has wrong base ring"))
        product *= factor
    end
    return product
end

function _apply_reduction_factors(factors, column::AbstractVector, R)
    n = length(column)
    return _factor_sequence_product(factors, R, n) * matrix(R, n, 1, collect(column))
end

function _target_reduced_column(R, n::Int)
    target = zero_matrix(R, n, 1)
    target[n, 1] = one(R)
    return target
end

function _checked_reduction_factors(factors, column::AbstractVector, R, stage::AbstractString)
    n = length(column)
    _apply_reduction_factors(factors, column, R) == _target_reduced_column(R, n) ||
        throw(ErrorException("internal error: $(stage) produced factors that do not reduce the column to e_n"))
    return factors
end
```

- [ ] **Step 3: Add ordinary exact stages**

Add the ordinary pipeline and small-stage helper:

```julia
function _reduce_polynomial_unimodular_column_exact(column::AbstractVector, R)
    factors = _reduce_exact_small_column(column, R)
    factors !== nothing && return factors

    block_factors = _reduce_via_supported_three_block(column, R)
    block_factors !== nothing && return block_factors

    return nothing
end

function _has_at_least_two_generators(R)::Bool
    try
        return ngens(R) >= 2
    catch err
        err isa MethodError || err isa ErrorException || rethrow()
        return false
    end
end

function _reduce_exact_small_column(column::AbstractVector, R)
    factors = _reduce_supported_unimodular_column(column, R)
    factors !== nothing && return _checked_reduction_factors(factors, column, R, "unit or witness reduction")

    _has_at_least_two_generators(R) || return nothing
    normalization_factors = _reduce_after_monicity_normalization(column, R)
    normalization_factors !== nothing &&
        return _checked_reduction_factors(normalization_factors, column, R, "monicity normalization reduction")

    return nothing
end
```

- [ ] **Step 4: Add the supported 3-entry block stage**

Add these helpers below the ordinary pipeline:

```julia
function _reduce_via_supported_three_block(column::AbstractVector, R)
    n = length(column)
    n > 3 || return nothing

    for i in 1:(n - 2), j in (i + 1):(n - 1), k in (j + 1):n
        indices = (i, j, k)
        subcolumn = [column[idx] for idx in indices]
        is_unimodular_column(subcolumn, R) || continue

        subfactors = _reduce_exact_small_column(subcolumn, R)
        subfactors === nothing && continue

        return _embedded_three_block_reduction(column, R, indices, subfactors)
    end

    return nothing
end

function _embedded_three_block_reduction(column::AbstractVector, R, indices, subfactors)
    n = length(column)
    pivot_idx = indices[end]
    embedded_factors = [block_embedding(factor, n, indices) for factor in subfactors]
    after_block = _apply_reduction_factors(embedded_factors, column, R)

    elimination_factors = typeof(identity_matrix(R, n))[]
    for row in 1:n
        row == pivot_idx && continue
        coeff = -after_block[row, 1]
        coeff == zero(R) && continue
        push!(elimination_factors, elementary_matrix(n, row, pivot_idx, coeff, R))
    end

    return _checked_reduction_factors(
        vcat(elimination_factors, embedded_factors),
        column,
        R,
        "embedded 3-entry reduction",
    )
end
```

- [ ] **Step 5: Add Laurent normalization and lift-back**

Add the Laurent-specific helpers:

```julia
function _reduce_laurent_unimodular_column(column::AbstractVector, R)
    unit_idx = findfirst(is_unit, column)
    unit_idx !== nothing && return _reduce_unit_witness_column(column, unit_idx, R)

    normalization = normalize_laurent_object(column)
    poly_column = normalization.normalized_object
    P = normalization.metadata.polynomial_ring
    is_unimodular_column(poly_column, P) || return nothing

    poly_factors = _reduce_polynomial_unimodular_column_exact(poly_column, P)
    poly_factors === nothing && return nothing

    lifted_factors = [_lift_polynomial_reduction_factor(factor, R) for factor in poly_factors]
    shift = only(normalization.metadata.shift_monomials)
    inverse_shift = only(normalization.metadata.inverse_shift_monomials)
    normalization_factors = _unit_normalization_factors(length(column), inverse_shift, shift, R)

    return _checked_reduction_factors(
        vcat(normalization_factors, lifted_factors),
        column,
        R,
        "Laurent normalization reduction",
    )
end

function _lift_polynomial_reduction_factor(factor, R)
    entries = [
        _coerce_into_ring(R, factor[row, col], "lifted reduction factor entry")
        for col in 1:ncols(factor), row in 1:nrows(factor)
    ]
    return matrix(R, nrows(factor), ncols(factor), vec(entries))
end
```

- [ ] **Step 6: Add staged unsupported failure**

Replace the old final unsupported error with:

```julia
function _throw_unsupported_unimodular_column_reduction(column::AbstractVector, R)
    n = length(column)
    profile = _is_laurent_polynomial_ring(R) ? "Laurent-normalized" : "ordinary polynomial"
    throw(ArgumentError("unsupported exact unimodular column reduction for $(profile) column of length $(n): no supported unit, witness-unit, monicity-normalized, or 3-entry block reduction stage applies"))
end
```

- [ ] **Step 7: Run the focused issue test to verify GREEN**

Run:

```bash
julia --project=. -e 'include("test/expert/unimodular_reduction_exact.jl")'
```

Expected: pass.

- [ ] **Step 8: Run old focused expert test**

Run:

```bash
julia --project=. -e 'include("test/expert/unimodular_columns.jl")'
```

Expected: pass.

- [ ] **Step 9: Commit the implementation**

Commit only the implementation file:

```bash
git add src/algorithm/column_reduction.jl
git commit -m "feat: add exact unimodular column reduction pipeline"
```

## Final Verification

After both tasks are complete, run:

```bash
julia --project=. -e 'include("test/expert/unimodular_reduction_exact.jl")'
julia --project=. -e 'using Pkg; Pkg.test()'
julia --project=. test/runtests.jl all
```

Then use `superpowers:requesting-code-review` through the subagent-driven final
review step, address Critical and Important findings, and use
`superpowers:finishing-a-development-branch` with option 2, "Push and create a
Pull Request".

## Self-Review

- The plan covers every requirement in the design spec and issue body.
- There are no placeholder steps.
- The test task fails before production changes and the implementation task has
  focused green checks plus full final verification.
