# Issue 57 Laurent Column-Peel Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement recursive Laurent `SL_n` last-column peeling for the Issue #38 row and column determinant-one cores.

**Architecture:** Add a focused internal reducer in `src/algorithm/laurent_column_peel.jl`, then use it as the Laurent determinant-one fallback in `elementary_factorization` after the existing block-local `reduce_sln_to_sl3` path fails. The reducer stores replay metadata for each peel step and uses the existing local `SL_3` unit-pivot solver for the final `2 x 2` block.

**Tech Stack:** Julia, Oscar Laurent polynomial rings and matrices, existing Suslin elementary matrices, `reduce_unimodular_column`, `realize_sl3_local`, and Test stdlib.

## Global Constraints

- Do not factorize the original Issue #38 `Q` matrix with determinant `u*v`.
- Do not add broad matrix search.
- Always peel the current last column first.
- Use `reduce_unimodular_column` for column reduction.
- Use right elementary factors `E(d, j, -B[d, j])` to clear the bottom row after left reduction.
- Finish the final `2 x 2` block by embedding it as a `3 x 3` local target and calling `realize_sl3_local(...; check_monic=false)`.
- The row core final block must be `[u u*v; 0 u^-1]`.
- The column core final block must be `[v^-1 u*v; 0 v]`.
- Negative controls must make replay verification and `verify_factorization` fail.
- Focused verification command is `julia --project=. -e 'include("test/expert/laurent_column_peel_issue38.jl")'`.
- Package verification command is `julia --project=. -e 'using Pkg; Pkg.test()'`.
- Do not commit `Manifest.toml`.

---

## File Structure

- Create `src/algorithm/laurent_column_peel.jl`: internal metadata structs, recursive peel implementation, replay verifier, and helper functions.
- Modify `src/Suslin.jl`: include the new algorithm file after `algorithm/sl3_local.jl` and before `algorithm/factorization.jl`.
- Modify `src/algorithm/factorization.jl`: call the column-peel fallback for Laurent determinant-one matrices when the block-local reducer fails.
- Create `test/expert/laurent_column_peel_issue38.jl`: focused Issue #57 acceptance and negative controls.
- Modify `test/runtests.jl`: register the expert test.
- Modify `test/fixtures/toricbuilder_issue38_cases.jl`: update normalized-core status metadata to `:supported_column_peel`.
- Modify `test/internal/toricbuilder_issue38_fixture.jl`: validate normalized cores now factorize, while the original `Q` remains unsupported.

---

### Task 1: Add the Failing Issue 57 Expert Test

**Files:**
- Create: `test/expert/laurent_column_peel_issue38.jl`

**Interfaces:**
- Consumes: planned internal `Suslin._factor_laurent_sl_column_peel(A)`, `Suslin._verify_laurent_column_peel_replay(certificate)`, and `elementary_factorization(A)`.
- Produces: focused RED coverage for row and column Issue #38 determinant-one cores.

- [ ] **Step 1: Write the failing test**

Create `test/expert/laurent_column_peel_issue38.jl`:

```julia
using Test
using Suslin
using Oscar

include("../fixtures/toricbuilder_issue38_cases.jl")

function _issue57_product(factors, R, n::Int)
    product = identity_matrix(R, n)
    for factor in factors
        product *= factor
    end
    return product
end

function _issue57_assert_step(step)
    R = base_ring(step.input_matrix)
    left_product = _issue57_product(step.left_factors, R, step.dimension)
    right_product = _issue57_product(step.right_factors, R, step.dimension)
    @test left_product * step.input_matrix * right_product == step.peeled_matrix
    @test step.peeled_matrix[step.dimension, step.dimension] == one(R)
    @test all(step.peeled_matrix[row, step.dimension] == zero(R) for row in 1:(step.dimension - 1))
    @test all(step.peeled_matrix[step.dimension, col] == zero(R) for col in 1:(step.dimension - 1))
    @test step.next_block == matrix(R, [
        step.peeled_matrix[row, col]
        for row in 1:(step.dimension - 1), col in 1:(step.dimension - 1)
    ])
end

function _issue57_corrupt_left_factor(certificate)
    corrupted = collect(certificate.peel_steps)
    first_step = first(corrupted)
    bad_left = copy(first_step.left_factors)
    R = base_ring(first_step.input_matrix)
    bad_left[1] = elementary_matrix(first_step.dimension, 1, 2, one(R), R) * bad_left[1]
    corrupted[1] = Suslin.LaurentColumnPeelStep(
        first_step.dimension,
        first_step.input_matrix,
        first_step.last_column,
        bad_left,
        first_step.after_left_matrix,
        first_step.right_factors,
        first_step.peeled_matrix,
        first_step.next_block,
    )
    return Suslin.LaurentColumnPeelFactorization(
        certificate.original_matrix,
        certificate.final_block,
        certificate.final_local_target,
        certificate.final_local_factors,
        certificate.final_factors,
        certificate.factors,
        certificate.product,
        corrupted,
        certificate.verification,
    )
end

function _issue57_delete_clearing_factor(certificate)
    corrupted = collect(certificate.peel_steps)
    first_step = first(corrupted)
    shortened = first_step.right_factors[1:(end - 1)]
    corrupted[1] = Suslin.LaurentColumnPeelStep(
        first_step.dimension,
        first_step.input_matrix,
        first_step.last_column,
        first_step.left_factors,
        first_step.after_left_matrix,
        shortened,
        first_step.peeled_matrix,
        first_step.next_block,
    )
    return Suslin.LaurentColumnPeelFactorization(
        certificate.original_matrix,
        certificate.final_block,
        certificate.final_local_target,
        certificate.final_local_factors,
        certificate.final_factors,
        certificate.factors[1:(end - 1)],
        certificate.product,
        corrupted,
        certificate.verification,
    )
end

function _issue57_assert_core(core, expected_final_block)
    certificate = Suslin._factor_laurent_sl_column_peel(core)

    @test certificate.final_block == expected_final_block
    @test length(certificate.factors) > 0
    @test verify_factorization(core, certificate.factors)
    @test verify_factorization(core, elementary_factorization(core))
    @test Suslin._verify_laurent_column_peel_replay(certificate)
    @test [step.dimension for step in certificate.peel_steps] == [6, 5, 4, 3]
    for step in certificate.peel_steps
        _issue57_assert_step(step)
    end

    corrupted_left = _issue57_corrupt_left_factor(certificate)
    @test !Suslin._verify_laurent_column_peel_replay(corrupted_left)
    @test !verify_factorization(core, corrupted_left.factors)

    deleted_clear = _issue57_delete_clearing_factor(certificate)
    @test !Suslin._verify_laurent_column_peel_replay(deleted_clear)
    @test !verify_factorization(core, deleted_clear.factors)
end

@testset "Issue 38 Laurent column peel" begin
    entry = only(ToricBuilderIssue38Cases.catalog().cases)
    R = entry.ring.object
    u, v = entry.ring.generators

    row_expected = matrix(R, [
        u      u * v;
        zero(R) u^-1
    ])
    column_expected = matrix(R, [
        v^-1   u * v;
        zero(R) v
    ])

    _issue57_assert_core(entry.normalizations.row.core, row_expected)
    _issue57_assert_core(entry.normalizations.column.core, column_expected)

    original_err = try
        elementary_factorization(entry.inputs.matrix)
        nothing
    catch err
        err
    end
    @test original_err isa ArgumentError
    @test occursin("Laurent GL_n normalization boundary", sprint(showerror, original_err))
end
```

- [ ] **Step 2: Run the focused test and confirm RED**

Run:

```bash
julia --project=. -e 'include("test/expert/laurent_column_peel_issue38.jl")'
```

Expected: FAIL because `Suslin._factor_laurent_sl_column_peel` is not defined.

- [ ] **Step 3: Commit the failing test**

Run:

```bash
git add test/expert/laurent_column_peel_issue38.jl docs/superpowers/plans/2026-06-21-issue-57-laurent-column-peel.md
git commit -m "test: add Issue 38 Laurent column peel coverage"
```

Expected: commit succeeds with a known failing focused test.

---

### Task 2: Implement Recursive Laurent Column Peel

**Files:**
- Create: `src/algorithm/laurent_column_peel.jl`
- Modify: `src/Suslin.jl`
- Modify: `src/algorithm/factorization.jl`

**Interfaces:**
- Consumes: `reduce_unimodular_column`, `elementary_matrix`, `block_embedding`, `realize_sl3_local`, `verify_factorization`, and determinant helpers.
- Produces: `LaurentColumnPeelStep`, `LaurentColumnPeelFactorization`, `_factor_laurent_sl_column_peel(A)`, `_verify_laurent_column_peel_replay(certificate)`, and `_laurent_column_peel_factors(A)`.

- [ ] **Step 1: Add the implementation file**

Create `src/algorithm/laurent_column_peel.jl` with the concrete implementation described by these signatures:

```julia
struct LaurentColumnPeelStep
    dimension::Int
    input_matrix
    last_column::Vector
    left_factors::Vector
    after_left_matrix
    right_factors::Vector
    peeled_matrix
    next_block
end

struct LaurentColumnPeelFactorization
    original_matrix
    final_block
    final_local_target
    final_local_factors::Vector
    final_factors::Vector
    factors::Vector
    product
    peel_steps::Vector{LaurentColumnPeelStep}
    verification
end

function _factor_laurent_sl_column_peel(A)
    _validate_laurent_column_peel_input(A)
    factors, steps, final_block, final_target, final_local, final_2x2 =
        _laurent_column_peel_recursive(A)
    R = base_ring(A)
    product = _factor_product(factors, R, nrows(A))
    certificate = LaurentColumnPeelFactorization(
        A,
        final_block,
        final_target,
        final_local,
        final_2x2,
        factors,
        product,
        steps,
        nothing,
    )
    verification = _laurent_column_peel_verification(certificate)
    verification.overall_ok || error("internal Laurent column-peel verification failed")
    return LaurentColumnPeelFactorization(
        certificate.original_matrix,
        certificate.final_block,
        certificate.final_local_target,
        certificate.final_local_factors,
        certificate.final_factors,
        certificate.factors,
        certificate.product,
        certificate.peel_steps,
        verification,
    )
end
```

The helper implementation must:

```julia
function _validate_laurent_column_peel_input(A)
    nrows(A) == ncols(A) || throw(ArgumentError("Laurent column-peel factorization requires a square matrix"))
    nrows(A) >= 2 || throw(ArgumentError("Laurent column-peel factorization requires size at least 2"))
    R = base_ring(A)
    _is_laurent_polynomial_ring(R) || throw(ArgumentError("Laurent column-peel factorization requires a Laurent polynomial ring"))
    profile = classify_laurent_determinant(A)
    profile.classification == :one || throw(ArgumentError("Laurent column-peel factorization requires determinant-one input"))
    return nrows(A)
end

function _laurent_column_peel_step(current)
    R = base_ring(current)
    d = nrows(current)
    last_column = [current[row, d] for row in 1:d]
    left_factors = reduce_unimodular_column(last_column, R)
    left_product = _factor_product(left_factors, R, d)
    after_left = left_product * current
    after_left * matrix(R, d, 1, last_column) == after_left[:, d:d]
    right_factors = typeof(identity_matrix(R, d))[]
    for j in 1:(d - 1)
        coeff = -after_left[d, j]
        coeff == zero(R) && continue
        push!(right_factors, elementary_matrix(d, d, j, coeff, R))
    end
    right_product = _factor_product(right_factors, R, d)
    peeled = after_left * right_product
    next_block = matrix(R, [peeled[row, col] for row in 1:(d - 1), col in 1:(d - 1)])
    _is_valid_laurent_column_peel_step_data(d, current, left_factors, after_left, right_factors, peeled, next_block) ||
        throw(ArgumentError("Laurent column-peel step failed exact replay"))
    return LaurentColumnPeelStep(d, current, last_column, left_factors, after_left, right_factors, peeled, next_block)
end
```

Use `_inverse_elementary_sequence(factors)` to reverse the sequence and negate
the unique off-diagonal entry in each elementary factor. Use
`_embed_upper_left_factors(factors, R, d)` to embed recursive `(d - 1) x (d - 1)`
factors into the upper-left block of the current `d x d` stage.

- [ ] **Step 2: Include the new file**

In `src/Suslin.jl`, include the file after `algorithm/sl3_local.jl`:

```julia
include("algorithm/laurent_column_peel.jl")
```

- [ ] **Step 3: Route Laurent determinant-one fallbacks**

In `src/algorithm/factorization.jl`, replace the Laurent determinant-one
branch with a helper:

```julia
function _laurent_sl_fallback_factorization(A)
    try
        reduction = reduce_sln_to_sl3(A)
        verify_factorization(A, reduction.factors) && return reduction.factors
    catch err
        err isa InterruptException && rethrow()
        err isa ArgumentError || rethrow()
    end

    certificate = _factor_laurent_sl_column_peel(A)
    verify_factorization(A, certificate.factors) && return certificate.factors
    error("internal Laurent column-peel factorization failed exact verification")
end
```

Use `_laurent_sl_fallback_factorization(A)` when
`ring_profile == :laurent && normalization.determinant_classification == :one`.

- [ ] **Step 4: Run the focused test and confirm GREEN**

Run:

```bash
julia --project=. -e 'include("test/expert/laurent_column_peel_issue38.jl")'
```

Expected: PASS.

- [ ] **Step 5: Commit the implementation**

Run:

```bash
git add src/Suslin.jl src/algorithm/factorization.jl src/algorithm/laurent_column_peel.jl
git commit -m "feat: add Laurent column peel factorization"
```

Expected: commit succeeds.

---

### Task 3: Update Fixture Expectations, Registration, and Verification

**Files:**
- Modify: `test/fixtures/toricbuilder_issue38_cases.jl`
- Modify: `test/internal/toricbuilder_issue38_fixture.jl`
- Modify: `test/runtests.jl`

**Interfaces:**
- Consumes: `elementary_factorization` and the Issue #38 fixture.
- Produces: internal fixture metadata that reflects supported normalized cores and expert-suite registration for Issue #57.

- [ ] **Step 1: Update fixture status metadata**

In `test/fixtures/toricbuilder_issue38_cases.jl`, replace the normalized core
failure status with:

```julia
_supported_column_peel_status() = (;
    kind = :supported_column_peel,
    factorization_path = :laurent_column_peel,
)
```

Use `_supported_column_peel_status()` for both
`normalizations.row.expected_current_status` and
`normalizations.column.expected_current_status`. Keep the original matrix entry
as the same Issue #38 `Q` block with determinant `u*v`.

- [ ] **Step 2: Update the internal fixture validator**

In `test/internal/toricbuilder_issue38_fixture.jl`, keep
`TORICBUILDER_ISSUE38_FAILURE_SUBSTRINGS` for the original `Q` boundary and add:

```julia
function _assert_supported_column_peel_core(matrix, status, label::AbstractString)
    hasproperty(status, :kind) || throw(ArgumentError("fixture $(label) status missing kind"))
    status.kind == :supported_column_peel || throw(ArgumentError("fixture $(label) expected supported column peel status"))
    hasproperty(status, :factorization_path) || throw(ArgumentError("fixture $(label) status missing factorization path"))
    status.factorization_path == :laurent_column_peel || throw(ArgumentError("fixture $(label) expected Laurent column peel path"))
    factors = elementary_factorization(matrix)
    isempty(factors) && throw(ArgumentError("fixture $(label) expected nonempty factors"))
    verify_factorization(matrix, factors) || throw(ArgumentError("fixture $(label) factorization did not verify"))
    return true
end
```

Call `_assert_supported_column_peel_core` for row and column cores instead of
`_assert_expected_factorization_failure`. Keep `_assert_expected_factorization_failure`
available for a new original-`Q` check:

```julia
function _assert_original_q_remains_unsupported(entry)
    err = try
        elementary_factorization(entry.inputs.matrix)
        nothing
    catch caught
        caught
    end
    err isa ArgumentError || throw(ArgumentError("fixture $(entry.id) original Q must remain unsupported"))
    occursin("Laurent GL_n normalization boundary", sprint(showerror, err)) ||
        throw(ArgumentError("fixture $(entry.id) original Q unsupported boundary changed"))
    return true
end
```

Invoke `_assert_original_q_remains_unsupported(entry)` from
`validate_toricbuilder_issue38_fixture(entry)`.

- [ ] **Step 3: Register the expert test**

In `test/runtests.jl`, add after `expert/sln_to_sl3_diagnostics.jl`:

```julia
"expert/laurent_column_peel_issue38.jl",
```

- [ ] **Step 4: Run focused, internal, expert, and package verification**

Run:

```bash
julia --project=. -e 'include("test/expert/laurent_column_peel_issue38.jl")'
julia --project=. -e 'include("test/internal/toricbuilder_issue38_fixture.jl")'
julia --project=. test/runtests.jl expert
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: all commands pass.

- [ ] **Step 5: Commit fixture and registration changes**

Run:

```bash
git add test/fixtures/toricbuilder_issue38_cases.jl test/internal/toricbuilder_issue38_fixture.jl test/runtests.jl
git commit -m "test: route Issue 38 cores through column peel"
```

Expected: commit succeeds.

## Self-Review

- Spec coverage: the plan covers the column-peel algorithm, final `2 x 2`
  blocks, replay metadata, row and column Issue #38 cores, negative controls,
  original `Q` out-of-scope behavior, test registration, and package
  verification.
- Placeholder scan: no unresolved placeholders remain.
- Type consistency: planned helper names match the test and implementation
  names: `_factor_laurent_sl_column_peel`, `_verify_laurent_column_peel_replay`,
  `LaurentColumnPeelStep`, and `LaurentColumnPeelFactorization`.
