# Issue 15 SL_n to Local SL_3 Reduction Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an inspectable reduction layer from supported larger `SL_n` matrices to local `SL_3` obligations, then assemble exact embedded factor sequences.

**Architecture:** The reduction layer lives in a focused algorithm file and returns typed records for local obligations and top-level reassembly metadata. The supported domain is block-local: full consecutive 3x3 coordinate blocks, optional caller-specified 3-index block locations, and identity outside those blocks. `elementary_factorization` delegates to the reduction layer only after the current direct local `SL_3` path fails for `n > 3`.

**Tech Stack:** Julia, Oscar matrices and polynomial/Laurent polynomial rings, existing Suslin helpers `block_embedding`, `embed_factor_sequence`, `realize_sl3_local`, `normalize_laurent_gl_matrix`, and `verify_factorization`.

## Global Constraints

- Preserve exact arithmetic and parent rings; never coerce between different matrix parents except through existing `_coerce_into_ring` helper paths.
- Unsupported but valid matrices must throw `ArgumentError` with `staged SL_n to local SL_3 reduction failure`.
- Ordinary polynomial reduction is limited to univariate base rings.
- Laurent local solving uses `check_monic=false` only after determinant normalization has produced the local core.
- The issue-specific verification command is `julia --project=. -e 'include("test/expert/sln_to_sl3_reduction.jl")'`.
- The Agent Desk package command is `julia --project=. -e 'using Pkg; Pkg.test()'`.
- The documented full suite from #21 is `julia --project=. test/runtests.jl all`.

---

## File Structure

- Create `src/algorithm/sln_to_sl3_reduction.jl`: public record structs, reduction constructor, exact product helpers, staged errors, and verification.
- Modify `src/Suslin.jl`: export the new records/functions and include the new algorithm file.
- Modify `src/algorithm/factorization.jl`: delegate supported `n > 3` inputs to `reduce_sln_to_sl3`.
- Create `test/expert/sln_to_sl3_reduction.jl`: focused issue #15 coverage and negative controls.
- Modify `test/runtests.jl`: register the new expert test file.
- Modify `test/public/api_surface.jl`: assert the new public API is exported.
- Modify `test/public/factorization_driver_shell.jl`: update the old `n > 3` rejection expectation so supported block-local inputs now pass, while unsupported larger inputs still fail with a staged reduction error.

---

### Task 1: Focused Failing Tests

**Files:**
- Create: `test/expert/sln_to_sl3_reduction.jl`
- Modify: `test/runtests.jl`
- Modify: `test/public/api_surface.jl`
- Modify: `test/public/factorization_driver_shell.jl`

**Interfaces:**
- Consumes: planned `SL3LocalObligation`, `SLNToSL3Reduction`, `reduce_sln_to_sl3(A; block_locations=nothing)`, `verify_sln_to_sl3_reduction(reduction)`.
- Produces: failing tests that define the required public behavior.

- [ ] **Step 1: Write the focused expert test**

Create `test/expert/sln_to_sl3_reduction.jl` with this content:

```julia
using Test
using Suslin
using Oscar

function _issue15_product(factors, R, n::Int)
    product = identity_matrix(R, n)
    for factor in factors
        product *= factor
    end
    return product
end

function _issue15_local_block(p, q, r, s, R)
    return matrix(R, [
        p q zero(R);
        r s zero(R);
        zero(R) zero(R) one(R)
    ])
end

function _issue15_supported_matrix(R, blocks, n::Int)
    product = identity_matrix(R, n)
    for (block, indices) in blocks
        product *= block_embedding(block, n, indices)
    end
    return product
end

function _issue15_captured_error(f)
    try
        f()
        return nothing
    catch err
        return err
    end
end

function _issue15_assert_reduction(A, expected_obligations::Int)
    reduction = reduce_sln_to_sl3(A)
    @test reduction isa SLNToSL3Reduction
    @test length(reduction.obligations) == expected_obligations
    @test verify_sln_to_sl3_reduction(reduction)
    @test verify_factorization(A, reduction.factors)
    @test _issue15_product(reduction.factors, base_ring(A), nrows(A)) == A
    @test all(obligation -> obligation isa SL3LocalObligation, reduction.obligations)
    @test all(obligation -> obligation.reassembly_data.embedded_product_ok, reduction.obligations)
    @test all(obligation -> obligation.reassembly_data.local_product_ok, reduction.obligations)
    return reduction
end

@testset "SL_n to local SL3 reduction supported examples" begin
    R, (X,) = Oscar.polynomial_ring(QQ, ["X"])

    block_a = _issue15_local_block(one(R) + X, one(R), X, one(R), R)
    block_b = _issue15_local_block(one(R), one(R) + X, X, one(R) + X + X^2, R)
    block_c = _issue15_local_block(one(R) + X, X, one(R), one(R), R)

    matrix6 = _issue15_supported_matrix(R, [(block_a, [1, 2, 3]), (block_b, [4, 5, 6])], 6)
    reduction6 = _issue15_assert_reduction(matrix6, 2)
    @test reduction6.obligations[1].block_location == [1, 2, 3]
    @test reduction6.obligations[2].block_location == [4, 5, 6]
    @test elementary_factorization(matrix6) == reduction6.factors

    dropped6 = reduction6.obligations[1].embedded_factors
    @test !verify_factorization(matrix6, dropped6)

    matrix8 = _issue15_supported_matrix(R, [(block_b, [1, 2, 3]), (block_c, [4, 5, 6])], 8)
    reduction8 = _issue15_assert_reduction(matrix8, 2)
    @test reduction8.obligations[1].block_location == [1, 2, 3]
    @test reduction8.obligations[2].block_location == [4, 5, 6]
    @test matrix8[7, 7] == one(R)
    @test matrix8[8, 8] == one(R)

    custom = _issue15_supported_matrix(R, [(block_a, [2, 4, 6])], 6)
    custom_reduction = reduce_sln_to_sl3(custom; block_locations = [[2, 4, 6]])
    @test verify_sln_to_sl3_reduction(custom_reduction)
    @test verify_factorization(custom, custom_reduction.factors)
end

@testset "SL_n to local SL3 reduction staged failures" begin
    R, (X,) = Oscar.polynomial_ring(QQ, ["X"])

    unsupported = identity_matrix(R, 6)
    unsupported[1, 4] = X
    unsupported[4, 1] = zero(R)
    unsupported_err = _issue15_captured_error(() -> reduce_sln_to_sl3(unsupported))
    @test unsupported_err isa ArgumentError
    @test occursin("staged SL_n to local SL_3 reduction failure", sprint(showerror, unsupported_err))
    @test !occursin("local SL_3 special-form recognition failed", sprint(showerror, unsupported_err))

    bad_locations_err = _issue15_captured_error(() -> reduce_sln_to_sl3(identity_matrix(R, 6); block_locations = [[1, 2, 2]]))
    @test bad_locations_err isa ArgumentError
    @test occursin("block locations", sprint(showerror, bad_locations_err))

    S, (X, Y) = Oscar.polynomial_ring(QQ, ["X", "Y"])
    multivariate = identity_matrix(S, 6)
    multivariate[1:3, 1:3] = _issue15_local_block(one(S) + X, one(S), X, one(S), S)
    multivariate_err = _issue15_captured_error(() -> reduce_sln_to_sl3(multivariate))
    @test multivariate_err isa ArgumentError
    @test occursin("ordinary polynomial reduction currently requires a univariate base ring", sprint(showerror, multivariate_err))
end
```

- [ ] **Step 2: Register the expert test**

In `test/runtests.jl`, add `"expert/sln_to_sl3_reduction.jl"` to the expert group after `"expert/sl3_local_extended.jl"`.

- [ ] **Step 3: Extend API surface tests**

In `test/public/api_surface.jl`, add:

```julia
@test isdefined(Suslin, :SL3LocalObligation)
@test isdefined(Suslin, :SLNToSL3Reduction)
@test isdefined(Suslin, :reduce_sln_to_sl3)
@test isdefined(Suslin, :verify_sln_to_sl3_reduction)
@test Suslin.SL3LocalObligation === SL3LocalObligation
@test Suslin.SLNToSL3Reduction === SLNToSL3Reduction
@test Suslin.reduce_sln_to_sl3 === reduce_sln_to_sl3
@test Suslin.verify_sln_to_sl3_reduction === verify_sln_to_sl3_reduction
```

- [ ] **Step 4: Update public driver shell expectation**

In `test/public/factorization_driver_shell.jl`, replace the old identity-4 rejection check with:

```julia
larger_sl = identity_matrix(R, 4)
larger_factors = elementary_factorization(larger_sl)
@test isempty(larger_factors)
@test verify_factorization(larger_sl, larger_factors)

unsupported_larger = identity_matrix(R, 4)
unsupported_larger[1, 4] = X
larger_err = _captured_error(() -> elementary_factorization(unsupported_larger))
@test larger_err isa ArgumentError
@test occursin("staged SL_n to local SL_3 reduction failure", sprint(showerror, larger_err))
@test !occursin("currently supports only 3x3 matrices", sprint(showerror, larger_err))
```

- [ ] **Step 5: Run focused test to verify failure**

Run: `julia --project=. -e 'include("test/expert/sln_to_sl3_reduction.jl")'`

Expected: FAIL with `UndefVarError: reduce_sln_to_sl3 not defined`.

- [ ] **Step 6: Commit failing tests**

Run:

```bash
git add test/expert/sln_to_sl3_reduction.jl test/runtests.jl test/public/api_surface.jl test/public/factorization_driver_shell.jl
git commit -m "test: cover sln to sl3 reduction"
```

---

### Task 2: Core Reduction Records and Constructor

**Files:**
- Create: `src/algorithm/sln_to_sl3_reduction.jl`
- Modify: `src/Suslin.jl`

**Interfaces:**
- Consumes: tests from Task 1 and existing helpers `block_embedding`, `embed_factor_sequence`, `realize_sl3_local`, `normalize_laurent_gl_matrix`.
- Produces: `SL3LocalObligation`, `SLNToSL3Reduction`, `reduce_sln_to_sl3`, and `verify_sln_to_sl3_reduction`.

- [ ] **Step 1: Export and include the new API**

In `src/Suslin.jl`, add exports near the factorization exports:

```julia
export SL3LocalObligation
export SLNToSL3Reduction
export reduce_sln_to_sl3
export verify_sln_to_sl3_reduction
```

Add the include after `include("algorithm/factorization.jl")`:

```julia
include("algorithm/sln_to_sl3_reduction.jl")
```

- [ ] **Step 2: Add the reduction implementation**

Create `src/algorithm/sln_to_sl3_reduction.jl` with public structs, staged errors, block-location normalization, local obligation solving, exact product verification, and the top-level reducer.

The core signatures must be:

```julia
struct SL3LocalObligation
    block_location::Vector{Int}
    ring
    target_local_matrix
    required_assumptions::Vector{Symbol}
    embedded_target
    local_factors::Vector
    embedded_factors::Vector
    reassembly_data
end

struct SLNToSL3Reduction
    ring
    size::Int
    original_matrix
    normalized_matrix
    normalization
    obligations::Vector{SL3LocalObligation}
    factors::Vector
    product
    verification
end

function reduce_sln_to_sl3(A; block_locations=nothing)
    return _construct_sln_to_sl3_reduction(A, block_locations)
end

function verify_sln_to_sl3_reduction(reduction::SLNToSL3Reduction)::Bool
    return _sln_to_sl3_reduction_verification(reduction).overall_ok
end
```

- [ ] **Step 3: Run focused test to verify implementation**

Run: `julia --project=. -e 'include("test/expert/sln_to_sl3_reduction.jl")'`

Expected: PASS.

- [ ] **Step 4: Run API surface test**

Run: `julia --project=. -e 'include("test/public/api_surface.jl")'`

Expected: PASS.

- [ ] **Step 5: Commit core implementation**

Run:

```bash
git add src/Suslin.jl src/algorithm/sln_to_sl3_reduction.jl
git commit -m "feat: add sln to sl3 reduction records"
```

---

### Task 3: Factorization Driver Delegation

**Files:**
- Modify: `src/algorithm/factorization.jl`
- Test: `test/public/factorization_driver_shell.jl`
- Test: `test/expert/sln_to_sl3_reduction.jl`

**Interfaces:**
- Consumes: `reduce_sln_to_sl3(A)` from Task 2.
- Produces: `elementary_factorization(A)` support for block-local `n > 3` matrices whose returned factors exactly multiply to the original input.

- [ ] **Step 1: Add driver fallback**

In `src/algorithm/factorization.jl`, keep the existing direct local `SL_3` path. Replace the unconditional staged failure after `X === nothing` with:

```julia
    if X !== nothing
        return realize_sl3_local(
            normalized_A[1, 1],
            normalized_A[1, 2],
            normalized_A[2, 1],
            normalized_A[2, 2],
            X,
        )
    end

    if nrows(normalized_A) > 3
        reduction = reduce_sln_to_sl3(A)
        verify_factorization(A, reduction.factors) && return reduction.factors
    end

    _throw_staged_factorization_failure(normalized_A, ring_profile, normalization)
```

- [ ] **Step 2: Run driver shell test**

Run: `julia --project=. -e 'include("test/public/factorization_driver_shell.jl")'`

Expected: PASS.

- [ ] **Step 3: Run focused reduction test**

Run: `julia --project=. -e 'include("test/expert/sln_to_sl3_reduction.jl")'`

Expected: PASS.

- [ ] **Step 4: Commit driver integration**

Run:

```bash
git add src/algorithm/factorization.jl test/public/factorization_driver_shell.jl test/expert/sln_to_sl3_reduction.jl
git commit -m "feat: route supported sln factorizations through sl3 reductions"
```

---

### Task 4: Verification Sweep

**Files:**
- Modify only if verification reveals an issue in files touched by Tasks 1-3.

**Interfaces:**
- Consumes: complete implementation.
- Produces: verified branch ready for PR.

- [ ] **Step 1: Run issue-specific verification**

Run: `julia --project=. -e 'include("test/expert/sln_to_sl3_reduction.jl")'`

Expected: PASS.

- [ ] **Step 2: Run Agent Desk package command**

Run: `julia --project=. -e 'using Pkg; Pkg.test()'`

Expected: PASS.

- [ ] **Step 3: Run documented full suite**

Run: `julia --project=. test/runtests.jl all`

Expected: PASS.

- [ ] **Step 4: Inspect final diff**

Run: `git status --short`

Expected: clean except intentional committed changes.

Run: `git log --oneline -5`

Expected: includes the design, test, implementation, and driver commits for issue #15.

---

## Self-Review

- Spec coverage: Tasks 1-3 cover inspectable local obligations, exact reassembly, 6x6 and 8x8 examples, driver integration, unsupported staged errors, API exports, and the omitted-obligation negative control.
- Placeholder scan: no task contains open-ended placeholder instructions.
- Type consistency: public names are consistently `SL3LocalObligation`, `SLNToSL3Reduction`, `reduce_sln_to_sl3`, and `verify_sln_to_sl3_reduction`.

## Execution Choice

Plan complete and saved to `docs/superpowers/plans/2026-06-19-issue-15-sln-to-sl3-reduction.md`. Two execution options:

1. Subagent-Driven (recommended) - dispatch a fresh subagent per task, review between tasks, fast iteration.
2. Inline Execution - execute tasks in this session using executing-plans, batch execution with checkpoints.

Under the standing Agent Desk answer policy, choose option 1 because it is marked recommended.
