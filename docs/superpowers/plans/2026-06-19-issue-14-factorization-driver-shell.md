# Issue 14 Factorization Driver Shell Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Refactor `elementary_factorization(A)` into a general `n >= 3` public driver shell that preserves the supported local `SL_3` path and sends unsupported inputs through staged validation and algorithm-layer failures.

**Architecture:** Keep the public API in `src/algorithm/factorization.jl`, split the driver into private validation, normalization, supported-case detection, and staged-failure helpers, and add a focused public driver-shell test file. Laurent determinant handling stays delegated to `normalize_laurent_gl_matrix(A)` from issue #20.

**Tech Stack:** Julia, Oscar, `Test`, existing Suslin private helpers.

## Global Constraints

- Do not add new public exports for this issue.
- Preserve the exact factor sequence output for the existing supported `3 x 3` univariate polynomial local `SL_3` path.
- Generic matrix validation must run before Laurent determinant normalization.
- Laurent input must call `normalize_laurent_gl_matrix(A)` before algorithm support dispatch.
- Determinant-not-one polynomial input must fail on the determinant/unit precondition before unsupported-shape dispatch.
- Valid determinant-one unsupported input must fail with a staged algorithm-layer message, not a shape-specific `currently supports only 3x3 matrices` hard stop.
- Register the focused driver-shell test in the public test group.
- Do not commit `Manifest.toml`.

---

## File Structure

- Modify `src/algorithm/factorization.jl`: replace the monolithic public driver with private helper layers and preserve `verify_factorization(A, factors)`.
- Create `test/public/factorization_driver_shell.jl`: focused public tests for supported behavior, generic validation, determinant preconditions, Laurent normalization boundary, and staged algorithm-layer failures.
- Modify `test/public/api_surface.jl`: keep API coverage and adjust only if old error assumptions conflict with the new staged shell.
- Modify `test/internal/gl_laurent_normalization.jl`: update the issue #20 boundary expectation that currently checks for the old `3 x 3` hard stop.
- Modify `test/runtests.jl`: add the new public test file to the public group.

## Task 1: Driver Shell TDD Refactor

**Files:**
- Modify: `src/algorithm/factorization.jl`
- Create: `test/public/factorization_driver_shell.jl`
- Modify: `test/internal/gl_laurent_normalization.jl`
- Modify: `test/runtests.jl`

**Interfaces:**
- Consumes: `normalize_laurent_gl_matrix(A)` and `_is_laurent_polynomial_ring(R)` from existing core modules.
- Produces: `elementary_factorization(A)` with layered validation and staged dispatch; existing callers still receive a factor sequence for the supported case.

- [ ] **Step 1: Write the failing focused public test**

Create `test/public/factorization_driver_shell.jl` with this content:

```julia
using Suslin
using Test
using Oscar

function _captured_error(f)
    try
        f()
        return nothing
    catch err
        return err
    end
end

@testset "public factorization driver shell" begin
    R, (X,) = Oscar.polynomial_ring(QQ, ["X"])

    supported = matrix(R, [
        one(R)      one(R) + X       zero(R);
        X           one(R) + X + X^2 zero(R);
        zero(R)     zero(R)          one(R)
    ])
    factors = elementary_factorization(supported)
    @test verify_factorization(supported, factors)

    larger_sl = identity_matrix(R, 4)
    larger_err = _captured_error(() -> elementary_factorization(larger_sl))
    @test larger_err isa ArgumentError
    @test occursin("SL_n reduction layer", sprint(showerror, larger_err))
    @test occursin("not yet implemented", sprint(showerror, larger_err))
    @test !occursin("currently supports only 3x3 matrices", sprint(showerror, larger_err))

    nonsquare = zero_matrix(R, 3, 4)
    nonsquare_err = _captured_error(() -> elementary_factorization(nonsquare))
    @test nonsquare_err isa ArgumentError
    @test occursin("A must be square", sprint(showerror, nonsquare_err))

    determinant_not_one = matrix(R, [
        X + one(R) zero(R) zero(R) zero(R);
        zero(R)    one(R) zero(R) zero(R);
        zero(R)    zero(R) one(R) zero(R);
        zero(R)    zero(R) zero(R) one(R)
    ])
    determinant_err = _captured_error(() -> elementary_factorization(determinant_not_one))
    @test determinant_err isa ArgumentError
    @test occursin("determinant/unit precondition", sprint(showerror, determinant_err))
    @test occursin("outside the staged SL_n factorization path", sprint(showerror, determinant_err))
    @test !occursin("SL_n reduction layer", sprint(showerror, determinant_err))
    @test !occursin("currently supports only 3x3 matrices", sprint(showerror, determinant_err))

    L, (x,) = suslin_laurent_polynomial_ring(GF(2), ["x"])
    normalizable_laurent = matrix(L, [
        x       zero(L) zero(L);
        zero(L) one(L)  zero(L);
        zero(L) zero(L) one(L)
    ])
    laurent_err = _captured_error(() -> elementary_factorization(normalizable_laurent))
    @test laurent_err isa ArgumentError
    @test occursin("Laurent GL_n normalization boundary", sprint(showerror, laurent_err))
    @test occursin("Laurent SL_n reduction layer", sprint(showerror, laurent_err))
    @test occursin("not yet implemented", sprint(showerror, laurent_err))

    non_normalizable_laurent = matrix(L, [
        x + one(L) zero(L) zero(L);
        zero(L)    one(L)  zero(L);
        zero(L)    zero(L) one(L)
    ])
    non_normalizable_err = _captured_error(() -> elementary_factorization(non_normalizable_laurent))
    @test non_normalizable_err isa ArgumentError
    @test occursin("unsupported Laurent GL_n determinant", sprint(showerror, non_normalizable_err))
    @test occursin("outside the staged SL_n factorization path", sprint(showerror, non_normalizable_err))
    @test !occursin("SL_n reduction layer", sprint(showerror, non_normalizable_err))
    @test !occursin("currently supports only 3x3 matrices", sprint(showerror, non_normalizable_err))
end
```

- [ ] **Step 2: Register the public test**

In `test/runtests.jl`, change the public group to:

```julia
const TEST_GROUP_FILES = Dict(
    "public" => [
        "public/api_surface.jl",
        "public/factorization_driver_shell.jl",
    ],
```

- [ ] **Step 3: Run the focused test to verify RED**

Run:

```bash
julia --project=. -e 'include("test/public/factorization_driver_shell.jl")'
```

Expected: fail because `elementary_factorization(identity_matrix(R, 4))` still throws the old `elementary_factorization currently supports only 3x3 matrices` message, and normalizable Laurent input does not yet report the new Laurent algorithm-layer message.

- [ ] **Step 4: Implement layered driver helpers**

Replace `elementary_factorization(A)` in `src/algorithm/factorization.jl` with this helper-based implementation and leave `verify_factorization(A, factors)` below it:

```julia
function _validate_factorization_matrix(A)
    nrows(A) == ncols(A) || throw(ArgumentError("A must be square"))
    n = nrows(A)
    n >= 3 || throw(ArgumentError("elementary_factorization requires matrices of size at least 3"))
    return n
end

function _factorization_ring_profile(R)
    _is_laurent_polynomial_ring(R) && return :laurent
    try
        collect(gens(R))
        return :polynomial
    catch err
        err isa MethodError || rethrow()
        throw(ArgumentError("A base ring is outside the supported exact polynomial or Laurent polynomial factorization path"))
    end
end

function _normalize_factorization_input(A, ring_profile::Symbol)
    if ring_profile == :laurent
        normalization = normalize_laurent_gl_matrix(A)
        return normalization.normalized_matrix, normalization
    end

    return A, nothing
end

function _require_polynomial_sl_determinant(A)
    R = base_ring(A)
    det(A) == one(R) || throw(ArgumentError("determinant/unit precondition failed: polynomial inputs must have determinant 1 and lie in the staged SL_n factorization path"))
    return nothing
end

function _supported_local_sl3_generator(A, R, ring_profile::Symbol)
    ring_profile == :polynomial || return nothing
    nrows(A) == 3 || return nothing

    ring_gens = collect(gens(R))
    length(ring_gens) == 1 || return nothing

    if A[1, 3] != zero(R) || A[2, 3] != zero(R) ||
            A[3, 1] != zero(R) || A[3, 2] != zero(R) ||
            A[3, 3] != one(R)
        return nothing
    end

    return ring_gens[1]
end

function _throw_staged_factorization_failure(A, ring_profile::Symbol, normalization)
    n = nrows(A)

    if ring_profile == :laurent
        classification = normalization.determinant_classification
        throw(ArgumentError("Laurent GL_n normalization boundary succeeded with determinant classification $(classification), but the Laurent SL_n reduction layer is not yet implemented"))
    end

    if n > 3
        throw(ArgumentError("SL_n reduction layer for n > 3 is not yet implemented in elementary_factorization"))
    end

    throw(ArgumentError("staged reduction to the supported univariate local SL_3 slice is not yet implemented in elementary_factorization"))
end

function elementary_factorization(A)
    _validate_factorization_matrix(A)

    R = base_ring(A)
    ring_profile = _factorization_ring_profile(R)
    normalized_A, normalization = _normalize_factorization_input(A, ring_profile)
    normalized_R = base_ring(normalized_A)

    if ring_profile == :polynomial
        _require_polynomial_sl_determinant(normalized_A)
    end

    X = _supported_local_sl3_generator(normalized_A, normalized_R, ring_profile)
    X === nothing && _throw_staged_factorization_failure(normalized_A, ring_profile, normalization)

    return realize_sl3_local(
        normalized_A[1, 1],
        normalized_A[1, 2],
        normalized_A[2, 1],
        normalized_A[2, 2],
        X,
    )
end
```

- [ ] **Step 5: Update the issue #20 boundary test expectation**

In `test/internal/gl_laurent_normalization.jl`, change the `normalized_then_rejected` matrix to size `3 x 3` and replace the old `currently supports only 3x3 matrices` assertion with the new Laurent staged dispatch assertion:

```julia
    normalized_then_rejected = matrix(R, [
        x zero(R) zero(R);
        zero(R) one(R) zero(R);
        zero(R) zero(R) one(R)
    ])
    err = try
        elementary_factorization(normalized_then_rejected)
        nothing
    catch caught
        caught
    end
    @test err isa ArgumentError
    @test occursin("Laurent GL_n normalization boundary", sprint(showerror, err))
    @test occursin("Laurent SL_n reduction layer", sprint(showerror, err))
```

- [ ] **Step 6: Run the focused test to verify GREEN**

Run:

```bash
julia --project=. -e 'include("test/public/factorization_driver_shell.jl")'
```

Expected: pass.

- [ ] **Step 7: Run default package tests**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: public and internal tests pass.

- [ ] **Step 8: Run full suite**

Run:

```bash
julia --project=. test/runtests.jl all
```

Expected: public, internal, and expert tests pass.

- [ ] **Step 9: Commit implementation**

Run:

```bash
git status --short
git add src/algorithm/factorization.jl test/public/factorization_driver_shell.jl test/internal/gl_laurent_normalization.jl test/runtests.jl docs/superpowers/plans/2026-06-19-issue-14-factorization-driver-shell.md
git commit -m "feat: generalize factorization driver shell"
```

Expected: commit succeeds and `Manifest.toml` is not staged.

## Plan Self-Review

- Spec coverage: the task covers generic validation, ring validation, Laurent normalization, supported-case detection, staged unsupported dispatch, focused public registration, and required verification commands.
- Placeholder scan: no task uses incomplete markers or references missing implementation details.
- Type consistency: all helpers are private Julia functions in `src/algorithm/factorization.jl`; the public interface remains `elementary_factorization(A)` returning a factor sequence for supported input.
