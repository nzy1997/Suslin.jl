# Issue 111 Polynomial Factorization Route Shell Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Route ordinary-polynomial `elementary_factorization(A)` through the verified route certificate shell while preserving the public factor-returning API.

**Architecture:** Add a small internal polynomial shell in `src/algorithm/factorization.jl`: build a route certificate, verify it, then return only the verified factors for supported certificates. Keep Laurent dispatch on the existing path, and use staged certificate evidence only to reproduce existing polynomial `ArgumentError` messages for unsupported inputs.

**Tech Stack:** Julia, Oscar polynomial matrices, Suslin internal route certificates, Test stdlib.

## Global Constraints

- Preserve `elementary_factorization(A)` as a public factor-returning API.
- No public route may return factors unless the route certificate replay and `verify_factorization(A, factors)` both succeed.
- Keep existing Laurent behavior untouched.
- Laurent determinant-one column-peel fallback and Laurent `GL_n` determinant-correction errors must not change.
- Do not broaden the supported polynomial family.
- Do not add Quillen routing.
- Preserve polynomial route ordering: local `SL_3`, then existing block-local `SL_n` reduction, then later staged failures.
- Unsupported polynomial inputs must throw clear staged `ArgumentError` messages rather than returning empty or unchecked factor lists.
- Required focused public command: `julia --project=. -e 'include("test/public/factorization_driver_shell.jl")'`.
- Required public suite command: `julia --project=. test/runtests.jl public`.
- Required package verification command: `julia --project=. -e 'using Pkg; Pkg.test()'`.

---

## File Structure

- Modify `src/algorithm/factorization.jl`: add `_polynomial_verified_route_factors`, `_polynomial_verified_certificate_factors`, staged evidence helpers that do not recurse through the public dispatcher, and route the ordinary-polynomial branch of `elementary_factorization(A)` through the shell.
- Modify `test/public/factorization_driver_shell.jl`: assert supported polynomial public factors match verified route certificate factors, and keep staged unsupported error assertions.
- Modify `test/expert/park_woodburn_route_certificate.jl`: add a public-dispatch negative control that temporarily forces the route constructor to return a corrupted supported certificate and checks that `elementary_factorization(A)` refuses it.
- Create `docs/superpowers/specs/2026-06-23-issue-111-polynomial-factorization-route-shell-design.md`: design record.
- Create `docs/superpowers/plans/2026-06-23-issue-111-polynomial-factorization-route-shell.md`: this implementation plan.

---

### Task 1: Add RED Coverage for Verified Public Shell Behavior

**Files:**
- Modify: `test/public/factorization_driver_shell.jl`
- Modify: `test/expert/park_woodburn_route_certificate.jl`

**Interfaces:**
- Consumes: existing `_polynomial_factorization_route_certificate(A; route=nothing)`.
- Consumes: existing `PolynomialFactorizationRouteCertificate`.
- Produces: failing coverage that requires public polynomial dispatch to consult route certificates before returning factors.

- [ ] **Step 1: Update public shell assertions**

In `test/public/factorization_driver_shell.jl`, after:

```julia
factors = elementary_factorization(supported)
@test verify_factorization(supported, factors)
```

add:

```julia
supported_cert = Suslin._polynomial_factorization_route_certificate(supported)
@test supported_cert.route == :fast_local_sl3
@test Suslin._verify_polynomial_factorization_route_certificate(supported_cert)
@test factors == supported_cert.factors
```

After:

```julia
larger_factors = elementary_factorization(larger_sl)
@test isempty(larger_factors)
@test verify_factorization(larger_sl, larger_factors)
```

add:

```julia
larger_cert = Suslin._polynomial_factorization_route_certificate(larger_sl)
@test larger_cert.route == :disjoint_local_blocks
@test Suslin._verify_polynomial_factorization_route_certificate(larger_cert)
@test larger_factors == larger_cert.factors
```

- [ ] **Step 2: Add expert public-dispatch corruption control**

In `test/expert/park_woodburn_route_certificate.jl`, add this helper near the other helpers:

```julia
function _pw_captured_error(f)
    try
        f()
        return nothing
    catch err
        return err
    end
end
```

At the end of the `"Park-Woodburn polynomial route certificates"` testset, add:

```julia
    public_bad_factors = copy(fast_cert.factors)
    public_bad_factors[1] = identity_matrix(R, n)
    public_bad_cert = _pw_replace_certificate(fast_cert; factors = public_bad_factors)
    matrix_type = typeof(fast_cert.matrix)
    @eval Suslin function _polynomial_factorization_route_certificate(
            A::$matrix_type;
            route = nothing)
        return $public_bad_cert
    end
    injected_method = which(
        Suslin._polynomial_factorization_route_certificate,
        (matrix_type,),
    )
    try
        public_err = _pw_captured_error(() -> elementary_factorization(fast_cert.matrix))
        @test public_err isa ErrorException
        @test occursin(
            "internal polynomial factorization route certificate verification failed",
            sprint(showerror, public_err),
        )
    finally
        Base.delete_method(injected_method)
    end
```

- [ ] **Step 3: Run public focused test**

Run:

```bash
julia --project=. -e 'include("test/public/factorization_driver_shell.jl")'
```

Expected: PASS. These assertions document existing supported factor outputs before adding the negative control.

- [ ] **Step 4: Run expert focused test and confirm RED**

Run:

```bash
julia --project=. -e 'include("test/expert/park_woodburn_route_certificate.jl")'
```

Expected: FAIL because `elementary_factorization(fast_cert.matrix)` does not yet consult the route certificate constructor and therefore returns unchecked factors instead of the injected verification error.

---

### Task 2: Implement the Verified Ordinary-Polynomial Route Shell

**Files:**
- Modify: `src/algorithm/factorization.jl`

**Interfaces:**
- Produces: `_polynomial_verified_route_factors(A)`.
- Produces: `_polynomial_verified_certificate_factors(cert)`.
- Updates: `elementary_factorization(A)` ordinary-polynomial dispatch.
- Updates: `_polynomial_staged_failure_evidence(A)` to avoid recursion through `elementary_factorization(A)`.

- [ ] **Step 1: Make staged evidence replay independent of public dispatch**

Replace `_polynomial_staged_failure_evidence(A)` with:

```julia
function _polynomial_staged_failure_evidence(A)
    R = base_ring(A)
    ring_profile = _factorization_ring_profile(R)
    X = _supported_local_sl3_generator(A, R, ring_profile)
    if X !== nothing
        return (; error_type = :none, message = "")
    end

    if nrows(A) > 3
        try
            reduce_sln_to_sl3(A)
            return (; error_type = :none, message = "")
        catch err
            err isa InterruptException && rethrow()
            err isa ArgumentError || rethrow()
            return (;
                error_type = Symbol(nameof(typeof(err))),
                message = sprint(showerror, err),
            )
        end
    end

    try
        _throw_staged_factorization_failure(A, :polynomial, nothing)
    catch err
        err isa InterruptException && rethrow()
        err isa ArgumentError || rethrow()
        return (;
            error_type = Symbol(nameof(typeof(err))),
            message = sprint(showerror, err),
        )
    end

    return (; error_type = :none, message = "")
end
```

- [ ] **Step 2: Add the verified certificate factor shell**

Add these helpers before `_laurent_sl_fallback_factorization(A)`:

```julia
function _polynomial_verified_route_factors(A)
    certificate = _polynomial_factorization_route_certificate(A)
    return _polynomial_verified_certificate_factors(certificate)
end

function _polynomial_verified_certificate_factors(certificate)
    if certificate.status == :supported
        if _verify_polynomial_factorization_route_certificate(certificate) &&
                verify_factorization(certificate.matrix, certificate.factors)
            return certificate.factors
        end
        error("internal polynomial factorization route certificate verification failed")
    elseif certificate.status == :staged
        _throw_polynomial_staged_certificate_failure(certificate)
    end

    throw(ArgumentError("unsupported polynomial factorization route certificate status $(certificate.status)"))
end

function _throw_polynomial_staged_certificate_failure(certificate)
    evidence = certificate.evidence
    if hasproperty(evidence, :message) &&
            evidence.message isa AbstractString &&
            !isempty(evidence.message)
        throw(ArgumentError(evidence.message))
    end

    _throw_staged_factorization_failure(certificate.matrix, :polynomial, nothing)
end
```

- [ ] **Step 3: Route polynomial public dispatch through the shell**

In `elementary_factorization(A)`, replace the polynomial direct route block:

```julia
if ring_profile == :polynomial
    _require_polynomial_sl_determinant(normalized_A)
end

X = _supported_local_sl3_generator(normalized_A, normalized_R, ring_profile)
if X !== nothing
    return realize_sl3_local(
        normalized_A[1, 1],
        normalized_A[1, 2],
        normalized_A[2, 1],
        normalized_A[2, 2],
        X,
    )
end
```

with:

```julia
if ring_profile == :polynomial
    _require_polynomial_sl_determinant(normalized_A)
    return _polynomial_verified_route_factors(normalized_A)
end

X = _supported_local_sl3_generator(normalized_A, normalized_R, ring_profile)
if X !== nothing
    return realize_sl3_local(
        normalized_A[1, 1],
        normalized_A[1, 2],
        normalized_A[2, 1],
        normalized_A[2, 2],
        X,
    )
end
```

Leave the Laurent fallback block below this unchanged.

- [ ] **Step 4: Run focused tests and confirm GREEN**

Run:

```bash
julia --project=. -e 'include("test/expert/park_woodburn_route_certificate.jl")'
julia --project=. -e 'include("test/public/factorization_driver_shell.jl")'
```

Expected: both PASS.

---

### Task 3: Final Verification and Pull Request Commit

**Files:**
- Modify: implementation, tests, design, and plan files from Tasks 1-2.

**Interfaces:**
- Consumes completed verified polynomial shell.
- Produces a pull request branch containing this issue's implementation.

- [ ] **Step 1: Run issue-required public commands**

Run:

```bash
julia --project=. -e 'include("test/public/factorization_driver_shell.jl")'
julia --project=. test/runtests.jl public
```

Expected: both exit 0.

- [ ] **Step 2: Run package verification**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: exits 0.

- [ ] **Step 3: Prepare the PR commit**

Because this Agent Desk sandbox may prevent local `.git` ref updates, prepare
the PR commit from the verified working tree contents. The commit must include:

```text
docs/superpowers/specs/2026-06-23-issue-111-polynomial-factorization-route-shell-design.md
docs/superpowers/plans/2026-06-23-issue-111-polynomial-factorization-route-shell.md
src/algorithm/factorization.jl
test/public/factorization_driver_shell.jl
test/expert/park_woodburn_route_certificate.jl
```

Commit message:

```text
feat: dispatch polynomial factorization through route shell
```

---

## Plan Self-Review

- The plan covers public polynomial routing through the verified certificate shell.
- The expert negative control fails before implementation and passes only when public dispatch consults the route certificate constructor.
- The implementation explicitly avoids recursive staged evidence through `elementary_factorization(A)`.
- Laurent code remains on the existing branch and is not modified.
- Verification includes the two issue commands and the required package command.
