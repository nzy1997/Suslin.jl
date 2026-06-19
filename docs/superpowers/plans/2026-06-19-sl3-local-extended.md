# Extended Local SL3 Solver Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Extend `realize_sl3_local` to recognize and solve explicit embedded `SL_2` unit-pivot special forms inside local `SL_3`.

**Architecture:** Keep the current parameter-tuple API and add a matrix-form convenience method. Split the local path into recognition, constructive solving, and exact verification helpers. Preserve the two existing unipotent open Gauss-cell factor sequences before trying the broader unit-pivot branches.

**Tech Stack:** Julia, Oscar, AbstractAlgebra, Test.

## Global Constraints

- Do not implement general `SL_n` reduction.
- Do not add a ToricBuilder acceptance harness.
- Do not add a general `SL_2` solver.
- Preserve the existing `realize_sl3_local(p, q, r, s, X; check_monic=true)` entry point.
- Keep the optional monicity assumption explicit; Laurent examples must use `check_monic=false`.
- Unsupported determinant-one embedded local `SL_3` inputs outside the implemented family must throw an `ArgumentError` whose message includes `staged local SL_3 solver failure`.
- Focused issue command: `julia --project=. -e 'include("test/expert/sl3_local_extended.jl")'`.
- Package verification command: `julia --project=. -e 'using Pkg; Pkg.test()'`.
- Full-suite verification command from issue #21: `julia --project=. test/runtests.jl all`.
- Do not commit `Manifest.toml`.

---

## File Structure

- Modify `src/algorithm/sl3_local.jl`: add matrix recognition, family recognition, unit-pivot constructive formulas, and local exact verification.
- Create `test/expert/sl3_local_extended.jl`: focused coverage for preserved open slices, new polynomial and Laurent unit-pivot families, exact products, matrix input, and staged unsupported failure.
- Modify `test/runtests.jl`: register the new expert file after `expert/sl3_local.jl`.

---

### Task 1: Extended Local SL3 Special Forms

**Files:**
- Modify: `src/algorithm/sl3_local.jl`
- Create: `test/expert/sl3_local_extended.jl`
- Modify: `test/runtests.jl`

**Interfaces:**
- Consumes: `elementary_matrix(n, i, j, a, R)`, `_same_base_ring(left, right)`, `_is_laurent_polynomial_ring(R)`, and `verify_factorization(A, factors)`.
- Produces: `realize_sl3_local(A, X; check_monic=true)` and internal helpers `_recognize_sl3_local_matrix`, `_recognize_sl3_local_parameters`, `_realize_sl3_local_form`, `_sl3_diagonal_unit_factors`, `_unit_inverse_or_nothing`, `_verify_sl3_local_factorization`, and `_throw_staged_sl3_local_failure`.

- [ ] **Step 1: Write the failing extended expert test**

Create `test/expert/sl3_local_extended.jl` with:

```julia
using Test
using Suslin
using Oscar

function _sl3_extended_target(p, q, r, s, R)
    return matrix(R, [
        p       q       zero(R);
        r       s       zero(R);
        zero(R) zero(R) one(R)
    ])
end

function _sl3_extended_product(factors, R)
    product = identity_matrix(R, 3)
    for factor in factors
        product *= factor
    end
    return product
end

function _sl3_extended_captured_error(f)
    try
        f()
        return nothing
    catch err
        return err
    end
end

function _test_exact_sl3_factorization(target, factors)
    R = base_ring(target)
    @test target == _sl3_extended_product(factors, R)
    @test verify_factorization(target, factors)
end

@testset "extended local SL3 special-form realization" begin
    R, (X,) = Oscar.polynomial_ring(QQ, ["X"])

    q_open = one(R)
    r_open = X
    p_open = one(R) + q_open * r_open
    s_open = one(R)
    open_target = _sl3_extended_target(p_open, q_open, r_open, s_open, R)
    open_factors = Suslin.realize_sl3_local(p_open, q_open, r_open, s_open, X)
    @test length(open_factors) == 2
    _test_exact_sl3_factorization(open_target, open_factors)

    p_dual = one(R)
    q_dual = X + one(R)
    r_dual = X
    s_dual = one(R) + q_dual * r_dual
    dual_target = _sl3_extended_target(p_dual, q_dual, r_dual, s_dual, R)
    dual_factors = Suslin.realize_sl3_local(p_dual, q_dual, r_dual, s_dual, X)
    @test length(dual_factors) == 2
    _test_exact_sl3_factorization(dual_target, dual_factors)

    s_unit = R(2)
    q_unit = one(R)
    r_unit = 2 * X
    p_unit = X + R(1 // 2)
    unit_s_target = _sl3_extended_target(p_unit, q_unit, r_unit, s_unit, R)
    unit_s_factors = Suslin.realize_sl3_local(unit_s_target, X)
    @test length(unit_s_factors) > 2
    _test_exact_sl3_factorization(unit_s_target, unit_s_factors)

    L, (x, y) = suslin_laurent_polynomial_ring(QQ, ["x", "y"])
    p_laurent = x
    q_laurent = x * y
    r_laurent = one(L)
    s_laurent = x^-1 + y
    laurent_target = _sl3_extended_target(p_laurent, q_laurent, r_laurent, s_laurent, L)
    laurent_factors = Suslin.realize_sl3_local(
        p_laurent,
        q_laurent,
        r_laurent,
        s_laurent,
        x;
        check_monic=false,
    )
    @test length(laurent_factors) > 2
    _test_exact_sl3_factorization(laurent_target, laurent_factors)

    unsupported = _sl3_extended_target(X, -one(R), one(R), zero(R), R)
    unsupported_err = _sl3_extended_captured_error(() -> Suslin.realize_sl3_local(unsupported, X))
    @test unsupported_err isa ArgumentError
    @test occursin("staged local SL_3 solver failure", sprint(showerror, unsupported_err))
end
```

Modify `test/runtests.jl` so the expert list contains the new file immediately after `expert/sl3_local.jl`:

```julia
        "expert/sl3_local.jl",
        "expert/sl3_local_extended.jl",
```

- [ ] **Step 2: Run the focused test to verify RED**

Run:

```bash
julia --project=. -e 'include("test/expert/sl3_local_extended.jl")'
```

Expected: FAIL before implementation because `Suslin.realize_sl3_local(unit_s_target, X)` has no matrix-method implementation yet. The failure should be a `MethodError` for `realize_sl3_local(::MatElem, ...)` or the local Oscar matrix type.

- [ ] **Step 3: Implement recognition, construction, and exact verification**

Replace `src/algorithm/sl3_local.jl` with:

```julia
function realize_sl3_local(A, X; check_monic::Bool=true)
    form = _recognize_sl3_local_matrix(A, X; check_monic)
    return _realize_sl3_local_form(form)
end

function realize_sl3_local(p, q, r, s, X; check_monic::Bool=true)
    form = _recognize_sl3_local_parameters(p, q, r, s, X; check_monic)
    return _realize_sl3_local_form(form)
end

function _recognize_sl3_local_matrix(A, X; check_monic::Bool=true)
    nrows(A) == 3 || throw(ArgumentError("local SL_3 special-form recognition failed: A must be 3x3"))
    ncols(A) == 3 || throw(ArgumentError("local SL_3 special-form recognition failed: A must be 3x3"))

    R = base_ring(A)
    parent(X) === R || throw(ArgumentError("A and X must lie in the same polynomial or Laurent polynomial ring"))
    if A[1, 3] != zero(R) || A[2, 3] != zero(R) ||
            A[3, 1] != zero(R) || A[3, 2] != zero(R) ||
            A[3, 3] != one(R)
        throw(ArgumentError("local SL_3 special-form recognition failed: A must be an embedded 2x2 block with trailing identity"))
    end

    return _recognize_sl3_local_parameters(A[1, 1], A[1, 2], A[2, 1], A[2, 2], X; check_monic)
end

function _recognize_sl3_local_parameters(p, q, r, s, X; check_monic::Bool=true)
    R = parent(p)
    parent(q) === R || throw(ArgumentError("p, q, r, s, and X must lie in the same polynomial or Laurent polynomial ring"))
    parent(r) === R || throw(ArgumentError("p, q, r, s, and X must lie in the same polynomial or Laurent polynomial ring"))
    parent(s) === R || throw(ArgumentError("p, q, r, s, and X must lie in the same polynomial or Laurent polynomial ring"))
    parent(X) === R || throw(ArgumentError("p, q, r, s, and X must lie in the same polynomial or Laurent polynomial ring"))

    ring_gens = collect(gens(R))
    var_idx = findfirst(isequal(X), ring_gens)
    var_idx === nothing && throw(ArgumentError("X must be one of the polynomial or Laurent polynomial ring generators"))

    if check_monic
        if _is_laurent_polynomial_ring(R)
            throw(ArgumentError("p monicity check is only supported for ordinary polynomial local SL_3 inputs; pass check_monic=false only when the caller has discharged the local monicity assumption"))
        end
        _is_monic_in_variable(p, var_idx, R) || throw(ArgumentError("p must be monic in X"))
    end

    target = matrix(R, [
        p q zero(R);
        r s zero(R);
        zero(R) zero(R) one(R)
    ])
    det(target) == one(R) || throw(ArgumentError("constructed matrix must have determinant 1"))

    if s == one(R) && p == one(R) + q * r
        return (; family = :open_s_one, R, p, q, r, s, target)
    end

    if p == one(R) && s == one(R) + q * r
        return (; family = :open_p_one, R, p, q, r, s, target)
    end

    s_inverse = _unit_inverse_or_nothing(s)
    if s_inverse !== nothing
        return (; family = :s_unit, R, p, q, r, s, target, pivot_inverse = s_inverse)
    end

    p_inverse = _unit_inverse_or_nothing(p)
    if p_inverse !== nothing
        return (; family = :p_unit, R, p, q, r, s, target, pivot_inverse = p_inverse)
    end

    _throw_staged_sl3_local_failure("supported families require one open unipotent slice or a unit diagonal pivot p or s")
end

function _realize_sl3_local_form(form)
    R = form.R

    factors = if form.family == :open_s_one
        [
            elementary_matrix(3, 1, 2, form.q, R),
            elementary_matrix(3, 2, 1, form.r, R),
        ]
    elseif form.family == :open_p_one
        [
            elementary_matrix(3, 2, 1, form.r, R),
            elementary_matrix(3, 1, 2, form.q, R),
        ]
    elseif form.family == :s_unit
        s_inverse = form.pivot_inverse
        vcat(
            [elementary_matrix(3, 1, 2, form.q * s_inverse, R)],
            _sl3_diagonal_unit_factors(s_inverse, R),
            [elementary_matrix(3, 2, 1, form.r * s_inverse, R)],
        )
    elseif form.family == :p_unit
        p_inverse = form.pivot_inverse
        vcat(
            [elementary_matrix(3, 2, 1, form.r * p_inverse, R)],
            _sl3_diagonal_unit_factors(form.p, R),
            [elementary_matrix(3, 1, 2, form.q * p_inverse, R)],
        )
    else
        _throw_staged_sl3_local_failure("unrecognized local solver family $(form.family)")
    end

    _verify_sl3_local_factorization(form.target, factors)
    return factors
end

function _sl3_diagonal_unit_factors(u, R)
    u_inverse = inv(u)
    return [
        elementary_matrix(3, 1, 2, u, R),
        elementary_matrix(3, 2, 1, -u_inverse, R),
        elementary_matrix(3, 1, 2, u, R),
        elementary_matrix(3, 1, 2, -one(R), R),
        elementary_matrix(3, 2, 1, one(R), R),
        elementary_matrix(3, 1, 2, -one(R), R),
    ]
end

function _unit_inverse_or_nothing(value)
    try
        is_unit(value) || return nothing
        return inv(value)
    catch err
        if err isa ArgumentError || err isa MethodError || err isa ErrorException
            return nothing
        end
        rethrow()
    end
end

function _verify_sl3_local_factorization(target, factors)
    R = base_ring(target)
    product = identity_matrix(R, 3)
    for factor in factors
        nrows(factor) == 3 || error("local SL_3 exact verification failed: factor has wrong row count")
        ncols(factor) == 3 || error("local SL_3 exact verification failed: factor has wrong column count")
        _same_base_ring(base_ring(factor), R) || error("local SL_3 exact verification failed: factor base ring mismatch")
        product *= factor
    end

    product == target || error("local SL_3 exact verification failed: factor product does not equal target")
    return nothing
end

function _throw_staged_sl3_local_failure(reason::AbstractString)
    throw(ArgumentError("staged local SL_3 solver failure: $(reason)"))
end

function _is_monic_in_variable(p, var_idx::Int, R)
    iszero(p) && return false

    target_degree = degree(p, var_idx)
    target_degree < 0 && return false

    ring_gens = collect(gens(R))
    total = zero(R)
    for (coeff, exponents) in zip(AbstractAlgebra.coefficients(p), AbstractAlgebra.exponent_vectors(p))
        exponents[var_idx] == target_degree || continue
        term = R(coeff)
        for idx in 1:length(ring_gens)
            idx == var_idx && continue
            exponent = exponents[idx]
            exponent == 0 && continue
            term *= ring_gens[idx]^exponent
        end
        total += term
    end

    return total == one(R)
end
```

- [ ] **Step 4: Run focused GREEN verification**

Run:

```bash
julia --project=. -e 'include("test/expert/sl3_local_extended.jl")'
```

Expected: PASS for the `extended local SL3 special-form realization` testset.

- [ ] **Step 5: Run local regression checks**

Run:

```bash
julia --project=. -e 'include("test/expert/sl3_local.jl")'
julia --project=. test/runtests.jl expert
```

Expected: both commands exit 0. The original `local SL3 special-form realization` testset should still report 7 passing tests.

- [ ] **Step 6: Commit implementation**

Run:

```bash
git add src/algorithm/sl3_local.jl test/expert/sl3_local_extended.jl test/runtests.jl
git commit -m "feat: extend local sl3 unit-pivot solver"
```

Expected: commit succeeds and `Manifest.toml` is not staged.

---

## Plan Self-Review

- Spec coverage: Task 1 covers recognition, constructive solving, exact verification, preserved open paths, two new supported families, matrix input, focused tests, expert registration, and staged negative control.
- Placeholder scan: no incomplete-marker patterns remain.
- Type consistency: the produced helper names match the interfaces and implementation snippets in the task.
