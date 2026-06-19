# Issue 9 Laurent Linear Solving Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a Suslin-owned exact Laurent linear solver based on the ToricBuilder reference behavior.

**Architecture:** The solver lives in one new core file and is exported as `solve_laurent_linear(A, B)`. It validates Laurent parent and row dimensions before solving, tries Oscar's native right-side solver first, and falls back to the ToricBuilder module-coordinate path only for known native Laurent capability errors.

**Tech Stack:** Julia, Oscar, AbstractAlgebra matrix APIs, existing Suslin Laurent validators, Test stdlib.

## Global Constraints

- Repository base branch is `main`; worker branch is `agent/issue-9-decide-ownership-of-laurent-linear-solving-for-s-run-1`.
- The implementation must stay scoped to Issue #9 and must not broaden into general factorization.
- Inputs are exact Laurent linear systems `A * U = B` over supported Laurent rings.
- `solve_laurent_linear(A, B)` returns an exact matrix `U` satisfying `A * U == B`.
- `base_ring(A)` must be a Laurent polynomial ring.
- `base_ring(B)` must be the same Laurent ring object as `base_ring(A)`.
- `nrows(B)` must equal `nrows(A)`.
- Invalid parent or dimension errors must happen before attempting a solve.
- No-solution failures must throw `ErrorException("No exact solution exists for A * U = B")`.
- Backend capability failures must not be confused with mathematical no-solution failures.
- Use the ToricBuilder implementation in `/Users/nzy/jcode/topological-code-decoupling/julia_code/ToricBuilder/src/core/laurent_linear_solve.jl` as reference behavior.
- Run `julia --project=. -e 'include("test/internal/laurent_linear_solve.jl")'`.
- Run `julia --project=. -e 'using Pkg; Pkg.test()'`.
- Run the documented full Suslin suite from #21: `julia --project=. test/runtests.jl all`.

---

### Task 1: Add Focused Failing Tests

**Files:**
- Create: `test/internal/laurent_linear_solve.jl`
- Modify: `test/runtests.jl`
- Modify: `test/public/api_surface.jl`

**Interfaces:**
- Consumes: existing `suslin_laurent_polynomial_ring`, `zero_matrix`, Oscar matrix multiplication.
- Produces: failing expectations for `solve_laurent_linear(A, B)` and its public export.

- [ ] **Step 1: Add the focused test file**

Create `test/internal/laurent_linear_solve.jl` with this content:

```julia
using Test
using Suslin
using Oscar

function _issue9_column_matrix(R, values)
    M = zero_matrix(R, length(values), 1)
    for (i, value) in enumerate(values)
        M[i, 1] = value
    end
    return M
end

function _issue9_error_message(err)
    return sprint(showerror, err)
end

@testset "Laurent linear solve" begin
    R, (x, y) = suslin_laurent_polynomial_ring(GF(2), ["x", "y"])

    A = matrix(R, [
        one(R) x;
        zero(R) one(R)
    ])

    single_solution = _issue9_column_matrix(R, [one(R), y^-1])
    single_rhs = A * single_solution
    computed_single = solve_laurent_linear(A, single_rhs)
    @test computed_single == single_solution
    @test A * computed_single == single_rhs

    multi_solution = matrix(R, [
        one(R)      y;
        x^-1        one(R) + y^-1
    ])
    multi_rhs = A * multi_solution
    computed_multi = solve_laurent_linear(A, multi_rhs)
    @test computed_multi == multi_solution
    @test A * computed_multi == multi_rhs

    unsolvable_A = zero_matrix(R, 1, 1)
    unsolvable_B = _issue9_column_matrix(R, [one(R)])
    try
        solve_laurent_linear(unsolvable_A, unsolvable_B)
        error("expected no-solution failure")
    catch err
        @test err isa ErrorException
        @test _issue9_error_message(err) == "No exact solution exists for A * U = B"
    end

    S, (u, v) = suslin_laurent_polynomial_ring(GF(2), ["u", "v"])
    wrong_parent_rhs = _issue9_column_matrix(S, [one(S), v])
    @test_throws ArgumentError solve_laurent_linear(A, wrong_parent_rhs)

    wrong_rows_rhs = _issue9_column_matrix(R, [one(R)])
    @test_throws DimensionMismatch solve_laurent_linear(A, wrong_rows_rhs)
end
```

- [ ] **Step 2: Wire the focused test into the internal group**

In `test/runtests.jl`, add `"internal/laurent_linear_solve.jl"` to the `"internal"` group after `"internal/laurent_normalization.jl"`:

```julia
    "internal" => [
        "internal/rings.jl",
        "internal/laurent_rings.jl",
        "internal/laurent_fixtures.jl",
        "internal/laurent_normalization.jl",
        "internal/laurent_linear_solve.jl",
        "internal/gl_laurent_normalization.jl",
        "internal/toricbuilder_contract.jl",
    ],
```

- [ ] **Step 3: Add public API expectations**

In `test/public/api_surface.jl`, add `solve_laurent_linear` to the `isdefined` and export-identity checks:

```julia
    @test isdefined(Suslin, :solve_laurent_linear)
```

and:

```julia
    @test Suslin.solve_laurent_linear === solve_laurent_linear
```

- [ ] **Step 4: Run the focused test and confirm the RED failure**

Run:

```bash
julia --project=. -e 'include("test/internal/laurent_linear_solve.jl")'
```

Expected: failure because `solve_laurent_linear` is not defined.

- [ ] **Step 5: Commit the failing tests**

```bash
git add test/internal/laurent_linear_solve.jl test/runtests.jl test/public/api_surface.jl
git commit -m "test: cover laurent linear solving contract"
```

---

### Task 2: Implement the Suslin-Owned Solver

**Files:**
- Create: `src/core/laurent_linear_solve.jl`
- Modify: `src/Suslin.jl`

**Interfaces:**
- Consumes: tests from Task 1, `_require_laurent_polynomial_ring`.
- Produces: exported `solve_laurent_linear(A, B)` with the Issue #9 contract.

- [ ] **Step 1: Add the solver implementation**

Create `src/core/laurent_linear_solve.jl` with this content:

```julia
function _validate_laurent_linear_inputs(A, B)
    R = _require_laurent_polynomial_ring(base_ring(A); label="A base ring")

    base_ring(B) === R || throw(ArgumentError("Matrices A and B must be over the same Laurent polynomial ring"))

    rows_A = nrows(A)
    rows_B = nrows(B)
    rows_B == rows_A || throw(DimensionMismatch("Number of rows in A ($(rows_A)) must match number of rows in B ($(rows_B))"))

    return R
end

function _solve_laurent_linear_native(A, B)
    solvable, solution = can_solve_with_solution(A, B; side=:right)
    solvable && return solution
    throw(ErrorException("No exact solution exists for A * U = B"))
end

function _is_native_laurent_solver_capability_error(err, backtrace)
    err isa MethodError || return false

    missing_name = string(err.f)
    (missing_name == "gcdxx" || missing_name == "annihilator") || return false

    frames = stacktrace(backtrace)
    return any(frame -> startswith(String(frame.func), "can_solve_with_solution"), frames)
end

function _laurent_linear_module_data(A, quotient_map)
    A_quo = map_entries(quotient_map, A)
    quotient_ring = base_ring(A_quo)
    polynomial_ring = base_ring(quotient_ring)
    quotient_ideal = modulus(quotient_ring)
    A_poly = map_entries(x -> polynomial_ring(Oscar.lift(x)), A_quo)

    free = free_module(polynomial_ring, nrows(A_poly))
    coefficient_count = ncols(A_poly)
    generators = [free(collect(A_poly[:, j])) for j in 1:coefficient_count]

    row_count = nrows(A_poly)
    for generator in gens(quotient_ideal)
        for row in 1:row_count
            relation = [zero(polynomial_ring) for _ in 1:row_count]
            relation[row] = generator
            push!(generators, free(relation))
        end
    end

    submodule, _ = sub(free, generators)
    return (;
        free,
        submodule,
        coefficient_count,
        polynomial_ring,
        quotient_ring,
    )
end

function _solve_laurent_linear_column(b_col, module_data, quotient_map, R)
    (; free, submodule, coefficient_count, polynomial_ring, quotient_ring) = module_data

    b_quo = map_entries(quotient_map, b_col)
    b_poly = map_entries(x -> polynomial_ring(Oscar.lift(x)), b_quo)
    b_vec = free(collect(b_poly[:, 1]))
    coords = Oscar.coordinates(b_vec, submodule)

    entries = [
        preimage(quotient_map, quotient_ring(coords[i]))
        for i in 1:coefficient_count
    ]
    return matrix(R, length(entries), 1, entries)
end

function _solve_laurent_linear_fallback(A, B)
    R = _validate_laurent_linear_inputs(A, B)
    quotient_map = Oscar._polyringquo(R)
    rhs_cols = ncols(B)
    module_data = _laurent_linear_module_data(A, quotient_map)

    try
        rhs_cols == 1 && return _solve_laurent_linear_column(B, module_data, quotient_map, R)

        solutions = [
            _solve_laurent_linear_column(B[:, j:j], module_data, quotient_map, R)
            for j in 1:rhs_cols
        ]
        return hcat(solutions...)
    catch err
        if err isa ErrorException && occursin("not liftable to the given generating system", sprint(showerror, err))
            throw(ErrorException("No exact solution exists for A * U = B"))
        end
        rethrow()
    end
end

"""
    solve_laurent_linear(A, B)

Solve the Laurent-polynomial linear system `A * U = B` exactly.

The implementation prefers Oscar's native right-side linear solver. If the
installed Oscar/AbstractAlgebra stack does not support Laurent-polynomial
matrices for that path, it falls back to a module-based exact solver while
preserving the same public API and no-solution error contract.

# Arguments
- `A`: `m x n` matrix over a Laurent polynomial ring
- `B`: `m x k` matrix over the same Laurent polynomial ring

# Returns
- `U`: `n x k` solution matrix such that `A * U == B`

# Throws
- `ErrorException` if no exact solution exists
- `ArgumentError` if the inputs are not over the same Laurent polynomial ring
- `DimensionMismatch` if `nrows(B) != nrows(A)`
"""
function solve_laurent_linear(A, B)
    _validate_laurent_linear_inputs(A, B)

    try
        return _solve_laurent_linear_native(A, B)
    catch err
        if _is_native_laurent_solver_capability_error(err, catch_backtrace())
            return _solve_laurent_linear_fallback(A, B)
        end
        rethrow()
    end
end
```

- [ ] **Step 2: Export and include the solver**

In `src/Suslin.jl`, add:

```julia
export solve_laurent_linear
```

near the other public Laurent exports, and add:

```julia
include("core/laurent_linear_solve.jl")
```

after `include("core/groebner_tools.jl")` and before later Laurent consumers.

- [ ] **Step 3: Run the focused test and confirm GREEN**

Run:

```bash
julia --project=. -e 'include("test/internal/laurent_linear_solve.jl")'
```

Expected: the `Laurent linear solve` test set passes.

- [ ] **Step 4: Run package tests**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: public and internal tests pass.

- [ ] **Step 5: Commit the implementation**

```bash
git add src/core/laurent_linear_solve.jl src/Suslin.jl
git commit -m "feat: add laurent linear solver"
```

---

### Task 3: Verify Full Suite and Review Scope

**Files:**
- No source files should be modified in this task unless verification exposes a defect.

**Interfaces:**
- Consumes: committed tests and implementation from Tasks 1 and 2.
- Produces: final verification evidence for PR creation.

- [ ] **Step 1: Run the focused Issue #9 test**

Run:

```bash
julia --project=. -e 'include("test/internal/laurent_linear_solve.jl")'
```

Expected: focused tests pass.

- [ ] **Step 2: Run the package test entry point**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: package tests pass.

- [ ] **Step 3: Run the documented full suite from #21**

Run:

```bash
julia --project=. test/runtests.jl all
```

Expected: public, internal, and expert groups pass.

- [ ] **Step 4: Review the final diff**

Run:

```bash
git diff --stat origin/main...HEAD
git diff --name-only origin/main...HEAD
```

Expected changed files are limited to the Issue #9 spec, plan, solver source,
exports/includes, and tests.
