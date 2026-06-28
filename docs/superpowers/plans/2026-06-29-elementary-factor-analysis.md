# Elementary Factor Analysis Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [x]`) syntax for tracking.

**Goal:** Add two public helpers that analyze already-computed elementary factor sequences.

**Architecture:** Keep the feature in `src/core/elementary_matrices.jl`, beside elementary matrix construction and sequence helpers. Export the public API from `src/Suslin.jl`; tests live with the existing elementary-matrix expert tests and public API surface tests.

**Tech Stack:** Julia, Oscar.jl, AbstractAlgebra polynomial/Laurent term metadata, Julia `Test`.

## Global Constraints

- Do not recompute a factorization from an input matrix.
- Accept factor sequences such as `elementary_factorization(A)`, `.factors`, and `.core_factors`.
- Scan only non-diagonal entries.
- Count `x^5*y^-7` as exponent weight `12`.
- Empty factor sequences return `0` for both helpers.
- Do not add an aggregate stats function.
- Do not change existing factorization or certificate logic.

---

## File Structure

- Modify `src/Suslin.jl`: export `max_elementary_factor_monomial_degree` and `total_elementary_factor_offdiagonal_monomials`.
- Modify `src/core/elementary_matrices.jl`: add private term helpers and the two public factor-sequence analysis functions.
- Modify `test/expert/elementary_matrices.jl`: add behavior tests for Laurent exponent weights, term counts, diagonal ignoring, and empty factor sequences.
- Modify `test/public/api_surface.jl`: assert the new helpers are exported and bound to `Suslin`.

### Task 1: Add Failing Tests

**Files:**
- Modify: `test/expert/elementary_matrices.jl`
- Modify: `test/public/api_surface.jl`

**Interfaces:**
- Consumes: desired public functions `max_elementary_factor_monomial_degree(factors)` and `total_elementary_factor_offdiagonal_monomials(factors)`.
- Produces: failing tests that define the behavior and export surface.

- [x] **Step 1: Add behavior tests**

Append this testset to `test/expert/elementary_matrices.jl`:

```julia
@testset "elementary factor sequence analysis" begin
    L, (x, y) = suslin_laurent_polynomial_ring(QQ, ["x", "y"])
    laurent_factor = elementary_matrix(3, 1, 2, x^5 * y^-7 + x^-2 + one(L), L)
    diagonal_noise = identity_matrix(L, 3)
    diagonal_noise[1, 1] = x^99 * y^-99

    @test max_elementary_factor_monomial_degree([laurent_factor, diagonal_noise]) == 12
    @test total_elementary_factor_offdiagonal_monomials([laurent_factor, diagonal_noise]) == 3

    R, (u, v) = suslin_polynomial_ring(QQ, ["u", "v"])
    polynomial_factor = elementary_matrix(3, 2, 3, u^2 * v + u + 3, R)
    constant_factor = elementary_matrix(3, 3, 1, R(7), R)
    polynomial_diagonal_noise = identity_matrix(R, 3)
    polynomial_diagonal_noise[2, 2] = u^100

    @test max_elementary_factor_monomial_degree([
        polynomial_factor,
        constant_factor,
        polynomial_diagonal_noise,
    ]) == 3
    @test total_elementary_factor_offdiagonal_monomials([
        polynomial_factor,
        constant_factor,
        polynomial_diagonal_noise,
    ]) == 4
    @test max_elementary_factor_monomial_degree(Matrix{Any}[]) == 0
    @test total_elementary_factor_offdiagonal_monomials(Matrix{Any}[]) == 0
end
```

- [x] **Step 2: Add API surface tests**

Add these assertions to `test/public/api_surface.jl` inside the existing API testset, near the other elementary exports:

```julia
@test isdefined(Suslin, :max_elementary_factor_monomial_degree)
@test isdefined(Suslin, :total_elementary_factor_offdiagonal_monomials)
```

Add these identity assertions near the other `Suslin.foo === foo` assertions:

```julia
@test Suslin.max_elementary_factor_monomial_degree === max_elementary_factor_monomial_degree
@test Suslin.total_elementary_factor_offdiagonal_monomials === total_elementary_factor_offdiagonal_monomials
```

- [x] **Step 3: Run tests to verify RED**

Run:

```bash
julia --project=. -e 'include("test/expert/elementary_matrices.jl")'
julia --project=. -e 'include("test/public/api_surface.jl")'
```

Expected: both commands fail because the new public functions are not defined/exported.

### Task 2: Implement Factor-Sequence Analysis Helpers

**Files:**
- Modify: `src/Suslin.jl`
- Modify: `src/core/elementary_matrices.jl`

**Interfaces:**
- Consumes: the tests from Task 1.
- Produces:
  - `max_elementary_factor_monomial_degree(factors)::Int`
  - `total_elementary_factor_offdiagonal_monomials(factors)::Int`

- [x] **Step 1: Export public helpers**

Add these exports to `src/Suslin.jl` beside the other elementary-matrix exports:

```julia
export max_elementary_factor_monomial_degree
export total_elementary_factor_offdiagonal_monomials
```

- [x] **Step 2: Add private term helpers and public functions**

Add this implementation to `src/core/elementary_matrices.jl` after `elementary_matrix`:

```julia
function _elementary_factor_term_exponents(value)
    try
        return collect(exponents(value))
    catch err
        err isa InterruptException && rethrow()
        err isa MethodError || err isa ErrorException || rethrow()
    end

    try
        return collect(AbstractAlgebra.exponent_vectors(value))
    catch err
        err isa InterruptException && rethrow()
        err isa MethodError || err isa ErrorException || rethrow()
        throw(ArgumentError("elementary factor analysis requires polynomial or Laurent polynomial entries"))
    end
end

function _elementary_factor_term_count(value)::Int
    iszero(value) && return 0
    return length(_elementary_factor_term_exponents(value))
end

function _elementary_factor_monomial_degree(raw_exponents)::Int
    return sum(abs, Int.(collect(raw_exponents)))
end

function _require_square_elementary_analysis_factor(factor)
    nrows(factor) == ncols(factor) || throw(ArgumentError("factor must be square"))
    return nrows(factor)
end

function max_elementary_factor_monomial_degree(factors)::Int
    max_degree = 0
    for factor in factors
        n = _require_square_elementary_analysis_factor(factor)
        for row in 1:n, col in 1:n
            row == col && continue
            value = factor[row, col]
            iszero(value) && continue
            for raw_exponents in _elementary_factor_term_exponents(value)
                max_degree = max(max_degree, _elementary_factor_monomial_degree(raw_exponents))
            end
        end
    end
    return max_degree
end

function total_elementary_factor_offdiagonal_monomials(factors)::Int
    total = 0
    for factor in factors
        n = _require_square_elementary_analysis_factor(factor)
        for row in 1:n, col in 1:n
            row == col && continue
            total += _elementary_factor_term_count(factor[row, col])
        end
    end
    return total
end
```

- [x] **Step 3: Run tests to verify GREEN**

Run:

```bash
julia --project=. -e 'include("test/expert/elementary_matrices.jl")'
julia --project=. -e 'include("test/public/api_surface.jl")'
```

Expected: both commands exit 0.

### Task 3: Verify Full Package and Commit

**Files:**
- Verify: `src/Suslin.jl`
- Verify: `src/core/elementary_matrices.jl`
- Verify: `test/expert/elementary_matrices.jl`
- Verify: `test/public/api_surface.jl`

**Interfaces:**
- Consumes: implemented public helpers from Task 2.
- Produces: committed implementation ready to push and open as a PR.

- [x] **Step 1: Run package tests**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: command exits 0.

- [x] **Step 2: Inspect diff**

Run:

```bash
git diff -- src/Suslin.jl src/core/elementary_matrices.jl test/expert/elementary_matrices.jl test/public/api_surface.jl
git status --short
```

Expected: only the implementation and test files are modified, plus this plan file if it has not yet been committed.

- [x] **Step 3: Commit implementation**

Run:

```bash
git add docs/superpowers/plans/2026-06-29-elementary-factor-analysis.md src/Suslin.jl src/core/elementary_matrices.jl test/expert/elementary_matrices.jl test/public/api_surface.jl
git commit -m "feat: add elementary factor analysis helpers"
```

Expected: commit succeeds.
