# Issue 12 Laurent Elementary Core Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make elementary, Cohn-type, conjugated elementary, and factorization verification routines work exactly over Laurent polynomial rings.

**Architecture:** Keep the public API unchanged and lift the existing generic algebra. Add focused Laurent expert coverage, route Laurent normality through a small adjugate-based inverse helper because Oscar's generic matrix `inv` is not implemented for Laurent matrix parents, and reuse the shared base-ring comparison helper for exact Laurent parent checks.

**Tech Stack:** Julia, Oscar, AbstractAlgebra, Test.

## Global Constraints

- Do not add a parallel Laurent-only public API.
- Do not enable the full Laurent `elementary_factorization` driver.
- Keep determinant normalization out of this issue.
- Preserve existing ordinary polynomial behavior.
- Focused issue command: `julia --project=. -e 'include("test/expert/laurent_elementary_core.jl")'`.
- Package verification command: `julia --project=. -e 'using Pkg; Pkg.test()'`.
- Full-suite verification command from issue #21: `julia --project=. test/runtests.jl all`.

---

## File Structure

- Create `test/expert/laurent_elementary_core.jl`: issue-specific expert tests for Laurent elementary construction, Cohn-type products, normality products, `verify_factorization`, and the mixed-parent negative control.
- Modify `test/runtests.jl`: register the new expert test file.
- Modify `src/algorithm/normality.jl`: add the Laurent inverse helper and use it in `realize_conjugate_elementary`.
- Modify `src/core/elementary_matrices.jl`: make `_same_base_ring` use identity comparison for Laurent parents.
- Modify `src/algorithm/factorization.jl`: use `_same_base_ring` inside `verify_factorization`.

### Task 1: Laurent Elementary Core Tests And Implementation

**Files:**
- Create: `test/expert/laurent_elementary_core.jl`
- Modify: `test/runtests.jl`
- Modify: `src/algorithm/normality.jl`
- Modify: `src/core/elementary_matrices.jl`
- Modify: `src/algorithm/factorization.jl`

**Interfaces:**
- Consumes: `suslin_laurent_polynomial_ring`, `elementary_matrix`, `realize_cohn_type`, `realize_conjugate_elementary`, and `verify_factorization`.
- Produces: no new public API. Adds internal helpers `_inverse_matrix_over_base_ring(M)`, `_laurent_inverse_matrix(M)`, `_adjugate_matrix(M)`, and `_matrix_without_row_col(M, skipped_row, skipped_col)`.

- [ ] **Step 1: Write the failing Laurent expert test**

Create `test/expert/laurent_elementary_core.jl` with:

```julia
using Test
using Suslin
using Oscar

const LAURENT_FIXTURE_CATALOG_PATH = joinpath(@__DIR__, "..", "fixtures", "laurent_cases.jl")
include(LAURENT_FIXTURE_CATALOG_PATH)

function laurent_product_of_factors(factors)
    isempty(factors) && throw(ArgumentError("factor list must be nonempty"))

    R = base_ring(first(factors))
    product = identity_matrix(R, size(first(factors), 1))
    for factor in factors
        product *= factor
    end
    return product
end

function laurent_cohn_type_target(n::Int, i::Int, j::Int, a, v::AbstractVector, R)
    target = identity_matrix(R, n)
    for row in 1:n
        target[row, i] += a * v[row] * v[j]
        target[row, j] -= a * v[row] * v[i]
    end
    return target
end

@testset "Laurent elementary core" begin
    catalog = LaurentFixtureCatalog.catalog()
    fixture_ring = catalog.ring.object
    fx, fy = catalog.ring.generators

    fixture_entry = only(filter(entry -> entry.id == "laurent-negative-exponent-normalization", catalog.cases))
    elementary_entry = fixture_entry.inputs.vector[1, 1]
    E = elementary_matrix(3, 1, 2, elementary_entry, fixture_ring)

    @test base_ring(E) === fixture_ring
    @test E[1, 2] == elementary_entry
    @test E[1, 1] == one(fixture_ring)
    @test E[2, 1] == zero(fixture_ring)
    @test parent(fx^-1 * fy) === fixture_ring

    R, (x, y) = suslin_laurent_polynomial_ring(QQ, ["x", "y"])
    a = x * y^-1 - y

    v3 = [one(R), x^-1, y]
    target3 = laurent_cohn_type_target(3, 1, 2, a, v3, R)
    factors3 = Suslin.realize_cohn_type(3, 1, 2, a, v3, R)

    @test target3 == laurent_product_of_factors(factors3)
    @test verify_factorization(target3, factors3)
    @test all(factor -> base_ring(factor) === R, factors3)

    v4 = [one(R), x^-1, y, x * y^-1]
    target4 = laurent_cohn_type_target(4, 2, 4, a, v4, R)
    factors4 = Suslin.realize_cohn_type(4, 2, 4, a, v4, R)

    @test target4 == laurent_product_of_factors(factors4)
    @test verify_factorization(target4, factors4)
    @test all(factor -> base_ring(factor) === R, factors4)

    B = matrix(R, [
        1      0  0;
        x^-1   1  0;
        y      x  1
    ])
    Binv = matrix(R, [
        1          0   0;
        -x^-1      1   0;
        1 - y     -x   1
    ])
    normality_entry = x + y^-1
    normality_elementary = elementary_matrix(3, 1, 3, normality_entry, R)
    normality_target = B * normality_elementary * Binv
    normality_factors = Suslin.realize_conjugate_elementary(B, 1, 3, normality_entry)

    @test B * Binv == identity_matrix(R, 3)
    @test Binv * B == identity_matrix(R, 3)
    @test normality_target == laurent_product_of_factors(normality_factors)
    @test verify_factorization(normality_target, normality_factors)
    @test all(factor -> base_ring(factor) === R, normality_factors)

    S, (u, _) = suslin_laurent_polynomial_ring(QQ, ["u", "v"])
    wrong_parent_factor = elementary_matrix(3, 1, 2, u, S)

    @test_throws ArgumentError verify_factorization(identity_matrix(R, 3), [wrong_parent_factor])
end
```

Append the new expert file to `test/runtests.jl` immediately after `expert/normality.jl`:

```julia
        "expert/normality.jl",
        "expert/laurent_elementary_core.jl",
```

- [ ] **Step 2: Run the focused test to verify RED**

Run:

```bash
julia --project=. -e 'include("test/expert/laurent_elementary_core.jl")'
```

Expected: FAIL with an `ErrorException("Not implemented")` stack trace from the Laurent `inv(B)` call inside `realize_conjugate_elementary`.

- [ ] **Step 3: Implement Laurent inverse support and exact parent comparison**

In `src/core/elementary_matrices.jl`, replace `_same_base_ring` with:

```julia
function _same_base_ring(left, right)::Bool
    if _is_laurent_polynomial_ring(left) || _is_laurent_polynomial_ring(right)
        return left === right
    end

    return left == right || left === right
end
```

In `src/algorithm/factorization.jl`, replace the direct base-ring comparison in `verify_factorization` with:

```julia
        _same_base_ring(base_ring(factor), R) || throw(ArgumentError("all factors must lie in the same base ring as A"))
```

In `src/algorithm/normality.jl`, replace:

```julia
    Binv = inv(B)
```

with:

```julia
    Binv = _inverse_matrix_over_base_ring(B)
```

Then add these helpers after `realize_conjugate_elementary` and before `_dot`:

```julia
function _inverse_matrix_over_base_ring(M)
    R = base_ring(M)
    _is_laurent_polynomial_ring(R) && return _laurent_inverse_matrix(M)

    return inv(M)
end

function _laurent_inverse_matrix(M)
    n = _require_square_matrix(M, "B")
    R = _require_laurent_polynomial_ring(base_ring(M); label="B base ring")
    determinant = det(M)
    determinant == zero(R) && throw(ArgumentError("B must be invertible over its Laurent polynomial ring"))

    determinant_inverse = try
        inv(determinant)
    catch err
        if err isa ArgumentError || err isa ErrorException || err isa MethodError
            throw(ArgumentError("B must have a Laurent-unit determinant"))
        end
        rethrow()
    end

    inverse = _adjugate_matrix(M)
    for row in 1:n, col in 1:n
        inverse[row, col] *= determinant_inverse
    end

    identity = identity_matrix(R, n)
    inverse * M == identity || throw(ArgumentError("failed to invert B over its Laurent polynomial ring"))
    M * inverse == identity || throw(ArgumentError("failed to invert B over its Laurent polynomial ring"))
    return inverse
end

function _adjugate_matrix(M)
    n = _require_square_matrix(M, "matrix")
    R = base_ring(M)

    if n == 1
        adjugate = zero_matrix(R, 1, 1)
        adjugate[1, 1] = one(R)
        return adjugate
    end

    adjugate = zero_matrix(R, n, n)
    for row in 1:n, col in 1:n
        sign = isodd(row + col) ? -one(R) : one(R)
        adjugate[col, row] = sign * det(_matrix_without_row_col(M, row, col))
    end
    return adjugate
end

function _matrix_without_row_col(M, skipped_row::Int, skipped_col::Int)
    n = _require_square_matrix(M, "matrix")
    R = base_ring(M)
    minor = zero_matrix(R, n - 1, n - 1)
    minor_row = 1

    for row in 1:n
        row == skipped_row && continue
        minor_col = 1
        for col in 1:n
            col == skipped_col && continue
            minor[minor_row, minor_col] = M[row, col]
            minor_col += 1
        end
        minor_row += 1
    end

    return minor
end
```

- [ ] **Step 4: Run focused GREEN verification**

Run:

```bash
julia --project=. -e 'include("test/expert/laurent_elementary_core.jl")'
```

Expected: PASS for the `Laurent elementary core` testset.

- [ ] **Step 5: Run regression checks**

Run:

```bash
julia --project=. test/runtests.jl expert
julia --project=. -e 'using Pkg; Pkg.test()'
julia --project=. test/runtests.jl all
```

Expected: all commands exit 0. The package command should run the default public and internal groups. The full-suite command should run public, internal, and expert groups including the new Laurent expert file.

- [ ] **Step 6: Commit implementation**

Run:

```bash
git add src/core/elementary_matrices.jl src/algorithm/factorization.jl src/algorithm/normality.jl test/expert/laurent_elementary_core.jl test/runtests.jl
git commit -m "feat: lift elementary core to Laurent rings"
```

## Self-Review

- Spec coverage: Task 1 covers the requested Laurent elementary matrix, Cohn-type realization for `n = 3` and `n = 4`, conjugated elementary realization, `verify_factorization` success cases, and mixed-parent negative control.
- Placeholder scan: no placeholders or deferred implementation steps remain.
- Type consistency: all helper names used in the task are produced in the same task and remain internal to existing source files.
