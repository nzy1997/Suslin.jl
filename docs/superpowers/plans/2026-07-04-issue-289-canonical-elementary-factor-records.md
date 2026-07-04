# Issue 289 Canonical Elementary Factor Records Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add private canonical records and a round-trip constructor for ordinary-polynomial elementary matrix factors used by Steinberg rewrites.

**Architecture:** Keep canonicalization in `src/core/elementary_matrices.jl` beside `elementary_matrix` and shape/ring helpers. Add a focused expert test file for the internal helpers and register it with the expert group; do not export new public API or alter factorization routes.

**Tech Stack:** Julia, Oscar matrices and polynomial rings, Suslin elementary matrix helpers, Julia `Test`.

## Global Constraints

- Keep the helpers internal at this stage.
- Reuse `elementary_matrix`, `_same_base_ring`, and existing shape checks.
- Accept identity matrices as identity records.
- Reject matrices with diagonal entries other than one.
- Reject matrices with more than one nonzero off-diagonal entry.
- Do not optimize sequences.
- Do not export a public API.
- Do not change current factorization routes.

---

## File Structure

- Modify `src/core/elementary_matrices.jl`: add `_canonical_elementary_factor_record(factor)` and `_elementary_factor_record_matrix(record)`.
- Create `test/expert/steinberg_factor_count_optimization.jl`: focused tests for accepted records, round trips, and negative controls.
- Modify `test/runtests.jl`: register the focused expert test.

### Task 1: Add Canonical Record Tests

**Files:**
- Create: `test/expert/steinberg_factor_count_optimization.jl`
- Modify: `test/runtests.jl`

**Interfaces:**
- Consumes: private functions `Suslin._canonical_elementary_factor_record(factor)` and `Suslin._elementary_factor_record_matrix(record)`.
- Produces: failing tests that define canonical identity and elementary records plus malformed-factor rejection.

- [x] **Step 1: Write the failing test file**

Create `test/expert/steinberg_factor_count_optimization.jl` with:

```julia
using Test
using Suslin
using Oscar

@testset "Steinberg canonical elementary factor records" begin
    R, (x, y) = Oscar.polynomial_ring(QQ, ["x", "y"])
    coefficient = x + y + one(R)

    elementary_factor = elementary_matrix(3, 1, 2, coefficient, R)
    elementary_record = Suslin._canonical_elementary_factor_record(elementary_factor)

    @test elementary_record.kind == :elementary
    @test elementary_record.n == 3
    @test Suslin._same_base_ring(elementary_record.ring, R)
    @test elementary_record.row == 1
    @test elementary_record.col == 2
    @test elementary_record.coefficient == coefficient
    @test Suslin._elementary_factor_record_matrix(elementary_record) == elementary_factor

    zero_elementary_factor = elementary_matrix(3, 1, 2, zero(R), R)
    identity_record = Suslin._canonical_elementary_factor_record(zero_elementary_factor)

    @test identity_record.kind == :identity
    @test identity_record.n == 3
    @test Suslin._same_base_ring(identity_record.ring, R)
    @test Suslin._elementary_factor_record_matrix(identity_record) == zero_elementary_factor

    nonsquare_factor = zero_matrix(R, 2, 3)
    bad_diagonal = identity_matrix(R, 3)
    bad_diagonal[2, 2] = x
    two_offdiagonal = identity_matrix(R, 3)
    two_offdiagonal[1, 2] = x
    two_offdiagonal[2, 3] = y

    @test_throws DimensionMismatch Suslin._canonical_elementary_factor_record(nonsquare_factor)
    @test_throws ArgumentError Suslin._canonical_elementary_factor_record(bad_diagonal)
    @test_throws ArgumentError Suslin._canonical_elementary_factor_record(two_offdiagonal)
end
```

- [x] **Step 2: Register the expert test**

Add `"expert/steinberg_factor_count_optimization.jl"` to the expert file list in `test/runtests.jl` immediately after `"expert/elementary_matrices.jl"`.

- [x] **Step 3: Run RED verification**

Run:

```bash
julia --project=. -e 'include("test/expert/steinberg_factor_count_optimization.jl")'
```

Expected: FAIL with `UndefVarError` for `_canonical_elementary_factor_record`.

- [x] **Step 4: Commit the red tests**

Run:

```bash
git add test/expert/steinberg_factor_count_optimization.jl test/runtests.jl
git commit -m "test: add Steinberg elementary factor record coverage"
```

### Task 2: Implement Canonical Record Helpers

**Files:**
- Modify: `src/core/elementary_matrices.jl`

**Interfaces:**
- Consumes: tests from Task 1 and existing helpers `elementary_matrix`, `_require_square_matrix`, and `_same_base_ring`.
- Produces:
  - `_canonical_elementary_factor_record(factor)`
  - `_elementary_factor_record_matrix(record)`

- [x] **Step 1: Add the private implementation**

Add this implementation near the existing elementary factor analysis helpers in `src/core/elementary_matrices.jl`:

```julia
function _canonical_elementary_factor_record(factor)
    n = _require_square_matrix(factor, "elementary factor")
    R = base_ring(factor)
    row = 0
    col = 0
    coefficient = zero(R)

    for i in 1:n, j in 1:n
        entry = factor[i, j]
        if i == j
            entry == one(R) || throw(ArgumentError("elementary factor diagonal must be one"))
            continue
        end

        entry == zero(R) && continue
        row == 0 || throw(ArgumentError("elementary factor must have at most one nonzero off-diagonal entry"))
        row = i
        col = j
        coefficient = entry
    end

    row == 0 && return (; kind = :identity, n, ring = R)
    return (; kind = :elementary, n, ring = R, row, col, coefficient)
end

function _elementary_factor_record_matrix(record)
    record.kind == :identity && return identity_matrix(record.ring, record.n)
    record.kind == :elementary &&
        return elementary_matrix(record.n, record.row, record.col, record.coefficient, record.ring)
    throw(ArgumentError("unknown elementary factor record kind"))
end
```

- [x] **Step 2: Run GREEN focused verification**

Run:

```bash
julia --project=. -e 'include("test/expert/steinberg_factor_count_optimization.jl")'
```

Expected: PASS.

- [x] **Step 3: Run expert registration verification**

Run:

```bash
julia --project=. test/runtests.jl expert
```

Expected: PASS, including the new expert file.

- [x] **Step 4: Commit the implementation**

Run:

```bash
git add src/core/elementary_matrices.jl
git commit -m "feat: canonicalize elementary factor records"
```

### Task 3: Final Verification And Review

**Files:**
- Read/review: branch diff and test output.

**Interfaces:**
- Consumes: Task 1 and Task 2 commits.
- Produces: verified branch ready for PR.

- [x] **Step 1: Run issue-required focused command**

Run:

```bash
julia --project=. -e 'include("test/expert/steinberg_factor_count_optimization.jl")'
```

Expected: PASS.

- [x] **Step 2: Run required package suite**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: PASS for the default public/internal suite.

- [x] **Step 3: Review branch diff**

Run:

```bash
git diff --stat origin/main..HEAD
git diff origin/main..HEAD
```

Expected: only docs, `src/core/elementary_matrices.jl`, `test/expert/steinberg_factor_count_optimization.jl`, and `test/runtests.jl` contain intentional changes.

- [x] **Step 4: Commit any verification-only plan updates**

If the plan checkboxes were updated, run:

```bash
git add docs/superpowers/plans/2026-07-04-issue-289-canonical-elementary-factor-records.md
git commit -m "docs: update issue 289 implementation plan"
```

## Self-Review

The plan covers all requested files, keeps helpers internal, includes TDD red and green steps, covers identity and elementary round trips, covers all three negative controls, and runs both the issue-required focused command and the required package suite.
