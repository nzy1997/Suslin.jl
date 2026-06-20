# Issue 43 Side-Aware Elementary Preconditioning Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add exact side-aware elementary row/column preconditioning primitives with ordered replay and verification.

**Architecture:** Keep the helpers in `src/core/elementary_matrices.jl` beside `elementary_matrix`, because they are matrix primitives rather than reduction algorithms. Represent operations as named tuples with explicit `side`, `factor`, and `transformed_matrix` fields, and make verification replay the recorded factor sides in order.

**Tech Stack:** Julia, Oscar exact matrices, Suslin polynomial/Laurent ring helpers, `Test`.

## Global Constraints

- Work in the existing linked worktree on branch `agent/issue-43-add-side-aware-elementary-preconditioning-primit-run-1`; do not create another worktree.
- Follow TDD: write the failing helper tests first, run them and confirm the expected failure, then write production code.
- Keep the API small: no optimizer, no sequence search, no determinant normalization, and no certificate type.
- Public helpers are `elementary_preconditioning_step(A, side, target, source, coefficient)`, `replay_elementary_preconditioning(A, steps)`, and `verify_elementary_preconditioning(A, steps, expected)`.
- For `:left`, target/source mean `row[target] += coefficient * row[source]` and the factor has entry `E[target, source] = coefficient`.
- For `:right`, target/source mean `column[target] += coefficient * column[source]` and the factor has entry `E[source, target] = coefficient`.
- Left factors must match `nrows(current_matrix)`; right factors must match `ncols(current_matrix)`.
- Coefficients must be coerced into `base_ring(A)` with the existing `_coerce_into_ring` behavior.
- Verification must return `false` rather than throw for invalid replay metadata, wrong dimensions, wrong base ring, or failed exact reconstruction.
- The issue-specific verification command is `julia --project=. -e 'include("test/expert/elementary_preconditioning.jl")'`.
- The Agent Desk package verification command is `julia --project=. -e 'using Pkg; Pkg.test()'`.
- The documented full-suite command is `julia --project=. test/runtests.jl all`.

---

## File Structure

- Modify `src/Suslin.jl` to export the three side-aware preconditioning helpers.
- Modify `src/core/elementary_matrices.jl` to add side validation, index validation, side-aware factor construction, replay, and verification.
- Create `test/expert/elementary_preconditioning.jl` for the issue-specific Laurent sequence and negative controls.
- Modify `test/runtests.jl` to register `expert/elementary_preconditioning.jl` in the expert group.
- Modify `test/public/api_surface.jl` to assert the new public helpers are exported and bound.

### Task 1: Side-Aware Elementary Preconditioning Helpers

**Files:**
- Modify: `src/Suslin.jl`
- Modify: `src/core/elementary_matrices.jl`
- Create: `test/expert/elementary_preconditioning.jl`
- Modify: `test/runtests.jl`
- Modify: `test/public/api_surface.jl`

**Interfaces:**
- Consumes: `elementary_matrix(n::Int, i::Int, j::Int, a, R)`, `_coerce_into_ring`, `_same_base_ring`, Oscar `identity_matrix`, `base_ring`, `nrows`, and `ncols`.
- Produces: `elementary_preconditioning_step(A, side, target, source, coefficient)`, `replay_elementary_preconditioning(A, steps)`, and `verify_elementary_preconditioning(A, steps, expected)`.

- [ ] **Step 1: Write the failing expert tests and public API assertions**

Create `test/expert/elementary_preconditioning.jl` with this content:

```julia
using Test
using Suslin
using Oscar

@testset "side-aware elementary preconditioning" begin
    R, (x, y) = suslin_laurent_polynomial_ring(GF(2), ["x", "y"])
    A = matrix(R, [
        one(R)      x                 y^-1;
        y           one(R) + x * y    x^-1;
        zero(R)     y^-1              one(R) + x
    ])

    left_step = elementary_preconditioning_step(A, :left, 2, 1, x^-1)
    expected_left_factor = elementary_matrix(3, 2, 1, x^-1, R)
    expected_left_matrix = expected_left_factor * A

    @test left_step.side == :left
    @test left_step.target == 2
    @test left_step.source == 1
    @test left_step.coefficient == x^-1
    @test left_step.factor == expected_left_factor
    @test left_step.transformed_matrix == expected_left_matrix

    right_step = elementary_preconditioning_step(left_step.transformed_matrix, :right, 3, 1, y)
    expected_right_factor = elementary_matrix(3, 1, 3, y, R)
    expected_right_matrix = left_step.transformed_matrix * expected_right_factor

    @test right_step.side == :right
    @test right_step.target == 3
    @test right_step.source == 1
    @test right_step.coefficient == y
    @test right_step.factor == expected_right_factor
    @test right_step.transformed_matrix == expected_right_matrix

    final_step = elementary_preconditioning_step(right_step.transformed_matrix, :left, 1, 3, x * y^-1)
    expected_final_factor = elementary_matrix(3, 1, 3, x * y^-1, R)
    expected_final_matrix = expected_final_factor * right_step.transformed_matrix
    steps = [left_step, right_step, final_step]

    @test final_step.factor == expected_final_factor
    @test final_step.transformed_matrix == expected_final_matrix
    @test replay_elementary_preconditioning(A, steps) == expected_final_matrix
    @test verify_elementary_preconditioning(A, steps, expected_final_matrix)

    swapped_steps = [
        (; step..., side = step.side == :left ? :right : :left)
        for step in steps
    ]
    @test !verify_elementary_preconditioning(A, swapped_steps, expected_final_matrix)
end

@testset "elementary preconditioning validation" begin
    R, (x, y) = suslin_laurent_polynomial_ring(GF(2), ["x", "y"])
    S, (u, v) = suslin_laurent_polynomial_ring(GF(2), ["u", "v"])
    A = matrix(R, [
        one(R)  x        zero(R);
        y       one(R)   x^-1;
        zero(R) y^-1     one(R)
    ])

    @test_throws ArgumentError elementary_preconditioning_step(A, :middle, 1, 2, x)
    @test_throws ArgumentError elementary_preconditioning_step(A, :left, 1, 1, x)
    @test_throws ArgumentError elementary_preconditioning_step(A, :left, 0, 1, x)
    @test_throws ArgumentError elementary_preconditioning_step(A, :right, 4, 1, x)
    @test_throws ArgumentError elementary_preconditioning_step(A, :right, 1, 2, u)

    wrong_ring_factor = elementary_matrix(3, 1, 2, u, S)
    wrong_size_factor = elementary_matrix(2, 1, 2, one(R), R)

    @test_throws ArgumentError replay_elementary_preconditioning(A, [(; side = :left, factor = wrong_ring_factor)])
    @test_throws ArgumentError replay_elementary_preconditioning(A, [(; side = :right, factor = wrong_size_factor)])
    @test_throws ArgumentError replay_elementary_preconditioning(A, [(; factor = identity_matrix(R, 3))])

    @test !verify_elementary_preconditioning(A, [(; side = :left, factor = wrong_ring_factor)], A)
    @test !verify_elementary_preconditioning(A, [(; side = :right, factor = wrong_size_factor)], A)
end
```

Modify `test/runtests.jl` by adding the new expert file immediately after
`expert/elementary_matrices.jl`:

```julia
"expert/elementary_preconditioning.jl",
```

Modify `test/public/api_surface.jl` in the public API testset by adding these
`isdefined` checks near the existing elementary helper checks:

```julia
@test isdefined(Suslin, :elementary_preconditioning_step)
@test isdefined(Suslin, :replay_elementary_preconditioning)
@test isdefined(Suslin, :verify_elementary_preconditioning)
```

Then add these binding checks near the existing elementary helper binding
checks:

```julia
@test Suslin.elementary_preconditioning_step === elementary_preconditioning_step
@test Suslin.replay_elementary_preconditioning === replay_elementary_preconditioning
@test Suslin.verify_elementary_preconditioning === verify_elementary_preconditioning
```

- [ ] **Step 2: Run the issue-specific test and verify RED**

Run:

```bash
julia --project=. -e 'include("test/expert/elementary_preconditioning.jl")'
```

Expected: FAIL with `UndefVarError` for `elementary_preconditioning_step`.

- [ ] **Step 3: Export the new helpers**

Modify `src/Suslin.jl` by adding these exports after `export elementary_matrix`:

```julia
export elementary_preconditioning_step
export replay_elementary_preconditioning
export verify_elementary_preconditioning
```

- [ ] **Step 4: Implement the minimal helper code**

Add this code to `src/core/elementary_matrices.jl` after `elementary_matrix`:

```julia
function _require_elementary_preconditioning_side(side)
    side isa Symbol || throw(ArgumentError("side must be :left or :right"))
    (side == :left || side == :right) || throw(ArgumentError("side must be :left or :right"))
    return side
end

function _require_preconditioning_index(index, limit::Int, label::AbstractString)
    index isa Integer || throw(ArgumentError("$label must be an integer"))
    idx = Int(index)
    1 <= idx <= limit || throw(ArgumentError("$label must be between 1 and $limit"))
    return idx
end

function _preconditioning_factor_size(A, side::Symbol)
    return side == :left ? nrows(A) : ncols(A)
end

function _preconditioning_factor_indices(side::Symbol, target::Int, source::Int)
    return side == :left ? (target, source) : (source, target)
end

function elementary_preconditioning_step(A, side, target, source, coefficient)
    checked_side = _require_elementary_preconditioning_side(side)
    R = base_ring(A)
    factor_size = _preconditioning_factor_size(A, checked_side)
    target_idx = _require_preconditioning_index(target, factor_size, "target")
    source_idx = _require_preconditioning_index(source, factor_size, "source")
    target_idx == source_idx && throw(ArgumentError("target and source must differ"))

    coerced_coefficient = _coerce_into_ring(R, coefficient, "coefficient")
    factor_row, factor_col = _preconditioning_factor_indices(checked_side, target_idx, source_idx)
    factor = elementary_matrix(factor_size, factor_row, factor_col, coerced_coefficient, R)
    transformed_matrix = checked_side == :left ? factor * A : A * factor

    return (;
        side = checked_side,
        target = target_idx,
        source = source_idx,
        coefficient = coerced_coefficient,
        factor,
        transformed_matrix,
    )
end

function _require_preconditioning_step_property(step, property::Symbol)
    property in propertynames(step) || throw(ArgumentError("preconditioning step must include $(property)"))
    return getproperty(step, property)
end

function _require_preconditioning_factor(current, side::Symbol, factor)
    expected_size = _preconditioning_factor_size(current, side)
    nrows(factor) == expected_size || throw(ArgumentError("preconditioning factor has wrong row count for $(side) side"))
    ncols(factor) == expected_size || throw(ArgumentError("preconditioning factor has wrong column count for $(side) side"))
    _same_base_ring(base_ring(factor), base_ring(current)) ||
        throw(ArgumentError("preconditioning factor must have the same base ring as the current matrix"))
    return factor
end

function _apply_preconditioning_factor(current, side, factor)
    checked_side = _require_elementary_preconditioning_side(side)
    checked_factor = _require_preconditioning_factor(current, checked_side, factor)
    return checked_side == :left ? checked_factor * current : current * checked_factor
end

function replay_elementary_preconditioning(A, steps)
    current = A
    for step in steps
        side = _require_preconditioning_step_property(step, :side)
        factor = _require_preconditioning_step_property(step, :factor)
        current = _apply_preconditioning_factor(current, side, factor)
    end
    return current
end

function verify_elementary_preconditioning(A, steps, expected)::Bool
    try
        replayed = replay_elementary_preconditioning(A, steps)
        nrows(replayed) == nrows(expected) || return false
        ncols(replayed) == ncols(expected) || return false
        _same_base_ring(base_ring(replayed), base_ring(expected)) || return false
        return replayed == expected
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end
```

- [ ] **Step 5: Run the issue-specific test and verify GREEN**

Run:

```bash
julia --project=. -e 'include("test/expert/elementary_preconditioning.jl")'
```

Expected: PASS. The test summary should report the side-aware elementary
preconditioning tests and validation tests with no failures.

- [ ] **Step 6: Run package verification**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: PASS for the public and internal default groups.

- [ ] **Step 7: Run full-suite verification**

Run:

```bash
julia --project=. test/runtests.jl all
```

Expected: PASS for public, internal, and expert groups, including
`expert/elementary_preconditioning.jl`.

- [ ] **Step 8: Commit the implementation**

Run:

```bash
git add src/Suslin.jl src/core/elementary_matrices.jl test/expert/elementary_preconditioning.jl test/runtests.jl test/public/api_surface.jl docs/superpowers/plans/2026-06-20-issue-43-side-aware-elementary-preconditioning.md
git commit -m "feat: add side-aware elementary preconditioning"
```

Expected: commit succeeds on branch
`agent/issue-43-add-side-aware-elementary-preconditioning-primit-run-1`.
