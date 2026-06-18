# Issue 7 Laurent Ring Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a Laurent polynomial ring constructor and strict Laurent parent validators for later Suslin Laurent workflows.

**Architecture:** Keep the feature in `src/core/rings.jl`, alongside the existing polynomial-ring wrapper. Export only the new constructor from `src/Suslin.jl`; keep validator helpers internal and reusable through module-qualified names. Cover the new behavior with focused internal tests and a minimal public API surface assertion.

**Tech Stack:** Julia, Oscar/AbstractAlgebra Laurent polynomial rings, Test stdlib, existing Suslin grouped test runner.

## Global Constraints

- Keep `suslin_polynomial_ring(F, names::Vector{String})` unchanged.
- `suslin_laurent_polynomial_ring(F, names::Vector{String})` must return an Oscar Laurent polynomial ring plus a concrete generator vector.
- Returned Laurent generators must belong to the reported parent ring.
- Validators must reject ordinary polynomial parents and ordinary polynomial elements where Laurent input is required.
- Validators must reject elements from a different Laurent parent with `ArgumentError`.
- Parent matching must use identity (`===`) rather than coercion or equality.
- Register the focused internal Laurent test file in the existing `internal` test group.
- Do not commit `Manifest.toml`.

---

## File Structure

- Modify `src/Suslin.jl`: export `suslin_laurent_polynomial_ring`.
- Modify `src/core/rings.jl`: add the Laurent constructor and internal validator helpers.
- Create `test/internal/laurent_rings.jl`: focused constructor and validator tests.
- Modify `test/runtests.jl`: include `internal/laurent_rings.jl` in the internal group.
- Modify `test/public/api_surface.jl`: assert the new constructor is defined and exported.

---

### Task 1: Laurent Ring Constructor and Validators

**Files:**
- Modify: `src/Suslin.jl`
- Modify: `src/core/rings.jl`
- Create: `test/internal/laurent_rings.jl`
- Modify: `test/runtests.jl`
- Modify: `test/public/api_surface.jl`

**Interfaces:**
- Consumes: `Oscar.laurent_polynomial_ring(F, names)` and existing Suslin export/include structure.
- Produces: `suslin_laurent_polynomial_ring(F, names::Vector{String})`,
  `Suslin._is_laurent_polynomial_ring(R)`,
  `Suslin._require_laurent_polynomial_ring(R; label="ring")`,
  `Suslin._require_laurent_element(value; label="value")`,
  `Suslin._require_laurent_element(value, R; label="value")`, and
  `Suslin._require_same_laurent_parent(values; label="values")`.

- [ ] **Step 1: Write the failing Laurent tests**

Create `test/internal/laurent_rings.jl` with:

```julia
using Test
using Suslin
using Oscar

@testset "suslin Laurent polynomial ring" begin
    R, vars = suslin_laurent_polynomial_ring(GF(2), ["x", "y"])
    x, y = vars

    @test length(vars) == 2
    @test ngens(R) == 2
    @test parent(x) === R
    @test parent(y) === R
    @test string(x) == "x"
    @test string(y) == "y"

    x_inv = x^-1
    mixed = x^-1 * y
    @test parent(x_inv) === R
    @test parent(mixed) === R
    @test string(x_inv) == "x^-1"
    @test string(mixed) == "x^-1*y"
end

@testset "Laurent validators" begin
    R, (x, y) = suslin_laurent_polynomial_ring(GF(2), ["x", "y"])
    S, (u, v) = suslin_laurent_polynomial_ring(GF(2), ["u", "v"])
    P, (px, py) = suslin_polynomial_ring(GF(2), ["x", "y"])

    @test Suslin._is_laurent_polynomial_ring(R)
    @test !Suslin._is_laurent_polynomial_ring(P)
    @test Suslin._require_laurent_polynomial_ring(R) === R
    @test_throws ArgumentError Suslin._require_laurent_polynomial_ring(P)

    value = x^-1 * y
    @test Suslin._require_laurent_element(value) === value
    @test Suslin._require_laurent_element(value, R) === value
    @test_throws ArgumentError Suslin._require_laurent_element(u, R)
    @test_throws ArgumentError Suslin._require_laurent_element(px)
    @test_throws ArgumentError Suslin._require_laurent_element(px, R)

    @test Suslin._require_same_laurent_parent([x, value, y^-1]) === R
    @test_throws ArgumentError Suslin._require_same_laurent_parent([x, u])
    @test_throws ArgumentError Suslin._require_same_laurent_parent([x, px])
    @test_throws ArgumentError Suslin._require_same_laurent_parent([])
end
```

Update `test/runtests.jl` so the `internal` group contains:

```julia
"internal" => [
    "internal/rings.jl",
    "internal/laurent_rings.jl",
],
```

Update `test/public/api_surface.jl` inside the existing `api surface` testset:

```julia
@test isdefined(Suslin, :suslin_laurent_polynomial_ring)
@test Suslin.suslin_laurent_polynomial_ring === suslin_laurent_polynomial_ring
```

- [ ] **Step 2: Run the new test to verify it fails for the missing constructor**

Run:

```bash
julia --project=. -e 'include(test/internal/laurent_rings.jl)'
```

Expected: FAIL with `UndefVarError: suslin_laurent_polynomial_ring not defined` or an equivalent failure showing the constructor/helper API is missing.

- [ ] **Step 3: Implement the minimal constructor and validators**

Update `src/Suslin.jl` near the existing ring export:

```julia
export suslin_polynomial_ring
export suslin_laurent_polynomial_ring
```

Replace `src/core/rings.jl` with:

```julia
function suslin_polynomial_ring(F, names::Vector{String})
    R, vars = Oscar.polynomial_ring(F, names)
    return R, collect(vars)
end

function suslin_laurent_polynomial_ring(F, names::Vector{String})
    R, vars = Oscar.laurent_polynomial_ring(F, names)
    return R, collect(vars)
end

function _is_laurent_polynomial_ring(R)
    return R isa LaurentPolyRing || R isa LaurentMPolyRing
end

function _require_laurent_polynomial_ring(R; label::AbstractString="ring")
    _is_laurent_polynomial_ring(R) && return R
    throw(ArgumentError("$label must be a Laurent polynomial ring"))
end

function _parent_for_validation(value, label::AbstractString)
    try
        return parent(value)
    catch err
        err isa MethodError || rethrow()
        throw(ArgumentError("$label must be an element of a Laurent polynomial ring"))
    end
end

function _require_laurent_element(value; label::AbstractString="value")
    _require_laurent_polynomial_ring(_parent_for_validation(value, label); label="$label parent")
    return value
end

function _require_laurent_element(value, R; label::AbstractString="value")
    _require_laurent_polynomial_ring(R; label="expected parent")
    parent_value = _parent_for_validation(value, label)
    parent_value === R || throw(ArgumentError("$label must belong to the expected Laurent polynomial ring"))
    return value
end

function _require_same_laurent_parent(values; label::AbstractString="values")
    state = iterate(values)
    state === nothing && throw(ArgumentError("$label must be nonempty"))

    first_value, next_state = state
    R = _require_laurent_polynomial_ring(_parent_for_validation(first_value, label); label="$label parent")

    while true
        state = iterate(values, next_state)
        state === nothing && return R
        value, next_state = state
        _require_laurent_element(value, R; label=label)
    end
end
```

- [ ] **Step 4: Run the focused Laurent test to verify it passes**

Run:

```bash
julia --project=. -e 'include(test/internal/laurent_rings.jl)'
```

Expected: PASS for both Laurent testsets.

- [ ] **Step 5: Run the default package test entry point**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: PASS for the default `public` and `internal` groups.

- [ ] **Step 6: Commit the implementation**

Run:

```bash
git status --short
git add src/Suslin.jl src/core/rings.jl test/internal/laurent_rings.jl test/runtests.jl test/public/api_surface.jl
git commit -m "feat: add laurent ring helpers"
```

Expected: commit succeeds and `Manifest.toml` is not staged.

---

## Plan Self-Review

- Spec coverage: The task covers the constructor, internal validators, strict
  parent identity checks, focused tests, public API assertion, default package
  test command, and no-Manifest constraint.
- Placeholder scan: No unresolved markers or unspecified edge handling remains.
- Type consistency: The produced function names match the design spec and the
  tests use only exported constructor names or `Suslin.`-qualified internal
  helpers.

## Execution Choice

Plan complete and saved to
`docs/superpowers/plans/2026-06-18-issue-7-laurent-ring-plan.md`.

Two execution options:

1. Subagent-Driven (recommended) - dispatch a fresh subagent per task and review
   between tasks.
2. Inline Execution - execute tasks in this session using executing-plans, batch
   execution with checkpoints.

Under the standing answer policy for this non-interactive Agent Desk run, choose
Option 1 because it is marked recommended.
