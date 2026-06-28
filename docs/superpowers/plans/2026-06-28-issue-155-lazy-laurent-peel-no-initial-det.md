# Issue 155 Lazy Laurent Peel No Initial Determinant Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a behavioral expert guard proving the lazy Laurent peel route completes one real peel step before determinant classification of the original input can occur.

**Architecture:** Add one internal lazy peel entry point that validates only shape/ring, performs one existing column-peel step, emits a completed progress event, then invokes an injectable determinant probe on the smaller block. Add an optional determinant probe keyword to the existing eager peel for the negative control while preserving its default eager behavior.

**Tech Stack:** Julia, Oscar Laurent polynomial matrices, Suslin internal Laurent column-peel helpers, Test stdlib.

## Global Constraints

- Do not implement the full lazy correction certificate.
- Do not add or change public exports.
- Lazy route name is `Suslin._factor_laurent_gl_lazy_determinant_peel(A; determinant_probe = classify_laurent_determinant, progress_callback = nothing)`.
- The lazy determinant probe must not be called at the original matrix dimension for inputs with at least one peel step.
- The lazy determinant probe must be called only after `completed_steps >= 1`.
- The negative control must prove the eager path invokes its determinant probe at the original matrix before any completed peel step.
- Expert test file is `test/expert/laurent_lazy_peel_no_initial_det.jl`.
- Register the test in the `expert` group in `test/runtests.jl`.
- Focused verification command is `julia --project=. -e 'include("test/expert/laurent_lazy_peel_no_initial_det.jl")'`.
- Package verification command is `julia --project=. -e 'using Pkg; Pkg.test()'`.
- Do not commit `Manifest.toml`.

---

## File Structure

- Create `test/expert/laurent_lazy_peel_no_initial_det.jl`: loads the #154 lazy determinant fixture catalog, defines sentinel determinant probes, asserts lazy deferral, and asserts eager-path negative control.
- Modify `test/runtests.jl`: add `"expert/laurent_lazy_peel_no_initial_det.jl"` after the existing Laurent column peel expert test.
- Modify `src/algorithm/laurent_column_peel.jl`: add the eager determinant probe keyword and the new lazy one-peel entry point.
- Keep `docs/superpowers/specs/2026-06-28-issue-155-lazy-laurent-peel-no-initial-det-design.md` and this plan as workflow artifacts.

---

### Task 1: Add the Failing Expert Guard

**Files:**
- Create: `test/expert/laurent_lazy_peel_no_initial_det.jl`
- Modify: `test/runtests.jl`

**Interfaces:**
- Consumes: `LaurentLazyDeterminantCases.catalog()` from `test/fixtures/laurent_lazy_determinant_cases.jl`.
- Produces: an expert test that calls planned internal functions
  `Suslin._factor_laurent_gl_lazy_determinant_peel(A; determinant_probe, progress_callback)` and
  `Suslin._factor_laurent_sl_column_peel(A; determinant_probe, progress_callback)`.

- [ ] **Step 1: Write the failing expert test**

Create `test/expert/laurent_lazy_peel_no_initial_det.jl` with:

```julia
using Test
using Suslin
using Oscar

include("../fixtures/laurent_lazy_determinant_cases.jl")

struct _Issue155InitialDeterminantProbeError <: Exception
    message::String
end

Base.showerror(io::IO, err::_Issue155InitialDeterminantProbeError) =
    print(io, err.message)

function _issue155_fixture(id::AbstractString)
    catalog = LaurentLazyDeterminantCases.catalog()
    return only(filter(entry -> entry.id == id, catalog.cases))
end

function _issue155_max_completed_steps(progress_records)
    isempty(progress_records) && return 0
    return maximum(record.completed_steps for record in progress_records)
end

function _issue155_lazy_probe(original_dimension::Int, progress_records, probe_records)
    return function (candidate)
        completed_before_probe = _issue155_max_completed_steps(progress_records)
        push!(probe_records, (;
            size = (nrows(candidate), ncols(candidate)),
            completed_before_probe,
        ))
        if nrows(candidate) == original_dimension || completed_before_probe < 1
            throw(_Issue155InitialDeterminantProbeError(
                "initial determinant classification invoked before lazy Laurent peel completed a step",
            ))
        end
        return Suslin.classify_laurent_determinant(candidate)
    end
end

function _issue155_eager_probe(probe_records)
    return function (candidate)
        push!(probe_records, (; size = (nrows(candidate), ncols(candidate))))
        throw(_Issue155InitialDeterminantProbeError(
            "eager determinant classification invoked at the original Laurent matrix",
        ))
    end
end

@testset "lazy Laurent peel defers initial determinant classification" begin
    entry = _issue155_fixture("monomial-unit-row-column-cores")
    A = entry.inputs.matrix
    original_size = (nrows(A), ncols(A))

    lazy_progress = Any[]
    lazy_probes = Any[]
    lazy_err = try
        Suslin._factor_laurent_gl_lazy_determinant_peel(
            A;
            progress_callback = record -> push!(lazy_progress, record),
            determinant_probe = _issue155_lazy_probe(original_size[1], lazy_progress, lazy_probes),
        )
        nothing
    catch err
        err
    end

    @test lazy_err isa ArgumentError
    @test occursin("lazy Laurent determinant correction", sprint(showerror, lazy_err))
    @test !isempty(lazy_progress)
    @test any(record -> record.completed_steps >= 1, lazy_progress)
    @test !isempty(lazy_probes)
    @test first(lazy_probes).size[1] < original_size[1]
    @test first(lazy_probes).size[2] < original_size[2]
    @test first(lazy_probes).completed_before_probe >= 1
    @test !(lazy_err isa _Issue155InitialDeterminantProbeError)

    eager_progress = Any[]
    eager_probes = Any[]
    eager_err = try
        Suslin._factor_laurent_sl_column_peel(
            A;
            progress_callback = record -> push!(eager_progress, record),
            determinant_probe = _issue155_eager_probe(eager_probes),
        )
        nothing
    catch err
        err
    end

    @test eager_err isa _Issue155InitialDeterminantProbeError
    @test occursin("eager determinant classification", sprint(showerror, eager_err))
    @test !isempty(eager_progress)
    @test _issue155_max_completed_steps(eager_progress) == 0
    @test !isempty(eager_probes)
    @test first(eager_probes).size == original_size
end
```

- [ ] **Step 2: Register the expert test**

In `test/runtests.jl`, add the new file immediately after
`"expert/laurent_column_peel_issue38.jl",`:

```julia
        "expert/laurent_lazy_peel_no_initial_det.jl",
```

- [ ] **Step 3: Run the focused red test**

Run:

```bash
julia --project=. -e 'include("test/expert/laurent_lazy_peel_no_initial_det.jl")'
```

Expected: FAIL because `Suslin._factor_laurent_gl_lazy_determinant_peel` is not defined, or because `_factor_laurent_sl_column_peel` does not yet accept `determinant_probe`.

---

### Task 2: Implement the Lazy and Eager Probe Surfaces

**Files:**
- Modify: `src/algorithm/laurent_column_peel.jl`

**Interfaces:**
- Consumes: existing `_validate_laurent_column_peel_input_shape_and_ring(A)`,
  `_emit_laurent_column_peel_progress`, `_laurent_column_peel_step`,
  `_laurent_column_peel_recursive`, `_inverse_elementary_sequence`,
  `_embed_upper_left_factors`, and `LaurentColumnPeelFactorization`.
- Produces:
  `Suslin._factor_laurent_gl_lazy_determinant_peel(A; determinant_probe, progress_callback)` and
  `Suslin._factor_laurent_sl_column_peel(A; progress_callback, determinant_probe)`.

- [ ] **Step 1: Add the determinant probe keyword to eager validation**

Change `_validate_laurent_column_peel_input_determinant_one` to:

```julia
function _validate_laurent_column_peel_input_determinant_one(
    A;
    determinant_probe = classify_laurent_determinant,
)
    profile = determinant_probe(A)
    profile.classification == :one || throw(ArgumentError("Laurent column-peel factorization requires determinant-one input"))
    return nrows(A)
end
```

- [ ] **Step 2: Thread the probe through the eager peel entry point**

Change the eager signature and validation call to:

```julia
function _factor_laurent_sl_column_peel(
    A;
    progress_callback = nothing,
    determinant_probe = classify_laurent_determinant,
)
    _validate_laurent_column_peel_input_shape_and_ring(A)
    _emit_laurent_column_peel_progress(
        progress_callback,
        A,
        0,
        _laurent_column_peel_empty_last_completed(),
    )
    _validate_laurent_column_peel_input_determinant_one(A; determinant_probe)
```

Keep the existing recursive body after that validation call unchanged.

- [ ] **Step 3: Add the lazy one-peel entry point**

Add this function after `_factor_laurent_sl_column_peel`:

```julia
function _factor_laurent_gl_lazy_determinant_peel(
    A;
    determinant_probe = classify_laurent_determinant,
    progress_callback = nothing,
)
    _validate_laurent_column_peel_input_shape_and_ring(A)
    n = nrows(A)
    n >= 3 || throw(ArgumentError("Lazy Laurent determinant peel requires size at least 3 so one peel step can complete before determinant classification"))
    empty_last_completed = _laurent_column_peel_empty_last_completed()
    _emit_laurent_column_peel_progress(progress_callback, A, 0, empty_last_completed)

    step_started_at = time()
    step = _laurent_column_peel_step(A)
    step_elapsed = round(time() - step_started_at; digits = 3)
    completed_steps = 1
    last_completed = (;
        dimension = n,
        elapsed_seconds = step_elapsed,
        left_factors = length(step.left_factors),
        right_factors = length(step.right_factors),
    )
    _emit_laurent_column_peel_progress(progress_callback, step.next_block, completed_steps, last_completed)

    deferred_profile = determinant_probe(step.next_block)
    deferred_profile.classification == :one || throw(ArgumentError(
        "lazy Laurent determinant correction after initial peel is not implemented for deferred determinant classification $(deferred_profile.classification)",
    ))

    next_factors, next_steps, final_block, final_target, final_local, final_2x2 =
        _laurent_column_peel_recursive(
            step.next_block;
            progress_callback,
            completed_steps,
            last_completed,
            emit_current = false,
        )
    R = base_ring(A)
    factors = vcat(
        _inverse_elementary_sequence(step.left_factors),
        _embed_upper_left_factors(next_factors, R, n),
        _inverse_elementary_sequence(step.right_factors),
    )
    product = _factor_product(factors, R, n)
    certificate = LaurentColumnPeelFactorization(
        A,
        final_block,
        final_target,
        final_local,
        final_2x2,
        factors,
        product,
        vcat(LaurentColumnPeelStep[step], next_steps),
        nothing,
    )
    verification = _laurent_column_peel_verification(certificate)
    verification.overall_ok || error("internal lazy Laurent column-peel verification failed")
    return LaurentColumnPeelFactorization(
        certificate.original_matrix,
        certificate.final_block,
        certificate.final_local_target,
        certificate.final_local_factors,
        certificate.final_factors,
        certificate.factors,
        certificate.product,
        certificate.peel_steps,
        verification,
    )
end
```

- [ ] **Step 4: Run the focused green test**

Run:

```bash
julia --project=. -e 'include("test/expert/laurent_lazy_peel_no_initial_det.jl")'
```

Expected: PASS.

---

### Task 3: Verify and Commit

**Files:**
- Verify: `test/expert/laurent_lazy_peel_no_initial_det.jl`
- Verify: `src/algorithm/laurent_column_peel.jl`
- Verify: `test/runtests.jl`
- Commit: implementation files and workflow plan/spec artifacts.

**Interfaces:**
- Consumes: Task 1 test and Task 2 implementation.
- Produces: a reviewed, committed branch ready for PR creation.

- [ ] **Step 1: Run package verification**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: PASS with public and internal test groups passing.

- [ ] **Step 2: Run the expert guard one more time**

Run:

```bash
julia --project=. -e 'include("test/expert/laurent_lazy_peel_no_initial_det.jl")'
```

Expected: PASS.

- [ ] **Step 3: Inspect the diff**

Run:

```bash
git status --short
git diff --stat
git diff -- src/algorithm/laurent_column_peel.jl test/expert/laurent_lazy_peel_no_initial_det.jl test/runtests.jl
```

Expected: only the intended source, test, test registry, and workflow artifact files changed.

- [ ] **Step 4: Commit**

Run:

```bash
git add src/algorithm/laurent_column_peel.jl test/expert/laurent_lazy_peel_no_initial_det.jl test/runtests.jl docs/superpowers/plans/2026-06-28-issue-155-lazy-laurent-peel-no-initial-det.md
git commit -m "test: guard lazy laurent peel determinant deferral"
```

Expected: commit succeeds.
