# Issue 156 Determinant-Deferred Laurent Peel Metadata Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an internal determinant-deferred Laurent peel certificate that records elementary peel steps before determinant classification and verifies replay to `blockdiag(deferred_submatrix, I)`.

**Architecture:** Extend the existing Laurent column-peel module with a sibling internal certificate and replay verifier. The certificate reuses `LaurentColumnPeelStep`, embeds accumulated left/right peel factors into the original dimension, and stops at the deferred submatrix instead of normalizing its determinant.

**Tech Stack:** Julia, Oscar matrices over Laurent polynomial rings, existing Suslin Laurent column-peel helpers, Test stdlib.

## Global Constraints

- Keep the public API unchanged.
- Modify `src/algorithm/laurent_column_peel.jl`; no new include is required unless the implementation is split.
- Add `test/expert/laurent_lazy_peel_certificate.jl`.
- Register the expert test in `test/runtests.jl`.
- Reuse `LaurentColumnPeelStep` and existing replay helpers where possible.
- The certificate must include the original matrix, peel steps, deferred submatrix, `determinant_source = :deferred_submatrix`, and replay data.
- The replay relation is `left_product * original_matrix * right_product == blockdiag(deferred_submatrix, I)`.
- Do not normalize or correct a non-one deferred determinant.
- Do not add row/column correction options or public certificate fields.
- Focused verification command is `julia --project=. -e 'include("test/expert/laurent_lazy_peel_certificate.jl")'`.
- Package verification command is `julia --project=. -e 'using Pkg; Pkg.test()'`.
- Do not commit `Manifest.toml`.

---

## File Structure

- Modify `src/algorithm/laurent_column_peel.jl`: add `LaurentDeterminantDeferredPeelCertificate`, replay-data helpers, `_laurent_determinant_deferred_peel_certificate`, and `_verify_laurent_determinant_deferred_peel_replay`.
- Create `test/expert/laurent_lazy_peel_certificate.jl`: focused tests for determinant-one replay parity, deferred metadata, and tamper rejection.
- Modify `test/runtests.jl`: add the new expert test file near `expert/laurent_lazy_peel_no_initial_det.jl`.
- Keep this plan and `docs/superpowers/specs/2026-06-28-issue-156-determinant-deferred-laurent-peel-design.md` as workflow artifacts.

---

### Task 1: Add the Failing Expert Test

**Files:**
- Create: `test/expert/laurent_lazy_peel_certificate.jl`
- Modify: `test/runtests.jl`

**Interfaces:**
- Consumes: planned `Suslin._laurent_determinant_deferred_peel_certificate(A)` and `Suslin._verify_laurent_determinant_deferred_peel_replay(certificate)`.
- Produces: a focused expert test that fails before production code exists.

- [ ] **Step 1: Write the failing test**

Create `test/expert/laurent_lazy_peel_certificate.jl`:

```julia
using Test
using Suslin
using Oscar

include("../fixtures/laurent_lazy_determinant_cases.jl")

function _issue156_fixture(id::AbstractString)
    catalog = LaurentLazyDeterminantCases.catalog()
    return only(filter(entry -> entry.id == id, catalog.cases))
end

function _issue156_blockdiag_deferred(deferred, original_dimension::Int)
    return block_embedding(deferred, original_dimension, collect(1:nrows(deferred)))
end

function _issue156_replace_first_factor_with_identity(factors, R, n::Int)
    tampered = copy(factors)
    tampered[1] = identity_matrix(R, n)
    return tampered
end

function _issue156_tamper_first_step_factor(certificate)
    steps = collect(certificate.peel_steps)
    first_step = first(steps)
    R = base_ring(first_step.input_matrix)

    if !isempty(first_step.left_factors)
        bad_left = _issue156_replace_first_factor_with_identity(
            first_step.left_factors,
            R,
            first_step.dimension,
        )
        steps[1] = Suslin.LaurentColumnPeelStep(
            first_step.dimension,
            first_step.input_matrix,
            first_step.last_column,
            bad_left,
            first_step.after_left_matrix,
            first_step.right_factors,
            first_step.peeled_matrix,
            first_step.next_block,
        )
    elseif !isempty(first_step.right_factors)
        bad_right = _issue156_replace_first_factor_with_identity(
            first_step.right_factors,
            R,
            first_step.dimension,
        )
        steps[1] = Suslin.LaurentColumnPeelStep(
            first_step.dimension,
            first_step.input_matrix,
            first_step.last_column,
            first_step.left_factors,
            first_step.after_left_matrix,
            bad_right,
            first_step.peeled_matrix,
            first_step.next_block,
        )
    else
        error("fixture did not produce a tamperable peel factor")
    end

    return Suslin.LaurentDeterminantDeferredPeelCertificate(
        certificate.original_matrix,
        steps,
        certificate.deferred_submatrix,
        certificate.determinant_source,
        certificate.verification,
    )
end

@testset "determinant-deferred lazy Laurent peel certificate" begin
    entry = _issue156_fixture("determinant-one-triangular")
    A = entry.inputs.matrix
    n = nrows(A)

    deferred = Suslin._laurent_determinant_deferred_peel_certificate(A)
    sl_certificate = Suslin._factor_laurent_sl_column_peel(A)

    @test deferred.original_matrix == A
    @test deferred.determinant_source == :deferred_submatrix
    @test !isempty(deferred.peel_steps)
    @test nrows(deferred.deferred_submatrix) < n
    @test ncols(deferred.deferred_submatrix) < ncols(A)
    @test deferred.deferred_submatrix == first(sl_certificate.peel_steps).next_block
    @test first(deferred.peel_steps).left_factors == first(sl_certificate.peel_steps).left_factors
    @test first(deferred.peel_steps).right_factors == first(sl_certificate.peel_steps).right_factors
    @test Suslin._verify_laurent_determinant_deferred_peel_replay(deferred)
    @test deferred.left_product * A * deferred.right_product == deferred.target_matrix
    @test deferred.target_matrix == _issue156_blockdiag_deferred(deferred.deferred_submatrix, n)
    @test deferred.verification.overall_ok
    @test deferred.verification.relation_ok
    @test deferred.verification.replay_metadata_ok

    first_step = first(deferred.peel_steps)
    @test !isempty(first_step.left_factors) || !isempty(first_step.right_factors)
    tampered = _issue156_tamper_first_step_factor(deferred)
    replay_ok = try
        Suslin._verify_laurent_determinant_deferred_peel_replay(tampered)
    catch err
        err isa InterruptException && rethrow()
        false
    end
    @test !replay_ok
end
```

Register the test in `test/runtests.jl` by adding:

```julia
"expert/laurent_lazy_peel_certificate.jl",
```

immediately after `"expert/laurent_lazy_peel_no_initial_det.jl",`.

- [ ] **Step 2: Run the focused test to verify RED**

Run:

```bash
julia --project=. -e 'include("test/expert/laurent_lazy_peel_certificate.jl")'
```

Expected: FAIL with `UndefVarError` for `_laurent_determinant_deferred_peel_certificate` or `LaurentDeterminantDeferredPeelCertificate`.

- [ ] **Step 3: Commit the failing test**

Run:

```bash
git add test/expert/laurent_lazy_peel_certificate.jl test/runtests.jl
git commit -m "test: add determinant-deferred laurent peel certificate coverage"
```

---

### Task 2: Implement the Internal Deferred Certificate

**Files:**
- Modify: `src/algorithm/laurent_column_peel.jl`

**Interfaces:**
- Consumes: `LaurentColumnPeelStep`, `_laurent_column_peel_step`, `_is_valid_laurent_column_peel_step_data`, `_factor_product`, `block_embedding`, and `_emit_laurent_column_peel_progress`.
- Produces: `LaurentDeterminantDeferredPeelCertificate`, `_laurent_determinant_deferred_peel_certificate`, `_verify_laurent_determinant_deferred_peel_replay`, and replay helpers.

- [ ] **Step 1: Add the certificate type**

Add this struct after `LaurentColumnPeelStep`:

```julia
struct LaurentDeterminantDeferredPeelCertificate
    original_matrix
    peel_steps::Vector{LaurentColumnPeelStep}
    deferred_submatrix
    determinant_source::Symbol
    left_factors::Vector
    right_factors::Vector
    left_product
    right_product
    target_matrix
    verification

    function LaurentDeterminantDeferredPeelCertificate(
        original_matrix,
        peel_steps::Vector{LaurentColumnPeelStep},
        deferred_submatrix,
        determinant_source::Symbol,
        verification,
    )
        replay = _laurent_determinant_deferred_peel_replay_data(
            original_matrix,
            peel_steps,
            deferred_submatrix,
        )
        return new(
            original_matrix,
            collect(peel_steps),
            deferred_submatrix,
            determinant_source,
            replay.left_factors,
            replay.right_factors,
            replay.left_product,
            replay.right_product,
            replay.target_matrix,
            verification,
        )
    end
end
```

- [ ] **Step 2: Add replay-data helpers**

Add helper functions near `_embed_upper_left_factors`:

```julia
function _embed_laurent_deferred_peel_factors(factors, R, original_dimension::Int, factor_dimension::Int)
    collected = collect(factors)
    if isempty(collected)
        return typeof(identity_matrix(R, original_dimension))[]
    end
    factor_dimension == original_dimension && return collected
    return [block_embedding(factor, original_dimension, collect(1:factor_dimension)) for factor in collected]
end

function _laurent_determinant_deferred_target(deferred_submatrix, original_dimension::Int)
    return block_embedding(deferred_submatrix, original_dimension, collect(1:nrows(deferred_submatrix)))
end

function _laurent_determinant_deferred_peel_replay_data(original_matrix, peel_steps, deferred_submatrix)
    R = base_ring(original_matrix)
    original_dimension = nrows(original_matrix)
    left_factors = typeof(identity_matrix(R, original_dimension))[]
    right_factors = typeof(identity_matrix(R, original_dimension))[]
    steps = collect(peel_steps)

    for step in Iterators.reverse(steps)
        append!(
            left_factors,
            _embed_laurent_deferred_peel_factors(
                step.left_factors,
                R,
                original_dimension,
                step.dimension,
            ),
        )
    end

    for step in steps
        append!(
            right_factors,
            _embed_laurent_deferred_peel_factors(
                step.right_factors,
                R,
                original_dimension,
                step.dimension,
            ),
        )
    end

    left_product = _factor_product(left_factors, R, original_dimension)
    right_product = _factor_product(right_factors, R, original_dimension)
    target_matrix = _laurent_determinant_deferred_target(deferred_submatrix, original_dimension)
    return (;
        left_factors,
        right_factors,
        left_product,
        right_product,
        target_matrix,
    )
end
```

- [ ] **Step 3: Add the entry point and verifier**

Add these functions after `_factor_laurent_gl_lazy_determinant_peel`:

```julia
function _laurent_determinant_deferred_peel_certificate(
    A;
    min_steps::Int = 1,
    progress_callback = nothing,
)
    _validate_laurent_column_peel_input_shape_and_ring(A)
    min_steps >= 1 || throw(ArgumentError("determinant-deferred Laurent peel requires at least one peel step"))
    n = nrows(A)
    n - min_steps >= 2 || throw(ArgumentError("determinant-deferred Laurent peel must leave a deferred submatrix of size at least 2"))

    steps = LaurentColumnPeelStep[]
    current = A
    completed_steps = 0
    last_completed = _laurent_column_peel_empty_last_completed()
    _emit_laurent_column_peel_progress(progress_callback, current, completed_steps, last_completed)

    for _ in 1:min_steps
        step_started_at = time()
        step = _laurent_column_peel_step(current)
        step_elapsed = round(time() - step_started_at; digits = 3)
        push!(steps, step)
        completed_steps += 1
        last_completed = (;
            dimension = step.dimension,
            elapsed_seconds = step_elapsed,
            left_factors = length(step.left_factors),
            right_factors = length(step.right_factors),
        )
        current = step.next_block
        _emit_laurent_column_peel_progress(progress_callback, current, completed_steps, last_completed)
    end

    certificate = LaurentDeterminantDeferredPeelCertificate(
        A,
        steps,
        current,
        :deferred_submatrix,
        nothing,
    )
    verification = _laurent_determinant_deferred_peel_verification(certificate)
    verification.overall_ok || error("internal determinant-deferred Laurent peel verification failed")
    return LaurentDeterminantDeferredPeelCertificate(
        certificate.original_matrix,
        certificate.peel_steps,
        certificate.deferred_submatrix,
        certificate.determinant_source,
        verification,
    )
end

function _verify_laurent_determinant_deferred_peel_replay(certificate)::Bool
    return _laurent_determinant_deferred_peel_verification(certificate).overall_ok
end
```

Add the verification helpers after `_laurent_column_peel_verification`:

```julia
function _laurent_determinant_deferred_peel_verification(certificate)
    R = base_ring(certificate.original_matrix)
    n = nrows(certificate.original_matrix)
    size_ok = nrows(certificate.original_matrix) == ncols(certificate.original_matrix) >= 3
    determinant_source_ok = certificate.determinant_source == :deferred_submatrix
    step_chain_ok = _laurent_determinant_deferred_step_chain_ok(
        certificate.peel_steps,
        certificate.original_matrix,
        certificate.deferred_submatrix,
    )
    steps_ok = try
        all(step -> _is_valid_laurent_column_peel_step_data(
                step.dimension,
                step.input_matrix,
                step.last_column,
                step.left_factors,
                step.after_left_matrix,
                step.right_factors,
                step.peeled_matrix,
                step.next_block,
            ),
            certificate.peel_steps,
        )
    catch err
        err isa InterruptException && rethrow()
        false
    end
    deferred_shape_ok = nrows(certificate.deferred_submatrix) == ncols(certificate.deferred_submatrix) &&
        nrows(certificate.deferred_submatrix) >= 2 &&
        nrows(certificate.deferred_submatrix) < n
    replay = try
        _laurent_determinant_deferred_peel_replay_data(
            certificate.original_matrix,
            certificate.peel_steps,
            certificate.deferred_submatrix,
        )
    catch err
        err isa InterruptException && rethrow()
        nothing
    end
    replay_metadata_ok = replay !== nothing &&
        _factor_sequences_equal(certificate.left_factors, replay.left_factors) &&
        _factor_sequences_equal(certificate.right_factors, replay.right_factors) &&
        certificate.left_product == replay.left_product &&
        certificate.right_product == replay.right_product &&
        certificate.target_matrix == replay.target_matrix
    target_ok = replay !== nothing &&
        certificate.target_matrix == _laurent_determinant_deferred_target(certificate.deferred_submatrix, n)
    relation_ok = try
        certificate.left_product * certificate.original_matrix * certificate.right_product ==
            certificate.target_matrix
    catch err
        err isa InterruptException && rethrow()
        false
    end
    overall_ok = size_ok && determinant_source_ok && step_chain_ok && steps_ok &&
        deferred_shape_ok && replay_metadata_ok && target_ok && relation_ok
    return (;
        overall_ok,
        size_ok,
        determinant_source_ok,
        step_chain_ok,
        steps_ok,
        deferred_shape_ok,
        replay_metadata_ok,
        target_ok,
        relation_ok,
    )
end

function _laurent_determinant_deferred_step_chain_ok(steps, original_matrix, deferred_submatrix)::Bool
    collected = collect(steps)
    isempty(collected) && return false
    current = original_matrix

    for step in collected
        step.input_matrix == current || return false
        step.dimension == nrows(current) || return false
        step.dimension == ncols(current) || return false
        current = step.next_block
    end

    return current == deferred_submatrix
end
```

- [ ] **Step 4: Run focused test to verify GREEN**

Run:

```bash
julia --project=. -e 'include("test/expert/laurent_lazy_peel_certificate.jl")'
```

Expected: PASS.

- [ ] **Step 5: Commit the implementation**

Run:

```bash
git add src/algorithm/laurent_column_peel.jl
git commit -m "feat: add determinant-deferred laurent peel metadata"
```

---

### Task 3: Verify and Review the Branch

**Files:**
- Review: `src/algorithm/laurent_column_peel.jl`
- Review: `test/expert/laurent_lazy_peel_certificate.jl`
- Review: `test/runtests.jl`

**Interfaces:**
- Consumes: completed Tasks 1 and 2.
- Produces: fresh verification evidence for the final branch.

- [ ] **Step 1: Run focused verification**

Run:

```bash
julia --project=. -e 'include("test/expert/laurent_lazy_peel_certificate.jl")'
```

Expected: PASS.

- [ ] **Step 2: Run package verification**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: PASS for the default public and internal groups.

- [ ] **Step 3: Inspect git diff**

Run:

```bash
git diff --stat main...HEAD
git diff --check
git status --short
```

Expected: no whitespace errors and no unintended untracked files.

- [ ] **Step 4: Request final code review**

Use `superpowers:requesting-code-review` with the merge base against `main` and
the current `HEAD`. Fix any Critical or Important findings, then rerun the
focused and package verification commands.

---

## Plan Self-Review

- Spec coverage: Task 1 covers the expert test and negative control; Task 2
  covers the internal certificate, deferred source metadata, replay data, and
  verifier; Task 3 covers required verification and review.
- Placeholder scan: no task contains placeholders or vague implementation
  instructions.
- Type consistency: the planned certificate, helper, and verifier names match
  across tasks and tests.
