# Issue 157 Deferred Laurent Submatrix Normalization Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Enrich determinant-deferred Laurent peel certificates with deferred-submatrix determinant classification, supported monomial-unit normalization, and structured unsupported boundaries.

**Architecture:** Add an internal metadata layer in `src/algorithm/laurent_column_peel.jl` that consumes `LaurentDeterminantDeferredPeelCertificate` without changing public exports. The existing lazy GL peel route keeps returning the current factorization for determinant-one deferred cores, but returns enriched metadata for non-one deferred determinants instead of throwing.

**Tech Stack:** Julia, Oscar matrices over Laurent polynomial rings, existing Suslin Laurent GL normalization helpers, Test stdlib.

## Global Constraints

- Keep the public API unchanged and do not export new names.
- Do not recompute the determinant of the original full matrix in the lazy deferred path.
- Classify only `certificate.deferred_submatrix` for the new metadata.
- Supported deferred determinant classes are exactly `:one` and `:laurent_monomial_unit`.
- Non-unit and unsupported unit classes must return structured staged boundary data.
- Do not hoist the deferred correction back to the original matrix-level certificate.
- Add `test/expert/laurent_lazy_submatrix_normalization.jl`.
- Register the new expert test in `test/runtests.jl`.
- Update the existing lazy no-initial-det expert test so the monomial-unit route returns metadata rather than an error.
- Focused verification command is `julia --project=. -e 'include("test/expert/laurent_lazy_submatrix_normalization.jl")'`.
- Existing lazy-route verification command is `julia --project=. -e 'include("test/expert/laurent_lazy_peel_no_initial_det.jl")'`.
- Package verification command is `julia --project=. -e 'using Pkg; Pkg.test()'`.
- Do not commit `Manifest.toml`.

---

## File Structure

- Modify `src/algorithm/laurent_column_peel.jl`: add deferred-submatrix normalization helpers, metadata verification, and update `_factor_laurent_gl_lazy_determinant_peel` at the existing deferred determinant boundary.
- Create `test/expert/laurent_lazy_submatrix_normalization.jl`: focused tests for determinant-one, monomial-unit normalization, deferred-only probing, and non-unit staged boundary behavior.
- Modify `test/expert/laurent_lazy_peel_no_initial_det.jl`: replace the old monomial-unit error expectation with metadata expectations while preserving the eager-determinant negative control.
- Modify `test/runtests.jl`: register the new expert test next to the existing lazy Laurent peel tests.

---

### Task 1: Add the Failing Expert Tests

**Files:**
- Create: `test/expert/laurent_lazy_submatrix_normalization.jl`
- Modify: `test/expert/laurent_lazy_peel_no_initial_det.jl`
- Modify: `test/runtests.jl`

**Interfaces:**
- Consumes planned `Suslin._normalize_laurent_determinant_deferred_submatrix(certificate; determinant_probe = ...)` and `Suslin._verify_laurent_determinant_deferred_submatrix_normalization(metadata)`.
- Produces failing tests that define the exact metadata fields and lazy route behavior.

- [ ] **Step 1: Write the failing test file**

Create `test/expert/laurent_lazy_submatrix_normalization.jl` with:

```julia
using Test
using Suslin
using Oscar

include("../fixtures/laurent_lazy_determinant_cases.jl")

function _issue157_fixture(id::AbstractString)
    catalog = LaurentLazyDeterminantCases.catalog()
    return only(filter(entry -> entry.id == id, catalog.cases))
end

function _issue157_deferred_probe(certificate, records)
    return function (candidate)
        push!(records, (; size = (nrows(candidate), ncols(candidate)), candidate))
        @test candidate == certificate.deferred_submatrix
        @test candidate != certificate.original_matrix
        return Suslin.classify_laurent_determinant(candidate)
    end
end

function _issue157_embedded_deferred(deferred)
    R = base_ring(deferred)
    return block_embedding(deferred, nrows(deferred) + 1, collect(1:nrows(deferred)))
end

@testset "deferred Laurent submatrix determinant-one normalization metadata" begin
    entry = _issue157_fixture("determinant-one-triangular")
    certificate = Suslin._laurent_determinant_deferred_peel_certificate(entry.inputs.matrix)
    metadata = Suslin._normalize_laurent_determinant_deferred_submatrix(certificate)
    R = base_ring(certificate.deferred_submatrix)

    @test metadata.peel_certificate == certificate
    @test metadata.determinant_source == :deferred_submatrix
    @test metadata.overall_determinant == one(R)
    @test metadata.determinant_classification == :one
    @test metadata.supported
    @test metadata.deferred_correction.kind == :identity
    @test metadata.deferred_diagonal_correction === nothing
    @test metadata.normalized_deferred_core == certificate.deferred_submatrix
    @test det(metadata.normalized_deferred_core) == one(R)
    @test metadata.staged_boundary === nothing
    @test metadata.verification.overall_ok
    @test Suslin._verify_laurent_determinant_deferred_submatrix_normalization(metadata)
end

@testset "deferred Laurent submatrix monomial-unit normalization metadata" begin
    R, (u, v) = suslin_laurent_polynomial_ring(GF(2), ["u", "v"])
    deferred = matrix(R, [
        u * v one(R);
        zero(R) one(R)
    ])
    A = _issue157_embedded_deferred(deferred)
    certificate = Suslin._laurent_determinant_deferred_peel_certificate(A)
    probe_records = Any[]
    metadata = Suslin._normalize_laurent_determinant_deferred_submatrix(
        certificate;
        determinant_probe = _issue157_deferred_probe(certificate, probe_records),
    )

    @test length(probe_records) == 1
    @test only(probe_records).size == (2, 2)
    @test metadata.overall_determinant == u * v
    @test metadata.determinant_classification == :laurent_monomial_unit
    @test metadata.supported
    @test metadata.deferred_correction.kind == :left_diagonal_determinant_correction
    @test metadata.deferred_diagonal_correction == metadata.deferred_correction
    @test metadata.deferred_correction.factor * metadata.normalized_deferred_core ==
        certificate.deferred_submatrix
    @test det(metadata.normalized_deferred_core) == one(R)
    @test metadata.staged_boundary === nothing
    @test metadata.verification.overall_ok
    @test metadata.verification.normalized_core_det_ok
    @test Suslin._verify_laurent_determinant_deferred_submatrix_normalization(metadata)
end

@testset "deferred Laurent submatrix non-unit staged boundary" begin
    entry = _issue157_fixture("non-unit-determinant-negative")
    A = _issue157_embedded_deferred(entry.inputs.matrix)
    certificate = Suslin._laurent_determinant_deferred_peel_certificate(A)
    metadata = Suslin._normalize_laurent_determinant_deferred_submatrix(certificate)

    @test certificate.deferred_submatrix == entry.inputs.matrix
    @test metadata.overall_determinant == entry.determinant_profile.expected_determinant
    @test metadata.determinant_classification == :non_unit
    @test !metadata.supported
    @test metadata.deferred_correction === nothing
    @test metadata.deferred_diagonal_correction === nothing
    @test metadata.normalized_deferred_core === nothing
    @test metadata.staged_boundary !== nothing
    @test metadata.staged_boundary.kind == :unsupported_deferred_laurent_determinant
    @test metadata.staged_boundary.reason == :non_unit_deferred_determinant
    @test metadata.staged_boundary.overall_determinant == metadata.overall_determinant
    @test metadata.verification.overall_ok
    @test metadata.verification.boundary_ok
    @test Suslin._verify_laurent_determinant_deferred_submatrix_normalization(metadata)
end

@testset "lazy GL peel returns enriched metadata for non-one deferred determinant" begin
    entry = _issue157_fixture("monomial-unit-row-column-cores")
    metadata = Suslin._factor_laurent_gl_lazy_determinant_peel(entry.inputs.matrix)

    @test metadata.determinant_source == :deferred_submatrix
    @test metadata.overall_determinant == entry.determinant_profile.expected_determinant
    @test metadata.determinant_classification == :laurent_monomial_unit
    @test metadata.supported
    @test metadata.normalized_deferred_core !== nothing
    @test det(metadata.normalized_deferred_core) == one(base_ring(metadata.normalized_deferred_core))
    @test metadata.staged_boundary === nothing
end
```

- [ ] **Step 2: Update the old lazy route test**

In `test/expert/laurent_lazy_peel_no_initial_det.jl`, replace the first testset's
`lazy_err = try ... catch` block and assertions for the monomial-unit fixture
with:

```julia
    lazy_metadata = Suslin._factor_laurent_gl_lazy_determinant_peel(
        A;
        progress_callback = record -> push!(lazy_progress, record),
        determinant_probe = _issue155_lazy_probe(original_size[1], lazy_progress, lazy_probes),
    )

    @test lazy_metadata.determinant_source == :deferred_submatrix
    @test lazy_metadata.determinant_classification == :laurent_monomial_unit
    @test lazy_metadata.supported
    @test lazy_metadata.normalized_deferred_core !== nothing
    @test det(lazy_metadata.normalized_deferred_core) ==
        one(base_ring(lazy_metadata.normalized_deferred_core))
    @test lazy_metadata.staged_boundary === nothing
```

Keep the existing progress, probe, and eager-route assertions after that block,
except remove `@test !(lazy_err isa _Issue155InitialDeterminantProbeError)`.

- [ ] **Step 3: Register the new expert test**

In `test/runtests.jl`, add:

```julia
"expert/laurent_lazy_submatrix_normalization.jl",
```

immediately after `"expert/laurent_lazy_peel_certificate.jl",`.

- [ ] **Step 4: Run focused tests to verify RED**

Run:

```bash
julia --project=. -e 'include("test/expert/laurent_lazy_submatrix_normalization.jl")'
```

Expected: FAIL with `UndefVarError` for `_normalize_laurent_determinant_deferred_submatrix`.

Run:

```bash
julia --project=. -e 'include("test/expert/laurent_lazy_peel_no_initial_det.jl")'
```

Expected: FAIL because `_factor_laurent_gl_lazy_determinant_peel` still throws on the monomial-unit deferred determinant.

- [ ] **Step 5: Commit the failing tests**

Run:

```bash
git add test/expert/laurent_lazy_submatrix_normalization.jl test/expert/laurent_lazy_peel_no_initial_det.jl test/runtests.jl
git commit -m "test: cover deferred laurent submatrix normalization"
```

---

### Task 2: Implement Deferred Submatrix Normalization Metadata

**Files:**
- Modify: `src/algorithm/laurent_column_peel.jl`

**Interfaces:**
- Consumes: `LaurentDeterminantDeferredPeelCertificate`, `_verify_laurent_determinant_deferred_peel_replay`, `_identity_correction`, `_left_diagonal_determinant_correction`, and `classify_laurent_determinant`.
- Produces: `_normalize_laurent_determinant_deferred_submatrix`, `_verify_laurent_determinant_deferred_submatrix_normalization`, and helper functions.

- [ ] **Step 1: Add supported-class and correction helpers**

Add these helpers after `_verify_laurent_determinant_deferred_peel_replay`:

```julia
function _is_supported_deferred_laurent_determinant_class(classification::Symbol)::Bool
    return classification == :one || classification == :laurent_monomial_unit
end

function _scoped_deferred_correction(correction)
    return merge(correction, (; scope = :deferred_submatrix))
end

function _deferred_laurent_determinant_boundary(certificate, determinant_profile)
    classification = determinant_profile.classification
    reason = classification == :non_unit ?
        :non_unit_deferred_determinant :
        :unsupported_deferred_unit_class
    deferred = certificate.deferred_submatrix
    return (;
        kind = :unsupported_deferred_laurent_determinant,
        determinant_source = certificate.determinant_source,
        overall_determinant = determinant_profile.determinant,
        determinant_classification = classification,
        deferred_submatrix_size = (nrows(deferred), ncols(deferred)),
        supported = false,
        reason,
    )
end
```

- [ ] **Step 2: Add the enrichment function**

Add:

```julia
function _normalize_laurent_determinant_deferred_submatrix(
    certificate;
    determinant_probe = classify_laurent_determinant,
)
    _verify_laurent_determinant_deferred_peel_replay(certificate) ||
        error("invalid determinant-deferred Laurent peel certificate replay")

    deferred = certificate.deferred_submatrix
    R = base_ring(deferred)
    n = nrows(deferred)
    determinant_profile = determinant_probe(deferred)
    classification = determinant_profile.classification
    determinant = determinant_profile.determinant

    deferred_correction = nothing
    deferred_diagonal_correction = nothing
    normalized_deferred_core = nothing
    staged_boundary = nothing

    if classification == :one
        deferred_correction = _scoped_deferred_correction(
            _identity_correction(R, n, determinant),
        )
        normalized_deferred_core = deferred
    elseif classification == :laurent_monomial_unit
        deferred_correction = _scoped_deferred_correction(
            _left_diagonal_determinant_correction(R, n, determinant),
        )
        deferred_diagonal_correction = deferred_correction
        normalized_deferred_core = deferred_correction.inverse_factor * deferred
    else
        staged_boundary = _deferred_laurent_determinant_boundary(
            certificate,
            determinant_profile,
        )
    end

    metadata = (;
        peel_certificate = certificate,
        deferred_submatrix = deferred,
        determinant_source = certificate.determinant_source,
        determinant_profile,
        overall_determinant = determinant,
        determinant_classification = classification,
        supported = _is_supported_deferred_laurent_determinant_class(classification),
        deferred_correction,
        deferred_diagonal_correction,
        normalized_deferred_core,
        staged_boundary,
        verification = nothing,
    )
    verification = _laurent_determinant_deferred_submatrix_normalization_verification(metadata)
    verification.overall_ok ||
        error("internal deferred Laurent submatrix normalization metadata verification failed")
    return merge(metadata, (; verification))
end
```

- [ ] **Step 3: Add metadata verification**

Add:

```julia
function _verify_laurent_determinant_deferred_submatrix_normalization(metadata)::Bool
    return _laurent_determinant_deferred_submatrix_normalization_verification(metadata).overall_ok
end

function _laurent_determinant_deferred_submatrix_normalization_verification(metadata)
    certificate_ok = try
        _verify_laurent_determinant_deferred_peel_replay(metadata.peel_certificate)
    catch err
        err isa InterruptException && rethrow()
        false
    end

    determinant_ok = false
    correction_ok = false
    normalized_core_det_ok = false
    boundary_ok = false

    try
        deferred = metadata.peel_certificate.deferred_submatrix
        R = base_ring(deferred)
        n = nrows(deferred)
        profile = classify_laurent_determinant(deferred)
        classification = profile.classification
        determinant_ok =
            metadata.deferred_submatrix == deferred &&
            metadata.determinant_source == :deferred_submatrix &&
            metadata.overall_determinant == profile.determinant &&
            metadata.determinant_classification == classification &&
            metadata.determinant_profile == profile &&
            metadata.supported == _is_supported_deferred_laurent_determinant_class(classification)

        if metadata.supported
            correction = metadata.deferred_correction
            identity = identity_matrix(R, n)
            correction_ok =
                metadata.staged_boundary === nothing &&
                correction !== nothing &&
                correction.scope == :deferred_submatrix &&
                correction.determinant == metadata.overall_determinant &&
                nrows(correction.factor) == n &&
                ncols(correction.factor) == n &&
                nrows(correction.inverse_factor) == n &&
                ncols(correction.inverse_factor) == n &&
                correction.factor * correction.inverse_factor == identity &&
                correction.inverse_factor * correction.factor == identity

            if classification == :one
                correction_ok = correction_ok &&
                    correction.kind == :identity &&
                    metadata.deferred_diagonal_correction === nothing &&
                    metadata.normalized_deferred_core == deferred
            elseif classification == :laurent_monomial_unit
                correction_ok = correction_ok &&
                    correction.kind == :left_diagonal_determinant_correction &&
                    metadata.deferred_diagonal_correction == correction &&
                    correction.factor * metadata.normalized_deferred_core == deferred
            end

            normalized_core_det_ok =
                metadata.normalized_deferred_core !== nothing &&
                det(metadata.normalized_deferred_core) == one(R)
        else
            boundary = metadata.staged_boundary
            boundary_ok =
                metadata.deferred_correction === nothing &&
                metadata.deferred_diagonal_correction === nothing &&
                metadata.normalized_deferred_core === nothing &&
                boundary !== nothing &&
                boundary.kind == :unsupported_deferred_laurent_determinant &&
                boundary.determinant_source == :deferred_submatrix &&
                boundary.overall_determinant == metadata.overall_determinant &&
                boundary.determinant_classification == metadata.determinant_classification &&
                boundary.deferred_submatrix_size == (n, n) &&
                boundary.supported == false
        end
    catch err
        err isa InterruptException && rethrow()
    end

    supported_ok = metadata.supported ? (correction_ok && normalized_core_det_ok) : boundary_ok
    overall_ok = certificate_ok && determinant_ok && supported_ok
    return (;
        overall_ok,
        certificate_ok,
        determinant_ok,
        correction_ok,
        normalized_core_det_ok,
        boundary_ok,
    )
end
```

- [ ] **Step 4: Update the lazy GL peel route**

In `_factor_laurent_gl_lazy_determinant_peel`, replace the direct
`deferred_profile = determinant_probe(step.next_block)` block with:

```julia
    deferred_certificate = LaurentDeterminantDeferredPeelCertificate(
        A,
        LaurentColumnPeelStep[step],
        step.next_block,
        :deferred_submatrix,
        nothing,
    )
    deferred_metadata = _normalize_laurent_determinant_deferred_submatrix(
        deferred_certificate;
        determinant_probe,
    )
    deferred_metadata.determinant_classification == :one || return deferred_metadata
```

Keep the existing recursive continuation for the determinant-one path.

- [ ] **Step 5: Run focused tests to verify GREEN**

Run:

```bash
julia --project=. -e 'include("test/expert/laurent_lazy_submatrix_normalization.jl")'
julia --project=. -e 'include("test/expert/laurent_lazy_peel_no_initial_det.jl")'
```

Expected: both commands PASS.

- [ ] **Step 6: Commit the implementation**

Run:

```bash
git add src/algorithm/laurent_column_peel.jl
git commit -m "feat: normalize deferred laurent submatrix determinants"
```

---

### Task 3: Final Verification and Review

**Files:**
- Review: `src/algorithm/laurent_column_peel.jl`
- Review: `test/expert/laurent_lazy_submatrix_normalization.jl`
- Review: `test/expert/laurent_lazy_peel_no_initial_det.jl`
- Review: `test/runtests.jl`

**Interfaces:**
- Consumes completed Tasks 1 and 2.
- Produces fresh verification evidence and a clean branch ready for PR.

- [ ] **Step 1: Run focused verification**

Run:

```bash
julia --project=. -e 'include("test/expert/laurent_lazy_submatrix_normalization.jl")'
julia --project=. -e 'include("test/expert/laurent_lazy_peel_no_initial_det.jl")'
```

Expected: both commands PASS.

- [ ] **Step 2: Run package verification**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: PASS for the default public and internal groups.

- [ ] **Step 3: Inspect git state**

Run:

```bash
git diff --check
git status --short --ignored
git log --oneline main..HEAD
```

Expected: no whitespace errors; no unintended tracked or untracked files;
`Manifest.toml` may appear only as ignored.

- [ ] **Step 4: Request final code review**

Use `superpowers:requesting-code-review` with the merge base against `main` and
the current `HEAD`. Fix any Critical or Important findings, then rerun the
focused and package verification commands.

---

## Plan Self-Review

- Spec coverage: Task 1 covers all required output fields, the `u*v` monomial
  fixture, deferred-only determinant probing, #154 non-unit negative control,
  and lazy route behavior. Task 2 implements the metadata and route behavior.
  Task 3 covers verification and review.
- Placeholder scan: no placeholders, incomplete steps, or vague test
  instructions remain.
- Type consistency: function names, metadata field names, and verification
  field names match across tests and implementation tasks.
