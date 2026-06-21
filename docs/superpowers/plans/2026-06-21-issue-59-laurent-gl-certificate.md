# Issue 59 Laurent GL Certificate Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an exact Laurent `GL_n` certificate for the original Issue #38 `Q` block.

**Architecture:** Add a focused certificate layer that composes `normalize_laurent_gl_matrix(Q)` with the existing Issue #57 determinant-one column-peel factorization of the normalized core. Keep `elementary_factorization(Q)` as a pure elementary-factor API that rejects the original `det(Q)=u*v` input.

**Tech Stack:** Julia, Oscar Laurent polynomial rings and matrices, existing Suslin Laurent normalization, existing Laurent column-peel factorization, and Test stdlib.

## Global Constraints

- Do not claim that the original Issue #38 `Q` matrix is a pure elementary product.
- Do not broaden support to arbitrary Laurent `GL_n` inputs beyond the existing normalization and Issue #38 core factorization path.
- Do not change the mathematical meaning of `verify_factorization(A, factors)`.
- Use the left-correction variant from `normalize_laurent_gl_matrix(A)`: `A == correction.factor * normalized_core`.
- The certificate verifier must recompute exact reconstruction and must not trust stored fields.
- Negative controls must make the certificate verifier return `false` after tampering with either the determinant correction factor or one core factor.
- Focused verification command is `julia --project=. -e 'include("test/expert/issue38_laurent_gl_certificate.jl")'`.
- Package verification command is `julia --project=. -e 'using Pkg; Pkg.test()'`.
- Do not commit `Manifest.toml`.

---

## File Structure

- Create `test/expert/issue38_laurent_gl_certificate.jl`: focused RED/GREEN acceptance and negative controls for the original Issue #38 `Q`.
- Modify `test/runtests.jl`: register the focused expert test.
- Create `src/algorithm/laurent_gl_certificate.jl`: certificate type, constructor API, exact verifier, and helper verification metadata.
- Modify `src/Suslin.jl`: export the certificate API and include the new algorithm file after `algorithm/laurent_column_peel.jl`.
- Modify `test/public/api_surface.jl`: assert the new exported API is visible.

---

### Task 1: Add the Failing Issue 59 Expert Test

**Files:**
- Create: `test/expert/issue38_laurent_gl_certificate.jl`
- Modify: `test/runtests.jl`

**Interfaces:**
- Consumes: planned `laurent_gl_factorization_certificate(A)`, `verify_laurent_gl_factorization_certificate(certificate)`, and `LaurentGLFactorizationCertificate`.
- Produces: focused RED coverage for the original Issue #38 Laurent `GL_n` certificate.

- [ ] **Step 1: Write the failing test**

Create `test/expert/issue38_laurent_gl_certificate.jl`:

```julia
using Test
using Suslin
using Oscar

include("../fixtures/toricbuilder_issue38_cases.jl")

function _issue59_product(factors, R, n::Int)
    product = identity_matrix(R, n)
    for factor in factors
        product *= factor
    end
    return product
end

function _issue59_tamper_correction(certificate)
    R = base_ring(certificate.original_matrix)
    bad_factor = copy(certificate.correction.factor)
    bad_factor[1, 1] = one(R)
    bad_correction = merge(certificate.correction, (; factor = bad_factor))
    return Suslin.LaurentGLFactorizationCertificate(
        certificate.original_matrix,
        certificate.determinant_profile,
        certificate.normalization,
        bad_correction,
        certificate.inverse_correction,
        certificate.normalized_core,
        certificate.core_factorization,
        certificate.core_factors,
        certificate.reconstructed_product,
        certificate.verification,
    )
end

function _issue59_tamper_core_factor(certificate)
    R = base_ring(certificate.original_matrix)
    n = nrows(certificate.original_matrix)
    bad_factors = copy(certificate.core_factors)
    bad_factors[1] = identity_matrix(R, n)
    return Suslin.LaurentGLFactorizationCertificate(
        certificate.original_matrix,
        certificate.determinant_profile,
        certificate.normalization,
        certificate.correction,
        certificate.inverse_correction,
        certificate.normalized_core,
        certificate.core_factorization,
        bad_factors,
        certificate.reconstructed_product,
        certificate.verification,
    )
end

@testset "Issue 38 Laurent GL certificate" begin
    entry = only(ToricBuilderIssue38Cases.catalog().cases)
    Q = entry.inputs.matrix
    R = entry.ring.object
    u, v = entry.ring.generators

    determinant_profile = classify_laurent_determinant(Q)
    @test determinant_profile.classification == :laurent_monomial_unit
    @test determinant_profile.determinant == u * v

    certificate = laurent_gl_factorization_certificate(Q)
    @test certificate.original_matrix == Q
    @test certificate.determinant_profile.classification == :laurent_monomial_unit
    @test certificate.determinant_profile.determinant == u * v
    @test certificate.correction == certificate.normalization.correction
    @test certificate.inverse_correction == certificate.correction.inverse_factor
    @test certificate.normalized_core == certificate.normalization.normalized_matrix
    @test det(certificate.normalized_core) == one(R)
    @test length(certificate.core_factors) > 0
    @test verify_factorization(certificate.normalized_core, certificate.core_factors)
    @test !verify_factorization(Q, certificate.core_factors)

    core_product = _issue59_product(certificate.core_factors, R, nrows(Q))
    @test core_product == certificate.normalized_core
    @test certificate.correction.factor * core_product == Q
    @test certificate.reconstructed_product == Q
    @test verify_laurent_gl_factorization_certificate(certificate)

    tampered_correction = _issue59_tamper_correction(certificate)
    @test !verify_laurent_gl_factorization_certificate(tampered_correction)

    tampered_factor = _issue59_tamper_core_factor(certificate)
    @test !verify_laurent_gl_factorization_certificate(tampered_factor)

    original_err = try
        elementary_factorization(Q)
        nothing
    catch err
        err
    end
    @test original_err isa ArgumentError
    @test occursin("Laurent GL_n normalization boundary", sprint(showerror, original_err))
end
```

In `test/runtests.jl`, add the new expert file immediately after
`expert/laurent_column_peel_issue38.jl`:

```julia
"expert/issue38_laurent_gl_certificate.jl",
```

- [ ] **Step 2: Run the focused test to verify RED**

Run:

```bash
julia --project=. -e 'include("test/expert/issue38_laurent_gl_certificate.jl")'
```

Expected: FAIL with `UndefVarError: laurent_gl_factorization_certificate not defined`.

- [ ] **Step 3: Commit the RED test**

```bash
git add test/expert/issue38_laurent_gl_certificate.jl test/runtests.jl
git commit -m "test: cover Issue 38 Laurent GL certificate"
```

---

### Task 2: Implement the Laurent GL Certificate API

**Files:**
- Create: `src/algorithm/laurent_gl_certificate.jl`
- Modify: `src/Suslin.jl`
- Modify: `test/public/api_surface.jl`

**Interfaces:**
- Consumes: `normalize_laurent_gl_matrix(A)`, `_factor_laurent_sl_column_peel(core)`, `_verify_laurent_column_peel_replay(certificate)`, `_factor_product(factors, R, n)`, `_factor_sequences_equal(left, right)`, and `verify_factorization(A, factors)`.
- Produces: `LaurentGLFactorizationCertificate`, `laurent_gl_factorization_certificate(A)`, and `verify_laurent_gl_factorization_certificate(certificate)::Bool`.

- [ ] **Step 1: Implement the certificate module**

Create `src/algorithm/laurent_gl_certificate.jl`:

```julia
struct LaurentGLFactorizationCertificate
    original_matrix
    determinant_profile
    normalization
    correction
    inverse_correction
    normalized_core
    core_factorization
    core_factors::Vector
    reconstructed_product
    verification

    function LaurentGLFactorizationCertificate(
        original_matrix,
        determinant_profile,
        normalization,
        correction,
        inverse_correction,
        normalized_core,
        core_factorization,
        core_factors::Vector,
        reconstructed_product,
        verification,
    )
        R = base_ring(original_matrix)
        product = correction.factor * _factor_product(core_factors, R, nrows(original_matrix))
        return new(
            original_matrix,
            determinant_profile,
            normalization,
            correction,
            inverse_correction,
            normalized_core,
            core_factorization,
            core_factors,
            product,
            verification,
        )
    end
end

function laurent_gl_factorization_certificate(A)
    normalization = normalize_laurent_gl_matrix(A)
    core = normalization.normalized_matrix
    core_factorization = _factor_laurent_sl_column_peel(core)
    core_factors = core_factorization.factors
    R = base_ring(A)
    reconstructed_product = normalization.correction.factor * _factor_product(core_factors, R, nrows(A))
    certificate = LaurentGLFactorizationCertificate(
        A,
        normalization.determinant_profile,
        normalization,
        normalization.correction,
        normalization.correction.inverse_factor,
        core,
        core_factorization,
        core_factors,
        reconstructed_product,
        nothing,
    )
    verification = _laurent_gl_factorization_certificate_verification(certificate)
    verification.overall_ok || error("internal Laurent GL factorization certificate verification failed")
    return LaurentGLFactorizationCertificate(
        certificate.original_matrix,
        certificate.determinant_profile,
        certificate.normalization,
        certificate.correction,
        certificate.inverse_correction,
        certificate.normalized_core,
        certificate.core_factorization,
        certificate.core_factors,
        certificate.reconstructed_product,
        verification,
    )
end

function verify_laurent_gl_factorization_certificate(certificate)::Bool
    return _laurent_gl_factorization_certificate_verification(certificate).overall_ok
end

function _laurent_gl_factorization_certificate_verification(certificate)
    size_ok = try
        nrows(certificate.original_matrix) == ncols(certificate.original_matrix)
    catch err
        err isa InterruptException && rethrow()
        false
    end

    normalization_ok = false
    determinant_profile_ok = false
    correction_ok = false
    inverse_correction_ok = false
    normalized_core_ok = false
    core_det_ok = false
    core_replay_ok = false
    core_factors_match_ok = false
    core_factors_ok = false
    reconstructed_product_ok = false

    if size_ok
        try
            A = certificate.original_matrix
            R = base_ring(A)
            n = nrows(A)
            recomputed = normalize_laurent_gl_matrix(A)
            normalization_ok = verify_laurent_gl_normalization(A, certificate.normalization)
            determinant_profile_ok =
                certificate.determinant_profile == recomputed.determinant_profile &&
                certificate.normalization.determinant_profile == recomputed.determinant_profile
            correction_ok =
                certificate.correction == recomputed.correction &&
                certificate.normalization.correction == recomputed.correction
            inverse_correction_ok =
                certificate.inverse_correction == recomputed.correction.inverse_factor &&
                certificate.correction.inverse_factor == recomputed.correction.inverse_factor
            normalized_core_ok =
                certificate.normalized_core == recomputed.normalized_matrix &&
                certificate.normalization.normalized_matrix == recomputed.normalized_matrix
            core_det_ok = det(certificate.normalized_core) == one(R)
            core_replay_ok =
                certificate.core_factorization.original_matrix == certificate.normalized_core &&
                _verify_laurent_column_peel_replay(certificate.core_factorization)
            core_factors_match_ok = _factor_sequences_equal(
                certificate.core_factors,
                certificate.core_factorization.factors,
            )
            core_factors_ok = verify_factorization(certificate.normalized_core, certificate.core_factors)
            rebuilt_core = _factor_product(certificate.core_factors, R, n)
            rebuilt = certificate.correction.factor * rebuilt_core
            reconstructed_product_ok =
                rebuilt_core == certificate.normalized_core &&
                rebuilt == A &&
                certificate.reconstructed_product == rebuilt
        catch err
            err isa InterruptException && rethrow()
        end
    end

    overall_ok = size_ok && normalization_ok && determinant_profile_ok &&
        correction_ok && inverse_correction_ok && normalized_core_ok &&
        core_det_ok && core_replay_ok && core_factors_match_ok &&
        core_factors_ok && reconstructed_product_ok

    return (
        overall_ok = overall_ok,
        size_ok = size_ok,
        normalization_ok = normalization_ok,
        determinant_profile_ok = determinant_profile_ok,
        correction_ok = correction_ok,
        inverse_correction_ok = inverse_correction_ok,
        normalized_core_ok = normalized_core_ok,
        core_det_ok = core_det_ok,
        core_replay_ok = core_replay_ok,
        core_factors_match_ok = core_factors_match_ok,
        core_factors_ok = core_factors_ok,
        reconstructed_product_ok = reconstructed_product_ok,
    )
end
```

- [ ] **Step 2: Wire exports and include**

In `src/Suslin.jl`, add exports near the other Laurent APIs:

```julia
export LaurentGLFactorizationCertificate
export laurent_gl_factorization_certificate
export verify_laurent_gl_factorization_certificate
```

Include the new file after `include("algorithm/laurent_column_peel.jl")`:

```julia
include("algorithm/laurent_gl_certificate.jl")
```

- [ ] **Step 3: Add API surface checks**

In `test/public/api_surface.jl`, add `isdefined` checks:

```julia
@test isdefined(Suslin, :LaurentGLFactorizationCertificate)
@test isdefined(Suslin, :laurent_gl_factorization_certificate)
@test isdefined(Suslin, :verify_laurent_gl_factorization_certificate)
```

Add exported binding checks:

```julia
@test Suslin.LaurentGLFactorizationCertificate === LaurentGLFactorizationCertificate
@test Suslin.laurent_gl_factorization_certificate === laurent_gl_factorization_certificate
@test Suslin.verify_laurent_gl_factorization_certificate === verify_laurent_gl_factorization_certificate
```

- [ ] **Step 4: Run focused GREEN verification**

Run:

```bash
julia --project=. -e 'include("test/expert/issue38_laurent_gl_certificate.jl")'
```

Expected: PASS.

- [ ] **Step 5: Run package verification**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: PASS for the default public and internal groups.

- [ ] **Step 6: Commit the implementation**

```bash
git add src/Suslin.jl src/algorithm/laurent_gl_certificate.jl test/public/api_surface.jl
git commit -m "feat: add Laurent GL certificate API"
```

---

## Plan Self-Review

- Spec coverage: Task 1 covers the exact Issue #38 fixture, determinant metadata, core factor verification, exact reconstruction, pure-factor negative control, and tampering controls. Task 2 implements the certificate object, verifier, exports, and API surface.
- Placeholder scan: no TODO, TBD, placeholder, or incomplete implementation markers remain.
- Type consistency: the plan consistently uses `LaurentGLFactorizationCertificate`, `laurent_gl_factorization_certificate(A)`, `verify_laurent_gl_factorization_certificate(certificate)::Bool`, and the existing `normalization.correction.factor` left-correction model.
