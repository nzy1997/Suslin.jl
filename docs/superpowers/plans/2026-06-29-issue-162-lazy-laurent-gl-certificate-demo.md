# Issue 162 Lazy Laurent GL Certificate Demo Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Demonstrate the public lazy determinant Laurent `GL_n` certificate on the original issue #38 `Q` block while preserving the staged `elementary_factorization(Q)` boundary.

**Architecture:** Extend the existing issue #38 expert certificate test with lazy original-input assertions and negative controls, then update the issue #38 example to print a concise lazy certificate PASS line. Refresh only the README/docs scope wording needed to describe the supported lazy monomial-unit path without claiming arbitrary Laurent `GL_n` support.

**Tech Stack:** Julia, Suslin.jl, Oscar, Test standard library, existing ToricBuilder issue #38 fixture.

## Global Constraints

- Do not change `elementary_factorization(Q)` to return original Laurent `GL_n` factors.
- Do not claim arbitrary Laurent `GL_n` support.
- Use the public lazy certificate API from #160: `laurent_gl_factorization_certificate(Q; determinant_strategy = :lazy, correction_side = :row)` and `verify_laurent_gl_factorization_certificate(certificate)`.
- The example PASS line must include determinant, correction side, determinant source, factor count, and `verified=true`.
- The expert test must assert exact reconstruction of the original issue #38 `Q` through the lazy certificate.
- The expert test must tamper with the lazy determinant correction, the reported correction side, or one hoisted elementary factor and assert that verification returns `false`.
- The expert test must keep asserting that `elementary_factorization(Q)` remains a staged boundary for the original `GL_n` input.

---

## File Structure

- Modify `test/expert/issue38_laurent_gl_certificate.jl`: add local lazy-certificate rebuild/tamper helpers and extend the issue #38 testset with original-input lazy certificate assertions and negative controls.
- Modify `example/toric_decoupling/issue38_unit_correction.jl`: add a helper that constructs and verifies the public lazy row certificate and prints the concise PASS line.
- Modify `README.md`: update current-scope language for the lazy supported monomial-unit correction path.
- Modify `docs/src/index.md`: mirror the README scope wording for generated docs.

### Task 1: Expert Lazy Original-Input Certificate Coverage

**Files:**
- Modify: `test/expert/issue38_laurent_gl_certificate.jl`

**Interfaces:**
- Consumes: `LaurentLazyGLHoistCertificate`, `laurent_gl_factorization_certificate(Q; determinant_strategy = :lazy, correction_side = :row)`, `verify_laurent_gl_factorization_certificate(certificate)`, existing `_issue59_product`.
- Produces: focused issue #38 assertions and helpers `_issue162_rebuild_lazy_certificate`, `_issue162_tamper_lazy_correction`, `_issue162_tamper_lazy_correction_side`, and `_issue162_tamper_lazy_hoisted_factor`.

- [ ] **Step 1: Write the lazy expert assertions and tamper helpers**

Add these helpers after `_issue59_malformed_normalization_certificate`:

```julia
function _issue162_rebuild_lazy_certificate(
    certificate;
    overall_determinant = certificate.overall_determinant,
    determinant_source = certificate.determinant_source,
    correction_side = certificate.correction_side,
    reconstruction_relation = certificate.reconstruction_relation,
    correction = certificate.correction,
    inverse_correction = certificate.inverse_correction,
    normalized_deferred_core = certificate.normalized_deferred_core,
    normalized_deferred_factorization = certificate.normalized_deferred_factorization,
    normalized_deferred_factors = certificate.normalized_deferred_factors,
    rewritten_left_factors = certificate.rewritten_left_factors,
    rewritten_right_factors = certificate.rewritten_right_factors,
    elementary_factors = certificate.elementary_factors,
    elementary_product = certificate.elementary_product,
    reconstructed_product = certificate.reconstructed_product,
    verification = certificate.verification,
)
    return Suslin.LaurentLazyGLHoistCertificate(
        certificate.original_matrix,
        certificate.deferred_metadata,
        overall_determinant,
        determinant_source,
        correction_side,
        reconstruction_relation,
        correction,
        inverse_correction,
        normalized_deferred_core,
        normalized_deferred_factorization,
        normalized_deferred_factors,
        rewritten_left_factors,
        rewritten_right_factors,
        elementary_factors,
        elementary_product,
        reconstructed_product,
        verification,
    )
end

function _issue162_tamper_lazy_correction(certificate)
    R = base_ring(certificate.original_matrix)
    bad_factor = copy(certificate.correction.factor)
    bad_factor[1, 1] = one(R)
    bad_correction = merge(certificate.correction, (; factor = bad_factor))
    return _issue162_rebuild_lazy_certificate(
        certificate;
        correction = bad_correction,
    )
end

function _issue162_tamper_lazy_correction_side(certificate)
    return _issue162_rebuild_lazy_certificate(
        certificate;
        correction_side = :column,
    )
end

function _issue162_tamper_lazy_hoisted_factor(certificate)
    R = base_ring(certificate.original_matrix)
    n = nrows(certificate.original_matrix)
    bad_factors = copy(certificate.elementary_factors)
    bad_factors[1] = identity_matrix(R, n)
    bad_product = _issue59_product(bad_factors, R, n)
    return _issue162_rebuild_lazy_certificate(
        certificate;
        elementary_factors = bad_factors,
        elementary_product = bad_product,
        reconstructed_product = certificate.correction.factor * bad_product,
    )
end
```

Then add this block before the existing `original_err` assertion:

```julia
    lazy_certificate = laurent_gl_factorization_certificate(
        Q;
        determinant_strategy = :lazy,
        correction_side = :row,
    )
    @test lazy_certificate isa LaurentLazyGLHoistCertificate
    @test lazy_certificate.original_matrix == Q
    @test lazy_certificate.overall_determinant == u * v
    @test lazy_certificate.determinant_source == :deferred_submatrix
    @test lazy_certificate.correction_side == :row
    @test lazy_certificate.reconstruction_relation ==
        :left_correction_times_elementary_product
    @test length(lazy_certificate.elementary_factors) > 0
    @test lazy_certificate.elementary_product ==
        _issue59_product(lazy_certificate.elementary_factors, R, nrows(Q))
    @test lazy_certificate.correction.factor * lazy_certificate.elementary_product == Q
    @test lazy_certificate.reconstructed_product == Q
    @test lazy_certificate.verification.overall_ok
    @test verify_laurent_gl_factorization_certificate(lazy_certificate)

    @test !verify_laurent_gl_factorization_certificate(
        _issue162_tamper_lazy_correction(lazy_certificate),
    )
    @test !verify_laurent_gl_factorization_certificate(
        _issue162_tamper_lazy_correction_side(lazy_certificate),
    )
    @test !verify_laurent_gl_factorization_certificate(
        _issue162_tamper_lazy_hoisted_factor(lazy_certificate),
    )
```

- [ ] **Step 2: Run the expert test**

Run:

```bash
julia --project=. -e 'include("test/expert/issue38_laurent_gl_certificate.jl")'
```

Expected: PASS. This task extends assertions against an already public API, so the failure evidence for the missing issue behavior comes from the baseline example command, which currently lacks the required `lazy_gl_certificate` line.

- [ ] **Step 3: Commit**

```bash
git add test/expert/issue38_laurent_gl_certificate.jl
git commit -m "test: cover issue 38 lazy gl certificate"
```

### Task 2: Example Lazy Certificate PASS Line

**Files:**
- Modify: `example/toric_decoupling/issue38_unit_correction.jl`

**Interfaces:**
- Consumes: public `laurent_gl_factorization_certificate` and `verify_laurent_gl_factorization_certificate`.
- Produces: helper `_run_lazy_gl_certificate(Q)` and example output line beginning `lazy_gl_certificate status=PASS`.

- [ ] **Step 1: Confirm the red example behavior**

Run:

```bash
julia --project=. example/toric_decoupling/issue38_unit_correction.jl | rg 'lazy_gl_certificate status=PASS'
```

Expected: FAIL with exit code 1 because the current example does not print the required lazy certificate line.

- [ ] **Step 2: Add the example helper**

Add this function after `_run_column_unit_correction(Q)`:

```julia
function _run_lazy_gl_certificate(Q)
    certificate = laurent_gl_factorization_certificate(
        Q;
        determinant_strategy = :lazy,
        correction_side = :row,
    )
    verified = verify_laurent_gl_factorization_certificate(certificate)

    certificate.reconstructed_product == Q ||
        error("lazy GL certificate does not reconstruct original Q")
    verified || error("lazy GL certificate did not verify")

    println(
        "lazy_gl_certificate status=PASS det=$(certificate.overall_determinant) correction_side=$(certificate.correction_side) determinant_source=$(certificate.determinant_source) factors=$(length(certificate.elementary_factors)) verified=$(verified)",
    )
    return certificate
end
```

In `main()`, call the helper after the column correction and before the final summary:

```julia
    lazy = _run_lazy_gl_certificate(Q)
```

Update the final summary to include the lazy factor count:

```julia
    println(
        "issue38_unit_correction_example status=PASS row_factors=$(length(row.factors)) column_factors=$(length(column.factors)) lazy_factors=$(length(lazy.elementary_factors))",
    )
```

- [ ] **Step 3: Run the example**

Run:

```bash
julia --project=. example/toric_decoupling/issue38_unit_correction.jl
```

Expected: PASS and output includes:

```text
lazy_gl_certificate status=PASS det=u*v correction_side=row determinant_source=deferred_submatrix
```

- [ ] **Step 4: Commit**

```bash
git add example/toric_decoupling/issue38_unit_correction.jl
git commit -m "example: show issue 38 lazy gl certificate"
```

### Task 3: Scope Wording Refresh

**Files:**
- Modify: `README.md`
- Modify: `docs/src/index.md`

**Interfaces:**
- Consumes: current README/docs scope bullets.
- Produces: conservative scope text that mentions the lazy supported monomial-unit certificate path and still excludes arbitrary Laurent `GL_n`.

- [ ] **Step 1: Update README scope wording**

Replace the two Laurent certificate/scope bullets in `README.md` with:

```markdown
- `laurent_gl_factorization_certificate(A)` defaults to the eager Laurent
  normalization/core certificate. With `determinant_strategy = :lazy`, it
  records the supported monomial-unit deferred determinant correction path for
  original Laurent `GL_n` inputs such as the issue #38 fixture.
- The implementation is not yet the full Park-Woodburn algorithm for arbitrary
  `SL_n(k[x_1, ..., x_m])`, `n >= 3`. General Quillen local realizability,
  coefficient-ring support beyond exact field-backed ordinary polynomial
  rings, arbitrary Laurent `GL_n` determinant correction, and factor-count
  optimization remain staged boundaries.
```

- [ ] **Step 2: Mirror docs index wording**

Apply the same replacement in `docs/src/index.md`.

- [ ] **Step 3: Run targeted verification**

Run:

```bash
julia --project=. example/toric_decoupling/issue38_unit_correction.jl
julia --project=. -e 'include("test/expert/issue38_laurent_gl_certificate.jl")'
```

Expected: both PASS.

- [ ] **Step 4: Commit**

```bash
git add README.md docs/src/index.md
git commit -m "docs: describe lazy laurent gl certificate scope"
```

### Task 4: Final Verification and PR Prep

**Files:**
- No direct source edits.

**Interfaces:**
- Consumes: tasks 1-3.
- Produces: verified branch ready for PR.

- [ ] **Step 1: Run issue verification commands**

Run:

```bash
julia --project=. example/toric_decoupling/issue38_unit_correction.jl
julia --project=. -e 'include("test/expert/issue38_laurent_gl_certificate.jl")'
```

Expected: both PASS; example prints `lazy_gl_certificate status=PASS det=u*v correction_side=row determinant_source=deferred_submatrix ... verified=true`.

- [ ] **Step 2: Run package test gate**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: PASS.

- [ ] **Step 3: Review git diff**

Run:

```bash
git diff --stat origin/main...HEAD
git diff --check
```

Expected: no whitespace errors; changed files are limited to the design/plan docs, issue #38 expert test, issue #38 example, README, and docs index.
