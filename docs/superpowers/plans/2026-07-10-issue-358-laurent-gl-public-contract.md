# Issue 358 Laurent GL Public Contract Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Finalize the public Laurent `GL_n` certificate contract while keeping `elementary_factorization(A)` an elementary-only determinant-one API.

**Architecture:** Reuse the existing eager `LaurentGLFactorizationCertificate` and strengthen its verifier rather than adding a mixed factor type. Update the Laurent determinant error so callers are directed to the certificate API, then align public docs, example output, and issue-38/public API tests with the final contract.

**Tech Stack:** Julia, Oscar, `Test`, repository-local Documenter markdown.

## Global Constraints

- Do not introduce heterogeneous factor types.
- Do not pretend a nontrivial monomial-unit determinant correction is elementary.
- Do not add non-monomial unit support or arbitrary coefficient-ring support.
- Do not bundle LaurentToPoly, d14 conversion, endpoint-reduction, or factor-count optimization work.
- Preserve `elementary_factorization(A)` as a factor-sequence API whose returned factors multiply exactly to `A`.
- Use exact reconstruction checks for every certificate acceptance path.

---

## File Structure

- Modify `src/core/elementary_matrices.jl`: add a private elementary-factor predicate shared by verifiers.
- Modify `src/algorithm/laurent_gl_certificate.jl`: add docstrings and strengthen eager certificate verification.
- Modify `src/algorithm/factorization.jl`: add a docstring and replace the Laurent monomial-unit error text.
- Modify `test/expert/issue38_laurent_gl_certificate.jl`: add issue-358 assertions and negative controls.
- Modify `test/public/laurent_gl_certificate_options.jl` and `test/public/api_surface.jl`: pin public error guidance and reconstruction relation.
- Modify `example/toric_decoupling/issue38_unit_correction.jl`: demonstrate the final certificate contract.
- Modify `README.md`, `docs/src/index.md`, and `docs/src/toricbuilder_contract.md`: document the split `SL_n`/Laurent `GL_n` public contract.

---

### Task 1: Pin The Public Laurent GL Contract In Tests

**Files:**
- Modify: `test/expert/issue38_laurent_gl_certificate.jl`
- Modify: `test/public/laurent_gl_certificate_options.jl`
- Modify: `test/public/api_surface.jl`

**Interfaces:**
- Consumes: `laurent_gl_factorization_certificate(A)`, `verify_laurent_gl_factorization_certificate(cert)`, `elementary_factorization(A)`, `classify_laurent_determinant(A)`.
- Produces: failing tests for `verification.core_factors_elementary_ok`, exact eager reconstruction relation, determinant-contract error guidance, correction-side tampering, reordered/modified core factors, copied reconstructed product rejection, and non-unit determinant rejection.

- [ ] **Step 1: Add issue-358 expert helpers**

In `test/expert/issue38_laurent_gl_certificate.jl`, add these helpers after `_issue59_tamper_correction`:

```julia
function _issue358_tamper_correction_side(certificate)
    bad_correction = merge(certificate.correction, (; side = :right))
    bad_normalization = merge(certificate.normalization, (; correction = bad_correction))
    return Suslin.LaurentGLFactorizationCertificate(
        certificate.original_matrix,
        certificate.determinant_profile,
        bad_normalization,
        bad_correction,
        certificate.inverse_correction,
        certificate.normalized_core,
        certificate.core_factorization,
        certificate.core_factors,
        certificate.reconstructed_product,
        certificate.verification,
    )
end

function _issue358_tamper_reordered_core_factor(certificate)
    length(certificate.core_factors) >= 2 ||
        error("issue-358 reorder control requires at least two core factors")
    bad_factors = copy(certificate.core_factors)
    bad_factors[1], bad_factors[2] = bad_factors[2], bad_factors[1]
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

function _issue358_non_unit_matrix(R, u)
    return diagonal_matrix(R, [u + one(R), one(R), one(R)])
end
```

- [ ] **Step 2: Add issue-358 expert assertions**

In the main issue-38 testset, after `@test det(certificate.normalized_core) == one(R)`, add:

```julia
    @test certificate.verification.overall_ok
    @test certificate.verification.core_factors_elementary_ok
    @test certificate.verification.reconstructed_product_ok
    @test certificate.correction.side == :left
    @test certificate.correction.factor * _issue59_product(
        certificate.core_factors,
        R,
        nrows(Q),
    ) == Q
```

After the existing `tampered_correction` assertion, add:

```julia
    tampered_correction_side = _issue358_tamper_correction_side(certificate)
    @test !verify_laurent_gl_factorization_certificate(tampered_correction_side)

    tampered_reordered_factor = _issue358_tamper_reordered_core_factor(certificate)
    @test !verify_laurent_gl_factorization_certificate(tampered_reordered_factor)
```

Before the `elementary_factorization(Q)` error block, add:

```julia
    non_unit = _issue358_non_unit_matrix(R, u)
    @test classify_laurent_determinant(non_unit).classification == :non_unit
    non_unit_err = try
        laurent_gl_factorization_certificate(non_unit)
        nothing
    catch err
        err
    end
    @test non_unit_err isa ArgumentError
    @test occursin("non-unit", sprint(showerror, non_unit_err))
```

Replace the old error-message assertion with:

```julia
    original_message = sprint(showerror, original_err)
    @test occursin("elementary_factorization(A) is an elementary-only SL_n API", original_message)
    @test occursin("requires determinant 1", original_message)
    @test occursin("laurent_gl_factorization_certificate(A)", original_message)
```

- [ ] **Step 3: Update public option test assertions**

In `test/public/laurent_gl_certificate_options.jl`, after the eager certificate reconstruction assertion, add:

```julia
    eager_core_product = _issue160_factor_product(
        eager_certificate.core_factors,
        base_ring(Q),
        nrows(Q),
    )
    @test eager_core_product == eager_certificate.normalized_core
    @test eager_certificate.correction.factor * eager_core_product == Q
    @test eager_certificate.verification.core_factors_elementary_ok
```

If `_issue160_factor_product` does not exist, add it after `_issue160_caught_error`:

```julia
function _issue160_factor_product(factors, R, n::Int)
    product = identity_matrix(R, n)
    for factor in factors
        product *= factor
    end
    return product
end
```

Replace the old original-error assertion with the same three `occursin`
checks from Step 2.

- [ ] **Step 4: Add public API surface monomial-unit error check**

In `test/public/api_surface.jl`, after the existing Laurent certificate option checks, add:

```julia
    monomial_unit_laurent_gl = diagonal_matrix(LR, [u, one(LR), one(LR)])
    monomial_err = try
        elementary_factorization(monomial_unit_laurent_gl)
        nothing
    catch err
        err
    end
    @test monomial_err isa ArgumentError
    monomial_message = sprint(showerror, monomial_err)
    @test occursin("elementary_factorization(A) is an elementary-only SL_n API", monomial_message)
    @test occursin("laurent_gl_factorization_certificate(A)", monomial_message)
```

- [ ] **Step 5: Run focused tests and confirm failure**

Run:

```bash
julia --project=. -e 'include("test/expert/issue38_laurent_gl_certificate.jl")'
julia --project=. -e 'include("test/public/laurent_gl_certificate_options.jl"); include("test/public/api_surface.jl")'
```

Expected: FAIL before implementation because `core_factors_elementary_ok` is missing and the old Laurent boundary message is still emitted.

- [ ] **Step 6: Commit failing tests**

```bash
git add test/expert/issue38_laurent_gl_certificate.jl test/public/laurent_gl_certificate_options.jl test/public/api_surface.jl
git commit -m "Test issue 358 Laurent GL public contract"
```

---

### Task 2: Implement The Verifier And Error Contract

**Files:**
- Modify: `src/core/elementary_matrices.jl`
- Modify: `src/algorithm/laurent_gl_certificate.jl`
- Modify: `src/algorithm/factorization.jl`

**Interfaces:**
- Consumes: `_canonical_elementary_factor_record(factor)`, `_same_base_ring(left, right)`, `normalize_laurent_gl_matrix(A)`, `_factor_product(factors, R, n)`.
- Produces: `_is_elementary_matrix_factor(factor, R, n)::Bool`, eager certificate verification field `core_factors_elementary_ok::Bool`, and determinant-contract error guidance.

- [ ] **Step 1: Add the private elementary-factor predicate**

In `src/core/elementary_matrices.jl`, after `_canonical_elementary_factor_record`, add:

```julia
function _is_elementary_matrix_factor(factor, R, n::Int)::Bool
    try
        nrows(factor) == n || return false
        ncols(factor) == n || return false
        _same_base_ring(base_ring(factor), R) || return false
        record = _canonical_elementary_factor_record(factor)
        return record.kind in (:identity, :elementary)
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end
```

- [ ] **Step 2: Add public docstrings for the certificate API**

At the top of `src/algorithm/laurent_gl_certificate.jl`, before `struct LaurentGLFactorizationCertificate`, add:

```julia
"""
    LaurentGLFactorizationCertificate

Verified decomposition certificate for supported Laurent `GL_n` inputs on the
eager determinant path. The certificate stores the original matrix, determinant
profile, diagonal monomial-unit correction, determinant-one normalized core,
elementary core factors, reconstructed product, and verification status.

For supported monomial-unit determinant inputs the public reconstruction
contract is:

```julia
certificate.original_matrix ==
    certificate.correction.factor * prod(certificate.core_factors)
```

The correction is not an elementary factor; it carries the original nontrivial
determinant while every elementary factor has determinant one.
"""
```

Before `function laurent_gl_factorization_certificate`, add:

```julia
"""
    laurent_gl_factorization_certificate(A; determinant_strategy = :eager, correction_side = nothing)

Return a verified Laurent `GL_n` decomposition certificate for supported square
Laurent matrices whose determinant is `1` or a supported Laurent monomial unit.
The eager certificate normalizes `A` to a determinant-one core, factors that
core into elementary matrices, and verifies exact reconstruction of the
original input by applying the stored diagonal correction.
"""
```

Before `function verify_laurent_gl_factorization_certificate`, add:

```julia
"""
    verify_laurent_gl_factorization_certificate(certificate)::Bool

Independently replay the stored Laurent `GL_n` certificate and return whether
the determinant classification, correction, determinant-one core,
elementary-factor product, and original-input reconstruction all verify.
"""
```

- [ ] **Step 3: Strengthen eager certificate verification**

In `_laurent_gl_factorization_certificate_verification`, add a local flag:

```julia
    core_factors_elementary_ok = false
```

After `core_factors_match_ok = ...`, add:

```julia
            core_factors_elementary_ok = all(
                factor -> _is_elementary_matrix_factor(factor, R, n),
                certificate.core_factors,
            )
```

Include it in `overall_ok`:

```julia
        core_det_ok && core_replay_ok && core_factors_match_ok &&
        core_factors_elementary_ok && core_factors_ok && reconstructed_product_ok
```

Include it in the returned named tuple:

```julia
        core_factors_elementary_ok = core_factors_elementary_ok,
```

- [ ] **Step 4: Strengthen lazy certificate elementary checks**

In `_laurent_gl_lazy_deferred_correction_certificate_verification`, add:

```julia
    elementary_factors_elementary_ok = false
```

After `elementary_factors_ok = ...`, add:

```julia
            elementary_factors_elementary_ok = all(
                factor -> _is_elementary_matrix_factor(factor, R, n),
                certificate.elementary_factors,
            )
```

Include it in `overall_ok` next to `elementary_factors_ok`, and return it in
the verification named tuple.

- [ ] **Step 5: Replace the Laurent determinant-contract error**

In `_throw_staged_factorization_failure`, replace the Laurent branch with:

```julia
    if ring_profile == :laurent
        profile = normalization.determinant_profile
        throw(ArgumentError(
            "elementary_factorization(A) is an elementary-only SL_n API and " *
            "requires determinant 1. This Laurent input has determinant " *
            "$(profile.determinant) with classification $(profile.classification); " *
            "use laurent_gl_factorization_certificate(A) for the supported " *
            "Laurent GL_n monomial-unit correction certificate.",
        ))
    end
```

- [ ] **Step 6: Add public docstring for `elementary_factorization`**

In `src/algorithm/factorization.jl`, immediately before `function elementary_factorization(A)`, add:

```julia
"""
    elementary_factorization(A)

Return an elementary factor sequence whose product is exactly `A` for supported
determinant-one (`SL_n`) inputs.

This API is intentionally elementary-only: every returned factor has
determinant one, so the input must have determinant one. For supported Laurent
`GL_n` inputs with a nontrivial monomial-unit determinant, use
`laurent_gl_factorization_certificate(A)` to obtain the diagonal correction
and the elementary factorization of the determinant-one core.
"""
```

- [ ] **Step 7: Run Task 1 tests and confirm pass**

Run:

```bash
julia --project=. -e 'include("test/expert/issue38_laurent_gl_certificate.jl")'
julia --project=. -e 'include("test/public/laurent_gl_certificate_options.jl"); include("test/public/api_surface.jl")'
```

Expected: PASS.

- [ ] **Step 8: Commit implementation**

```bash
git add src/core/elementary_matrices.jl src/algorithm/laurent_gl_certificate.jl src/algorithm/factorization.jl
git commit -m "Implement Laurent GL determinant contract"
```

---

### Task 3: Update The Issue-38 Example And Public Documentation

**Files:**
- Modify: `example/toric_decoupling/issue38_unit_correction.jl`
- Modify: `README.md`
- Modify: `docs/src/index.md`
- Modify: `docs/src/toricbuilder_contract.md`

**Interfaces:**
- Consumes: final eager `LaurentGLFactorizationCertificate` fields and verifier.
- Produces: user-facing wording that explains the `SL_n`/Laurent `GL_n` split and a runnable issue-38 example showing exact original reconstruction.

- [ ] **Step 1: Update the issue-38 original error reporter**

In `example/toric_decoupling/issue38_unit_correction.jl`, rename
`_report_original_boundary` to `_report_elementary_contract_error` and replace
the message checks and prints with:

```julia
    occursin("elementary_factorization(A) is an elementary-only SL_n API", message) ||
        error("unexpected original Q factorization error: $(message)")
    occursin("laurent_gl_factorization_certificate(A)", message) ||
        error("original Q error did not point to the Laurent GL certificate API: $(message)")

    println("original_q_elementary_factorization status=DETERMINANT_CONTRACT det=$(det(Q))")
    println("original_q_elementary_factorization error=\"$(message)\"")
```

- [ ] **Step 2: Add the eager public certificate example function**

In the same example file, add this function after `_run_column_unit_correction`:

```julia
function _run_public_gl_certificate(Q)
    certificate = laurent_gl_factorization_certificate(Q)
    verified = verify_laurent_gl_factorization_certificate(certificate)
    R = base_ring(Q)
    n = nrows(Q)
    core_product = _factor_product(certificate.core_factors, R, n)

    certificate.determinant_profile.classification == :laurent_monomial_unit ||
        error("issue-38 certificate did not record a Laurent monomial-unit determinant")
    det(certificate.normalized_core) == one(R) ||
        error("certificate normalized core is not determinant one")
    core_product == certificate.normalized_core ||
        error("certificate core factors do not multiply to the normalized core")
    certificate.correction.factor * core_product == Q ||
        error("certificate correction and core product do not reconstruct Q")
    certificate.reconstructed_product == Q ||
        error("certificate reconstructed product does not equal Q")
    certificate.verification.core_factors_elementary_ok ||
        error("certificate core factors were not verified elementary")
    verified || error("eager Laurent GL certificate did not verify")

    println(
        "public_gl_certificate status=PASS det=$(certificate.determinant_profile.determinant) class=$(certificate.determinant_profile.classification) det_core=$(det(certificate.normalized_core)) factors=$(length(certificate.core_factors)) reconstructed=true verified=$(verified)",
    )
    return certificate
end
```

In `main()`, call `_report_elementary_contract_error(Q)` and
`certificate = _run_public_gl_certificate(Q)`, then change the final print to:

```julia
    println(
        "issue38_unit_correction_example status=PASS certificate_factors=$(length(certificate.core_factors)) row_factors=$(length(row.factors)) column_factors=$(length(column.factors)) lazy_factors=$(length(lazy.elementary_factors))",
    )
    return (; Q, certificate, row, column, lazy)
```

- [ ] **Step 3: Update README and Documenter scope text**

In both `README.md` and `docs/src/index.md`, replace the Laurent certificate
bullet with:

```markdown
- `elementary_factorization(A)` is an elementary-only determinant-one
  (`SL_n`) factor sequence API: every returned factor has determinant one and
  the product must equal the original input. Supported Laurent `GL_n` inputs
  with nontrivial monomial-unit determinant are therefore intentionally routed
  to `laurent_gl_factorization_certificate(A)` instead of returning a
  misleading core-only factor sequence.
- `laurent_gl_factorization_certificate(A)` is the public decomposition API
  for supported Laurent `GL_n` inputs whose determinant is `1` or a supported
  Laurent monomial unit. The eager certificate records the original matrix,
  determinant metadata, diagonal correction, inverse correction,
  determinant-one normalized core, ordered elementary core factors,
  reconstructed product, and verification status. Its public reconstruction
  relation is `A == certificate.correction.factor * prod(certificate.core_factors)`.
  The diagonal correction is not part of the elementary sequence because it
  carries the nontrivial determinant.
```

Keep the existing ordinary-polynomial bullets intact; do not claim arbitrary
Laurent `GL_n` support.

- [ ] **Step 4: Update ToricBuilder contract text**

In `docs/src/toricbuilder_contract.md`, replace the final paragraph under
`GL_n Laurent Normalization Boundary` and the `Suslin Output Contract` staged
wording with:

```markdown
For supported Laurent monomial-unit `GL_n` inputs, the public decomposition is
the verified `laurent_gl_factorization_certificate(A)` relation
`A == certificate.correction.factor * prod(certificate.core_factors)`. The
core factors are elementary factors of the determinant-one normalized core; the
diagonal correction is recorded separately because it is not elementary when
its determinant is nontrivial.
```

- [ ] **Step 5: Run example and documentation smoke targets**

Run:

```bash
julia --project=. example/toric_decoupling/issue38_unit_correction.jl
julia --project=. -e 'include("test/expert/issue38_laurent_gl_certificate.jl")'
```

Expected: PASS, with example output containing `public_gl_certificate status=PASS` and no `STAGED_BOUNDARY` line.

- [ ] **Step 6: Commit docs and example**

```bash
git add example/toric_decoupling/issue38_unit_correction.jl README.md docs/src/index.md docs/src/toricbuilder_contract.md
git commit -m "Document Laurent GL certificate contract"
```

---

### Task 4: Verification, Review, And PR

**Files:**
- Modify if needed: tests or docs touched by failures only.

**Interfaces:**
- Consumes: all previous task outputs.
- Produces: verified branch pushed to GitHub and one focused pull request with `Closes #358`.

- [ ] **Step 1: Run required issue verification**

```bash
julia --project=. example/toric_decoupling/issue38_unit_correction.jl
julia --project=. -e 'include("test/expert/issue38_laurent_gl_certificate.jl")'
```

Expected: both PASS.

- [ ] **Step 2: Run required package verification**

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: PASS for the default public/internal package tests.

- [ ] **Step 3: Inspect diff for scope**

```bash
git status --short
git diff --stat origin/main...HEAD
git diff --check origin/main...HEAD
```

Expected: only Laurent GL contract, docs, example, test, and Superpowers workflow files are changed; no whitespace errors.

- [ ] **Step 4: Commit any verification fixes**

If Step 1 or Step 2 required fixes, commit them with:

```bash
git add src/core/elementary_matrices.jl src/algorithm/laurent_gl_certificate.jl src/algorithm/factorization.jl test/expert/issue38_laurent_gl_certificate.jl test/public/laurent_gl_certificate_options.jl test/public/api_surface.jl example/toric_decoupling/issue38_unit_correction.jl README.md docs/src/index.md docs/src/toricbuilder_contract.md
git commit -m "Finalize issue 358 verification"
```

- [ ] **Step 5: Push and create PR**

```bash
git push -u origin agent/issue-358-finalize-the-public-laurent-gl-decomposition-con-run-2
gh pr create --repo nzy1997/Suslin.jl --base main --head agent/issue-358-finalize-the-public-laurent-gl-decomposition-con-run-2 --title "Finalize Laurent GL decomposition contract" --body "Closes #358"
```

Expected: PR URL is returned. Do not merge.
