# Issue 239 SL3 Closeout Gate Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add final #184 closeout acceptance coverage and public scope wording for the evidence-backed ordinary-polynomial `SL_3` driver path.

**Architecture:** Keep the existing driver implementation narrow. Add public tests that prove `elementary_factorization` only returns factors for the replayed #235/#236/#237/#220 `SL_3` evidence path, add an expert tamper control for #236 witness replay, and update README/docs scope boundaries.

**Tech Stack:** Julia, Oscar.jl, `Test`, existing Suslin internal route-certificate helpers.

## Global Constraints

- Do not implement ECP, recursive `SL_n`, broad local-form discovery, Laurent/ToricBuilder mainline support, or factor-count optimization.
- Supported wording must be limited to exact ordinary-polynomial `3 x 3` determinant-one inputs with replayed #235 context, #236 local-form/variable-change witness, #237 ordinary Quillen local evidence, and verified #183/#220 global patch evidence.
- Staged wording must cover determinant-one `SL_3` inputs with no supported local-form, variable-change, normality/conjugation, Murthy, or Quillen evidence.
- Public `elementary_factorization(A)` must not return factors for determinant-not-one inputs or unsupported coefficient rings.
- Follow TDD: write the failing test or docs check first, run it, then implement the minimal change.

---

### Task 1: Public #184 Acceptance Gate

**Files:**
- Modify: `test/public/park_woodburn_polynomial_factorization.jl`

**Interfaces:**
- Consumes: `elementary_factorization(A)`, `verify_factorization(A, factors)`, `Suslin._polynomial_factorization_route_certificate(A)`, and existing route evidence verifiers.
- Produces: Public acceptance coverage showing the non-fixture multivariate `SL_3` path replays #235/#236/#237/#220 evidence and corrupted returned factors fail exact verification.

- [ ] **Step 1: Write the failing public acceptance assertions**

Add these assertions immediately after:

```julia
    @test issue238_factors == issue238_cert.factors
    @test Suslin._verify_polynomial_factorization_route_certificate(issue238_cert)
```

```julia
    @test issue238_cert.status == :supported
    issue238_evidence = issue238_cert.evidence
    @test issue238_evidence.context.catalog_metadata.context_issue_id == "#235"
    @test issue238_evidence.context.catalog_metadata.driver_issue_id == "#184"
    @test Suslin._verify_sl3_realization_input_context(issue238_evidence.context)
    @test issue238_evidence.witness_selection.local_form_witness.witness_issue_id == "#236"
    @test Suslin._verify_sl3_local_form_witness_selection(
        issue238_evidence.witness_selection,
    )
    @test issue238_evidence.local_evidence_provider.metadata.provider_issue_id == "#237"
    @test Suslin._verify_sl3_murthy_quillen_local_evidence_provider(
        issue238_evidence.local_evidence_provider,
    )
    @test Suslin.verify_quillen_murthy_adapter_consumption(
        issue238_evidence.quillen_consumption,
    )
    @test Suslin._verify_polynomial_quillen_patch_route_adapter(
        issue238_evidence.quillen_route_adapter,
    )
    issue238_route_metadata =
        issue238_evidence.quillen_route_adapter.quillen_patch.replay_metadata.metadata
    @test issue238_route_metadata.context_issue_id == "#235"
    @test issue238_route_metadata.witness_issue_id == "#236"
    @test issue238_route_metadata.provider_issue_id == "#237"
    @test issue238_route_metadata.patch_issue_id == "#220"

    corrupted_issue238_factors = copy(issue238_factors)
    corrupted_issue238_factors[1] =
        corrupted_issue238_factors[1] *
        elementary_matrix(3, 1, 2, one(issue238_R), issue238_R)
    @test !verify_factorization(issue238_A, corrupted_issue238_factors)

    missing_issue236_witness =
        elementary_matrix(3, 1, 3, issue238_X, issue238_R) *
        elementary_matrix(3, 2, 1, issue238_r, issue238_R)
    missing_issue236_factors, missing_issue236_err =
        _pw_acceptance_result_or_error(missing_issue236_witness)
    @test missing_issue236_factors === nothing
    @test missing_issue236_err isa ArgumentError
    @test occursin(
        "evidence-backed SL_3 polynomial route",
        sprint(showerror, missing_issue236_err),
    )
    @test occursin("#236 local-form witness", sprint(showerror, missing_issue236_err))
```

- [ ] **Step 2: Run the focused public test and verify it fails before production changes if any assertion exposes a gap**

Run:

```bash
julia --project=. -e 'include("test/public/park_woodburn_polynomial_factorization.jl")'
```

Expected before implementation: either PASS if the implementation already satisfies the closeout assertions, or FAIL on a specific missing assertion. If it passes, record that Task 1 is an acceptance-only test hardening task.

- [ ] **Step 3: Implement the minimal public-route fix only if Step 2 fails**

If Step 2 fails because the code lacks one of the required replay metadata fields or returns factors for `missing_issue236_witness`, fix the narrow route in `src/algorithm/factorization.jl` without broadening supported inputs. The expected minimal fix is to preserve the existing `PolynomialSL3QuillenMurthyRouteEvidence` metadata and ensure `_polynomial_sl3_quillen_murthy_route_error(A)` is used for multivariate `SL_3` staged diagnostics.

- [ ] **Step 4: Run the focused public test to verify green**

Run:

```bash
julia --project=. -e 'include("test/public/park_woodburn_polynomial_factorization.jl")'
```

Expected: PASS.

- [ ] **Step 5: Commit Task 1**

```bash
git add test/public/park_woodburn_polynomial_factorization.jl src/algorithm/factorization.jl
git commit -m "test: gate public sl3 closeout acceptance"
```

If `src/algorithm/factorization.jl` was not modified, omit it from `git add`.

### Task 2: Expert #236 Witness Tamper Control

**Files:**
- Modify: `test/expert/park_woodburn_route_certificate.jl`

**Interfaces:**
- Consumes: `_pw_rebuild`, `_pw_replace_certificate`, `PolynomialSL3QuillenMurthyRouteEvidence`, and `_verify_polynomial_factorization_route_certificate`.
- Produces: Expert evidence tamper coverage proving a corrupted #236 local-form witness invalidates the route evidence and enclosing route certificate.

- [ ] **Step 1: Write the failing expert tamper test**

Add this block immediately after:

```julia
    @test !Suslin._verify_polynomial_factorization_route_certificate(bad_provider_cert)
```

```julia
    tampered_issue236_witness = merge(
        sl3_cert.evidence.context.local_form_witness,
        (; monic_entry_position = (1, 2)),
    )
    tampered_issue236_context = _pw_rebuild(
        sl3_cert.evidence.context;
        local_form_witness = tampered_issue236_witness,
    )
    tampered_issue236_evidence = _pw_rebuild(
        sl3_cert.evidence;
        context = tampered_issue236_context,
    )
    tampered_issue236_cert =
        _pw_replace_certificate(sl3_cert; evidence = tampered_issue236_evidence)
    @test !Suslin._verify_sl3_realization_input_context(tampered_issue236_context)
    @test !Suslin._verify_polynomial_sl3_quillen_murthy_route_evidence(
        tampered_issue236_evidence,
    )
    @test !Suslin._verify_polynomial_factorization_route_certificate(
        tampered_issue236_cert,
    )
```

- [ ] **Step 2: Run the focused expert test and verify the new assertion catches the boundary**

Run:

```bash
julia --project=. -e 'include("test/expert/park_woodburn_route_certificate.jl")'
```

Expected before implementation: either PASS if the verifier already rejects the tamper, or FAIL if a verifier accepts corrupted witness evidence.

- [ ] **Step 3: Implement the minimal verifier fix only if Step 2 fails**

If the tamper is accepted, update `src/algorithm/factorization.jl` so `SL3RealizationInputContext` creation identity and local-form witness replay are required by `_verify_sl3_realization_input_context` and by `_polynomial_sl3_quillen_murthy_route_core_verification`. Do not add new supported witness discovery.

- [ ] **Step 4: Run the focused expert test to verify green**

Run:

```bash
julia --project=. -e 'include("test/expert/park_woodburn_route_certificate.jl")'
```

Expected: PASS.

- [ ] **Step 5: Commit Task 2**

```bash
git add test/expert/park_woodburn_route_certificate.jl src/algorithm/factorization.jl
git commit -m "test: reject tampered sl3 local-form evidence"
```

If `src/algorithm/factorization.jl` was not modified, omit it from `git add`.

### Task 3: Public Support Boundary Documentation

**Files:**
- Modify: `README.md`
- Modify: `docs/src/index.md`

**Interfaces:**
- Consumes: Existing "Current scope" / "Scope" wording.
- Produces: Public documentation that moves #184 from staged to supported only for the evidence-backed ordinary-polynomial `SL_3` path and keeps #185-#187 plus other excluded work staged/out of scope.

- [ ] **Step 1: Write the failing documentation checks**

Run these checks before editing docs:

```bash
rg -n 'general `SL_3` \\(#184\\)|#184.*staged|evidence-backed ordinary-polynomial `SL_3`' README.md docs/src/index.md
```

Expected before implementation: current docs still describe general `SL_3` (#184) as staged and do not use the final evidence-backed support wording.

- [ ] **Step 2: Replace the first `elementary_factorization(A)` scope bullet in both docs**

In `README.md`, replace the first and final `elementary_factorization(A)` scope bullets with wording equivalent to:

```markdown
- `elementary_factorization(A)` is staged but now supports the evidence-backed
  ordinary-polynomial `SL_3` driver path for exact field-backed `3 x 3`
  determinant-one inputs whose #235 checked context, #236 local-form or
  variable-change witness, #237 ordinary Quillen local evidence, and #183/#220
  global patch evidence all replay before factors are returned. This includes
  the implemented Park-Woodburn `SL_3` special-form route and preserves the
  existing univariate local `SL_3`, selected `n > 3` ordinary-polynomial
  block-local/column-peel fixtures, and #183 Quillen patch gates. The #184
  closeout does not claim arbitrary determinant-one multivariate `SL_3`
  realization unless the corresponding local-form, variable-change, or
  normality/conjugation witness is implemented and replayed.
```

Then add a separate boundary bullet:

```markdown
- Staged ordinary-polynomial `SL_3` inputs include determinant-one matrices with
  no supported local-form, variable-change, normality/conjugation, Murthy, or
  Quillen evidence path. ECP (#185), recursive `SL_n` (#186), full public
  Park-Woodburn acceptance (#187), coefficient-ring support beyond exact
  field-backed ordinary polynomial rings, arbitrary Laurent `GL_n` determinant
  correction, Laurent/ToricBuilder mainline acceptance, and Steinberg
  factor-count optimization remain out of scope.
```

Make the same conceptual replacement in `docs/src/index.md`, matching that file's "Scope" section.

- [ ] **Step 3: Run documentation wording checks**

Run:

```bash
rg -n 'evidence-backed ordinary-polynomial `SL_3`|#235 checked context|#236 local-form|#237 ordinary Quillen|#183/#220 global patch|ECP \\(#185\\)|recursive `SL_n` \\(#186\\)|full public Park-Woodburn acceptance \\(#187\\)|Laurent/ToricBuilder mainline acceptance|factor-count optimization|arbitrary determinant-one multivariate `SL_3`' README.md docs/src/index.md
```

Expected: both files contain the supported, staged, and out-of-scope wording.

- [ ] **Step 4: Run the public test after docs edits**

Run:

```bash
julia --project=. -e 'include("test/public/park_woodburn_polynomial_factorization.jl")'
```

Expected: PASS.

- [ ] **Step 5: Commit Task 3**

```bash
git add README.md docs/src/index.md
git commit -m "docs: document sl3 closeout boundary"
```

### Task 4: Full Verification and Final Review

**Files:**
- No direct edits expected unless verification finds a defect.

**Interfaces:**
- Consumes: All changes from Tasks 1-3.
- Produces: Required verification evidence for issue #239 and a reviewed branch ready for PR.

- [ ] **Step 1: Run required focused expert verification**

```bash
julia --project=. -e 'include("test/expert/park_woodburn_route_certificate.jl")'
```

Expected: PASS.

- [ ] **Step 2: Run required focused public verification**

```bash
julia --project=. -e 'include("test/public/park_woodburn_polynomial_factorization.jl")'
```

Expected: PASS.

- [ ] **Step 3: Run required full-suite verification**

```bash
julia --project=. test/runtests.jl all
```

Expected: PASS.

- [ ] **Step 4: Run required package verification**

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: PASS.

- [ ] **Step 5: Request final code review**

Use `superpowers:requesting-code-review` through the subagent-driven final review flow. Fix Critical and Important findings before finishing.

- [ ] **Step 6: Finish branch**

Use `superpowers:verification-before-completion` and `superpowers:finishing-a-development-branch`. When prompted, choose "Push and create a Pull Request" per the Agent Desk standing instruction.

## Plan Self-Review

- Spec coverage: Tasks 1 and 2 cover public/expert acceptance and tamper controls; Task 3 covers README/docs support boundary; Task 4 covers required verification and PR readiness.
- Red-flag scan: no incomplete-marker wording remains in actionable steps.
- Type consistency: helper names and evidence types match existing code in `test/public/park_woodburn_polynomial_factorization.jl`, `test/expert/park_woodburn_route_certificate.jl`, and `src/algorithm/factorization.jl`.
