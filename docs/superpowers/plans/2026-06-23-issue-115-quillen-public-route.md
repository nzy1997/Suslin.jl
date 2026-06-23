# Issue 115 Quillen Public Route Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Route the deterministic multivariate Park-Woodburn Quillen fixture through public `elementary_factorization(A)`.

**Architecture:** Add narrow fixture recognition in `src/algorithm/factorization.jl`, rebuild verified Quillen patch evidence with existing certificate constructors, and pass that patch through the issue 114 adapter. Update fixture metadata and focused expert/public tests to prove success and fail-closed behavior.

**Tech Stack:** Julia, Oscar exact polynomial matrices, existing Suslin Quillen patch verifiers, Test stdlib.

## Global Constraints

- No `AGENTS.md`, `CLAUDE.md`, `GEMINI.md`, or `CODEX.md` file is present in this checkout.
- The worker branch is `agent/issue-115-wire-multivariate-quillen-routing-into-elementar-run-1`.
- Accept only deterministic fixture-backed Quillen routing for this issue.
- Reuse issue 99 fixture shape and issue 105 constructive Quillen certificate constructors; do not implement new denominator-cover logic.
- Reuse the issue 114 `PolynomialQuillenPatchRouteAdapter`; do not bypass `verify_quillen_patch(patch)`.
- Determinant/unit precondition failures must remain distinct from determinant-one missing-witness failures.
- Do not add the final broad public acceptance file from issue 116.
- Do not broaden to arbitrary multivariate local realizability.
- Do not handle Laurent `GL_n`.
- Required focused commands: `julia --project=. -e 'include("test/expert/park_woodburn_quillen_route_adapter.jl")'`, `julia --project=. -e 'include("test/expert/park_woodburn_route_certificate.jl")'`, and `julia --project=. -e 'include("test/public/factorization_driver_shell.jl")'`.
- Required package command: `julia --project=. -e 'using Pkg; Pkg.test()'`.

---

## File Structure

- Modify `src/algorithm/factorization.jl`: add the automatic fixture-backed Quillen route, fixture patch construction helpers, and multivariate missing-witness staged error.
- Modify `test/fixtures/park_woodburn_polynomial_cases.jl`: mark the deterministic Quillen fixture as `:quillen_patch` and `:supported`.
- Modify `test/internal/park_woodburn_polynomial_fixtures.jl`: update expected Quillen catalog metadata.
- Modify `test/expert/park_woodburn_quillen_route_adapter.jl`: update automatic route expectations.
- Modify `test/expert/park_woodburn_route_certificate.jl`: add Quillen route tag and nested adapter tamper coverage.
- Modify `test/public/factorization_driver_shell.jl`: add public Quillen success and close negative control.

---

### Task 1: Add RED Coverage for Public Quillen Routing

**Files:**
- Modify: `test/fixtures/park_woodburn_polynomial_cases.jl`
- Modify: `test/internal/park_woodburn_polynomial_fixtures.jl`
- Modify: `test/expert/park_woodburn_quillen_route_adapter.jl`
- Modify: `test/expert/park_woodburn_route_certificate.jl`
- Modify: `test/public/factorization_driver_shell.jl`

**Interfaces:**
- Consumes existing `ParkWoodburnPolynomialFixtureCatalog`.
- Consumes existing `PolynomialQuillenPatchRouteAdapter`.
- Produces failing expectations for automatic `:quillen_patch` routing before production code is changed.

- [ ] **Step 1: Update Quillen fixture metadata**

Change the `quillen-patched-substitution-witness-qq` Park-Woodburn fixture to:

```julia
route = :quillen_patch,
status = :supported,
```

In the internal validator, update the `:multivariate_quillen` role and expected
entry metadata to require `(:quillen_patch, :supported)`. Keep the referenced
issue 99 target-matrix equality check.

- [ ] **Step 2: Update adapter route expectation**

In `test/expert/park_woodburn_quillen_route_adapter.jl`, change the automatic
route assertion for the Quillen fixture from staged failure to:

```julia
auto_cert = Suslin._polynomial_factorization_route_certificate(route_entry.matrix)
@test auto_cert.route == :quillen_patch
@test auto_cert.status == :supported
@test auto_cert.evidence isa Suslin.PolynomialQuillenPatchRouteAdapter
@test verify_factorization(route_entry.matrix, auto_cert.factors)
@test Suslin._verify_polynomial_factorization_route_certificate(auto_cert)
```

- [ ] **Step 3: Add route-certificate Quillen assertions**

In `test/expert/park_woodburn_route_certificate.jl`, add a helper that replaces
Quillen adapter fields:

```julia
function _pw_replace_quillen_adapter(
        adapter;
        target = adapter.target,
        route = adapter.route,
        quillen_patch = adapter.quillen_patch,
        global_elementary_factors = adapter.global_elementary_factors,
        product = adapter.product,
        target_matrix = adapter.target_matrix,
        replay_metadata = adapter.replay_metadata,
        verification = adapter.verification)
    return Suslin.PolynomialQuillenPatchRouteAdapter(
        target,
        route,
        quillen_patch,
        global_elementary_factors,
        product,
        target_matrix,
        replay_metadata,
        verification,
    )
end
```

Add a route test for the Quillen entry:

```julia
quillen_entry = entries["quillen-patched-substitution-witness-qq"]
quillen_cert = Suslin._polynomial_factorization_route_certificate(quillen_entry.matrix)
@test quillen_cert.route == :quillen_patch
@test quillen_cert.status == :supported
@test quillen_cert.evidence isa Suslin.PolynomialQuillenPatchRouteAdapter
@test Suslin._verify_polynomial_quillen_patch_route_adapter(quillen_cert.evidence)
@test verify_factorization(quillen_entry.matrix, quillen_cert.factors)
@test Suslin._verify_polynomial_factorization_route_certificate(quillen_cert)

bad_quillen_factors = copy(quillen_cert.evidence.global_elementary_factors)
bad_quillen_factors[1] =
    bad_quillen_factors[1] *
    elementary_matrix(nrows(quillen_entry.matrix), 1, 3, one(base_ring(quillen_entry.matrix)), base_ring(quillen_entry.matrix))
bad_quillen_evidence = _pw_replace_quillen_adapter(
    quillen_cert.evidence;
    global_elementary_factors = bad_quillen_factors,
)
bad_quillen_cert = _pw_replace_certificate(quillen_cert; evidence = bad_quillen_evidence)
@test verify_factorization(bad_quillen_cert.matrix, bad_quillen_cert.factors)
@test !Suslin._verify_polynomial_factorization_route_certificate(bad_quillen_cert)
```

- [ ] **Step 4: Add public success and negative control**

In `test/public/factorization_driver_shell.jl`, after the recursive column-peel
catalog checks, add:

```julia
quillen_supported = entries["quillen-patched-substitution-witness-qq"].matrix
quillen_factors = elementary_factorization(quillen_supported)
@test verify_factorization(quillen_supported, quillen_factors)
quillen_cert = Suslin._polynomial_factorization_route_certificate(quillen_supported)
@test quillen_cert.route == :quillen_patch
@test quillen_factors == quillen_cert.factors
@test quillen_cert.evidence isa Suslin.PolynomialQuillenPatchRouteAdapter
@test Suslin._verify_polynomial_quillen_patch_route_adapter(quillen_cert.evidence)

S = base_ring(quillen_supported)
SX, Sr, Sg = collect(gens(S))
quillen_unsupported = elementary_matrix(
    3,
    1,
    2,
    SX + Sr^2 * Sg + Sg + Sr + one(S),
    S,
)
quillen_unsupported_err =
    _captured_error(() -> elementary_factorization(quillen_unsupported))
@test quillen_unsupported_err isa ArgumentError
@test occursin("missing Quillen/local realizability witness", sprint(showerror, quillen_unsupported_err))
```

Update the existing multivariate determinant-one staged assertion to expect the
same missing-witness message.

- [ ] **Step 5: Run focused RED verification**

Run:

```bash
julia --project=. -e 'include("test/expert/park_woodburn_quillen_route_adapter.jl")'
julia --project=. -e 'include("test/expert/park_woodburn_route_certificate.jl")'
julia --project=. -e 'include("test/public/factorization_driver_shell.jl")'
```

Expected: FAIL because automatic route selection still returns staged failure
for the Quillen fixture.

---

### Task 2: Implement Narrow Quillen Fixture Routing

**Files:**
- Modify: `src/algorithm/factorization.jl`

**Interfaces:**
- Produces automatic `:quillen_patch` route certificates for the supported fixture.
- Produces a staged missing-witness `ArgumentError` for unsupported determinant-one multivariate matrices.

- [ ] **Step 1: Add fixture route selection**

In `_polynomial_factorization_route_certificate(A)`, after the disjoint
local-block attempt and before recursive column peel, add:

```julia
quillen_certificate = _polynomial_quillen_fixture_route_certificate(A)
quillen_certificate !== nothing && return quillen_certificate
```

In `_polynomial_staged_failure_evidence(A; ...)`, return `error_type = :none`
when `_polynomial_quillen_fixture_route_certificate(A)` succeeds.

- [ ] **Step 2: Add fixture patch construction helpers**

Add helpers that recognize exactly `E_12(X + r^2*g + g + 1)` over `QQ[X, r, g]`,
construct two local certificates over denominators `r` and `1 - r`, build the
cover with multipliers `1, 1`, normalize the local contributions, assemble the
deterministic patch, and verify it through `_polynomial_quillen_patch_route_certificate`.

- [ ] **Step 3: Add staged missing-witness message**

Update `_throw_staged_factorization_failure(A, :polynomial, nothing)` for
multivariate ordinary-polynomial inputs to throw:

```julia
ArgumentError("determinant-one polynomial input is outside the implemented fixture-backed Quillen route: missing Quillen/local realizability witness")
```

Keep the determinant/unit precondition message unchanged in
`_require_polynomial_sl_determinant(A)`.

- [ ] **Step 4: Run focused GREEN verification**

Run:

```bash
julia --project=. -e 'include("test/expert/park_woodburn_quillen_route_adapter.jl")'
julia --project=. -e 'include("test/expert/park_woodburn_route_certificate.jl")'
julia --project=. -e 'include("test/public/factorization_driver_shell.jl")'
```

Expected: all commands exit 0.

---

### Task 3: Full Verification and PR Preparation

**Files:**
- Modify: files from Tasks 1-2.

**Interfaces:**
- Produces a committed branch ready to push and open as a PR.

- [ ] **Step 1: Run issue-required focused commands**

Run:

```bash
julia --project=. -e 'include("test/expert/park_woodburn_quillen_route_adapter.jl")'
julia --project=. -e 'include("test/expert/park_woodburn_route_certificate.jl")'
julia --project=. -e 'include("test/public/factorization_driver_shell.jl")'
```

Expected: all exit 0.

- [ ] **Step 2: Run package verification**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: exits 0.

- [ ] **Step 3: Run diff hygiene check**

Run:

```bash
git diff --check origin/main..HEAD
```

Expected: no output, exit 0.

- [ ] **Step 4: Commit implementation**

Run:

```bash
git add docs/superpowers/specs/2026-06-23-issue-115-quillen-public-route-design.md docs/superpowers/plans/2026-06-23-issue-115-quillen-public-route.md src/algorithm/factorization.jl test/fixtures/park_woodburn_polynomial_cases.jl test/internal/park_woodburn_polynomial_fixtures.jl test/expert/park_woodburn_quillen_route_adapter.jl test/expert/park_woodburn_route_certificate.jl test/public/factorization_driver_shell.jl
git commit -m "feat: route quillen fixture through factorization"
```

Expected: commit succeeds.

---

## Plan Self-Review

- Spec coverage: the plan covers automatic public Quillen dispatch, fixture metadata, nested route evidence replay, public success, negative control, and staged missing-witness errors.
- Placeholder scan: no incomplete markers remain.
- Type consistency: route tag `:quillen_patch`, adapter type names, and certificate helper names match existing issue 114 code.
- Scope check: the plan accepts only one deterministic multivariate fixture and leaves broad public acceptance to issue 116.

## Execution Choice

Under the Standing Answer Policy, choose **Subagent-Driven (recommended)** from
the writing-plans execution options. Use `superpowers:subagent-driven-development`
to execute this plan task-by-task.
