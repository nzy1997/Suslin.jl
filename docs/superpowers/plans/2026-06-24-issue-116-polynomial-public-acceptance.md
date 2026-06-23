# Issue 116 Polynomial Public Acceptance Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add the final public Park-Woodburn ordinary-polynomial factorization acceptance command and document the supported driver scope.

**Architecture:** Add a catalog-backed public acceptance file that calls `elementary_factorization(A)` for the supported route families and asserts exact verification. Register that file in the public test group and update README/docs with conservative route support and staged boundaries.

**Tech Stack:** Julia, Oscar exact polynomial matrices, Suslin route certificates for test evidence, Test stdlib, Documenter markdown.

## Global Constraints

- No `AGENTS.md` file is present in this checkout.
- The worker branch is `agent/issue-116-add-final-public-park-woodburn-polynomial-factor-run-1`.
- Keep the public API unchanged; tests call `elementary_factorization(A)` and `verify_factorization(A, factors)`.
- Include one univariate `SL_3` fast-local example, one `n > 3` recursive polynomial column-peel example, and one deterministic multivariate Quillen example.
- Negative controls must include determinant not equal to `1` and determinant-one outside implemented witness families; both must throw staged `ArgumentError`s and return no factors.
- Do not handle Laurent `GL_n` determinant corrections.
- Do not optimize factor count.
- Do not claim arbitrary coefficient-ring support beyond supported exact field-backed ordinary polynomial rings.
- Required focused command: `julia --project=. -e 'include("test/public/park_woodburn_polynomial_factorization.jl")'`.
- Required full command: `julia --project=. test/runtests.jl all`.
- Required package command: `julia --project=. -e 'using Pkg; Pkg.test()'`.

---

## File Structure

- Create `test/public/park_woodburn_polynomial_factorization.jl`: final public acceptance tests for supported routes and negative controls.
- Modify `test/runtests.jl`: register the acceptance file in the `public` group.
- Modify `README.md`: update public scope text.
- Modify `docs/src/index.md`: update documentation scope text.
- Keep this design and plan in `docs/superpowers/` because prior issue branches in this repository commit Superpowers specs and plans.

---

### Task 1: Public Acceptance Test and Registration

**Files:**
- Create: `test/public/park_woodburn_polynomial_factorization.jl`
- Modify: `test/runtests.jl`

**Interfaces:**
- Consumes `ParkWoodburnPolynomialFixtureCatalog.cases_by_id()` and `.catalog().negative_controls`.
- Produces the final public acceptance command required by issue 64 and issue 116.

- [ ] **Step 1: Write the acceptance test**

Create `test/public/park_woodburn_polynomial_factorization.jl` with helper
functions equivalent to:

```julia
using Suslin
using Test
using Oscar

const PARK_WOODBURN_ACCEPTANCE_CATALOG_PATH =
    joinpath(@__DIR__, "..", "fixtures", "park_woodburn_polynomial_cases.jl")

function _pw_acceptance_result_or_error(A)
    factors = nothing
    try
        factors = elementary_factorization(A)
        return factors, nothing
    catch err
        return factors, err
    end
end
```

Load the catalog if needed. Test:

```julia
fast_local = entries["pw-poly-univariate-sl3-fast-local-qq"].matrix
fast_factors, fast_err = _pw_acceptance_result_or_error(fast_local)
@test fast_err === nothing
@test fast_factors !== nothing
@test verify_factorization(fast_local, fast_factors)
fast_cert = Suslin._polynomial_factorization_route_certificate(fast_local)
@test fast_cert.route == :fast_local_sl3
@test fast_factors == fast_cert.factors

recursive = entries["pw-poly-recursive-column-peel-sln-block-qq"].matrix
recursive_factors, recursive_err = _pw_acceptance_result_or_error(recursive)
@test recursive_err === nothing
@test recursive_factors !== nothing
@test verify_factorization(recursive, recursive_factors)
recursive_cert = Suslin._polynomial_factorization_route_certificate(recursive)
@test recursive_cert.route == :polynomial_column_peel
@test recursive_cert.evidence isa Suslin.PolynomialColumnPeelCertificate
@test recursive_cert.evidence.final_certificate.route == :disjoint_local_blocks

quillen = entries["quillen-patched-substitution-witness-qq"].matrix
quillen_factors, quillen_err = _pw_acceptance_result_or_error(quillen)
@test quillen_err === nothing
@test quillen_factors !== nothing
@test verify_factorization(quillen, quillen_factors)
quillen_cert = Suslin._polynomial_factorization_route_certificate(quillen)
@test quillen_cert.route == :quillen_patch
@test quillen_cert.evidence isa Suslin.PolynomialQuillenPatchRouteAdapter
@test Suslin._verify_polynomial_quillen_patch_route_adapter(quillen_cert.evidence)
```

Add negative controls:

```julia
negative_entries = Dict(entry.id => entry for entry in catalog.negative_controls)
det_factors, det_err =
    _pw_acceptance_result_or_error(negative_entries["pw-poly-det-not-one-control"].matrix)
@test det_factors === nothing
@test det_err isa ArgumentError
@test occursin("determinant/unit precondition", sprint(showerror, det_err))

outside_factors, outside_err =
    _pw_acceptance_result_or_error(negative_entries["pw-poly-det-one-outside-witness-control"].matrix)
@test outside_factors === nothing
@test outside_err isa ArgumentError
@test occursin("staged reduction to the supported univariate local SL_3 slice", sprint(showerror, outside_err))
```

- [ ] **Step 2: Run RED check before registration**

Run:

```bash
julia --project=. -e 'include("test/public/park_woodburn_polynomial_factorization.jl")'
```

Expected before the file exists: FAIL with `No such file or directory`. This
has already been observed in this Agent Desk run and proves the public
acceptance command was absent.

- [ ] **Step 3: Register the test**

Add `"public/park_woodburn_polynomial_factorization.jl"` to the `public` group
in `test/runtests.jl`, after `public/factorization_driver_shell.jl`.

- [ ] **Step 4: Run focused GREEN check**

Run:

```bash
julia --project=. -e 'include("test/public/park_woodburn_polynomial_factorization.jl")'
```

Expected: exit 0.

---

### Task 2: Documentation Scope Update

**Files:**
- Modify: `README.md`
- Modify: `docs/src/index.md`

**Interfaces:**
- Consumes the issue 64/116 scope language.
- Produces user-facing documentation for supported and staged polynomial driver routes.

- [ ] **Step 1: Update README scope**

Replace the first and fourth `Current scope` bullets with language that says
`elementary_factorization(A)` supports:

- univariate local `SL_3` ordinary-polynomial matrices;
- selected `n > 3` ordinary-polynomial matrices through block-local reduction
  and recursive polynomial column peel;
- deterministic multivariate Quillen/local-to-global fixture-backed matrices;
- determinant-one Laurent inputs only through the existing Laurent SL path.

Keep the staged boundary language for arbitrary `SL_n(k[x_1, ..., x_m])`,
general Quillen local realizability, coefficient-ring support beyond exact
field-backed ordinary polynomial rings, Laurent `GL_n` determinant correction,
and factor-count optimization.

- [ ] **Step 2: Update docs index scope**

Mirror the README scope in `docs/src/index.md` with concise bullets under
`## Scope`.

- [ ] **Step 3: Run documentation-adjacent public command**

Run:

```bash
julia --project=. -e 'include("test/public/park_woodburn_polynomial_factorization.jl")'
```

Expected: exit 0, proving the docs do not advertise unsupported behavior beyond
the tested acceptance routes.

---

### Task 3: Verification, Review, and PR Preparation

**Files:**
- All files from Tasks 1 and 2.

**Interfaces:**
- Produces a verified branch ready for a draft PR.

- [ ] **Step 1: Run issue-required focused command**

Run:

```bash
julia --project=. -e 'include("test/public/park_woodburn_polynomial_factorization.jl")'
```

Expected: exit 0.

- [ ] **Step 2: Run issue-required full suite**

Run:

```bash
julia --project=. test/runtests.jl all
```

Expected: exit 0.

- [ ] **Step 3: Run package verification**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: exit 0.

- [ ] **Step 4: Run diff hygiene check**

Run:

```bash
git diff --check origin/main..HEAD
```

Expected: no output, exit 0.

- [ ] **Step 5: Commit implementation**

Run:

```bash
git add docs/superpowers/specs/2026-06-24-issue-116-polynomial-public-acceptance-design.md docs/superpowers/plans/2026-06-24-issue-116-polynomial-public-acceptance.md test/public/park_woodburn_polynomial_factorization.jl test/runtests.jl README.md docs/src/index.md
git commit -m "test: add polynomial factorization public acceptance"
```

Expected: commit succeeds.

---

## Plan Self-Review

- Spec coverage: the plan covers the focused acceptance command, public group
  registration, README/docs scope, focused verification, full suite, package
  verification, and PR preparation.
- Placeholder scan: no incomplete markers remain.
- Type consistency: route tags and certificate type names match existing code.
- Scope check: the plan does not implement new algorithm routes; it only proves
  and documents the conservative routes already wired by prior issues.

## Execution Choice

Under the Standing Answer Policy, choose **Subagent-Driven (recommended)** from
the writing-plans execution options. Use `superpowers:subagent-driven-development`
to execute this plan task-by-task.
