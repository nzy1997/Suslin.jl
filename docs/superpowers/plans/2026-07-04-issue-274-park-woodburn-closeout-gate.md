# Issue 274 Park-Woodburn Closeout Gate Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add the final parent #187 closeout gate tying public tests, expert documentation smoke coverage, documentation, and a coverage audit into one accepted ordinary-polynomial Park-Woodburn claim.

**Architecture:** Keep factorization behavior unchanged. Add direct public-test inventory assertions over the existing final mainline acceptance catalog, add a short Markdown audit note, and extend documentation smoke coverage so README, Documenter, and the audit note all agree about the narrow #187 support boundary and #188/non-polynomial separation.

**Tech Stack:** Julia, Oscar, `Test`, Markdown documentation, existing Suslin fixture catalogs.

## Global Constraints

- Do not change factorization algorithms, fixture matrices, public APIs, or route certificate semantics.
- Supported #187 scope is exact field-backed ordinary polynomial rings only.
- Supported #187 inputs are determinant-one `SL_3` and `SL_n`, `n > 3`, through the implemented evidence-backed route only.
- The public suite must name the final #187 success cases for multivariate `SL_3`, recursive `SL_n`, and README-style public examples.
- The public suite must name the final #187 negative controls for determinant-not-one, unsupported coefficient ring, missing `SL_3` local-form evidence, missing `SL_3` Quillen evidence, missing ECP evidence, missing final `SL_3` evidence, and Laurent boundary.
- The audit note must map #181/#195, #182/#212, #183/#220, #184/#239, #185/#249, and #186/#266 to the #270/#271/#272 acceptance cases that consume each layer.
- Unsupported coefficient rings, arbitrary Laurent `GL_n`, ToricBuilder mainline acceptance, and Steinberg factor-count optimization (#188) remain separate from #187.
- README and `docs/src/index.md` must not claim Laurent/ToricBuilder support, unsupported coefficient-ring support, or Steinberg optimization as part of #187.
- In this Agent Desk sandbox, use `JULIAUP_DEPOT_PATH=/private/tmp/agentdesk-juliaup JULIA_DEPOT_PATH=/private/tmp/agentdesk-julia-depot JULIA_PKG_PRECOMPILE_AUTO=0` before Julia commands so juliaup and compiled-cache writes use writable paths.

---

### Task 1: Public Final Acceptance Inventory Gate

**Files:**
- Modify: `test/public/park_woodburn_polynomial_factorization.jl`

**Interfaces:**
- Consumes: `ParkWoodburnMainlineAcceptanceFixtureCatalog.cases_by_id()` and `ParkWoodburnMainlineAcceptanceFixtureCatalog.catalog().negative_controls`.
- Produces: public assertions that fail if a required #187 final success case or negative-control case id is removed from the catalog.

- [ ] **Step 1: Write the failing public inventory assertion call**

In `test/public/park_woodburn_polynomial_factorization.jl`, after `mainline_entries` is assigned, create `mainline_negative_entries` and call a helper that does not exist yet:

```julia
    mainline_negative_entries = Dict(
        entry.id => entry for
        entry in ParkWoodburnMainlineAcceptanceFixtureCatalog.catalog().negative_controls
    )
    _pw_assert_required_mainline_catalog_inventory(
        mainline_entries,
        mainline_negative_entries,
    )
    for entry in values(mainline_negative_entries)
        _pw_poly_assert_mainline_negative_public_failure(entry)
    end
```

Replace the existing loop:

```julia
    for entry in ParkWoodburnMainlineAcceptanceFixtureCatalog.catalog().negative_controls
        _pw_poly_assert_mainline_negative_public_failure(entry)
    end
```

- [ ] **Step 2: Run the red public command**

Run:

```bash
JULIAUP_DEPOT_PATH=/private/tmp/agentdesk-juliaup JULIA_DEPOT_PATH=/private/tmp/agentdesk-julia-depot JULIA_PKG_PRECOMPILE_AUTO=0 julia --project=. -e 'include("test/public/park_woodburn_polynomial_factorization.jl")'
```

Expected: FAIL with `UndefVarError: _pw_assert_required_mainline_catalog_inventory not defined`.

- [ ] **Step 3: Implement the inventory helper**

In `test/public/park_woodburn_polynomial_factorization.jl`, after `_pw_assert_issue187_recursive_catalog_acceptance`, add:

```julia
const REQUIRED_MAINLINE_ACCEPTANCE_CASE_IDS = (
    "pw-mainline-sl3-multivariate-issue184-qq",
    "pw-mainline-sln-recursive-issue185-186-gf2",
    "pw-mainline-readme-ordinary-polynomial-qq",
)

const REQUIRED_MAINLINE_NEGATIVE_CONTROL_IDS = (
    "pw-mainline-negative-det-not-one",
    "pw-mainline-negative-unsupported-coefficient-ring",
    "pw-mainline-negative-missing-sl3-local-form-evidence",
    "pw-mainline-negative-missing-sl3-quillen-evidence",
    "pw-mainline-negative-missing-ecp-evidence",
    "pw-mainline-negative-missing-evidence",
    "pw-mainline-negative-missing-final-sl3-evidence",
    "pw-mainline-negative-laurent-boundary",
)

function _pw_assert_required_mainline_catalog_inventory(
        mainline_entries,
        mainline_negative_entries)
    @test all(id -> haskey(mainline_entries, id), REQUIRED_MAINLINE_ACCEPTANCE_CASE_IDS)
    @test all(
        id -> haskey(mainline_negative_entries, id),
        REQUIRED_MAINLINE_NEGATIVE_CONTROL_IDS,
    )
    @test mainline_entries["pw-mainline-sl3-multivariate-issue184-qq"].entry_class ==
          :issue184_sl3_multivariate
    @test mainline_entries["pw-mainline-sln-recursive-issue185-186-gf2"].entry_class ==
          :issue185_186_sln_recursive
    @test mainline_entries["pw-mainline-readme-ordinary-polynomial-qq"].entry_class ==
          :readme_public_example
    @test mainline_negative_entries[
        "pw-mainline-negative-unsupported-coefficient-ring"
    ].negative_kind == :unsupported_coefficient_ring
    @test mainline_negative_entries["pw-mainline-negative-laurent-boundary"].negative_kind ==
          :laurent_boundary

    pruned_acceptance = filter(
        pair -> pair.first != "pw-mainline-sln-recursive-issue185-186-gf2",
        mainline_entries,
    )
    @test !all(
        id -> haskey(pruned_acceptance, id),
        REQUIRED_MAINLINE_ACCEPTANCE_CASE_IDS,
    )
    pruned_negative = filter(
        pair -> pair.first != "pw-mainline-negative-unsupported-coefficient-ring",
        mainline_negative_entries,
    )
    @test !all(
        id -> haskey(pruned_negative, id),
        REQUIRED_MAINLINE_NEGATIVE_CONTROL_IDS,
    )
    return nothing
end
```

- [ ] **Step 4: Run the green public command**

Run:

```bash
JULIAUP_DEPOT_PATH=/private/tmp/agentdesk-juliaup JULIA_DEPOT_PATH=/private/tmp/agentdesk-julia-depot JULIA_PKG_PRECOMPILE_AUTO=0 julia --project=. -e 'include("test/public/park_woodburn_polynomial_factorization.jl")'
```

Expected: PASS.

- [ ] **Step 5: Commit Task 1**

Run:

```bash
git add test/public/park_woodburn_polynomial_factorization.jl
git commit -m "test: gate Park-Woodburn final public inventory"
```

---

### Task 2: Audit Note And Documentation Smoke Gate

**Files:**
- Modify: `test/expert/documentation_smoke.jl`
- Modify: `README.md`
- Modify: `docs/src/index.md`
- Create: `docs/audits/2026-07-04-issue-187-park-woodburn-mainline-acceptance.md`

**Interfaces:**
- Consumes: README scope text, Documenter scope text, and the new audit note.
- Produces: expert smoke coverage that fails if the audit note omits an upstream gate row, omits a final acceptance/negative case id, or if public docs fold #188 or Laurent/ToricBuilder into #187.

- [ ] **Step 1: Write failing documentation smoke assertions**

In `test/expert/documentation_smoke.jl`, after `SUPPORT_BOUNDARY_EVIDENCE_PAGE`, add:

```julia
const ISSUE187_ACCEPTANCE_AUDIT_PAGE = joinpath(
    @__DIR__,
    "..",
    "..",
    "docs",
    "audits",
    "2026-07-04-issue-187-park-woodburn-mainline-acceptance.md",
)
```

After `_read_support_boundary_evidence`, add:

```julia
function _read_issue187_acceptance_audit()
    @test isfile(ISSUE187_ACCEPTANCE_AUDIT_PAGE)
    return read(ISSUE187_ACCEPTANCE_AUDIT_PAGE, String)
end
```

After `_assert_issue187_public_contract`, add:

```julia
function _assert_issue187_acceptance_audit(text)
    squashed = _squash_whitespace(text)
    for pair in (
        ("#181", "#195"),
        ("#182", "#212"),
        ("#183", "#220"),
        ("#184", "#239"),
        ("#185", "#249"),
        ("#186", "#266"),
    )
        @test occursin(pair[1], squashed)
        @test occursin(pair[2], squashed)
    end
    for issue_id in ("#270", "#271", "#272")
        @test occursin(issue_id, squashed)
    end
    for case_id in (
        "pw-mainline-sl3-multivariate-issue184-qq",
        "pw-mainline-sln-recursive-issue185-186-gf2",
        "pw-mainline-readme-ordinary-polynomial-qq",
        "pw-mainline-negative-unsupported-coefficient-ring",
        "pw-mainline-negative-missing-ecp-evidence",
        "pw-mainline-negative-missing-final-sl3-evidence",
        "pw-mainline-negative-laurent-boundary",
    )
        @test occursin(case_id, squashed)
    end
    @test occursin("Steinberg factor-count optimization (#188) remains separate", squashed)
    @test occursin("Laurent/ToricBuilder mainline support remains separate", squashed)
    @test occursin("Unsupported coefficient rings remain negative controls", squashed)
end
```

Inside the `"support boundary evidence page"` nested testset, after the existing assertions, add:

```julia
        audit = _read_issue187_acceptance_audit()
        _assert_issue187_acceptance_audit(audit)
```

Inside `_assert_issue187_public_contract`, add:

```julia
    @test occursin("#187 closeout coverage audit", squashed)
```

- [ ] **Step 2: Run the red documentation command**

Run:

```bash
JULIAUP_DEPOT_PATH=/private/tmp/agentdesk-juliaup JULIA_DEPOT_PATH=/private/tmp/agentdesk-julia-depot JULIA_PKG_PRECOMPILE_AUTO=0 julia --project=. -e 'include("test/expert/documentation_smoke.jl")'
```

Expected: FAIL because `docs/audits/2026-07-04-issue-187-park-woodburn-mainline-acceptance.md` does not exist and README/docs do not yet contain the audit pointer phrase.

- [ ] **Step 3: Create the audit note**

Create `docs/audits/2026-07-04-issue-187-park-woodburn-mainline-acceptance.md`:

````markdown
# Issue 187 Park-Woodburn Mainline Acceptance Audit

Date: 2026-07-04

Issue #274 is the parent-level closeout gate for #187. It records that the
final ordinary-polynomial Park-Woodburn public contract is accepted only for
exact field-backed ordinary-polynomial determinant-one `SL_3` and `SL_n`,
`n > 3`, inputs through the implemented evidence-backed route.

## Upstream Layer Map

| Parent layer | Upstream closeout | Evidence boundary | #187 consumers |
| --- | --- | --- | --- |
| #181 normality/conjugation | #195 | Cohn-type, rank-one, conjugated-elementary, and ECP nested-normality certificates | #270 catalog metadata for `pw-mainline-sl3-multivariate-issue184-qq`, `pw-mainline-readme-ordinary-polynomial-qq`, and `pw-mainline-sln-recursive-issue185-186-gf2` through the #184/#185/#186 route chain |
| #182 Murthy local `SL_3` | #212 | Local denominator and ordinary factor replay for supported `SL_3` slices | #271 public `SL_3` success coverage for `pw-mainline-sl3-multivariate-issue184-qq` and README-style public coverage for `pw-mainline-readme-ordinary-polynomial-qq` |
| #183 Quillen patching | #220 | Supplied-evidence global patch replay and base-term controls | #271 `SL_3` and README-style success cases; #272 `pw-mainline-negative-missing-sl3-quillen-evidence` |
| #184 `SL_3` mainline | #239 | Evidence-backed final `SL_3` route with #235/#236/#237/#183/#220 replay | #271 `pw-mainline-sl3-multivariate-issue184-qq`, `pw-mainline-readme-ordinary-polynomial-qq`, and recursive final-block evidence consumed by `pw-mainline-sln-recursive-issue185-186-gf2`; #272 `pw-mainline-negative-missing-sl3-local-form-evidence`, `pw-mainline-negative-missing-sl3-quillen-evidence`, `pw-mainline-negative-missing-evidence`, and `pw-mainline-negative-missing-final-sl3-evidence` |
| #185 ECP reducer | #249 | Verified ordinary-polynomial ECP peel certificates and public reducer dispatch | #271 recursive `pw-mainline-sln-recursive-issue185-186-gf2`; #272 `pw-mainline-negative-missing-ecp-evidence` |
| #186 recursive `SL_n` | #266 | Recursive column-peel route with #185 ECP steps, #184 final `SL_3` route, and #186 provenance | #271 recursive `pw-mainline-sln-recursive-issue185-186-gf2`; #272 recursive staged controls `pw-mainline-negative-missing-ecp-evidence`, `pw-mainline-negative-missing-evidence`, and `pw-mainline-negative-missing-final-sl3-evidence` |

## Final Acceptance Inventory

- #270 catalog entries accepted by #187:
  `pw-mainline-sl3-multivariate-issue184-qq`,
  `pw-mainline-sln-recursive-issue185-186-gf2`, and
  `pw-mainline-readme-ordinary-polynomial-qq`.
- #271 public success coverage proves these accepted entries through
  `test/public/park_woodburn_polynomial_factorization.jl` and
  `test/public/factorization_driver_shell.jl`.
- #272 negative controls remain part of the closeout gate:
  `pw-mainline-negative-det-not-one`,
  `pw-mainline-negative-unsupported-coefficient-ring`,
  `pw-mainline-negative-missing-sl3-local-form-evidence`,
  `pw-mainline-negative-missing-sl3-quillen-evidence`,
  `pw-mainline-negative-missing-ecp-evidence`,
  `pw-mainline-negative-missing-evidence`,
  `pw-mainline-negative-missing-final-sl3-evidence`, and
  `pw-mainline-negative-laurent-boundary`.

## Non-Claims

- Laurent/ToricBuilder mainline support remains separate from #187.
- Unsupported coefficient rings remain negative controls and out of scope.
- Steinberg factor-count optimization (#188) remains separate from #187.
- This audit does not add algorithmic support beyond the landed #181-#186
  ordinary-polynomial chain.

## Verification

The closeout gate is the full suite plus package test:

```bash
julia --project=. test/runtests.jl all
julia --project=. -e 'using Pkg; Pkg.test()'
```
````

- [ ] **Step 4: Add the audit pointer to README**

In the #187 bullet in `README.md`, after the sentence ending `covered by the #187 acceptance tests.`, add:

```markdown
  The #187 closeout coverage audit maps the accepted public cases and upstream
  gates in `docs/audits/2026-07-04-issue-187-park-woodburn-mainline-acceptance.md`.
```

- [ ] **Step 5: Add the audit pointer to Documenter docs**

Make the same sentence addition in the #187 bullet in `docs/src/index.md`.

- [ ] **Step 6: Run the green documentation command**

Run:

```bash
JULIAUP_DEPOT_PATH=/private/tmp/agentdesk-juliaup JULIA_DEPOT_PATH=/private/tmp/agentdesk-julia-depot JULIA_PKG_PRECOMPILE_AUTO=0 julia --project=. -e 'include("test/expert/documentation_smoke.jl")'
```

Expected: PASS.

- [ ] **Step 7: Commit Task 2**

Run:

```bash
git add README.md docs/src/index.md test/expert/documentation_smoke.jl docs/audits/2026-07-04-issue-187-park-woodburn-mainline-acceptance.md
git commit -m "docs: audit Park-Woodburn issue 187 closeout"
```

---

### Task 3: Full Verification And Pull Request

**Files:**
- No additional source edits expected.

**Interfaces:**
- Consumes: committed Task 1 and Task 2 changes.
- Produces: full-suite verification evidence, package-test evidence, final review, pushed branch, and PR against `main`.

- [ ] **Step 1: Run the issue-required full suite**

Run:

```bash
JULIAUP_DEPOT_PATH=/private/tmp/agentdesk-juliaup JULIA_DEPOT_PATH=/private/tmp/agentdesk-julia-depot JULIA_PKG_PRECOMPILE_AUTO=0 julia --project=. test/runtests.jl all
```

Expected: command exits 0.

- [ ] **Step 2: Run the issue-required package test**

Run:

```bash
JULIAUP_DEPOT_PATH=/private/tmp/agentdesk-juliaup JULIA_DEPOT_PATH=/private/tmp/agentdesk-julia-depot JULIA_PKG_PRECOMPILE_AUTO=0 julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: command exits 0.

- [ ] **Step 3: Run diff hygiene**

Run:

```bash
git diff --check HEAD~2..HEAD
```

Expected: command exits 0.

- [ ] **Step 4: Request final code review**

Dispatch a final reviewer with the branch diff and these requirements:

```text
Issue #274 closeout must add only audit/test/docs coverage. It must not change
factorization behavior, fixture matrices, public APIs, or route certificate
semantics. The public tests must name final #187 success cases and negative
controls. The audit note must map #181/#195, #182/#212, #183/#220, #184/#239,
#185/#249, and #186/#266 to #270/#271/#272 acceptance cases. README and docs
must keep unsupported coefficient rings, arbitrary Laurent GL_n, ToricBuilder
mainline acceptance, and Steinberg factor-count optimization (#188) separate
from #187.
```

Fix any Critical or Important findings and re-run the affected tests.

- [ ] **Step 5: Finish branch**

Use `superpowers:finishing-a-development-branch`, choose option 2
`Push and create a Pull Request`, push
`agent/issue-274-gate-parent-187-with-full-suite-park-woodburn-ac-run-1`, and
create a PR against `main` with `Closes #274` in the body.
