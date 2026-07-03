# Issue 271 Final Public Mainline Success Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Prove the final public Park-Woodburn ordinary-polynomial success cases through public APIs and route certificates.

**Architecture:** This is a test-only gate. Public tests consume the #270 final acceptance catalog and call `elementary_factorization(A)` plus `verify_factorization(A, factors)` on the catalog matrices; expert tests add nested certificate corruption checks that preserve stored products but invalidate evidence.

**Tech Stack:** Julia, Oscar, Suslin.jl test fixtures, `Test`.

## Global Constraints

- Do not add new mathematical support beyond the landed #184/#185/#186 routes.
- Do not claim Laurent/ToricBuilder support, unsupported coefficient rings, or Steinberg factor-count optimization.
- Use catalog ids only to select examples; the public route must be derived from the input matrix and replayed evidence.
- Keep older fast-local, legacy Quillen, and disjoint-block examples as regression coverage, but they must not count as the new #187 recursive mainline success.
- Required verification commands:
  `julia --project=. -e 'include("test/public/park_woodburn_polynomial_factorization.jl")'`,
  `julia --project=. -e 'include("test/expert/park_woodburn_route_certificate.jl")'`,
  and `julia --project=. -e 'using Pkg; Pkg.test()'`.

---

## File Structure

- Modify `test/public/park_woodburn_polynomial_factorization.jl` to include the #270 catalog and assert final public success for the three positive catalog entries.
- Modify `test/public/factorization_driver_shell.jl` to include one compact #270 catalog shell proof for the recursive `SL_n` public route.
- Modify `test/expert/park_woodburn_route_certificate.jl` to include the #270 catalog and reject nested evidence corruption with stored product restored.
- No production source files are expected to change.

### Task 1: Public #270 Catalog Success Cases

**Files:**
- Modify: `test/public/park_woodburn_polynomial_factorization.jl`

**Interfaces:**
- Consumes: `ParkWoodburnMainlineAcceptanceFixtureCatalog.cases_by_id()`
- Produces: public tests proving `elementary_factorization` and route certificates for #270 positives

- [ ] **Step 1: Write the failing public catalog assertions**

Add the #270 catalog path near the existing catalog constants:

```julia
const PARK_WOODBURN_MAINLINE_ACCEPTANCE_CATALOG_PATH =
    joinpath(@__DIR__, "..", "fixtures", "park_woodburn_mainline_acceptance_cases.jl")
```

Add helpers before the testset:

```julia
function _pw_assert_issue184_sl3_public_acceptance(entry)
    A = entry.matrix
    @test nrows(A) == 3
    @test length(collect(gens(base_ring(A)))) > 1
    factors, err = _pw_acceptance_result_or_error(A)
    @test err === nothing
    @test factors !== nothing
    @test verify_factorization(A, factors)
    cert = Suslin._polynomial_factorization_route_certificate(A)
    @test cert.route == :quillen_patch
    @test cert.status == :supported
    @test factors == cert.factors
    @test Suslin._verify_polynomial_factorization_route_certificate(cert)
    @test cert.evidence isa Suslin.PolynomialSL3QuillenMurthyRouteEvidence ||
          cert.evidence isa Suslin.PolynomialSL3SuppliedQuillenRouteEvidence
    @test entry.public_route.issue_id == "#187"
    @test "#184" in entry.upstream_issue_ids
    return cert
end

function _pw_assert_readme_public_acceptance(entry)
    A = entry.matrix
    factors, err = _pw_acceptance_result_or_error(A)
    @test err === nothing
    @test factors !== nothing
    @test verify_factorization(A, factors)
    cert = Suslin._polynomial_factorization_route_certificate(A)
    @test cert.status == :supported
    @test factors == cert.factors
    @test Suslin._verify_polynomial_factorization_route_certificate(cert)
    @test entry.entry_class == :readme_public_example
    @test entry.public_route.issue_id == "#187"
    return cert
end

function _pw_assert_issue187_recursive_catalog_acceptance(entry, expected_step_count::Int)
    cert = _pw_assert_public_issue186_recursive_acceptance(entry.matrix, expected_step_count)
    @test entry.entry_class == :issue185_186_sln_recursive
    @test entry.public_route.route_marker == :issue186_recursive_mainline
    @test entry.public_route.issue_id == "#187"
    @test "#184" in entry.upstream_issue_ids
    @test "#185" in entry.upstream_issue_ids
    @test "#186" in entry.upstream_issue_ids
    @test hasproperty(entry.upstream_evidence, :ecp_case_id)
    @test hasproperty(entry.upstream_evidence, :final_sl3_case_id)
    return cert
end
```

Inside the testset, include the #270 catalog and select the three supported entries:

```julia
if !isdefined(Main, :ParkWoodburnMainlineAcceptanceFixtureCatalog)
    include(PARK_WOODBURN_MAINLINE_ACCEPTANCE_CATALOG_PATH)
end
mainline_entries = ParkWoodburnMainlineAcceptanceFixtureCatalog.cases_by_id()

sl3_mainline =
    mainline_entries["pw-mainline-sl3-multivariate-issue184-qq"]
_pw_assert_issue184_sl3_public_acceptance(sl3_mainline)

recursive_mainline =
    mainline_entries["pw-mainline-sln-recursive-issue185-186-gf2"]
_pw_assert_issue187_recursive_catalog_acceptance(recursive_mainline, 1)

readme_mainline =
    mainline_entries["pw-mainline-readme-ordinary-polynomial-qq"]
_pw_assert_readme_public_acceptance(readme_mainline)
```

- [ ] **Step 2: Run the focused public test**

Run:

```bash
julia --project=. -e 'include("test/public/park_woodburn_polynomial_factorization.jl")'
```

Expected: exit 0 after the helper code is integrated.

- [ ] **Step 3: Confirm the #270 cases are not replaced by legacy regressions**

Search the changed public test and confirm that the new success assertions use
`ParkWoodburnMainlineAcceptanceFixtureCatalog.cases_by_id()` entries named
`pw-mainline-sl3-multivariate-issue184-qq`,
`pw-mainline-sln-recursive-issue185-186-gf2`, and
`pw-mainline-readme-ordinary-polynomial-qq`.

Run:

```bash
rg -n "pw-mainline-(sl3-multivariate|sln-recursive|readme)" test/public/park_woodburn_polynomial_factorization.jl
```

Expected: all three ids are present.

- [ ] **Step 4: Run the focused public test again**

Run:

```bash
julia --project=. -e 'include("test/public/park_woodburn_polynomial_factorization.jl")'
```

Expected: exit 0.

- [ ] **Step 5: Commit Task 1**

```bash
git add test/public/park_woodburn_polynomial_factorization.jl
git commit -m "test: prove final public mainline catalog cases"
```

### Task 2: Public Driver Shell Catalog Proof

**Files:**
- Modify: `test/public/factorization_driver_shell.jl`

**Interfaces:**
- Consumes: `ParkWoodburnMainlineAcceptanceFixtureCatalog.cases_by_id()`
- Produces: a public driver-shell check that the #270 recursive case reaches the #186 route

- [ ] **Step 1: Add the #270 catalog include**

Add near the existing constants:

```julia
const PARK_WOODBURN_MAINLINE_ACCEPTANCE_CATALOG_PATH =
    joinpath(@__DIR__, "..", "fixtures", "park_woodburn_mainline_acceptance_cases.jl")
```

Inside the testset after the existing catalog includes:

```julia
if !isdefined(Main, :ParkWoodburnMainlineAcceptanceFixtureCatalog)
    include(PARK_WOODBURN_MAINLINE_ACCEPTANCE_CATALOG_PATH)
end
mainline_entries = ParkWoodburnMainlineAcceptanceFixtureCatalog.cases_by_id()
```

Add the recursive catalog shell assertion:

```julia
recursive_catalog = mainline_entries["pw-mainline-sln-recursive-issue185-186-gf2"]
recursive_catalog_cert = _assert_public_issue186_recursive_route(recursive_catalog.matrix, 1)
@test recursive_catalog.public_route.issue_id == "#187"
@test recursive_catalog.public_route.route_marker == :issue186_recursive_mainline
@test recursive_catalog_cert.evidence.mainline_support_metadata.issue_id == "#186"
```

- [ ] **Step 2: Run the public shell test**

Run:

```bash
julia --project=. -e 'include("test/public/factorization_driver_shell.jl")'
```

Expected: exit 0.

- [ ] **Step 3: Commit Task 2**

```bash
git add test/public/factorization_driver_shell.jl
git commit -m "test: cover final catalog in public driver shell"
```

### Task 3: Expert Nested Evidence Corruption Gate

**Files:**
- Modify: `test/expert/park_woodburn_route_certificate.jl`

**Interfaces:**
- Consumes: `ParkWoodburnMainlineAcceptanceFixtureCatalog.cases_by_id()`
- Produces: expert tests proving nested evidence corruption is rejected even with stored product restored

- [ ] **Step 1: Include the #270 catalog**

Add a constant near the other fixture constants:

```julia
const PARK_WOODBURN_MAINLINE_ACCEPTANCE_CATALOG_PATH =
    joinpath(@__DIR__, "..", "fixtures", "park_woodburn_mainline_acceptance_cases.jl")
```

Inside the testset:

```julia
if !isdefined(Main, :ParkWoodburnMainlineAcceptanceFixtureCatalog)
    include(PARK_WOODBURN_MAINLINE_ACCEPTANCE_CATALOG_PATH)
end
mainline_entries = ParkWoodburnMainlineAcceptanceFixtureCatalog.cases_by_id()
```

- [ ] **Step 2: Add corruption assertions**

After the existing `auto_peel_cert` checks, add:

```julia
mainline_acceptance_entry =
    mainline_entries["pw-mainline-sln-recursive-issue185-186-gf2"]
mainline_acceptance_cert = Suslin._polynomial_factorization_route_certificate(
    mainline_acceptance_entry.matrix,
)
@test mainline_acceptance_cert.route == :polynomial_column_peel
@test verify_factorization(
    mainline_acceptance_cert.matrix,
    mainline_acceptance_cert.factors,
)
@test Suslin._verify_polynomial_factorization_route_certificate(
    mainline_acceptance_cert,
)

bad_mainline_metadata = merge(
    mainline_acceptance_cert.evidence.mainline_support_metadata,
    (; final_route_issue184_ok = false),
)
bad_mainline_evidence = _pw_replace_column_peel_certificate(
    mainline_acceptance_cert.evidence;
    mainline_support_metadata = bad_mainline_metadata,
    product = mainline_acceptance_cert.evidence.product,
)
bad_mainline_cert = _pw_replace_certificate(
    mainline_acceptance_cert;
    evidence = bad_mainline_evidence,
    product = mainline_acceptance_cert.product,
)
@test verify_factorization(bad_mainline_cert.matrix, bad_mainline_cert.factors)
@test bad_mainline_cert.product == mainline_acceptance_cert.matrix
@test !Suslin._verify_polynomial_column_peel_certificate(bad_mainline_cert.evidence)
@test !Suslin._verify_polynomial_factorization_route_certificate(bad_mainline_cert)
```

This mutates nested evidence while preserving the stored product and factor product.

- [ ] **Step 3: Run the focused expert test**

Run:

```bash
julia --project=. -e 'include("test/expert/park_woodburn_route_certificate.jl")'
```

Expected: exit 0.

- [ ] **Step 4: Commit Task 3**

```bash
git add test/expert/park_woodburn_route_certificate.jl
git commit -m "test: reject corrupted final mainline route evidence"
```

### Task 4: Verification, Full Test, and PR

**Files:**
- Verify all modified files
- No source edits expected

**Interfaces:**
- Consumes: commits from Tasks 1-3
- Produces: pushed worker branch and pull request

- [ ] **Step 1: Run issue-required commands**

Run:

```bash
julia --project=. -e 'include("test/public/park_woodburn_polynomial_factorization.jl")'
julia --project=. -e 'include("test/expert/park_woodburn_route_certificate.jl")'
```

Expected: both exit 0.

- [ ] **Step 2: Run package tests**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: exit 0.

- [ ] **Step 3: Inspect final diff**

Run:

```bash
git status --short
git log --oneline --max-count=6
```

Expected: clean worktree after commits, recent commits include design and test commits.

- [ ] **Step 4: Push and create PR**

Run:

```bash
git push -u origin agent/issue-271-prove-final-public-success-cases-for-the-park-wo-run-1
gh pr create --base main --head agent/issue-271-prove-final-public-success-cases-for-the-park-wo-run-1 --title "Prove final Park-Woodburn public success cases" --body "Closes #271"
```

Expected: branch pushed and PR URL printed.

## Self-Review

Spec coverage: Task 1 covers the three #270 positive public cases; Task 2
covers the public driver shell route; Task 3 covers nested evidence corruption
with stored product restored; Task 4 covers all required verification and PR
creation.

Red-flag scan: no unfinished markers or open-ended validation steps remain.

Type consistency: all helpers consume existing `NamedTuple` fixture entries and
existing Suslin certificate verifier APIs; no new production interfaces are
introduced.
