# Issue 270 Mainline Acceptance Catalog Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a validated final #187 Park-Woodburn ordinary-polynomial mainline acceptance catalog without changing public factorization behavior.

**Architecture:** Add one focused fixture module for final acceptance entries and one internal validator that imports sibling evidence catalogs, validates metadata/schema/ring/determinant/source refs, and proves the negative controls fail. Register the validator in `test/runtests.jl` so `Pkg.test()` covers it.

**Tech Stack:** Julia, Oscar, `Test`, existing Suslin fixture catalogs for #184 `SL_3`, #185 ECP, and #186 recursive `SL_n`.

## Global Constraints

- Do not change public factorization behavior, implement ECP, implement recursive `SL_n`, broaden Laurent/ToricBuilder support, or optimize factor counts.
- The final catalog must be new: `test/fixtures/park_woodburn_mainline_acceptance_cases.jl`.
- The validator must be new: `test/internal/park_woodburn_mainline_acceptance_fixtures.jl`.
- Positive supported cases must use exact field-backed ordinary polynomial rings and determinant-one matrices.
- Required entries: evidence-backed multivariate `SL_3` consuming #184, `SL_n`, `n > 3`, consuming #185/#186, README-style ordinary-polynomial example, determinant-not-one negative control, unsupported coefficient-ring negative control, and missing-evidence negative control.
- Source refs must cover `refs/arXiv-alg-geom9405003v1` Sections 2, 3, 4, and 5 across the catalog, and each case must include the section for the layer it exercises.
- A route that claims #187 support must include required upstream evidence ids: #184 for `SL_3`, #185 and #186 for recursive `SL_n`, and parent coverage metadata for #181/#182/#183/#184/#185/#186.
- The validator must reject determinant-one metadata set to false for a supported entry, #187 support without required #184/#185/#186 evidence ids, unsupported coefficient rings, and missing Park-Woodburn layer source refs.

---

### Task 1: Final Acceptance Catalog And Validator

**Files:**
- Create: `test/fixtures/park_woodburn_mainline_acceptance_cases.jl`
- Create: `test/internal/park_woodburn_mainline_acceptance_fixtures.jl`
- Modify: `test/runtests.jl`

**Interfaces:**
- Consumes: `ParkWoodburnSL3DriverFixtureCatalog.cases_by_id()`, `ECPMainlineFixtureCatalog.cases_by_id()`, `ParkWoodburnSLnDriverFixtureCatalog.cases_by_id()`, Oscar ring/matrix predicates.
- Produces: `ParkWoodburnMainlineAcceptanceFixtureCatalog.catalog()`, `ParkWoodburnMainlineAcceptanceFixtureCatalog.cases_by_id()`, `validate_park_woodburn_mainline_acceptance_fixture(entry)`, and `validate_park_woodburn_mainline_acceptance_fixture_catalog(catalog)`.

- [ ] **Step 1: Write the failing validator and registration**

Create `test/internal/park_woodburn_mainline_acceptance_fixtures.jl` with:

```julia
using Test
using Oscar
using Suslin

const PARK_WOODBURN_MAINLINE_ACCEPTANCE_CATALOG_PATH =
    joinpath(@__DIR__, "..", "fixtures", "park_woodburn_mainline_acceptance_cases.jl")

const PW_MAINLINE_SECTION_REFS = Dict(
    :section2 => "refs/arXiv-alg-geom9405003v1 Section 2",
    :section3 => "refs/arXiv-alg-geom9405003v1 Section 3",
    :section4 => "refs/arXiv-alg-geom9405003v1 Section 4",
    :section5 => "refs/arXiv-alg-geom9405003v1 Section 5",
)
```

Then implement validator helpers with these exact exported behaviors:

```julia
function validate_park_woodburn_mainline_acceptance_fixture(entry)
    # Must throw ArgumentError for invalid entries and return true otherwise.
end

function validate_park_woodburn_mainline_acceptance_fixture_catalog(catalog)
    # Must validate schema, required ids/classes, unique ids, source refs,
    # supported/staged boundaries, and negative controls.
end
```

Add the new internal validator to `TEST_GROUP_FILES["internal"]` in `test/runtests.jl` immediately after `internal/park_woodburn_sln_driver_fixtures.jl`.

- [ ] **Step 2: Run the red validator command**

Run:

```bash
julia --project=. -e 'include("test/internal/park_woodburn_mainline_acceptance_fixtures.jl")'
```

Expected before the catalog exists: FAIL with an include/file-not-found error for `park_woodburn_mainline_acceptance_cases.jl`.

- [ ] **Step 3: Implement the fixture catalog**

Create `test/fixtures/park_woodburn_mainline_acceptance_cases.jl` as module `ParkWoodburnMainlineAcceptanceFixtureCatalog`. Include sibling catalogs:

```julia
const PARK_WOODBURN_SL3_DRIVER_CATALOG_PATH =
    joinpath(@__DIR__, "park_woodburn_sl3_driver_cases.jl")
const ECP_MAINLINE_CATALOG_PATH =
    joinpath(@__DIR__, "ecp_mainline_cases.jl")
const PARK_WOODBURN_SLN_DRIVER_CATALOG_PATH =
    joinpath(@__DIR__, "park_woodburn_sln_driver_cases.jl")
```

The catalog must return:

```julia
return (;
    cases = [
        sl3_multivariate_case,
        sln_recursive_case,
        readme_style_case,
        staged_missing_evidence_case,
    ],
    negative_controls = [
        det_not_one_control,
        unsupported_ring_control,
        missing_evidence_control,
    ],
)
```

Required ids:

```julia
"pw-mainline-sl3-multivariate-issue184-qq"
"pw-mainline-sln-recursive-issue185-186-gf2"
"pw-mainline-readme-ordinary-polynomial-qq"
"pw-mainline-staged-missing-evidence-qq"
"pw-mainline-negative-det-not-one"
"pw-mainline-negative-unsupported-coefficient-ring"
"pw-mainline-negative-missing-evidence"
```

- [ ] **Step 4: Fill validator rules**

Implement rules that require:

```julia
const REQUIRED_PW_MAINLINE_CASE_IDS = Set([
    "pw-mainline-sl3-multivariate-issue184-qq",
    "pw-mainline-sln-recursive-issue185-186-gf2",
    "pw-mainline-readme-ordinary-polynomial-qq",
    "pw-mainline-staged-missing-evidence-qq",
])

const REQUIRED_PW_MAINLINE_NEGATIVE_IDS = Set([
    "pw-mainline-negative-det-not-one",
    "pw-mainline-negative-unsupported-coefficient-ring",
    "pw-mainline-negative-missing-evidence",
])
```

Supported `:mainline_accepted` entries must validate determinant one, ordinary exact field-backed polynomial rings, required section refs, required upstream issue ids, and required upstream evidence ids. Staged entries must include `missing_evidence` and `staged_reason`. Negative controls must record `base_case_id` and `reason` and must throw `ArgumentError` when validated.

- [ ] **Step 5: Run focused green tests**

Run:

```bash
julia --project=. -e 'include("test/internal/park_woodburn_mainline_acceptance_fixtures.jl")'
julia --project=. test/runtests.jl internal
```

Expected: both commands exit 0. The focused validator must prove negative controls by using `@test_throws ArgumentError`.

- [ ] **Step 6: Commit Task 1**

Run:

```bash
git add test/fixtures/park_woodburn_mainline_acceptance_cases.jl test/internal/park_woodburn_mainline_acceptance_fixtures.jl test/runtests.jl docs/superpowers/plans/2026-07-03-issue-270-mainline-acceptance-catalog.md
git commit -m "test: add Park-Woodburn mainline acceptance catalog"
```
