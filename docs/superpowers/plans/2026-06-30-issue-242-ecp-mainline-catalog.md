# Issue 242 ECP Mainline Catalog Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a validated Park-Woodburn Section 4 ECP mainline fixture catalog for parent #185.

**Architecture:** Add a test-only fixture module with ordinary-polynomial column cases and a focused internal validator that checks exact rings, selected-variable and monicity metadata, unimodularity witnesses, source refs, #185 consumers, staged support boundaries, link/lower-variable expectations, and required negative controls. Keep the schema separate from `ecp_column_cases.jl` while reusing existing fixture families as provenance where useful.

**Tech Stack:** Julia, Oscar exact polynomial rings, existing `Suslin` ECP link witness/link step helpers, Julia `Test`.

## Global Constraints

- Do not implement monicity search, link extraction, `SL_3` realization, public dispatch, recursive matrix factorization, Laurent/ToricBuilder support, or factor-count optimization.
- Keep all changes test/catalog scoped; do not export new public API names from `src/Suslin.jl`.
- Every positive or staged-positive case must be an ordinary polynomial column over an exact field-backed polynomial ring.
- Every positive or staged-positive case must have length at least three and a replayable unimodularity witness summing exactly to one.
- Every positive or staged-positive case must have `source_refs` containing `refs/arXiv-alg-geom9405003v1 Section 4`.
- Every positive or staged-positive case must have `consumer_issue_ids` containing `#185`.
- Supported status must require replayable link-step evidence and lower-variable evidence; fixture-only or missing evidence remains staged.
- Negative controls must reject a non-unimodular column, a corrupted Bezout/resultant/link witness, a selected variable that is not a ring generator, and a supported-without-required-evidence entry.
- Focused verification command is `julia --project=. -e 'include("test/internal/ecp_mainline_fixtures.jl")'`.
- Package verification command is `julia --project=. -e 'using Pkg; Pkg.test()'`.

---

## File Structure

- Create `test/internal/ecp_mainline_fixtures.jl`: validator, required id checks, negative-control assertions, and focused testset.
- Create `test/fixtures/ecp_mainline_cases.jl`: `ECPMainlineFixtureCatalog` with five staged-positive cases and four negative controls.
- Modify `test/runtests.jl`: register `internal/ecp_mainline_fixtures.jl` in the internal test group after `internal/ecp_column_fixtures.jl`.

### Task 1: Add Red ECP Mainline Validator Contract

**Files:**
- Create: `test/internal/ecp_mainline_fixtures.jl`
- Modify: `test/runtests.jl`

**Interfaces:**
- Consumes future `ECPMainlineFixtureCatalog.catalog()` and `.cases_by_id()`.
- Produces `validate_ecp_mainline_fixture(entry)` and `validate_ecp_mainline_fixture_catalog(catalog)` for later #185 children.

- [ ] **Step 1: Write the failing internal validator**

Create `test/internal/ecp_mainline_fixtures.jl` with helpers for:

```julia
const ECP_MAINLINE_CATALOG_PATH = joinpath(@__DIR__, "..", "fixtures", "ecp_mainline_cases.jl")
const PARK_WOODBURN_SECTION_4_REF = "refs/arXiv-alg-geom9405003v1 Section 4"
const REQUIRED_ECP_MAINLINE_CASE_IDS = Set([
    "ecp-mainline-gf2-hard-slice",
    "ecp-mainline-qq-link-bezout",
    "ecp-mainline-length4-coupled-qq",
    "ecp-mainline-monicity-change-gf2",
    "ecp-mainline-sl3-route-qq",
])
const REQUIRED_ECP_MAINLINE_NEGATIVE_IDS = Set([
    "ecp-mainline-negative-non-unimodular",
    "ecp-mainline-negative-corrupt-link-witness",
    "ecp-mainline-negative-selected-variable-not-generator",
    "ecp-mainline-negative-supported-without-evidence",
])
```

The validator must define checks for required fields, exact field-backed ordinary polynomial rings, length at least three, selected generator metadata, monicity replay when claimed, unimodularity and witness replay, optional all-entry support evidence, source refs, #185 consumers, stage evidence fields, link witness replay, lower-variable expectation replay, supported/staged evidence boundaries, unique ids, and negative control rejection.

- [ ] **Step 2: Register the internal test**

Add `"internal/ecp_mainline_fixtures.jl"` to the `internal` list in `test/runtests.jl` immediately after `"internal/ecp_column_fixtures.jl"`.

- [ ] **Step 3: Run RED verification**

Run:

```bash
julia --project=. -e 'include("test/internal/ecp_mainline_fixtures.jl")'
```

Expected: FAIL with a missing catalog file/module error for `ecp_mainline_cases.jl`.

- [ ] **Step 4: Commit the red validator**

```bash
git add test/internal/ecp_mainline_fixtures.jl test/runtests.jl docs/superpowers/plans/2026-06-30-issue-242-ecp-mainline-catalog.md
git commit -m "test: add ecp mainline fixture validator contract"
```

### Task 2: Add ECP Mainline Fixture Catalog And Negative Controls

**Files:**
- Create: `test/fixtures/ecp_mainline_cases.jl`
- Test: `test/internal/ecp_mainline_fixtures.jl`

**Interfaces:**
- Consumes the validator from Task 1 and existing ECP column fixture ids as provenance.
- Produces `ECPMainlineFixtureCatalog.catalog()` and `.cases_by_id()`.

- [ ] **Step 1: Implement the fixture module**

Create `test/fixtures/ecp_mainline_cases.jl` with module `ECPMainlineFixtureCatalog`:

```julia
module ECPMainlineFixtureCatalog

using Oscar
using Suslin

function catalog()
    return (;
        cases = [
            gf2_hard_slice_case,
            qq_link_bezout_case,
            length4_coupled_case,
            monicity_change_case,
            sl3_route_case,
        ],
        negative_controls = [
            non_unimodular_control,
            corrupt_link_witness_control,
            bad_selected_variable_control,
            supported_without_evidence_control,
        ],
    )
end

function cases_by_id()
    return Dict(entry.id => entry for entry in catalog().cases)
end

end
```

The positive or staged-positive cases must include:

- `ecp-mainline-gf2-hard-slice`: GF(2) variable-change hard-slice case with supplied link witness metadata and lower-variable expectation.
- `ecp-mainline-qq-link-bezout`: QQ link-Bezout case reused from the old ECP catalog as provenance, with exact Bezout/resultant expectation and staged first-entry/link-step boundary.
- `ecp-mainline-length4-coupled-qq`: length-four Lagrange-style column whose witness uses all four entries and whose three-entry subcolumns are not unimodular.
- `ecp-mainline-monicity-change-gf2`: GF(2) monicity-changing selected-variable case with replayed substitution metadata.
- `ecp-mainline-sl3-route-qq`: QQ monic-first case with supplied link witness metadata and an explicit staged #184 `SL_3` link-realization expectation.

The negative controls must include the four ids listed in the validator.

- [ ] **Step 2: Run GREEN focused verification**

Run:

```bash
julia --project=. -e 'include("test/internal/ecp_mainline_fixtures.jl")'
```

Expected: PASS.

- [ ] **Step 3: Run package verification**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: PASS.

- [ ] **Step 4: Commit the catalog**

```bash
git add test/fixtures/ecp_mainline_cases.jl test/internal/ecp_mainline_fixtures.jl test/runtests.jl
git commit -m "test: add park woodburn ecp mainline catalog"
```

## Plan Self-Review

- The plan maps every issue requirement to a validator or catalog task.
- The task boundary preserves TDD: the validator fails before the catalog exists, then passes after the catalog is added.
- No public API or production route changes are included.
- The required Agent Desk package verification command is included.
- No placeholder markers remain.
