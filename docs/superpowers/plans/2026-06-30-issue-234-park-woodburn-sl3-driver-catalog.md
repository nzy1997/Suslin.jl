# Issue 234 Park-Woodburn SL3 Driver Catalog Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a validated source-grounded `SL_3` driver fixture catalog for parent #184.

**Architecture:** Add a new test-only fixture module for `3 x 3` ordinary-polynomial driver inputs and a focused internal validator that checks ring exactness, determinant-one matrices, selected-variable metadata, source references, unique ids, explicit #184 consumers, supported evidence replay, staged boundary reasons, and negative-control rejection. Keep the new schema separate from the older broad Park-Woodburn polynomial route catalog.

**Tech Stack:** Julia, Oscar exact polynomial rings and matrices, existing Quillen patch and Quillen mainline fixture catalogs, `Test`.

## Global Constraints

- Do not implement local solving, Quillen patch assembly, public `elementary_factorization` dispatch, ECP, recursive `SL_n`, Laurent/ToricBuilder support, coordinate-change search, or Steinberg optimization.
- Keep all changes test/catalog scoped; do not export new public API names from `src/Suslin.jl`.
- Every positive or staged-positive case must be a square `3 x 3` matrix over an exact field-backed ordinary polynomial ring.
- Every positive or staged-positive case must have determinant one, nonempty `source_refs`, and `consumer_issue_ids` containing `#184`.
- Supported status must require replayable explicit local-form monicity evidence or replayable upstream Quillen mainline evidence.
- Determinant-one matrices with no local-form, variable-change, normality/conjugation, or upstream evidence must be staged, not supported.
- Negative controls must reject determinant not one, unsupported coefficient ring, selected variable not a generator, claimed-but-missing local evidence, and supported-without-witness metadata.
- Focused verification command is `julia --project=. -e 'include("test/internal/park_woodburn_sl3_driver_fixtures.jl")'`.
- Package verification command is `julia --project=. -e 'using Pkg; Pkg.test()'`.

---

## File Structure

- Create `test/internal/park_woodburn_sl3_driver_fixtures.jl`: validator, required id checks, negative-control assertions, and focused testset.
- Create `test/fixtures/park_woodburn_sl3_driver_cases.jl`: `ParkWoodburnSL3DriverFixtureCatalog` with positive/staged cases and negative controls.
- Modify `test/runtests.jl`: register `internal/park_woodburn_sl3_driver_fixtures.jl` in the internal test group after the existing Park-Woodburn polynomial fixture validator.

### Task 1: Add Red Driver Validator Contract

**Files:**
- Create: `test/internal/park_woodburn_sl3_driver_fixtures.jl`
- Modify: `test/runtests.jl`

**Interfaces:**
- Consumes future `ParkWoodburnSL3DriverFixtureCatalog.catalog()` and `.cases_by_id()`.
- Produces `validate_park_woodburn_sl3_driver_fixture(entry)` and `validate_park_woodburn_sl3_driver_fixture_catalog(catalog)` for later #184 children.

- [ ] **Step 1: Write the failing internal validator**

Create `test/internal/park_woodburn_sl3_driver_fixtures.jl` with these concrete requirements:

```julia
using Test
using Oscar
using Suslin

const PARK_WOODBURN_SL3_DRIVER_CATALOG_PATH =
    joinpath(@__DIR__, "..", "fixtures", "park_woodburn_sl3_driver_cases.jl")
const QUILLEN_MAINLINE_CATALOG_PATH =
    joinpath(@__DIR__, "..", "fixtures", "quillen_mainline_cases.jl")
const QUILLEN_PATCH_CATALOG_PATH =
    joinpath(@__DIR__, "..", "fixtures", "quillen_patch_cases.jl")
const PARK_WOODBURN_SECTION_3_REF = "refs/arXiv-alg-geom9405003v1 Section 3"
const PARK_WOODBURN_SECTION_5_REF = "refs/arXiv-alg-geom9405003v1 Section 5"
```

The validator must define helper checks for required fields, exact
field-backed ordinary polynomial rings, `3 x 3` matrices, determinant one,
selected generator metadata, monic-in-selected-variable replay, upstream
Quillen mainline replay, supported/staged evidence boundaries, unique ids, and
negative control rejection.

- [ ] **Step 2: Register the internal test**

Add `"internal/park_woodburn_sl3_driver_fixtures.jl"` to the `internal` list in
`test/runtests.jl` immediately after
`"internal/park_woodburn_polynomial_fixtures.jl"`.

- [ ] **Step 3: Run RED verification**

Run:

```bash
julia --project=. -e 'include("test/internal/park_woodburn_sl3_driver_fixtures.jl")'
```

Expected: FAIL with a missing catalog file/module error for
`park_woodburn_sl3_driver_cases.jl`.

- [ ] **Step 4: Commit the red validator**

```bash
git add test/internal/park_woodburn_sl3_driver_fixtures.jl test/runtests.jl
git commit -m "test: add sl3 driver fixture validator contract"
```

### Task 2: Add Driver Fixture Catalog And Negative Controls

**Files:**
- Create: `test/fixtures/park_woodburn_sl3_driver_cases.jl`
- Test: `test/internal/park_woodburn_sl3_driver_fixtures.jl`

**Interfaces:**
- Consumes the validator from Task 1 plus existing `QuillenMainlineFixtureCatalog` and `QuillenPatchFixtureCatalog`.
- Produces `ParkWoodburnSL3DriverFixtureCatalog.catalog()` and `.cases_by_id()`.

- [ ] **Step 1: Implement the fixture module**

Create `test/fixtures/park_woodburn_sl3_driver_cases.jl` with module
`ParkWoodburnSL3DriverFixtureCatalog`. The module must provide:

```julia
function catalog()
    return (; cases = [...], negative_controls = [...])
end

function cases_by_id()
    return Dict(entry.id => entry for entry in catalog().cases)
end
```

The cases must include the five ids listed in the design:
`sl3-driver-univariate-fast-local-qq`,
`sl3-driver-multivariate-monic-special-form-qq`,
`sl3-driver-quillen-mainline-evidence-gf2`,
`sl3-driver-legacy-quillen-patched-substitution-qq`, and
`sl3-driver-det-one-no-witness-staged-qq`.

The negative controls must include the five ids listed in the validator:
`sl3-driver-negative-det-not-one`,
`sl3-driver-negative-unsupported-coefficient-ring`,
`sl3-driver-negative-selected-variable-not-generator`,
`sl3-driver-negative-claimed-local-evidence-missing`, and
`sl3-driver-negative-supported-without-witness`.

- [ ] **Step 2: Run GREEN focused verification**

Run:

```bash
julia --project=. -e 'include("test/internal/park_woodburn_sl3_driver_fixtures.jl")'
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
git add test/fixtures/park_woodburn_sl3_driver_cases.jl test/internal/park_woodburn_sl3_driver_fixtures.jl test/runtests.jl
git commit -m "test: add park woodburn sl3 driver catalog"
```

## Plan Self-Review

- The plan maps every issue requirement to a validator or catalog task.
- The task boundary preserves TDD: the validator fails before the catalog
  exists, then passes after the catalog is added.
- No public API or production route changes are included.
- The required Agent Desk package verification command is included.
