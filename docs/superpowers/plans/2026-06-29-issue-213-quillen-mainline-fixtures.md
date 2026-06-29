# Issue 213 Quillen Mainline Fixtures Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a Park-Woodburn Section 3 mainline Quillen fixture wrapper and validator that reuses existing #99 fixture ids.

**Architecture:** Keep `test/fixtures/quillen_patch_cases.jl` unchanged as the algebraic source of truth. Add `test/fixtures/quillen_mainline_cases.jl` as a thin metadata wrapper over selected #99 entries, then add `test/internal/quillen_mainline_fixtures.jl` to validate the wrapper. Register the validator in the internal test group.

**Tech Stack:** Julia, Oscar, Suslin test fixtures, `Test`.

## Global Constraints

- Preserve existing #99 fixture ids.
- Reuse or alias existing `test/fixtures/quillen_patch_cases.jl` entries; do not duplicate their matrices as a divergent schema.
- Positive cases must include the existing two-open QQ shape, a nontrivial multiplier case, and one case shaped for later Murthy-adapter consumption.
- Negative controls must corrupt denominator coverage and local evidence replay.
- Every positive mainline case must cite `refs/arXiv-alg-geom9405003v1 Section 3`.
- Every positive mainline case must include `#183` in `consumer_issue_ids`.
- The schema must expose raw local denominator provenance separately from exponent `l`.
- The schema must expose denominator-cover terms and multipliers.
- The schema must expose local factor records or placeholders for later #214 sequence records.
- The schema must expose patched-substitution chain metadata using the Park-Woodburn sign convention.
- The schema must expose whether `A(0) = I` is assumed, base-term factors are supplied, or base-term evidence is staged.
- Do not implement automatic cover solving, Murthy local solving, global route selection, ECP, recursive `SL_n`, Laurent/ToricBuilder support, or factor-count optimization.

---

## File Structure

- Create `test/fixtures/quillen_mainline_cases.jl`: `QuillenMainlineFixtureCatalog` module, helper constructors, positive entries, negative controls, and `cases_by_id()`.
- Create `test/internal/quillen_mainline_fixtures.jl`: focused validator and tests for positive entries and negative controls.
- Modify `test/runtests.jl`: include `internal/quillen_mainline_fixtures.jl` after `internal/quillen_patch_fixtures.jl`.

### Task 1: Add the Mainline Validator Red Test

**Files:**
- Create: `test/internal/quillen_mainline_fixtures.jl`

**Interfaces:**
- Consumes: `test/fixtures/quillen_mainline_cases.jl`, which does not exist yet in this red step.
- Produces: `validate_quillen_mainline_fixture(entry)` and `validate_quillen_mainline_fixture_catalog(catalog)` for Task 2.

- [ ] **Step 1: Write the failing validator test**

Create `test/internal/quillen_mainline_fixtures.jl` with validators that include the future fixture file and check the required metadata. The file should define:

```julia
using Test
using Oscar
using Suslin

const QUILLEN_MAINLINE_CATALOG_PATH = joinpath(@__DIR__, "..", "fixtures", "quillen_mainline_cases.jl")
const REQUIRED_QUILLEN_MAINLINE_IDS = Set([
    "quillen-two-open-cover-qq",
    "quillen-nontrivial-multipliers-qq",
    "quillen-patched-substitution-witness-qq",
    "quillen-constructive-acceptance-gf2",
])
const REQUIRED_QUILLEN_MAINLINE_NEGATIVE_IDS = Set([
    "quillen-mainline-uncovered-denominator-control",
    "quillen-mainline-tampered-local-evidence-control",
])
const PARK_WOODBURN_SECTION_3_REF = "refs/arXiv-alg-geom9405003v1 Section 3"

function _qml_field(entry, field::Symbol)
    hasproperty(entry, field) || throw(ArgumentError("mainline fixture entry missing field $(field)"))
    return getproperty(entry, field)
end

function _qml_product(factors, R, n::Int)
    product = identity_matrix(R, n)
    for factor in factors
        product *= factor
    end
    return product
end
```

Continue the file with validation helpers for ring metadata, denominator cover, raw denominator provenance, local evidence, patched chain, base term evidence, source refs, consumer ids, and catalog uniqueness. End with a `@testset "Quillen mainline fixture catalog"` that includes `QUILLEN_MAINLINE_CATALOG_PATH`, validates the catalog, checks required positive/negative ids, checks at least four positives and two negatives, and asserts both negative controls throw `ArgumentError`.

- [ ] **Step 2: Run focused red test**

Run:

```bash
julia --project=. -e 'include("test/internal/quillen_mainline_fixtures.jl")'
```

Expected: FAIL because `test/fixtures/quillen_mainline_cases.jl` is missing.

- [ ] **Step 3: Commit the red validator**

```bash
git add test/internal/quillen_mainline_fixtures.jl
git commit -m "test: add quillen mainline fixture validator"
```

### Task 2: Add the Mainline Fixture Wrapper

**Files:**
- Create: `test/fixtures/quillen_mainline_cases.jl`
- Modify: `test/internal/quillen_mainline_fixtures.jl` if Task 1 validator needs minor alignment with concrete field names.

**Interfaces:**
- Consumes: `QuillenPatchFixtureCatalog.cases_by_id()` from `test/fixtures/quillen_patch_cases.jl`.
- Produces: `QuillenMainlineFixtureCatalog.catalog()` and `QuillenMainlineFixtureCatalog.cases_by_id()`.

- [ ] **Step 1: Create the wrapper module**

Create `test/fixtures/quillen_mainline_cases.jl` with:

```julia
module QuillenMainlineFixtureCatalog

using Oscar
using Suslin

const QUILLEN_PATCH_CATALOG_PATH = joinpath(@__DIR__, "quillen_patch_cases.jl")
const PARK_WOODBURN_SECTION_3_REF = "refs/arXiv-alg-geom9405003v1 Section 3"

function _patch_catalog_module()
    if !isdefined(Main, :QuillenPatchFixtureCatalog)
        Base.include(Main, QUILLEN_PATCH_CATALOG_PATH)
    end
    return Main.QuillenPatchFixtureCatalog
end
```

Add helper constructors for `raw_denominator_provenance`, `denominator_cover`, `local_evidence`, `patched_substitution_chain`, `base_term_evidence`, and `_mainline_case`. The patched chain helper must compute each step as `next = current - variable * denominator^exponent_l * coverage_multiplier`.

- [ ] **Step 2: Add positive cases**

In `catalog()`, fetch `patch_entries = _patch_catalog_module().cases_by_id()` and build positive mainline cases from:

```julia
two_open = patch_entries["quillen-two-open-cover-qq"]
nontrivial = patch_entries["quillen-nontrivial-multipliers-qq"]
patched = patch_entries["quillen-patched-substitution-witness-qq"]
constructive = patch_entries["quillen-constructive-acceptance-gf2"]
```

Use `_mainline_case` for each. Add `#183` to every `consumer_issue_ids`; add `#214` where local sequence placeholders are present; add `#211` to the Murthy-shaped constructive case. For the Murthy-shaped case, set `murthy_adapter_handoff = (; status = :staged_until_adapter, issue_id = "#211", accepts_denominator_one_factors = true)`.

- [ ] **Step 3: Add negative controls**

Add:

```julia
bad_cover = _negative_control(
    "quillen-mainline-uncovered-denominator-control",
    two_open_mainline.id,
    "mutated coverage multiplier makes denominator cover sum differ from one",
    merge(two_open_mainline, (;
        denominator_cover = _denominator_cover_from_data((
            two_open.denominator_data[1],
            merge(two_open.denominator_data[2], (; coverage_multiplier = two_open.ring.generators[2])),
        ), two_open.ring.object),
    )),
)
```

Add a second negative control that mutates the first local evidence factor by multiplying it by a nontrivial elementary matrix, while leaving the recorded expected product unchanged.

- [ ] **Step 4: Run focused green test**

Run:

```bash
julia --project=. -e 'include("test/internal/quillen_mainline_fixtures.jl")'
```

Expected: PASS.

- [ ] **Step 5: Commit the wrapper**

```bash
git add test/fixtures/quillen_mainline_cases.jl test/internal/quillen_mainline_fixtures.jl
git commit -m "test: add quillen mainline fixture catalog"
```

### Task 3: Register the Internal Validator

**Files:**
- Modify: `test/runtests.jl`

**Interfaces:**
- Consumes: `test/internal/quillen_mainline_fixtures.jl`.
- Produces: internal-group coverage for the new validator.

- [ ] **Step 1: Register the test file**

In the `"internal"` list in `test/runtests.jl`, add:

```julia
"internal/quillen_mainline_fixtures.jl",
```

immediately after:

```julia
"internal/quillen_patch_fixtures.jl",
```

- [ ] **Step 2: Run internal group**

Run:

```bash
julia --project=. test/runtests.jl internal
```

Expected: PASS.

- [ ] **Step 3: Commit registration**

```bash
git add test/runtests.jl
git commit -m "test: register quillen mainline fixtures"
```

## Final Verification

Run:

```bash
julia --project=. -e 'include("test/internal/quillen_mainline_fixtures.jl")'
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: both commands exit 0.

## Self-Review

- Spec coverage: covered by wrapper, validator, negative controls, and test registration tasks.
- Placeholder scan: no `TBD`, `TODO`, or unresolved fields are present.
- Type consistency: all public test helper names use the `quillen_mainline` prefix and catalog module name `QuillenMainlineFixtureCatalog`.

## Execution Choice

Plan complete and saved to `docs/superpowers/plans/2026-06-29-issue-213-quillen-mainline-fixtures.md`. Two execution options:

1. Subagent-Driven (recommended) - dispatch a fresh subagent per task, review between tasks, fast iteration.
2. Inline Execution - execute tasks in this session using executing-plans, batch execution with checkpoints.

Automatic choice for this non-interactive run: Subagent-Driven, because it is marked recommended.
