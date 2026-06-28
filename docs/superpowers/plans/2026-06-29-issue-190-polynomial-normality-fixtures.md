# Issue 190 Polynomial Normality Fixture Catalog Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a small Park-Woodburn Section 2 ordinary-polynomial normality fixture catalog and expert verifier.

**Architecture:** Keep the catalog as test support under `test/fixtures`, mirroring existing fixture modules. The expert verifier includes the catalog, checks exact algebraic conventions for the three Section 2 layers, and rejects catalog negative controls.

**Tech Stack:** Julia, Oscar polynomial rings and matrices, Suslin elementary matrix helpers, Test stdlib.

## Global Constraints

- Do not implement the normality certificate API.
- Do not add Murthy, Quillen, ECP, Laurent, ToricBuilder, Steinberg factor-count, or factorization-driver logic.
- Use exact field-backed ordinary polynomial rings; this plan uses `QQ[x, y]`.
- Add exactly one representative positive fixture for each Park-Woodburn Section 2 layer requested by issue 190.
- Every positive case must check exact equality to its stored target matrix.
- Include negative controls and reject them through the same verifier.
- Register the expert test in `test/runtests.jl`.
- Required focused verification command: `julia --project=. -e 'include("test/expert/polynomial_normality_fixtures.jl")'`.
- Required package verification command: `julia --project=. -e 'using Pkg; Pkg.test()'`.

---

## File Structure

- Create `test/expert/polynomial_normality_fixtures.jl`: failing validator/test first, then final validator for all positive and negative cases.
- Create `test/fixtures/polynomial_normality_cases.jl`: catalog module with `catalog()`, `cases_by_id()`, and `negative_controls_by_id()`.
- Modify `test/runtests.jl`: add `expert/polynomial_normality_fixtures.jl` to the expert group near the existing normality tests.

---

### Task 1: Write the Expert Verifier Test First

**Files:**
- Create: `test/expert/polynomial_normality_fixtures.jl`

**Interfaces:**
- Consumes: `ParkWoodburnPolynomialNormalityFixtureCatalog.catalog()` and lookup helpers that Task 2 will add.
- Produces: `validate_polynomial_normality_fixture(entry)` and `validate_polynomial_normality_fixture_catalog(catalog)` test helpers.

- [ ] **Step 1: Create the verifier test file**

Write `test/expert/polynomial_normality_fixtures.jl` with helpers that include `test/fixtures/polynomial_normality_cases.jl`, define required ids, compute Cohn/rank-one/conjugated targets independently, validate metadata, and assert all negative controls throw `ArgumentError`.

- [ ] **Step 2: Run the focused command to confirm RED**

Run:

```bash
julia --project=. -e 'include("test/expert/polynomial_normality_fixtures.jl")'
```

Expected: FAIL because `test/fixtures/polynomial_normality_cases.jl` does not exist yet.

---

### Task 2: Add the Fixture Catalog

**Files:**
- Create: `test/fixtures/polynomial_normality_cases.jl`

**Interfaces:**
- Consumes: Oscar matrix/ring constructors and `Suslin.elementary_matrix`.
- Produces: `ParkWoodburnPolynomialNormalityFixtureCatalog.catalog()`, `cases_by_id()`, and `negative_controls_by_id()`.

- [ ] **Step 1: Implement the catalog**

Create one positive case each for:

- `pw-section2-cohn-type-qq`
- `pw-section2-orthogonal-rank-one-qq`
- `pw-section2-conjugated-elementary-qq`

Also create negative controls:

- `pw-section2-cohn-type-tampered-target-control`
- `pw-section2-rank-one-bad-orthogonality-control`
- `pw-section2-conjugated-elementary-tampered-target-control`

- [ ] **Step 2: Run the focused command to confirm GREEN**

Run:

```bash
julia --project=. -e 'include("test/expert/polynomial_normality_fixtures.jl")'
```

Expected: PASS with a Park-Woodburn Section 2 fixture-catalog testset.

---

### Task 3: Register and Verify the Expert Test

**Files:**
- Modify: `test/runtests.jl`

**Interfaces:**
- Consumes: `test/expert/polynomial_normality_fixtures.jl`.
- Produces: expert-suite coverage for the new catalog.

- [ ] **Step 1: Register the expert file**

Add `"expert/polynomial_normality_fixtures.jl"` to the `TEST_GROUP_FILES["expert"]` list near `expert/normality.jl`.

- [ ] **Step 2: Run focused and package verification**

Run:

```bash
julia --project=. -e 'include("test/expert/polynomial_normality_fixtures.jl")'
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: both commands exit 0.

## Plan Self-Review

- The plan starts with a failing expert test before adding the catalog.
- It covers every requested output file and the optional registration.
- It keeps the fixture count intentionally small: one positive case per Section
  2 layer.
- It excludes all certificate and unrelated algorithm work.
