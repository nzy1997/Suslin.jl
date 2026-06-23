# Issue 109 Park-Woodburn Polynomial Fixture Catalog Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a shared ordinary-polynomial Park-Woodburn driver fixture catalog and validator for issue 109.

**Architecture:** Keep the fixture catalog under `test/fixtures` as test support, mirroring existing Laurent/ECP/Quillen catalogs. The internal validator includes the fixture, checks algebraic and metadata invariants, and registers with the repository's `internal` test group without changing `elementary_factorization`.

**Tech Stack:** Julia, Oscar polynomial rings, Suslin test helpers, Test stdlib.

## Global Constraints

- Do not change `elementary_factorization`.
- Do not implement recursive polynomial column peeling or Quillen patching.
- Do not add Laurent or ToricBuilder fixtures.
- Use ordinary polynomial rings over exact coefficient fields: `QQ[X]`, `GF(2)[x, y]`, and `QQ[X, r, g]`.
- Include at least four positive fixture entries and three negative controls.
- Positive entries must have determinant exactly `one(R)`.
- The multivariate Quillen case must reuse a #99 Quillen fixture id and be blocked by #105, not marked currently supported.
- Register the focused validator in the `internal` group of `test/runtests.jl`.
- Required focused verification command: `julia --project=. -e 'include("test/internal/park_woodburn_polynomial_fixtures.jl")'`.
- Required package verification command: `julia --project=. -e 'using Pkg; Pkg.test()'`.

---

## File Structure

- Create `test/fixtures/park_woodburn_polynomial_cases.jl`: catalog module with ordinary polynomial matrix entries, metadata helpers, and negative controls.
- Create `test/internal/park_woodburn_polynomial_fixtures.jl`: validator functions and testset proving both positive entries and negative controls are enforced.
- Modify `test/runtests.jl`: add the validator to `TEST_GROUP_FILES["internal"]`.

---

### Task 1: Add the Catalog Validator Test First

**Files:**
- Create: `test/internal/park_woodburn_polynomial_fixtures.jl`
- Modify: `test/runtests.jl`

**Interfaces:**
- Consumes: `ParkWoodburnPolynomialFixtureCatalog.catalog()` and `cases_by_id()` from the fixture file that Task 2 will add.
- Produces: `validate_park_woodburn_polynomial_fixture(entry)` and `validate_park_woodburn_polynomial_fixture_catalog(catalog)` test helpers.

- [ ] **Step 1: Write the failing internal test**

Create `test/internal/park_woodburn_polynomial_fixtures.jl` with validators that include the not-yet-existing fixture file, require the issue 109 positive and negative ids, validate matrix/ring/route/status metadata, and assert route/status mutations are rejected.

- [ ] **Step 2: Register the internal test**

Add `"internal/park_woodburn_polynomial_fixtures.jl"` to the `internal` list in `test/runtests.jl` after the Quillen patch fixture validator.

- [ ] **Step 3: Run the focused test and confirm RED**

Run:

```bash
julia --project=. -e 'include("test/internal/park_woodburn_polynomial_fixtures.jl")'
```

Expected: FAIL because `test/fixtures/park_woodburn_polynomial_cases.jl` does not exist yet.

---

### Task 2: Add the Fixture Catalog

**Files:**
- Create: `test/fixtures/park_woodburn_polynomial_cases.jl`

**Interfaces:**
- Consumes: Oscar's `polynomial_ring`, `matrix`, `identity_matrix`, and Suslin's `elementary_matrix`.
- Produces: `ParkWoodburnPolynomialFixtureCatalog.catalog()` returning `(; cases, negative_controls)` and `cases_by_id()`.

- [ ] **Step 1: Implement literal catalog entries**

Add a fixture module with helpers for ring metadata, expected status metadata, and negative controls. Define positive entries for:

- `pw-poly-univariate-sl3-fast-local-qq`
- `pw-poly-univariate-sln-disjoint-blocks-qq`
- `pw-poly-recursive-column-peel-gf2`
- `quillen-patched-substitution-witness-qq`

Define negative controls for:

- `pw-poly-det-not-one-control`
- `pw-poly-det-one-outside-witness-control`
- `pw-poly-wrong-route-control`

- [ ] **Step 2: Run the focused test and confirm GREEN**

Run:

```bash
julia --project=. -e 'include("test/internal/park_woodburn_polynomial_fixtures.jl")'
```

Expected: PASS.

- [ ] **Step 3: Run package tests**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: PASS.

---

## Plan Self-Review

- The plan covers every requested file from issue 109.
- The validator is written before the catalog implementation, satisfying TDD.
- Route/status negative mutations are explicit, not decorative.
- The plan avoids public API changes and unrelated implementation work.
