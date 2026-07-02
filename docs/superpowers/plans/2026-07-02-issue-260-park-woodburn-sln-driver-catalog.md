# Issue 260 Park-Woodburn SLn Driver Catalog Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a validated ordinary-polynomial Park-Woodburn `SL_n` recursive driver fixture catalog for #186.

**Architecture:** Follow the repository's existing fixture pattern: a catalog module in `test/fixtures/`, a validator/test harness in `test/internal/`, and test-runner registration in `test/runtests.jl`. The catalog stores source-grounded matrices and replay metadata but does not implement recursive `SL_n` factorization.

**Tech Stack:** Julia, Oscar, Suslin internal certificate helpers, Test stdlib.

## Global Constraints

- Repository has no `AGENTS.md`; follow README test commands.
- Add `test/fixtures/park_woodburn_sln_driver_cases.jl`.
- Add `test/internal/park_woodburn_sln_driver_fixtures.jl`.
- Register the internal validator in `test/runtests.jl`.
- Keep the new catalog separate from `park_woodburn_polynomial_cases.jl`.
- Positive and staged-positive catalog matrices must be exact ordinary-polynomial determinant-one matrices of size `n >= 4` over field-backed polynomial rings.
- Record `consumer_issue_ids = ("#186",)` for every positive and staged-positive case.
- Record Park-Woodburn source refs for Reduction to `SL_3(k[x_1,...,x_m])`, Section 4 ECP, and Section 5 final `SL_3` evidence where applicable.
- Support roles must include `:issue186_mainline`, `:staged_issue186_candidate`, and `:legacy_regression`.
- Staged entries must use machine-readable reason codes from `:missing_ecp_evidence`, `:missing_final_sl3_route`, `:unsupported_coefficient_ring`, and `:legacy_regression_only`.
- Legacy fast-local/disjoint-block examples must not count as #186 mainline support.
- Do not implement ECP, `SL_3` realization, recursive peeling, public dispatch, Laurent/ToricBuilder support, #187 final acceptance, or Steinberg factor-count optimization.

---

## File Structure

- Create `test/internal/park_woodburn_sln_driver_fixtures.jl`: validator, required-id assertions, positive/staged-positive checks, negative-control checks.
- Create `test/fixtures/park_woodburn_sln_driver_cases.jl`: catalog module, matrix construction helpers, five positive/staged-positive cases, five negative controls, `cases_by_id()`.
- Modify `test/runtests.jl`: include the new internal validator in the `internal` group after the SL3 driver fixture validator.

### Task 1: Validator Red Test

**Files:**
- Create: `test/internal/park_woodburn_sln_driver_fixtures.jl`
- Modify: `test/runtests.jl`

**Interfaces:**
- Consumes: a future `ParkWoodburnSLnDriverFixtureCatalog.catalog()` and `.cases_by_id()`.
- Produces: `validate_park_woodburn_sln_driver_fixture(entry)` and `validate_park_woodburn_sln_driver_fixture_catalog(catalog)`.

- [ ] **Step 1: Write the failing validator test**

Create `test/internal/park_woodburn_sln_driver_fixtures.jl` with this structure:

```julia
using Test
using Oscar
using Suslin

const PARK_WOODBURN_SLN_DRIVER_CATALOG_PATH =
    joinpath(@__DIR__, "..", "fixtures", "park_woodburn_sln_driver_cases.jl")
const PARK_WOODBURN_SL3_DRIVER_CATALOG_PATH =
    joinpath(@__DIR__, "..", "fixtures", "park_woodburn_sl3_driver_cases.jl")
const PARK_WOODBURN_SLN_REDUCTION_REF =
    "refs/arXiv-alg-geom9405003v1 Section \"Reduction to SL_3(k[x_1,...,x_m])\""
const PARK_WOODBURN_SECTION_4_REF = "refs/arXiv-alg-geom9405003v1 Section 4"
const PARK_WOODBURN_SECTION_5_REF = "refs/arXiv-alg-geom9405003v1 Section 5"

const REQUIRED_SLN_DRIVER_CASE_IDS = Set([
    "sln-driver-legacy-recursive-column-peel-qq",
    "sln-driver-sl4-gf2-ecp-mainline",
    "sln-driver-sl5-gf2-two-step",
    "sln-driver-sl4-final-sl3-evidence-qq",
    "sln-driver-staged-missing-final-sl3-qq",
])

const REQUIRED_SLN_DRIVER_NEGATIVE_IDS = Set([
    "sln-driver-negative-det-not-one",
    "sln-driver-negative-unsupported-coefficient-ring",
    "sln-driver-negative-corrupt-peel-expectation",
    "sln-driver-negative-unknown-staged-reason",
    "sln-driver-negative-false-mainline-support",
])
```

The file must then define helper functions for field access, ring/matrix checks,
factor products, SL3 fixture lookup, per-step replay, final-route checks, and
catalog validation. The validator must throw `ArgumentError` for invalid
entries and return `true` for valid entries.

- [ ] **Step 2: Register the validator**

In `test/runtests.jl`, add:

```julia
        "internal/park_woodburn_sln_driver_fixtures.jl",
```

after:

```julia
        "internal/park_woodburn_sl3_driver_fixtures.jl",
```

- [ ] **Step 3: Run RED verification**

Run:

```bash
julia --project=. -e 'include("test/internal/park_woodburn_sln_driver_fixtures.jl")'
```

Expected: failure because `test/fixtures/park_woodburn_sln_driver_cases.jl` does not exist yet. The failure proves the new validator is wired to the intended catalog path.

- [ ] **Step 4: Commit the red test**

```bash
git add test/internal/park_woodburn_sln_driver_fixtures.jl test/runtests.jl
git commit -m "test: add issue 260 SLn driver validator"
```

### Task 2: Fixture Catalog Green Implementation

**Files:**
- Create: `test/fixtures/park_woodburn_sln_driver_cases.jl`
- Modify: `test/internal/park_woodburn_sln_driver_fixtures.jl` only if the red test exposes a validator defect.

**Interfaces:**
- Consumes: `validate_park_woodburn_sln_driver_fixture_catalog`.
- Produces: `ParkWoodburnSLnDriverFixtureCatalog.catalog()` and `.cases_by_id()`.

- [ ] **Step 1: Add the catalog module skeleton**

Create `test/fixtures/park_woodburn_sln_driver_cases.jl`:

```julia
module ParkWoodburnSLnDriverFixtureCatalog

using Oscar
using Suslin

const PARK_WOODBURN_SL3_DRIVER_CATALOG_PATH =
    joinpath(@__DIR__, "park_woodburn_sl3_driver_cases.jl")
const PARK_WOODBURN_POLYNOMIAL_CATALOG_PATH =
    joinpath(@__DIR__, "park_woodburn_polynomial_cases.jl")
const PARK_WOODBURN_SLN_REDUCTION_REF =
    "refs/arXiv-alg-geom9405003v1 Section \"Reduction to SL_3(k[x_1,...,x_m])\""
const PARK_WOODBURN_SECTION_4_REF = "refs/arXiv-alg-geom9405003v1 Section 4"
const PARK_WOODBURN_SECTION_5_REF = "refs/arXiv-alg-geom9405003v1 Section 5"

if !isdefined(@__MODULE__, :ParkWoodburnSL3DriverFixtureCatalog)
    include(PARK_WOODBURN_SL3_DRIVER_CATALOG_PATH)
end
if !isdefined(@__MODULE__, :ParkWoodburnPolynomialFixtureCatalog)
    include(PARK_WOODBURN_POLYNOMIAL_CATALOG_PATH)
end
```

Add helper functions `_ring_metadata`, `_ordinary_ring_constructor`, `_case`,
`_negative_control`, `_factor_product`, `_inverse_factor_product`,
`_matrix_column`, `_build_peel_step`, `_sl3_cases_by_id`, and
`_polynomial_cases_by_id`.

- [ ] **Step 2: Build replayed peel cases**

Use `_build_peel_step(next_block, last_column; bottom_row_entries, ecp_status,
ecp_source_case_id)` to construct determinant-one matrices:

```julia
function _build_peel_step(next_block, last_column; bottom_row_entries, ecp_status = :replayed, ecp_source_case_id = nothing)
    R = base_ring(next_block)
    d = nrows(next_block) + 1
    certificate = ecp_status == :replayed ?
        Suslin.ecp_column_reduction_certificate(collect(last_column), R) :
        nothing
    left_factors = certificate === nothing ? typeof(identity_matrix(R, d))[] : certificate.factors
    left_product = _factor_product(left_factors, R, d)
    after_left = block_embedding(next_block, d, collect(1:(d - 1)))
    for col in 1:(d - 1)
        after_left[d, col] = bottom_row_entries[col]
    end
    input_matrix = _inverse_factor_product(left_factors, R, d) * after_left
    right_factors = typeof(identity_matrix(R, d))[]
    for col in 1:(d - 1)
        coeff = -after_left[d, col]
        coeff == zero(R) || push!(right_factors, elementary_matrix(d, d, col, coeff, R))
    end
    right_product = _factor_product(right_factors, R, d)
    peeled_matrix = after_left * right_product
    return (;
        dimension = d,
        input_matrix = input_matrix,
        last_column = collect(last_column),
        last_column_ecp = (;
            status = ecp_status,
            certificate = certificate,
            source_case_id = ecp_source_case_id,
            target_column = Suslin._column_peel_target_column(R, d),
        ),
        right_clearing = (;
            status = :replayed,
            after_left_matrix = after_left,
            right_factors = right_factors,
            peeled_matrix = peeled_matrix,
            next_block = next_block,
        ),
        next_block = next_block,
    )
end
```

- [ ] **Step 3: Add five positive or staged-positive cases**

Create:

1. `sln-driver-legacy-recursive-column-peel-qq` from
   `ParkWoodburnPolynomialFixtureCatalog.cases_by_id()["pw-poly-recursive-column-peel-sl3-qq"]`,
   with `support_role = :legacy_regression`, `expected_status = :staged`, and
   `staged_reason_codes = (:legacy_regression_only,)`.
2. `sln-driver-sl4-gf2-ecp-mainline`, a multivariate `SL_4` case over the
   supported #184 GF(2) final block with a replayed embedded-three-block ECP
   last-column certificate.
3. `sln-driver-sl5-gf2-two-step`, an `SL_5` case over the same ring whose first
   peel descends to the previous `SL_4` case and whose second peel descends to
   the supported #184 `SL_3` final block.
4. `sln-driver-sl4-final-sl3-evidence-qq`, a QQ multivariate `SL_4` case whose
   final block references `sl3-driver-multivariate-monic-special-form-qq`.
5. `sln-driver-staged-missing-final-sl3-qq`, a determinant-one staged candidate
   with replayed ECP metadata but `final_route.status = :missing` and
   `staged_reason_codes = (:missing_final_sl3_route,)`.

- [ ] **Step 4: Add five negative controls**

Create controls by mutating valid entries:

1. `sln-driver-negative-det-not-one`: change a diagonal entry so determinant is
   not one.
2. `sln-driver-negative-unsupported-coefficient-ring`: use `ZZ["X"]`.
3. `sln-driver-negative-corrupt-peel-expectation`: change `expected_peel_count`
   or `descent_dimensions` to disagree with `peel_steps`.
4. `sln-driver-negative-unknown-staged-reason`: use
   `staged_reason_codes = (:unknown_reason_code,)`.
5. `sln-driver-negative-false-mainline-support`: mark a staged/missing-final
   entry as `:issue186_mainline` and `:supported` without final `SL_3`
   evidence.

- [ ] **Step 5: Run GREEN verification**

Run:

```bash
julia --project=. -e 'include("test/internal/park_woodburn_sln_driver_fixtures.jl")'
```

Expected: exit 0.

- [ ] **Step 6: Commit the catalog**

```bash
git add test/fixtures/park_woodburn_sln_driver_cases.jl test/internal/park_woodburn_sln_driver_fixtures.jl test/runtests.jl
git commit -m "test: add issue 260 SLn driver catalog"
```

### Task 3: Full Verification and Review

**Files:**
- No planned file edits.

**Interfaces:**
- Consumes: all changes from Tasks 1 and 2.
- Produces: verified branch ready for PR creation.

- [ ] **Step 1: Run focused verification**

Run:

```bash
julia --project=. -e 'include("test/internal/park_woodburn_sln_driver_fixtures.jl")'
```

Expected: exit 0.

- [ ] **Step 2: Run default package verification**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: exit 0.

- [ ] **Step 3: Review the diff**

Run:

```bash
git diff --stat origin/main...HEAD
git diff --check
```

Expected: no whitespace errors and the diff is scoped to the spec/plan,
fixture catalog, validator, and `test/runtests.jl`.

- [ ] **Step 4: Commit any final fixes**

If verification or review exposes a defect, fix it narrowly, rerun the affected
verification command, and commit with:

```bash
git add <changed-files>
git commit -m "test: fix issue 260 SLn driver validation"
```
