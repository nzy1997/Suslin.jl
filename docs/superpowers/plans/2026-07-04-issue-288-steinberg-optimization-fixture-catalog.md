# Issue 288 Steinberg Optimization Fixture Catalog Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a validated test fixture catalog for ordinary-polynomial Steinberg rewrite opportunities consumed by #188.

**Architecture:** Add one fixture module with positive and negative Steinberg optimization cases, one internal validator that proves schema/ring/product/rule correctness, and register the validator in the internal test group. Keep all behavior test-only and avoid production optimizer or public API changes.

**Tech Stack:** Julia, Oscar, Suslin elementary matrices, Julia `Test`, existing fixture-catalog patterns.

## Global Constraints

- Do not implement a Steinberg optimizer.
- Do not export a public API.
- Do not change `elementary_factorization`.
- Do not broaden Laurent or ToricBuilder support.
- Positive cases must use exact field-backed ordinary polynomial rings.
- Every positive case must include `refs/arXiv-alg-geom9405003v1 Section 6`.
- Every positive case must include `consumer_issue_ids = ("#188",)`.
- The validator must reject mismatched factor rings, stale expected products, and invalid commutator indices.

---

## File Structure

- Create `test/fixtures/steinberg_optimization_cases.jl`: private fixture catalog module `SteinbergOptimizationFixtureCatalog`.
- Create `test/internal/steinberg_optimization_fixtures.jl`: validator functions and focused testset.
- Modify `test/runtests.jl`: include the new validator in the internal group.
- Modify this plan file as tasks are completed.

### Task 1: Fixture Catalog And Validator

**Files:**
- Create: `test/fixtures/steinberg_optimization_cases.jl`
- Create: `test/internal/steinberg_optimization_fixtures.jl`
- Modify: `test/runtests.jl`

**Interfaces:**
- Consumes: `elementary_matrix(n, row, col, coefficient, R)`, Oscar `polynomial_ring`, `identity_matrix`, and matrix multiplication.
- Produces: `SteinbergOptimizationFixtureCatalog.catalog()`, `SteinbergOptimizationFixtureCatalog.cases_by_id()`, `validate_steinberg_optimization_fixture(entry)`, and `validate_steinberg_optimization_fixture_catalog(catalog)`.

- [ ] **Step 1: Write the failing validator and test registration**

Create `test/internal/steinberg_optimization_fixtures.jl` with constants and empty validator functions first:

```julia
using Test
using Oscar
using Suslin

const STEINBERG_OPTIMIZATION_CATALOG_PATH =
    joinpath(@__DIR__, "..", "fixtures", "steinberg_optimization_cases.jl")
const STEINBERG_SECTION_6_REF = "refs/arXiv-alg-geom9405003v1 Section 6"

const REQUIRED_STEINBERG_POSITIVE_IDS = Set([
    "steinberg-identity-removal-qq",
    "steinberg-same-position-merge-qq",
    "steinberg-inverse-cancellation-qq",
    "steinberg-commutator-forward-qq",
    "steinberg-commutator-reverse-qq",
    "steinberg-disjoint-commutator-identity-qq",
])

const REQUIRED_STEINBERG_NEGATIVE_IDS = Set([
    "steinberg-negative-mismatched-factor-rings",
    "steinberg-negative-stale-expected-product",
    "steinberg-negative-invalid-commutator-indices",
])

const REQUIRED_STEINBERG_RULE_NAMES = Set([
    :identity_removal,
    :same_position_merge,
    :inverse_cancellation,
    :commutator_forward,
    :commutator_reverse,
    :disjoint_commutator_identity,
])

function validate_steinberg_optimization_fixture(entry)
    throw(ArgumentError("fixture catalog is not implemented yet"))
end

function validate_steinberg_optimization_fixture_catalog(catalog)
    throw(ArgumentError("fixture catalog is not implemented yet"))
end

@testset "Steinberg optimization fixture catalog" begin
    include(STEINBERG_OPTIMIZATION_CATALOG_PATH)
    catalog = SteinbergOptimizationFixtureCatalog.catalog()
    validate_steinberg_optimization_fixture_catalog(catalog)
end
```

Add `"internal/steinberg_optimization_fixtures.jl"` to `TEST_GROUP_FILES["internal"]` in `test/runtests.jl` immediately after `"internal/park_woodburn_mainline_acceptance_fixtures.jl"`.

- [ ] **Step 2: Run the red validator command**

Run:

```bash
julia --project=. -e 'include("test/internal/steinberg_optimization_fixtures.jl")'
```

Expected: FAIL because `test/fixtures/steinberg_optimization_cases.jl` does not exist or because the initial validator throws `ArgumentError`.

- [ ] **Step 3: Implement the fixture module**

Create `test/fixtures/steinberg_optimization_cases.jl` as module `SteinbergOptimizationFixtureCatalog`. The module must define helpers:

```julia
function _ring_metadata(description, R, generator_names, generators)
function _ordinary_ring_constructor(coefficient, variables)
function _factor_metadata(row, col, coefficient)
function _product(factors, R, n::Int)
function _case(; id, rule_name, description, ring_constructor, ring, matrix_size,
    factor_metadata, factors, expected_rewrite_factors, original_product,
    rewritten_product, rewrite_span, rule_metadata, source_refs, consumer_issue_ids)
function _negative_control(id, base_case_id, reason, entry)
function catalog()
function cases_by_id()
```

The positive catalog must contain exactly these rule examples:

```julia
identity_removal_case = _case(
    id = "steinberg-identity-removal-qq",
    rule_name = :identity_removal,
    factors = (
        elementary_matrix(3, 1, 2, zero(R), R),
        elementary_matrix(3, 2, 3, x + one(R), R),
    ),
    expected_rewrite_factors = (elementary_matrix(3, 2, 3, x + one(R), R),),
)

same_position_merge_case = _case(
    id = "steinberg-same-position-merge-qq",
    rule_name = :same_position_merge,
    factors = (
        elementary_matrix(3, 1, 2, x, R),
        elementary_matrix(3, 1, 2, y + one(R), R),
    ),
    expected_rewrite_factors = (elementary_matrix(3, 1, 2, x + y + one(R), R),),
)

inverse_cancellation_case = _case(
    id = "steinberg-inverse-cancellation-qq",
    rule_name = :inverse_cancellation,
    factors = (
        elementary_matrix(3, 2, 3, x * y + one(R), R),
        elementary_matrix(3, 2, 3, -(x * y + one(R)), R),
    ),
    expected_rewrite_factors = (),
)

commutator_forward_case = _case(
    id = "steinberg-commutator-forward-qq",
    rule_name = :commutator_forward,
    factors = (
        elementary_matrix(3, 1, 2, x + one(R), R),
        elementary_matrix(3, 2, 3, y, R),
        elementary_matrix(3, 1, 2, -(x + one(R)), R),
        elementary_matrix(3, 2, 3, -y, R),
    ),
    expected_rewrite_factors = (elementary_matrix(3, 1, 3, (x + one(R)) * y, R),),
)

commutator_reverse_case = _case(
    id = "steinberg-commutator-reverse-qq",
    rule_name = :commutator_reverse,
    factors = (
        elementary_matrix(3, 2, 3, x, R),
        elementary_matrix(3, 1, 2, y + one(R), R),
        elementary_matrix(3, 2, 3, -x, R),
        elementary_matrix(3, 1, 2, -(y + one(R)), R),
    ),
    expected_rewrite_factors = (elementary_matrix(3, 1, 3, -x * (y + one(R)), R),),
)

disjoint_commutator_identity_case = _case(
    id = "steinberg-disjoint-commutator-identity-qq",
    rule_name = :disjoint_commutator_identity,
    factors = (
        elementary_matrix(4, 1, 2, x, R),
        elementary_matrix(4, 3, 4, y + one(R), R),
        elementary_matrix(4, 1, 2, -x, R),
        elementary_matrix(4, 3, 4, -(y + one(R)), R),
    ),
    expected_rewrite_factors = (),
)
```

For each case, fill the remaining fields with matching factor metadata, product values computed through `_product`, `rewrite_span`, Section 6 source refs, and `consumer_issue_ids = ("#188",)`.

Negative controls must be:

```julia
"steinberg-negative-mismatched-factor-rings"
"steinberg-negative-stale-expected-product"
"steinberg-negative-invalid-commutator-indices"
```

The mismatched-ring control should replace one factor with a factor over a separate polynomial ring. The stale-product control should keep valid factors but set `original_product = identity_matrix(R, 3)` for a nonidentity product. The invalid-commutator control should set `rule_metadata.indices` so the claimed commutator relation violates the Section 6 index condition.

- [ ] **Step 4: Implement validator rules**

Replace the initial validator with helpers that enforce:

```julia
function _steinberg_field(entry, field::Symbol)
function _steinberg_require_matrix_over(matrix_value, R, n::Int, label)
function _steinberg_factor_product(factors, R, n::Int)
function _steinberg_assert_ring(entry)
function _steinberg_assert_factor_sequence(entry, R, n::Int, field::Symbol, metadata_field::Symbol)
function _steinberg_assert_rule_metadata(entry, R)
function validate_steinberg_optimization_fixture(entry)
function validate_steinberg_optimization_fixture_catalog(catalog)
```

`validate_steinberg_optimization_fixture(entry)` must throw `ArgumentError` unless the entry has all required fields, a known `rule_name`, an exact field-backed `MPolyRing`, square matrix factors over one ring, factor metadata matching `elementary_matrix`, Section 6 source refs, `#188` consumer metadata, exact stored original and rewritten products, and exact equality between original and rewritten products.

`_steinberg_assert_rule_metadata(entry, R)` must check the rule-specific index conditions:

```julia
:identity_removal => removed factor coefficient is zero
:same_position_merge => two factors share row/col and rewritten coefficient is their sum
:inverse_cancellation => two factors share row/col, coefficients sum to zero, rewrite is empty
:commutator_forward => E_ij(a) E_jl(b) E_ij(-a) E_jl(-b), i != l, rewrite E_il(ab)
:commutator_reverse => E_ij(a) E_li(b) E_ij(-a) E_li(-b), j != l, rewrite E_lj(-ab)
:disjoint_commutator_identity => E_ij(a) E_lp(b) E_ij(-a) E_lp(-b), i != p, j != l, rewrite is empty
```

`validate_steinberg_optimization_fixture_catalog(catalog)` must require positive and negative ids, unique ids across both lists, one positive case for every required rule name, and negative controls that each throw `ArgumentError` when passed to `validate_steinberg_optimization_fixture`.

- [ ] **Step 5: Run focused green tests**

Run:

```bash
julia --project=. -e 'include("test/internal/steinberg_optimization_fixtures.jl")'
julia --project=. test/runtests.jl internal
```

Expected: both commands exit 0.

- [ ] **Step 6: Run package verification**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: exits 0.

- [ ] **Step 7: Inspect and commit implementation**

Run:

```bash
git diff -- test/fixtures/steinberg_optimization_cases.jl test/internal/steinberg_optimization_fixtures.jl test/runtests.jl docs/superpowers/plans/2026-07-04-issue-288-steinberg-optimization-fixture-catalog.md
git status --short
```

Expected: only the new fixture catalog, internal validator, test runner registration, and this plan are changed.

Commit:

```bash
git add test/fixtures/steinberg_optimization_cases.jl test/internal/steinberg_optimization_fixtures.jl test/runtests.jl docs/superpowers/plans/2026-07-04-issue-288-steinberg-optimization-fixture-catalog.md
git commit -m "test: add Steinberg optimization fixtures"
```

## Self-Review

Spec coverage: Task 1 creates the requested fixture file, internal validator,
and test registration; it records Section 6 refs and #188 consumer ids; it
covers all requested positive relations and required negative controls.

Completeness scan: no deferred implementation markers are
left in the plan.

Type consistency: fixture functions return NamedTuples and validator functions
return `true` or throw `ArgumentError`, matching existing repository validator
style.
