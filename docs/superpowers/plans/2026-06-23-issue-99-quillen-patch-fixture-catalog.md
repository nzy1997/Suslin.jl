# Issue 99 Quillen Patch Fixture Catalog Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a shared deterministic Quillen local-to-global patching fixture catalog and validator for ordinary polynomial rings.

**Architecture:** Keep the catalog under `test/fixtures/` as test support, not public API. The fixture module builds exact `QQ[X,r,g]` and `GF(2)[X,r,s]` cases with denominator coverage, local certificates, local factor matrices, expected global corrections, and optional patched-substitution witness data; the internal validator reconstructs and checks those claims exactly.

**Tech Stack:** Julia, Oscar ordinary polynomial rings, Suslin Quillen patching helpers, Test stdlib.

## Global Constraints

- Do not implement constructive patch assembly.
- Do not change the public `elementary_factorization` driver.
- Do not solve arbitrary local realizability.
- Do not add Laurent `GL_n` determinant correction or ToricBuilder-specific fixtures.
- Do not add new exported package APIs.
- Catalog file is `test/fixtures/quillen_patch_cases.jl`.
- Validator file is `test/internal/quillen_patch_fixtures.jl`.
- Register the validator in the `internal` group in `test/runtests.jl`.
- Catalog must validate at least five named fixture entries.
- Required valid fixture ids are `quillen-two-open-cover-qq`, `quillen-nontrivial-multipliers-qq`, `quillen-supplied-local-certificate-gf2`, `quillen-patched-substitution-witness-qq`, and `quillen-constructive-acceptance-gf2`.
- Required negative-control fixture ids are `quillen-uncovered-denominator-control` and `quillen-tampered-local-factor-control`.
- Each fixture must include a fixture id, target matrix or correction, substitution variable, local denominators, coverage multipliers, local factor data, expected status, and patched-substitution witness fields where present.
- Positive entries must check exact equations such as `sum(c_i * d_i) == 1` and local factor products against expected local corrections.
- Negative controls must prove the validator rejects a corrupt coverage multiplier and a corrupt local factor.
- Keep Laurent examples from `test/expert/quillen_patching_exact.jl` out of this catalog.
- Focused validator command is `julia --project=. -e 'include("test/internal/quillen_patch_fixtures.jl")'`.
- Internal group command is `julia --project=. test/runtests.jl internal`.
- Full package command is `julia --project=. -e 'using Pkg; Pkg.test()'`.

---

## File Structure

- Create `test/fixtures/quillen_patch_cases.jl`: owns catalog construction, helper constructors, positive entries, negative controls, and `cases_by_id()`.
- Create `test/internal/quillen_patch_fixtures.jl`: owns validation helpers and focused tests. This is test-only support and may define validator functions without exporting package APIs.
- Modify `test/runtests.jl`: add `internal/quillen_patch_fixtures.jl` to the `internal` test group after `internal/ecp_column_fixtures.jl`.

---

### Task 1: Quillen Fixture Validator Skeleton

**Files:**
- Create: `test/internal/quillen_patch_fixtures.jl`

**Interfaces:**
- Consumes: fixture entries shaped by Task 2.
- Produces: `validate_quillen_patch_fixture(entry)` and `validate_quillen_patch_fixture_catalog(catalog)`.

- [ ] **Step 1: Write the failing validator**

Create `test/internal/quillen_patch_fixtures.jl` with this header and required-id constants:

```julia
using Test
using Suslin
using Oscar

const QUILLEN_PATCH_CATALOG_PATH = joinpath(@__DIR__, "..", "fixtures", "quillen_patch_cases.jl")
const REQUIRED_QUILLEN_PATCH_FIELDS = (
    :id,
    :kind,
    :stage_coverage,
    :ring_constructor,
    :ring,
    :size,
    :substitution_variable,
    :target_matrix,
    :base_matrix,
    :denominator_data,
    :local_factors,
    :expected,
    :patched_substitution_witness,
    :source_refs,
    :consumer_issue_ids,
)
const REQUIRED_QUILLEN_PATCH_IDS = Set([
    "quillen-two-open-cover-qq",
    "quillen-nontrivial-multipliers-qq",
    "quillen-supplied-local-certificate-gf2",
    "quillen-patched-substitution-witness-qq",
    "quillen-constructive-acceptance-gf2",
])
const REQUIRED_QUILLEN_PATCH_NEGATIVE_IDS = Set([
    "quillen-uncovered-denominator-control",
    "quillen-tampered-local-factor-control",
])
```

Implement helpers with these names:

```julia
_quillen_field(entry, field::Symbol)
_quillen_factor_product(factors, R, n::Int)
_quillen_assert_metadata(entry)
_quillen_assert_status(entry)
_quillen_assert_denominator_cover(entry)
_quillen_assert_local_factor(entry, local_factor)
_quillen_assert_local_factors(entry)
_quillen_assert_patched_substitution_witness(entry)
validate_quillen_patch_fixture(entry)
validate_quillen_patch_fixture_catalog(catalog)
```

The validator must throw `ArgumentError` on invalid fixture metadata or invalid algebraic claims.

- [ ] **Step 2: Add a RED testset**

At the bottom of `test/internal/quillen_patch_fixtures.jl`, add a testset named `"Quillen patch fixture catalog"` that includes `QUILLEN_PATCH_CATALOG_PATH`, calls `QuillenPatchFixtureCatalog.catalog()`, calls `validate_quillen_patch_fixture_catalog(catalog)`, checks required ids, checks at least five valid entries, and checks each negative control throws.

- [ ] **Step 3: Run RED verification**

Run:

```bash
julia --project=. -e 'include("test/internal/quillen_patch_fixtures.jl")'
```

Expected: FAIL because `test/fixtures/quillen_patch_cases.jl` does not exist yet.

### Task 2: Quillen Fixture Catalog

**Files:**
- Create: `test/fixtures/quillen_patch_cases.jl`

**Interfaces:**
- Produces: `QuillenPatchFixtureCatalog.catalog()` returning `(; cases, negative_controls)`.
- Produces: `QuillenPatchFixtureCatalog.cases_by_id()` returning `Dict(entry.id => entry for entry in catalog().cases)`.

- [ ] **Step 1: Create the catalog module**

Create `test/fixtures/quillen_patch_cases.jl` with:

```julia
module QuillenPatchFixtureCatalog

using Oscar
using Suslin

function _ring_metadata(description, R, generator_names, generators)
    return (; description, object = R, generator_names, generators)
end

function _ordinary_ring_constructor(coefficient, variables)
    return (; function_name = :polynomial_ring, coefficient, variables)
end

function _case(; id, kind, stage_coverage, ring_constructor, ring, size, substitution_variable, target_matrix, base_matrix, denominator_data, local_factors, expected, patched_substitution_witness, source_refs, consumer_issue_ids)
    return (; id, kind, stage_coverage, ring_constructor, ring, size, substitution_variable, target_matrix, base_matrix, denominator_data, local_factors, expected, patched_substitution_witness, source_refs, consumer_issue_ids)
end

function _negative_control(id, base_case_id, reason, entry)
    return merge(entry, (; id, base_case_id, reason))
end
```

- [ ] **Step 2: Add local factor helpers**

Add helper constructors:

```julia
function _local_factor(; certificate, denominator, coverage_multiplier, correction, factor, expected_correction)
    return (; certificate, denominator, coverage_multiplier, correction, factor, expected_correction)
end

function _correction(row, col, entry)
    return (; row, col, entry)
end

function _certificate(indices, denominators)
    return (; indices, denominators)
end
```

- [ ] **Step 3: Add five valid ordinary polynomial cases**

Inside `catalog()`, construct:

```julia
RQ, (X, r, g) = Oscar.polynomial_ring(QQ, ["X", "r", "g"])
R2, (Y, u, s) = Oscar.polynomial_ring(GF(2), ["X", "r", "s"])
```

Use only ordinary polynomial rings. For elementary correction products, use
single-row/column repeated elementary factors so the product is exactly
`elementary_matrix(n, row, col, sum(weighted_entries), R)`.

The two-open `QQ` seed case uses denominators `r` and `1-r` with multipliers
`1, 1`. The nontrivial multiplier case also uses denominators `r` and `1-r`,
but with nonconstant Bezout multipliers `1 + g * (1-r)` and `1 - g * r`.
The `GF(2)` supplied-certificate case uses denominators `r` and `1+r` with
multipliers `1, 1`. The patched-substitution case records `matrix`, `variable`,
`denominator`, `exponent`, `shift`, and `expected_matrix`. The second `GF(2)`
constructive acceptance case uses a different elementary row/column and target
entry from the supplied-certificate case.

- [ ] **Step 4: Add negative controls**

Create `negative_controls` by merging valid entries:

```julia
bad_cover = _negative_control(
    "quillen-uncovered-denominator-control",
    "quillen-two-open-cover-qq",
    "mutated coverage multiplier",
    merge(two_open_case, (;
        denominator_data = (
            two_open_case.denominator_data[1],
            merge(two_open_case.denominator_data[2], (; coverage_multiplier = r)),
        ),
    )),
)

bad_factor = _negative_control(
    "quillen-tampered-local-factor-control",
    "quillen-supplied-local-certificate-gf2",
    "mutated local factor matrix",
    merge(supplied_certificate_case, (;
        local_factors = (
            merge(supplied_certificate_case.local_factors[1], (;
                factor = supplied_certificate_case.local_factors[1].factor *
                    elementary_matrix(supplied_certificate_case.size, 1, 3, one(R2), R2),
            )),
            supplied_certificate_case.local_factors[2],
        ),
    )),
)
```

### Task 3: Complete Validator Checks And Registration

**Files:**
- Modify: `test/internal/quillen_patch_fixtures.jl`
- Modify: `test/runtests.jl`

**Interfaces:**
- Consumes: `QuillenPatchFixtureCatalog.catalog()`.
- Produces: registered internal test coverage.

- [ ] **Step 1: Validate metadata and status**

In `_quillen_assert_metadata(entry)`, require all fields, ordinary
`polynomial_ring` constructor metadata, ring generator metadata, exact ordinary
polynomial rings, nonempty source refs, and consumer issue ids beginning with
`"#"`.

In `_quillen_assert_status(entry)`, support `:passes` and `:staged_fail`.
`:passes` entries must have determinant-one expected corrections. `:staged_fail`
entries must include `expected.missing` as a nonempty tuple.

- [ ] **Step 2: Validate denominator cover and factors**

In `_quillen_assert_denominator_cover(entry)`, compute:

```julia
sum(data.coverage_multiplier * data.denominator for data in entry.denominator_data; init = zero(R)) == one(R)
```

In `_quillen_assert_local_factor(entry, local_factor)`, check certificate
lengths, certificate denominators for correction row and column, parent rings,
and:

```julia
weighted_entry = local_factor.coverage_multiplier * local_factor.denominator * local_factor.correction.entry
local_factor.expected_correction == elementary_matrix(entry.size, local_factor.correction.row, local_factor.correction.col, weighted_entry, R)
local_factor.factor == local_factor.expected_correction
```

In `_quillen_assert_local_factors(entry)`, check all local factors and require
their product equals `entry.expected.global_correction` and `entry.target_matrix`.

- [ ] **Step 3: Validate patched-substitution witnesses**

If `entry.patched_substitution_witness !== nothing`, require fields
`matrix`, `variable`, `denominator`, `exponent`, `shift`, and
`expected_matrix`, then check:

```julia
Suslin.patched_substitution(witness.matrix, witness.variable, witness.denominator, witness.exponent, witness.shift) == witness.expected_matrix
```

- [ ] **Step 4: Add runtime mutation checks**

In the testset, mutate one valid entry's coverage multiplier and one valid
entry's first local factor, and assert both throw `ArgumentError` through
`validate_quillen_patch_fixture`.

- [ ] **Step 5: Register internal test**

Modify `test/runtests.jl` by adding:

```julia
"internal/quillen_patch_fixtures.jl",
```

after `"internal/ecp_column_fixtures.jl",`.

- [ ] **Step 6: Run focused and internal tests**

Run:

```bash
julia --project=. -e 'include("test/internal/quillen_patch_fixtures.jl")'
julia --project=. test/runtests.jl internal
```

Expected: both commands exit 0.

### Task 4: Full Verification And Publish Prep

**Files:**
- Verify all changed files.

**Interfaces:**
- Produces: branch ready for pull request.

- [ ] **Step 1: Run full package tests**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: exits 0.

- [ ] **Step 2: Inspect branch diff**

Run:

```bash
git status --short
git diff --stat origin/main
git diff --name-only origin/main
```

Expected: only the design doc, plan, catalog, internal validator, and
`test/runtests.jl` changed. Ignored local `Manifest.toml` from dependency
instantiation must not be committed.

- [ ] **Step 3: Publish**

Because Agent Desk's sandbox denies local writes to the shared external
`.git` directory, publish with the GitHub connector by creating a tree and
commit on branch `agent/issue-99-add-a-quillen-local-to-global-fixture-catalog-run-1`, then open a draft PR against `main` linking Issue 99.
