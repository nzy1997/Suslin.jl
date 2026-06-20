# Issue 39 ToricBuilder Issue 38 Fixture Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an offline, reusable test fixture for the Issue #38 ToricBuilder `Q` block and validate its current staged failure behavior.

**Architecture:** Keep the fixture under `test/fixtures/` as test support, not public API. The fixture module constructs the exact Oscar Laurent matrix and normalization metadata; the internal test owns validation helpers, failure assertions, and the negative control.

**Tech Stack:** Julia, Oscar Laurent polynomial rings and matrices, Suslin Laurent normalization/factorization APIs, Test stdlib.

## Global Constraints

- Do not make Issue #38 factorize in this issue.
- Do not depend on a local ToricBuilder checkout at test runtime.
- Fixture file is `test/fixtures/toricbuilder_issue38_cases.jl`.
- Internal validator file is `test/internal/toricbuilder_issue38_fixture.jl`.
- Register the validator in the `internal` group in `test/runtests.jl`.
- Fixture ID is `toricbuilder-issue-38-q-block`.
- Original determinant must be exactly `u*v`.
- Both row and column normalized cores must have determinant `1`.
- Both normalized cores must currently fail with `staged SL_n to local SL_3 reduction failure` and `failed to solve local SL_3 obligation`.
- Focused verification command is `julia --project=. -e 'include("test/internal/toricbuilder_issue38_fixture.jl")'`.
- Package verification command is `julia --project=. -e 'using Pkg; Pkg.test()'`.
- Do not commit `Manifest.toml`.

---

## File Structure

- Create `test/fixtures/toricbuilder_issue38_cases.jl`: module `ToricBuilderIssue38Cases`, fixture construction, determinant metadata, row normalization, column normalization, provenance, and expected current status.
- Create `test/internal/toricbuilder_issue38_fixture.jl`: validator functions and focused tests, including a corrupted-entry negative control.
- Modify `test/runtests.jl`: add `internal/toricbuilder_issue38_fixture.jl` to the default internal group.

---

### Task 1: Add the Failing Validator

**Files:**
- Create: `test/internal/toricbuilder_issue38_fixture.jl`

**Interfaces:**
- Consumes: planned module `ToricBuilderIssue38Cases.catalog()` from `test/fixtures/toricbuilder_issue38_cases.jl`.
- Produces: `validate_toricbuilder_issue38_fixture(entry)` and `validate_toricbuilder_issue38_catalog(catalog)` for internal validation.

- [ ] **Step 1: Create the validator with an expected missing-fixture failure**

Create `test/internal/toricbuilder_issue38_fixture.jl` with validation helpers that:

```julia
using Test
using Suslin
using Oscar

const TORICBUILDER_ISSUE38_FIXTURE_PATH =
    joinpath(@__DIR__, "..", "fixtures", "toricbuilder_issue38_cases.jl")

const REQUIRED_TORICBUILDER_ISSUE38_FIELDS = (
    :id,
    :kind,
    :ring,
    :dimensions,
    :inputs,
    :determinant_profile,
    :normalizations,
    :expected_current_status,
    :provenance,
    :consumer_test_ids,
)

function _fixture_matrix_size(A)
    return (nrows(A), ncols(A))
end

function _require_field(entry, field::Symbol)
    hasproperty(entry, field) || throw(ArgumentError("fixture entry missing field $(field)"))
    return getproperty(entry, field)
end
```

Add focused checks for metadata, determinant profile, row normalization, column normalization, and expected current failures. The final testset must begin with:

```julia
@testset "ToricBuilder Issue 38 Q block fixture" begin
    @test isfile(TORICBUILDER_ISSUE38_FIXTURE_PATH)

    include(TORICBUILDER_ISSUE38_FIXTURE_PATH)
    catalog = ToricBuilderIssue38Cases.catalog()
    @test validate_toricbuilder_issue38_catalog(catalog)
end
```

- [ ] **Step 2: Run the focused validator and verify it fails for the missing fixture**

Run:

```bash
julia --project=. -e 'include("test/internal/toricbuilder_issue38_fixture.jl")'
```

Expected: FAIL because `test/fixtures/toricbuilder_issue38_cases.jl` does not exist yet.

---

### Task 2: Add the Fixture Module

**Files:**
- Create: `test/fixtures/toricbuilder_issue38_cases.jl`
- Modify: `test/internal/toricbuilder_issue38_fixture.jl`

**Interfaces:**
- Consumes: `Suslin.normalize_laurent_gl_matrix`, `Suslin.classify_laurent_determinant`, Oscar `laurent_polynomial_ring`, `matrix`, and `diagonal_matrix`.
- Produces: `ToricBuilderIssue38Cases.catalog()` returning `(; cases = [entry])`.

- [ ] **Step 1: Implement fixture construction**

Create module `ToricBuilderIssue38Cases` with:

```julia
module ToricBuilderIssue38Cases

using Oscar
using Suslin

function _issue38_q_block(R, u, v)
    return matrix(R, [
        1 + v^-1              1      1  0  1 + v^-1              1 + v^-1;
        u*v^-1 + 1 + v^-1     u      1  1  u*v^-1 + 1 + v^-1     u*v^-1 + v^-1;
        u*v^-1                u*v    0  0  u*v^-1                u*v^-1;
        1                     1      1  0  1                     1;
        1 + v^-1              v      0  0  v^-1                  1 + v^-1;
        u + v                 u*v    v  v  u + v                 u + v
    ])
end
```

Then add helpers for:

```julia
_ring_metadata(R, u, v) = (;
    description = "GF(2)[u^+/-1, v^+/-1]",
    object = R,
    generators = (u, v),
    variables = ("u", "v"),
)

_issue38_provenance() = (;
    source = :toricbuilder_issue_38_mwe,
    issue = "#38",
    issue_url = "https://github.com/nzy1997/Suslin.jl/issues/38",
    source_description = "upper-left Q block of transfer_result.column_transformation for the 1 + x + x*y color-code example",
    reported_main_commit = "c985b1aac9fc9152d860e4e90d012964867bb27d",
)
```

Compute:

```julia
determinant = det(Q)
row_normalization = Suslin.normalize_laurent_gl_matrix(Q)
Dcol = diagonal_matrix(R, [inv(determinant), one(R), one(R), one(R), one(R), one(R)])
col_core = Q * Dcol
```

Return one fixture entry with `inputs.matrix = Q`, `determinant_profile.expected_determinant = determinant`, row and column normalization metadata, expected failure strings, provenance, and consumer test ID `issue-39-toricbuilder-issue38-fixture`.

- [ ] **Step 2: Complete validator logic**

Update `test/internal/toricbuilder_issue38_fixture.jl` so `validate_toricbuilder_issue38_fixture(entry)` checks:

```julia
_fixture_matrix_size(entry.inputs.matrix) == (6, 6)
det(entry.inputs.matrix) == entry.determinant_profile.expected_determinant
entry.determinant_profile.expected_determinant == u * v
classify_laurent_determinant(entry.inputs.matrix).classification ==
    entry.determinant_profile.expected_class
det(entry.normalizations.row.core) == one(R)
verify_laurent_gl_normalization(entry.inputs.matrix, entry.normalizations.row.normalization)
det(entry.normalizations.column.core) == one(R)
entry.inputs.matrix * entry.normalizations.column.correction_factor ==
    entry.normalizations.column.core
```

Add `_assert_expected_factorization_failure(matrix, status)` that catches
`elementary_factorization(matrix)` and verifies the two required message
substrings.

- [ ] **Step 3: Run the focused validator and verify it passes**

Run:

```bash
julia --project=. -e 'include("test/internal/toricbuilder_issue38_fixture.jl")'
```

Expected: PASS.

---

### Task 3: Register the Internal Test and Verify Negative Control

**Files:**
- Modify: `test/runtests.jl`
- Modify: `test/internal/toricbuilder_issue38_fixture.jl`

**Interfaces:**
- Consumes: `validate_toricbuilder_issue38_fixture(entry)`.
- Produces: internal suite coverage for the Issue #38 fixture.

- [ ] **Step 1: Add the internal test to the runner**

Add this file to `TEST_GROUP_FILES["internal"]` in `test/runtests.jl`:

```julia
"internal/toricbuilder_issue38_fixture.jl",
```

- [ ] **Step 2: Add the corrupted-entry negative control**

In the focused testset, add:

```julia
entry = only(catalog.cases)
corrupted_matrix = copy(entry.inputs.matrix)
corrupted_matrix[1, 1] += one(base_ring(corrupted_matrix))
bad_entry = merge(entry, (; inputs = merge(entry.inputs, (; matrix = corrupted_matrix))))
@test_throws ArgumentError validate_toricbuilder_issue38_fixture(bad_entry)
```

The expected failure is an `ArgumentError` from determinant or normalization
verification because the matrix no longer matches the recorded Issue #38
determinant/normalization metadata.

- [ ] **Step 3: Run focused and internal verification**

Run:

```bash
julia --project=. -e 'include("test/internal/toricbuilder_issue38_fixture.jl")'
julia --project=. test/runtests.jl internal
```

Expected: both PASS.

---

### Task 4: Package Verification and Commit

**Files:**
- Verify all changed files from Tasks 1-3.

**Interfaces:**
- Consumes: registered internal test suite.
- Produces: committed implementation ready for PR.

- [ ] **Step 1: Run required package verification**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: PASS.

- [ ] **Step 2: Check repository diff**

Run:

```bash
git status --short
git diff -- test/fixtures/toricbuilder_issue38_cases.jl test/internal/toricbuilder_issue38_fixture.jl test/runtests.jl
```

Expected: only the fixture, internal test, runner registration, and Superpowers docs/plans are changed.

- [ ] **Step 3: Commit implementation**

Run:

```bash
git add test/fixtures/toricbuilder_issue38_cases.jl test/internal/toricbuilder_issue38_fixture.jl test/runtests.jl docs/superpowers/plans/2026-06-20-issue-39-toricbuilder-issue38-fixture.md
git commit -m "Add Issue 38 ToricBuilder fixture"
```

Expected: commit succeeds without adding `Manifest.toml`.

## Self-Review

- The plan covers fixture construction, metadata validation, row and column normalization, expected current staged failures, negative control, test registration, and package verification.
- No unresolved placeholders remain.
- Property names are consistent between fixture and validator: `inputs.matrix`, `determinant_profile.expected_determinant`, `normalizations.row`, `normalizations.column`, and `expected_current_status`.
