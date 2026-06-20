# Issue 40 ToricBuilder Laurent Problem Catalog Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a shared ToricBuilder Laurent problem catalog that references the existing Issue #38, ToricBuilder contract, and synthetic large Laurent fixtures.

**Architecture:** Keep the manifest under `test/fixtures/` as test support, not public API. The catalog composes existing fixture modules and normalizes their metadata into one named-tuple schema; the internal test validates the schema, live matrix facts, required IDs, and negative controls.

**Tech Stack:** Julia, Oscar Laurent polynomial rings and matrices, Suslin determinant helpers, Test stdlib.

## Global Constraints

- Do not implement new factorization algorithms in this issue.
- Do not require ToricBuilder at runtime.
- Fixture file is `test/fixtures/toricbuilder_laurent_problem_catalog.jl`.
- Internal validator file is `test/internal/toricbuilder_problem_catalog.jl`.
- Register the validator in the `internal` group in `test/runtests.jl`.
- Catalog entries must include stable ID, provenance, determinant class, expected current status, verifier path, and consuming milestone/issue metadata.
- Catalog must include `toricbuilder-issue-38-q-block`, `toricbuilder-factor-toric-block-3-qinv`, `toricbuilder-factor-toric-block-3-pinv`, `laurent-block-local-40x40`, and `laurent-block-local-48x48`.
- Expected status symbols are limited to `:unsupported_now`, `:verified_contract`, `:supported_block_local`, and `:target_acceptance`.
- Focused verification command is `julia --project=. -e 'include("test/internal/toricbuilder_problem_catalog.jl")'`.
- Package verification command is `julia --project=. -e 'using Pkg; Pkg.test()'`.
- Do not commit `Manifest.toml`.

---

## File Structure

- Create `test/fixtures/toricbuilder_laurent_problem_catalog.jl`: module `ToricBuilderLaurentProblemCatalog`, includes existing fixture modules, builds normalized catalog entries, and exposes `catalog()`.
- Create `test/internal/toricbuilder_problem_catalog.jl`: schema and live-data validator, required-ID checks, duplicate-ID negative control, and missing-provenance negative control.
- Modify `test/runtests.jl`: add `internal/toricbuilder_problem_catalog.jl` to the default internal group.

---

### Task 1: Add the Failing Catalog Validator

**Files:**
- Create: `test/internal/toricbuilder_problem_catalog.jl`

**Interfaces:**
- Consumes: planned module `ToricBuilderLaurentProblemCatalog.catalog()` from `test/fixtures/toricbuilder_laurent_problem_catalog.jl`.
- Produces: `validate_toricbuilder_problem_entry(entry)` and `validate_toricbuilder_problem_catalog(catalog)`.

- [x] **Step 1: Write the failing validator test**

Create `test/internal/toricbuilder_problem_catalog.jl` with this initial test:

```julia
using Test
using Suslin
using Oscar

const TORICBUILDER_PROBLEM_CATALOG_PATH =
    joinpath(@__DIR__, "..", "fixtures", "toricbuilder_laurent_problem_catalog.jl")

const REQUIRED_TORICBUILDER_PROBLEM_FIELDS = (
    :id,
    :kind,
    :source_fixture,
    :ring,
    :dimensions,
    :matrix,
    :determinant_profile,
    :expected_current_status,
    :verifier,
    :provenance,
    :consumers,
)

const REQUIRED_TORICBUILDER_PROBLEM_IDS = Set([
    "toricbuilder-issue-38-q-block",
    "toricbuilder-factor-toric-block-3-qinv",
    "toricbuilder-factor-toric-block-3-pinv",
    "laurent-block-local-40x40",
    "laurent-block-local-48x48",
])

const ALLOWED_TORICBUILDER_PROBLEM_STATUSES = Set([
    :unsupported_now,
    :verified_contract,
    :supported_block_local,
    :target_acceptance,
])

function _problem_field(entry, field::Symbol)
    hasproperty(entry, field) || throw(ArgumentError("ToricBuilder problem entry missing field $(field)"))
    return getproperty(entry, field)
end

function _problem_matrix_size(A)
    return (nrows(A), ncols(A))
end

function _problem_determinant_classification(A)
    profile = Suslin.classify_laurent_determinant(A)
    return profile.classification
end

function _assert_nonempty_metadata(value, label::AbstractString)
    value === nothing && throw(ArgumentError("$(label) must not be nothing"))
    if value isa AbstractString
        isempty(value) && throw(ArgumentError("$(label) must not be empty"))
    elseif value isa Tuple || value isa AbstractVector
        isempty(value) && throw(ArgumentError("$(label) must not be empty"))
    elseif value isa NamedTuple
        isempty(keys(value)) && throw(ArgumentError("$(label) must not be empty"))
    end
    return true
end

function validate_toricbuilder_problem_entry(entry)
    for field in REQUIRED_TORICBUILDER_PROBLEM_FIELDS
        _problem_field(entry, field)
    end

    _assert_nonempty_metadata(entry.id, "problem id")
    _assert_nonempty_metadata(entry.provenance, "problem $(entry.id) provenance")
    _assert_nonempty_metadata(entry.consumers, "problem $(entry.id) consumers")
    entry.expected_current_status in ALLOWED_TORICBUILDER_PROBLEM_STATUSES ||
        throw(ArgumentError("problem $(entry.id) has unsupported expected status $(entry.expected_current_status)"))

    hasproperty(entry.verifier, :path) ||
        throw(ArgumentError("problem $(entry.id) missing verifier path"))
    isfile(joinpath(@__DIR__, "..", entry.verifier.path)) ||
        throw(ArgumentError("problem $(entry.id) verifier path $(entry.verifier.path) does not exist"))

    hasproperty(entry.dimensions, :matrix) ||
        throw(ArgumentError("problem $(entry.id) missing matrix dimensions"))
    _problem_matrix_size(entry.matrix) == entry.dimensions.matrix ||
        throw(ArgumentError("problem $(entry.id) matrix dimensions do not match metadata"))

    hasproperty(entry.determinant_profile, :expected_class) ||
        throw(ArgumentError("problem $(entry.id) missing determinant class"))
    _problem_determinant_classification(entry.matrix) == entry.determinant_profile.expected_class ||
        throw(ArgumentError("problem $(entry.id) determinant class does not match metadata"))

    return true
end

function validate_toricbuilder_problem_catalog(catalog)
    hasproperty(catalog, :cases) ||
        throw(ArgumentError("ToricBuilder problem catalog missing cases"))
    isempty(catalog.cases) &&
        throw(ArgumentError("ToricBuilder problem catalog must not be empty"))

    ids = [entry.id for entry in catalog.cases]
    length(ids) == length(unique(ids)) ||
        throw(ArgumentError("ToricBuilder problem catalog ids must be unique"))
    REQUIRED_TORICBUILDER_PROBLEM_IDS ⊆ Set(ids) ||
        throw(ArgumentError("ToricBuilder problem catalog missing required IDs"))

    for entry in catalog.cases
        validate_toricbuilder_problem_entry(entry)
    end

    return true
end

@testset "ToricBuilder Laurent problem catalog" begin
    @test isfile(TORICBUILDER_PROBLEM_CATALOG_PATH)

    include(TORICBUILDER_PROBLEM_CATALOG_PATH)
    catalog = ToricBuilderLaurentProblemCatalog.catalog()
    @test validate_toricbuilder_problem_catalog(catalog)
end
```

- [ ] **Step 2: Run the focused validator and verify it fails for the missing catalog**

Run:

```bash
julia --project=. -e 'include("test/internal/toricbuilder_problem_catalog.jl")'
```

Expected: FAIL because `test/fixtures/toricbuilder_laurent_problem_catalog.jl` does not exist yet.

- [ ] **Step 3: Commit the failing test**

Run:

```bash
git add test/internal/toricbuilder_problem_catalog.jl docs/superpowers/plans/2026-06-21-issue-40-toricbuilder-laurent-problem-catalog.md
git commit -m "test: add ToricBuilder problem catalog validator"
```

Expected: commit succeeds with a known failing focused validator.

---

### Task 2: Add the Shared Catalog Module

**Files:**
- Create: `test/fixtures/toricbuilder_laurent_problem_catalog.jl`
- Modify: `test/internal/toricbuilder_problem_catalog.jl`

**Interfaces:**
- Consumes: `ToricBuilderIssue38Cases.catalog()`, `ToricBuilderFactorToricBlock3Fixture.fixture()`, and `LaurentLargeAcceptanceCases.acceptance_catalog()`.
- Produces: `ToricBuilderLaurentProblemCatalog.catalog()` returning `(; cases = [...])`.

- [ ] **Step 1: Write the catalog module**

Create `test/fixtures/toricbuilder_laurent_problem_catalog.jl` with:

```julia
module ToricBuilderLaurentProblemCatalog

using Oscar
using Suslin

include("toricbuilder_issue38_cases.jl")
include("toricbuilder_factor_toric_block_3.jl")
include("laurent_large_acceptance_cases.jl")

_matrix_dimensions(A) = (nrows(A), ncols(A))

function _ring_metadata_from_matrix(A; description = nothing, generators = nothing)
    R = base_ring(A)
    return (;
        description = description === nothing ? string(R) : description,
        object = R,
        generators,
    )
end

function _determinant_profile(A)
    profile = Suslin.classify_laurent_determinant(A)
    return (;
        expected_class = profile.classification,
        expected_determinant = profile.determinant,
        monomial_exponents = hasproperty(profile, :monomial_exponents) ? profile.monomial_exponents : nothing,
        monomial_coefficient = hasproperty(profile, :monomial_coefficient) ? profile.monomial_coefficient : nothing,
    )
end

function _verifier(path::String, scenario::Symbol)
    return (;
        path,
        scenario,
    )
end

function _consumers(; milestone::Int, issues, tests)
    return (;
        milestone,
        issues = Tuple(issues),
        tests = Tuple(tests),
    )
end

function _issue38_problem()
    issue38 = only(ToricBuilderIssue38Cases.catalog().cases)
    return (;
        id = issue38.id,
        kind = :issue38_q_block,
        source_fixture = :toricbuilder_issue38_cases,
        ring = issue38.ring,
        dimensions = issue38.dimensions,
        matrix = issue38.inputs.matrix,
        determinant_profile = (;
            expected_class = issue38.determinant_profile.expected_class,
            expected_determinant = issue38.determinant_profile.expected_determinant,
            monomial_exponents = issue38.determinant_profile.monomial_exponents,
            monomial_coefficient = issue38.determinant_profile.monomial_coefficient,
        ),
        expected_current_status = :unsupported_now,
        verifier = _verifier("internal/toricbuilder_issue38_fixture.jl", :issue38_fixture_validator),
        provenance = merge(issue38.provenance, (;
            fixture_id = issue38.id,
            source_fixture = "test/fixtures/toricbuilder_issue38_cases.jl",
        )),
        consumers = _consumers(
            milestone = 4,
            issues = ("#38", "#39", "#40"),
            tests = ("test/internal/toricbuilder_issue38_fixture.jl", "test/internal/toricbuilder_problem_catalog.jl"),
        ),
    )
end

function _contract_problem(entry)
    id_suffix = lowercase(entry.toricbuilder_role)
    return (;
        id = "toricbuilder-factor-toric-block-3-$(id_suffix)",
        kind = :toricbuilder_contract,
        source_fixture = :toricbuilder_factor_toric_block_3,
        ring = _ring_metadata_from_matrix(entry.matrix; description = entry.ring),
        dimensions = (; matrix = entry.size),
        matrix = entry.matrix,
        determinant_profile = (;
            expected_class = Symbol(entry.determinant_classification),
            expected_determinant = det(entry.matrix),
            monomial_exponents = nothing,
            monomial_coefficient = nothing,
        ),
        expected_current_status = :verified_contract,
        verifier = _verifier("internal/toricbuilder_contract.jl", :toricbuilder_contract_validator),
        provenance = merge(entry.provenance, (;
            fixture_id = entry.name,
            toricbuilder_role = entry.toricbuilder_role,
            source_fixture = "test/fixtures/toricbuilder_factor_toric_block_3.jl",
        )),
        consumers = _consumers(
            milestone = 4,
            issues = ("#19", "#40"),
            tests = ("test/internal/toricbuilder_contract.jl", "test/internal/toricbuilder_problem_catalog.jl"),
        ),
    )
end

function _large_acceptance_problem(entry, status::Symbol)
    return (;
        id = entry.id,
        kind = :synthetic_block_local_acceptance,
        source_fixture = :laurent_large_acceptance_cases,
        ring = entry.ring,
        dimensions = (; matrix = entry.size),
        matrix = entry.matrix,
        determinant_profile = _determinant_profile(entry.matrix),
        expected_current_status = status,
        verifier = _verifier("public/laurent_large_acceptance.jl", :large_laurent_acceptance),
        provenance = merge(entry.provenance, (;
            fixture_id = entry.id,
            source_fixture = "test/fixtures/laurent_large_acceptance_cases.jl",
        )),
        consumers = _consumers(
            milestone = 4,
            issues = ("#17", "#40"),
            tests = ("test/public/laurent_large_acceptance.jl", "test/internal/toricbuilder_problem_catalog.jl"),
        ),
    )
end

function catalog()
    contract_fixture = ToricBuilderFactorToricBlock3Fixture.fixture()
    contract_cases = Dict(entry.name => entry for entry in contract_fixture.cases)
    large_cases = Dict(entry.id => entry for entry in LaurentLargeAcceptanceCases.acceptance_catalog().cases)

    return (;
        cases = [
            _issue38_problem(),
            _contract_problem(contract_cases["factor_toric_block_3_qinv"]),
            _contract_problem(contract_cases["factor_toric_block_3_pinv"]),
            _large_acceptance_problem(large_cases["laurent-block-local-40x40"], :supported_block_local),
            _large_acceptance_problem(large_cases["laurent-block-local-48x48"], :target_acceptance),
        ],
    )
end

end
```

- [ ] **Step 2: Add live metadata checks and negative controls to the validator**

Update `test/internal/toricbuilder_problem_catalog.jl` so `validate_toricbuilder_problem_entry(entry)` also checks:

```julia
hasproperty(entry.consumers, :milestone) ||
    throw(ArgumentError("problem $(entry.id) missing consumer milestone"))
hasproperty(entry.consumers, :issues) ||
    throw(ArgumentError("problem $(entry.id) missing consumer issues"))
hasproperty(entry.consumers, :tests) ||
    throw(ArgumentError("problem $(entry.id) missing consumer tests"))
isempty(entry.consumers.issues) &&
    throw(ArgumentError("problem $(entry.id) must record at least one consumer issue"))
isempty(entry.consumers.tests) &&
    throw(ArgumentError("problem $(entry.id) must record at least one consumer test"))
```

Add assertions in the testset:

```julia
by_id = Dict(entry.id => entry for entry in catalog.cases)
@test haskey(by_id, "toricbuilder-issue-38-q-block")
@test haskey(by_id, "toricbuilder-factor-toric-block-3-qinv")
@test haskey(by_id, "toricbuilder-factor-toric-block-3-pinv")
@test haskey(by_id, "laurent-block-local-40x40")
@test haskey(by_id, "laurent-block-local-48x48")
@test by_id["toricbuilder-issue-38-q-block"].expected_current_status == :unsupported_now
@test by_id["toricbuilder-factor-toric-block-3-qinv"].expected_current_status == :verified_contract
@test by_id["toricbuilder-factor-toric-block-3-pinv"].expected_current_status == :verified_contract
@test by_id["laurent-block-local-40x40"].expected_current_status == :supported_block_local
@test by_id["laurent-block-local-48x48"].expected_current_status == :target_acceptance

duplicate_cases = [catalog.cases; (merge(first(catalog.cases), (; id = catalog.cases[2].id,)))]
@test_throws ArgumentError validate_toricbuilder_problem_catalog((; cases = duplicate_cases))

missing_provenance = merge(first(catalog.cases), (; provenance = (;),))
@test_throws ArgumentError validate_toricbuilder_problem_entry(missing_provenance)
```

- [ ] **Step 3: Run the focused validator and verify it passes**

Run:

```bash
julia --project=. -e 'include("test/internal/toricbuilder_problem_catalog.jl")'
```

Expected: PASS.

- [ ] **Step 4: Commit the catalog implementation**

Run:

```bash
git add test/fixtures/toricbuilder_laurent_problem_catalog.jl test/internal/toricbuilder_problem_catalog.jl
git commit -m "test: add ToricBuilder Laurent problem catalog"
```

Expected: commit succeeds.

---

### Task 3: Register and Verify the Catalog Test

**Files:**
- Modify: `test/runtests.jl`

**Interfaces:**
- Consumes: `test/internal/toricbuilder_problem_catalog.jl`.
- Produces: default internal test-suite coverage for the shared catalog.

- [ ] **Step 1: Register the internal test**

Add this file to `TEST_GROUP_FILES["internal"]` in `test/runtests.jl` after the Issue #38 fixture:

```julia
"internal/toricbuilder_problem_catalog.jl",
```

- [ ] **Step 2: Run focused, internal, and package verification**

Run:

```bash
julia --project=. -e 'include("test/internal/toricbuilder_problem_catalog.jl")'
julia --project=. test/runtests.jl internal
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: all commands pass.

- [ ] **Step 3: Check repository diff**

Run:

```bash
git status --short
git diff -- test/fixtures/toricbuilder_laurent_problem_catalog.jl test/internal/toricbuilder_problem_catalog.jl test/runtests.jl
```

Expected: only the shared catalog fixture, internal test, runner registration, and Superpowers plan/spec artifacts are changed. `Manifest.toml` is not staged.

- [ ] **Step 4: Commit registration**

Run:

```bash
git add test/runtests.jl
git commit -m "test: register ToricBuilder problem catalog"
```

Expected: commit succeeds.

## Self-Review

- Spec coverage: the plan covers the shared manifest file, internal validator, registration, required IDs, required metadata, negative controls, focused verification, and package verification.
- Placeholder scan: no unresolved placeholders remain.
- Type consistency: `ToricBuilderLaurentProblemCatalog.catalog()` returns `cases`, and validator names match the produced fields.
