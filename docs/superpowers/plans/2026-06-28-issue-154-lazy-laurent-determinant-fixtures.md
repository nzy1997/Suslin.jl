# Issue 154 Lazy Laurent Determinant Fixtures Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an internal lazy Laurent determinant fixture catalog with source-grounded Issue #38 coverage, compact synthetic supported cases, and a non-unit negative case.

**Architecture:** Keep all deliverables in test support. The fixture module constructs named-tuple catalog entries; the internal validator owns exact algebraic checks, drift checks against the existing Issue #38 fixture, negative controls, and test registration.

**Tech Stack:** Julia, Oscar Laurent polynomial rings and matrices, Suslin Laurent determinant classification and normalization helpers, Test stdlib.

## Global Constraints

- Do not implement lazy determinant peeling.
- Do not change public APIs.
- Fixture file is `test/fixtures/laurent_lazy_determinant_cases.jl`.
- Internal validator file is `test/internal/laurent_lazy_determinant_fixtures.jl`.
- Register the validator in the `internal` group in `test/runtests.jl`.
- Reuse or wrap `test/fixtures/toricbuilder_issue38_cases.jl`; do not copy the Issue #38 matrix into the new fixture.
- The Issue #38 fixture entry must report determinant `u*v`.
- The Issue #38 row and column normalized cores must verify as determinant-one cores.
- The catalog must include a determinant-one case that needs no correction.
- The catalog must include a Laurent monomial-unit case whose row and column normalized cores both verify.
- The catalog must include a Laurent non-unit determinant negative input rejected by the supported-correction predicate.
- Focused verification command is `julia --project=. -e 'include("test/internal/laurent_lazy_determinant_fixtures.jl")'`.
- Package verification command is `julia --project=. -e 'using Pkg; Pkg.test()'`.
- Do not commit `Manifest.toml`.

---

## File Structure

- Create `test/fixtures/laurent_lazy_determinant_cases.jl`: module `LaurentLazyDeterminantCases`; wraps `ToricBuilderIssue38Cases.catalog()`; creates compact determinant-one, monomial-unit, and non-unit cases; returns `(; cases = [...])`.
- Create `test/internal/laurent_lazy_determinant_fixtures.jl`: validator functions, focused testset, drift checks, and negative controls.
- Modify `test/runtests.jl`: add `"internal/laurent_lazy_determinant_fixtures.jl"` to the internal test list.
- Keep `docs/superpowers/specs/2026-06-28-issue-154-lazy-laurent-determinant-fixtures-design.md` and this plan as workflow artifacts.

---

### Task 1: Add the Failing Validator

**Files:**
- Create: `test/internal/laurent_lazy_determinant_fixtures.jl`

**Interfaces:**
- Consumes: planned `LaurentLazyDeterminantCases.catalog()`.
- Produces: `validate_laurent_lazy_determinant_fixture(entry)`, `validate_laurent_lazy_determinant_catalog(catalog)`, and `_fixture_supports_lazy_determinant_correction(entry)`.

- [ ] **Step 1: Write the validator test first**

Create `test/internal/laurent_lazy_determinant_fixtures.jl` with this initial test shell:

```julia
using Test
using Suslin
using Oscar

const LAURENT_LAZY_DETERMINANT_FIXTURE_PATH =
    joinpath(@__DIR__, "..", "fixtures", "laurent_lazy_determinant_cases.jl")

const REQUIRED_LAZY_DETERMINANT_FIELDS = (
    :id,
    :kind,
    :ring,
    :dimensions,
    :inputs,
    :determinant_profile,
    :expected_correction,
    :normalizations,
    :negative_control,
    :provenance,
    :consumer_test_ids,
)

function _lazy_field(entry, field::Symbol)
    hasproperty(entry, field) || throw(ArgumentError("lazy determinant fixture $(get(entry, :id, "<unknown>")) missing field $(field)"))
    return getproperty(entry, field)
end

function _matrix_size(A)
    return (nrows(A), ncols(A))
end

function _fixture_supports_lazy_determinant_correction(entry)::Bool
    expected = _lazy_field(entry, :expected_correction)
    hasproperty(expected, :supported) || throw(ArgumentError("fixture $(entry.id) missing correction support flag"))
    hasproperty(entry, :inputs) && hasproperty(entry.inputs, :matrix) ||
        throw(ArgumentError("fixture $(entry.id) missing matrix input"))
    classification = classify_laurent_determinant(entry.inputs.matrix).classification
    actual_supported = classification in (:one, :laurent_monomial_unit, :permutation_sign_unit)
    return expected.supported && actual_supported
end

function validate_laurent_lazy_determinant_fixture(entry)
    for field in REQUIRED_LAZY_DETERMINANT_FIELDS
        _lazy_field(entry, field)
    end
    isempty(entry.consumer_test_ids) && throw(ArgumentError("fixture $(entry.id) must record at least one consumer test id"))
    _matrix_size(entry.inputs.matrix) == entry.dimensions.matrix ||
        throw(ArgumentError("fixture $(entry.id) matrix dimensions do not match metadata"))
    return true
end

function validate_laurent_lazy_determinant_catalog(catalog)
    hasproperty(catalog, :cases) || throw(ArgumentError("lazy determinant catalog missing cases"))
    ids = [entry.id for entry in catalog.cases]
    length(ids) == length(unique(ids)) || throw(ArgumentError("lazy determinant fixture ids must be unique"))
    for entry in catalog.cases
        validate_laurent_lazy_determinant_fixture(entry)
    end
    return true
end

@testset "lazy Laurent determinant fixture catalog" begin
    @test isfile(LAURENT_LAZY_DETERMINANT_FIXTURE_PATH)

    include(LAURENT_LAZY_DETERMINANT_FIXTURE_PATH)
    catalog = LaurentLazyDeterminantCases.catalog()
    @test validate_laurent_lazy_determinant_catalog(catalog)
end
```

- [ ] **Step 2: Run the focused validator and verify the RED result**

Run:

```bash
julia --project=. -e 'include("test/internal/laurent_lazy_determinant_fixtures.jl")'
```

Expected: the command exits nonzero because `test/fixtures/laurent_lazy_determinant_cases.jl` does not exist yet.

---

### Task 2: Add the Fixture Catalog

**Files:**
- Create: `test/fixtures/laurent_lazy_determinant_cases.jl`
- Modify: `test/internal/laurent_lazy_determinant_fixtures.jl`

**Interfaces:**
- Consumes: `ToricBuilderIssue38Cases.catalog()`, `Suslin.classify_laurent_determinant`, and `Suslin.normalize_laurent_gl_matrix`.
- Produces: `LaurentLazyDeterminantCases.catalog()` returning four entries.

- [ ] **Step 1: Implement the fixture module**

Create `test/fixtures/laurent_lazy_determinant_cases.jl` with module `LaurentLazyDeterminantCases`. It must include the existing Issue #38 fixture file:

```julia
module LaurentLazyDeterminantCases

using Oscar
using Suslin

include("toricbuilder_issue38_cases.jl")

_matrix_dimensions(A) = (nrows(A), ncols(A))

function _ring_metadata(R, generators, variables, description)
    return (;
        description,
        object = R,
        generators,
        variables,
    )
end
```

Use these helper contracts:

```julia
_synthetic_provenance(description) = (;
    source = :synthetic,
    issue = "#154",
    description,
)

function _determinant_profile(A)
    profile = Suslin.classify_laurent_determinant(A)
    return (;
        expected_determinant = profile.determinant,
        expected_class = profile.classification,
        monomial_exponents = profile.monomial_exponents,
        monomial_coefficient = profile.monomial_coefficient,
    )
end

function _expected_correction(profile)
    if profile.expected_class == :one
        return (; supported = true, kind = :identity, supports = (:row_core, :column_core), unsupported_reason = nothing)
    elseif profile.expected_class == :laurent_monomial_unit
        return (; supported = true, kind = :monomial_unit_diagonal, supports = (:row_core, :column_core), unsupported_reason = nothing)
    else
        return (; supported = false, kind = :unsupported, supports = (), unsupported_reason = profile.expected_class)
    end
end
```

Construct supported normalizations with exact row and column cores:

```julia
function _supported_normalizations(A)
    R = base_ring(A)
    n = nrows(A)
    row = Suslin.normalize_laurent_gl_matrix(A)
    determinant = det(A)
    column_correction = diagonal_matrix(R, [i == 1 ? inv(determinant) : one(R) for i in 1:n])
    return (;
        row = (; core = row.normalized_matrix, normalization = row),
        column = (; core = A * column_correction, correction_factor = column_correction),
    )
end
```

The catalog must create:

```julia
R, (x, y) = Suslin.suslin_laurent_polynomial_ring(GF(2), ["x", "y"])
determinant_one = matrix(R, [one(R) x; zero(R) one(R)])
monomial_unit = matrix(R, [x^-1 * y one(R); zero(R) one(R)])
non_unit = matrix(R, [x + one(R) zero(R); zero(R) one(R)])
```

Wrap Issue #38 like this:

```julia
issue38 = only(ToricBuilderIssue38Cases.catalog().cases)
```

The Issue #38 entry must use `issue38.inputs.matrix`, `issue38.ring`, `issue38.determinant_profile`, `issue38.normalizations`, and provenance that records `source_fixture_id = issue38.id`.

- [ ] **Step 2: Complete exact validator checks**

Extend `test/internal/laurent_lazy_determinant_fixtures.jl` so each entry validates:

```julia
actual_profile = classify_laurent_determinant(entry.inputs.matrix)
actual_profile.determinant == entry.determinant_profile.expected_determinant
actual_profile.classification == entry.determinant_profile.expected_class
det(entry.inputs.matrix) == entry.determinant_profile.expected_determinant
entry.expected_correction.supported == (actual_profile.classification in (:one, :laurent_monomial_unit, :permutation_sign_unit))
```

For supported entries, validate:

```julia
normalizations = entry.normalizations
row = normalizations.row
column = normalizations.column
verify_laurent_gl_normalization(entry.inputs.matrix, row.normalization)
row.core == row.normalization.normalized_matrix
entry.inputs.matrix * column.correction_factor == column.core
det(row.core) == one(base_ring(row.core))
det(column.core) == one(base_ring(column.core))
verify_factorization(row.core, elementary_factorization(row.core))
verify_factorization(column.core, elementary_factorization(column.core))
```

For unsupported entries, validate:

```julia
entry.normalizations === nothing
!_fixture_supports_lazy_determinant_correction(entry)
```

Add an Issue #38 drift helper that recomputes `only(LaurentLazyDeterminantCases.ToricBuilderIssue38Cases.catalog().cases)` and checks:

```julia
entry.inputs.matrix == source.inputs.matrix
entry.determinant_profile.expected_determinant == source.determinant_profile.expected_determinant
entry.determinant_profile.expected_determinant == u * v
entry.normalizations.row.core == source.normalizations.row.core
entry.normalizations.column.core == source.normalizations.column.core
```

- [ ] **Step 3: Run the focused validator and verify the GREEN result**

Run:

```bash
julia --project=. -e 'include("test/internal/laurent_lazy_determinant_fixtures.jl")'
```

Expected: the command exits 0 and prints a testset named `lazy Laurent determinant fixture catalog`.

---

### Task 3: Register the Test and Add Negative Controls

**Files:**
- Modify: `test/internal/laurent_lazy_determinant_fixtures.jl`
- Modify: `test/runtests.jl`

**Interfaces:**
- Consumes: `validate_laurent_lazy_determinant_fixture(entry)` and `_fixture_supports_lazy_determinant_correction(entry)`.
- Produces: registered internal test coverage and explicit negative controls.

- [ ] **Step 1: Add required negative controls**

In the focused testset, build `fixture_by_id = Dict(entry.id => entry for entry in catalog.cases)` and assert all required IDs exist:

```julia
@test haskey(fixture_by_id, "issue-38-q-block-lazy-determinant")
@test haskey(fixture_by_id, "determinant-one-triangular")
@test haskey(fixture_by_id, "monomial-unit-row-column-cores")
@test haskey(fixture_by_id, "non-unit-determinant-negative")
```

Add the deliberate bad class mutation:

```julia
monomial = fixture_by_id["monomial-unit-row-column-cores"]
bad_monomial = merge(
    monomial,
    (; determinant_profile = merge(monomial.determinant_profile, (; expected_class = :one))),
)
@test_throws ArgumentError validate_laurent_lazy_determinant_fixture(bad_monomial)
```

Assert the non-unit entry is not supported:

```julia
non_unit = fixture_by_id["non-unit-determinant-negative"]
@test !non_unit.expected_correction.supported
@test !_fixture_supports_lazy_determinant_correction(non_unit)
```

- [ ] **Step 2: Register the internal test**

Add this line to the `internal` group in `test/runtests.jl` after the existing Laurent normalization checks:

```julia
"internal/laurent_lazy_determinant_fixtures.jl",
```

- [ ] **Step 3: Run focused and internal verification**

Run:

```bash
julia --project=. -e 'include("test/internal/laurent_lazy_determinant_fixtures.jl")'
julia --project=. test/runtests.jl internal
```

Expected: both commands exit 0.

---

### Task 4: Package Verification, Review, and Commit

**Files:**
- Verify all files changed by Tasks 1-3.

**Interfaces:**
- Consumes: focused and internal test results.
- Produces: committed implementation ready for PR creation.

- [ ] **Step 1: Run required package verification**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: the command exits 0.

- [ ] **Step 2: Inspect repository state**

Run:

```bash
git status --short
git diff -- test/fixtures/laurent_lazy_determinant_cases.jl test/internal/laurent_lazy_determinant_fixtures.jl test/runtests.jl docs/superpowers/plans/2026-06-28-issue-154-lazy-laurent-determinant-fixtures.md
```

Expected: only the new fixture, new internal test, test runner registration, and plan file are uncommitted.

- [ ] **Step 3: Commit the implementation**

Run:

```bash
git add test/fixtures/laurent_lazy_determinant_cases.jl test/internal/laurent_lazy_determinant_fixtures.jl test/runtests.jl docs/superpowers/plans/2026-06-28-issue-154-lazy-laurent-determinant-fixtures.md
git commit -m "test: add lazy Laurent determinant fixtures"
```

Expected: commit succeeds without adding `Manifest.toml`.

## Self-Review

- The plan covers all four required fixture classes.
- The plan records exact file paths, case IDs, negative controls, and verification commands.
- The plan preserves the issue scope: fixture catalog and internal validator only.
- The Issue #38 matrix is wrapped through the existing fixture module rather than copied.
