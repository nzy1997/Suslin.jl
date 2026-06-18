# Issue 6 Laurent Fixture Catalog Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a small shared exact Laurent fixture catalog and validator for later Laurent algebra tests.

**Architecture:** Keep the catalog under `test/fixtures/` as test support, not public API. The catalog builds structured named-tuple fixtures with exact Oscar Laurent objects, while `test/internal/laurent_fixtures.jl` validates metadata, exact relations, determinant profiles, and negative controls.

**Tech Stack:** Julia, Oscar Laurent polynomial rings, Suslin test helpers, Test stdlib.

## Global Constraints

- Do not add ToricBuilder as a Suslin dependency.
- Do not add a public Suslin fixture API.
- Catalog file is `test/fixtures/laurent_cases.jl`.
- Validator file is `test/internal/laurent_fixtures.jl`.
- Required cases are one solvable linear system, one unsolvable linear system, one negative-exponent normalization case, and one ToricBuilder-motivated case from Issue 19.
- Use the smaller Issue 19 `Pinv` case as the ToricBuilder-motivated fixture.
- The ToricBuilder provenance commit is `fa7f82252d42fdc0b2726bc48af24ac4c70a8d73`.
- Focused validator command is `julia --project=. -e 'include("test/internal/laurent_fixtures.jl")'`.
- Package verification command is `julia --project=. -e 'using Pkg; Pkg.test()'`.
- Full-suite verification command is `julia --project=. test/runtests.jl all`.

---

## File Structure

- Create `test/internal/laurent_fixtures.jl`: validator functions, focused tests, and negative controls. This file owns validation logic for the shared catalog.
- Create `test/fixtures/laurent_cases.jl`: exact fixture catalog module. This file owns fixture construction and reuses the existing Issue 19 ToricBuilder fixture.
- Modify `test/runtests.jl`: add `internal/laurent_fixtures.jl` to the default internal test group.

---

### Task 1: Write the Failing Validator Test

**Files:**
- Create: `test/internal/laurent_fixtures.jl`

**Interfaces:**
- Consumes: planned catalog module `LaurentFixtureCatalog.catalog()` from `test/fixtures/laurent_cases.jl`.
- Produces: `validate_laurent_fixture(entry)` and `validate_laurent_fixture_catalog(catalog)` for internal test validation.

- [ ] **Step 1: Write the failing validator file**

Create `test/internal/laurent_fixtures.jl` with this content:

```julia
using Test
using Suslin
using Oscar

const LAURENT_FIXTURE_CATALOG_PATH = joinpath(@__DIR__, "..", "fixtures", "laurent_cases.jl")
const REQUIRED_LAURENT_FIXTURE_FIELDS = (
    :id,
    :kind,
    :ring_constructor,
    :ring,
    :dimensions,
    :inputs,
    :expected_relation,
    :provenance,
    :determinant_profile,
    :consumer_test_ids,
)

function _fixture_field(entry, field::Symbol)
    hasproperty(entry, field) || throw(ArgumentError("fixture entry missing field $(field)"))
    return getproperty(entry, field)
end

function _fixture_matrix_size(A)
    return (nrows(A), ncols(A))
end

function _fixture_determinant_classification(A)
    R = base_ring(A)
    d = det(A)
    d == one(R) && return "one"
    is_unit(d) && return "laurent_monomial_unit"
    return "non-unit"
end

function _assert_fixture_metadata(entry)
    for field in REQUIRED_LAURENT_FIXTURE_FIELDS
        _fixture_field(entry, field)
    end

    isempty(entry.consumer_test_ids) && throw(ArgumentError("fixture $(entry.id) must record at least one consumer test id"))
    hasproperty(entry.ring_constructor, :function_name) || throw(ArgumentError("fixture $(entry.id) missing ring constructor function name"))
    hasproperty(entry.ring_constructor, :coefficient) || throw(ArgumentError("fixture $(entry.id) missing ring constructor coefficient"))
    hasproperty(entry.ring_constructor, :variables) || throw(ArgumentError("fixture $(entry.id) missing ring constructor variables"))
    hasproperty(entry.ring, :description) || throw(ArgumentError("fixture $(entry.id) missing ring description"))
    hasproperty(entry.ring, :object) || throw(ArgumentError("fixture $(entry.id) missing ring object"))
    Suslin._require_laurent_polynomial_ring(entry.ring.object; label="fixture $(entry.id) ring")
    return true
end

function _assert_fixture_determinant_profile(entry)
    profile = entry.determinant_profile
    hasproperty(profile, :relevant) || throw(ArgumentError("fixture $(entry.id) missing determinant relevance"))
    profile.relevant || return true
    hasproperty(profile, :expected_class) || throw(ArgumentError("fixture $(entry.id) missing determinant class"))
    hasproperty(entry.inputs, :matrix) || throw(ArgumentError("fixture $(entry.id) determinant profile requires a matrix input"))

    actual_class = _fixture_determinant_classification(entry.inputs.matrix)
    actual_class == profile.expected_class || throw(ArgumentError("fixture $(entry.id) determinant classification $(profile.expected_class) does not match $(actual_class)"))
    return true
end

function _assert_solvable_linear_system_fixture(entry)
    A = entry.inputs.matrix
    rhs = entry.inputs.rhs
    solution = entry.inputs.expected_solution

    _fixture_matrix_size(A) == entry.dimensions.matrix || throw(ArgumentError("fixture $(entry.id) matrix dimensions do not match metadata"))
    _fixture_matrix_size(rhs) == entry.dimensions.rhs || throw(ArgumentError("fixture $(entry.id) rhs dimensions do not match metadata"))
    _fixture_matrix_size(solution) == entry.dimensions.solution || throw(ArgumentError("fixture $(entry.id) solution dimensions do not match metadata"))
    ncols(rhs) == 1 || throw(ArgumentError("fixture $(entry.id) rhs must be a column matrix"))
    ncols(solution) == 1 || throw(ArgumentError("fixture $(entry.id) solution must be a column matrix"))
    base_ring(A) == base_ring(rhs) || throw(ArgumentError("fixture $(entry.id) matrix and rhs use different rings"))
    base_ring(A) == base_ring(solution) || throw(ArgumentError("fixture $(entry.id) matrix and solution use different rings"))
    entry.expected_relation.kind == :linear_system_solution || throw(ArgumentError("fixture $(entry.id) has wrong expected relation kind"))
    A * solution == rhs || throw(ArgumentError("fixture $(entry.id) claimed solution does not satisfy matrix * solution == rhs"))
    return true
end

function _assert_unsolvable_linear_system_fixture(entry)
    A = entry.inputs.matrix
    rhs = entry.inputs.rhs
    certificate = entry.inputs.unsolvability_certificate

    _fixture_matrix_size(A) == entry.dimensions.matrix || throw(ArgumentError("fixture $(entry.id) matrix dimensions do not match metadata"))
    _fixture_matrix_size(rhs) == entry.dimensions.rhs || throw(ArgumentError("fixture $(entry.id) rhs dimensions do not match metadata"))
    base_ring(A) == base_ring(rhs) || throw(ArgumentError("fixture $(entry.id) matrix and rhs use different rings"))
    entry.expected_relation.kind == :linear_system_no_solution || throw(ArgumentError("fixture $(entry.id) has wrong expected relation kind"))
    certificate.kind == :zero_matrix_nonzero_rhs || throw(ArgumentError("fixture $(entry.id) uses an unsupported unsolvability certificate"))

    R = base_ring(A)
    A == zero_matrix(R, nrows(A), ncols(A)) || throw(ArgumentError("fixture $(entry.id) certificate requires a zero matrix"))
    i, j = certificate.rhs_index
    rhs[i, j] != zero(R) || throw(ArgumentError("fixture $(entry.id) certificate requires a nonzero rhs entry"))
    return true
end

function _assert_negative_normalization_fixture(entry)
    input_vector = entry.inputs.vector
    normalized_vector = entry.inputs.normalized_vector
    normalization_unit = entry.inputs.normalization_unit

    _fixture_matrix_size(input_vector) == entry.dimensions.vector || throw(ArgumentError("fixture $(entry.id) vector dimensions do not match metadata"))
    _fixture_matrix_size(normalized_vector) == entry.dimensions.normalized_vector || throw(ArgumentError("fixture $(entry.id) normalized dimensions do not match metadata"))
    base_ring(input_vector) == base_ring(normalized_vector) || throw(ArgumentError("fixture $(entry.id) vectors use different rings"))
    parent(normalization_unit) == base_ring(input_vector) || throw(ArgumentError("fixture $(entry.id) normalization unit uses the wrong ring"))
    entry.expected_relation.kind == :negative_exponent_normalization || throw(ArgumentError("fixture $(entry.id) has wrong expected relation kind"))
    normalization_unit * input_vector == normalized_vector || throw(ArgumentError("fixture $(entry.id) normalization relation failed"))
    return true
end

function _assert_toricbuilder_relation_fixture(entry)
    A = entry.inputs.matrix
    source = entry.inputs.source_matrix

    _fixture_matrix_size(A) == entry.dimensions.matrix || throw(ArgumentError("fixture $(entry.id) matrix dimensions do not match metadata"))
    _fixture_matrix_size(source) == entry.dimensions.source_matrix || throw(ArgumentError("fixture $(entry.id) source dimensions do not match metadata"))
    base_ring(A) == base_ring(source) || throw(ArgumentError("fixture $(entry.id) relation matrices use different rings"))
    entry.expected_relation.kind == :toricbuilder_left_inverse || throw(ArgumentError("fixture $(entry.id) has wrong expected relation kind"))
    entry.provenance.issue == "#19" || throw(ArgumentError("fixture $(entry.id) must trace to Issue 19"))
    entry.provenance.toricbuilder_commit == "fa7f82252d42fdc0b2726bc48af24ac4c70a8d73" || throw(ArgumentError("fixture $(entry.id) has unexpected ToricBuilder commit"))

    R = base_ring(A)
    source * A == identity_matrix(R, nrows(A)) || throw(ArgumentError("fixture $(entry.id) ToricBuilder relation failed"))
    return true
end

function validate_laurent_fixture(entry)
    _assert_fixture_metadata(entry)
    _assert_fixture_determinant_profile(entry)

    if entry.kind == :solvable_linear_system
        _assert_solvable_linear_system_fixture(entry)
    elseif entry.kind == :unsolvable_linear_system
        _assert_unsolvable_linear_system_fixture(entry)
    elseif entry.kind == :negative_exponent_normalization
        _assert_negative_normalization_fixture(entry)
    elseif entry.kind == :toricbuilder_relation
        _assert_toricbuilder_relation_fixture(entry)
    else
        throw(ArgumentError("fixture $(entry.id) has unsupported kind $(entry.kind)"))
    end

    return true
end

function validate_laurent_fixture_catalog(catalog)
    hasproperty(catalog, :cases) || throw(ArgumentError("Laurent fixture catalog missing cases"))
    isempty(catalog.cases) && throw(ArgumentError("Laurent fixture catalog must not be empty"))
    ids = [entry.id for entry in catalog.cases]
    length(ids) == length(unique(ids)) || throw(ArgumentError("Laurent fixture ids must be unique"))

    for entry in catalog.cases
        validate_laurent_fixture(entry)
    end
    return true
end

@testset "shared Laurent fixture catalog" begin
    @test isfile(LAURENT_FIXTURE_CATALOG_PATH)

    include(LAURENT_FIXTURE_CATALOG_PATH)
    catalog = LaurentFixtureCatalog.catalog()

    @test validate_laurent_fixture_catalog(catalog)

    fixture_by_id = Dict(entry.id => entry for entry in catalog.cases)
    @test haskey(fixture_by_id, "laurent-linear-system-solvable")
    @test haskey(fixture_by_id, "laurent-linear-system-unsolvable")
    @test haskey(fixture_by_id, "laurent-negative-exponent-normalization")
    @test haskey(fixture_by_id, "toricbuilder-factor-toric-block-3-pinv")

    solvable = fixture_by_id["laurent-linear-system-solvable"]
    bad_solution = copy(solvable.inputs.expected_solution)
    bad_solution[1, 1] += one(base_ring(bad_solution))
    bad_solvable = merge(solvable, (; inputs = merge(solvable.inputs, (; expected_solution = bad_solution))))
    @test_throws ArgumentError validate_laurent_fixture(bad_solvable)

    bad_det = merge(solvable, (; determinant_profile = merge(solvable.determinant_profile, (; expected_class = "non-unit"))))
    @test_throws ArgumentError validate_laurent_fixture(bad_det)

    toricbuilder = fixture_by_id["toricbuilder-factor-toric-block-3-pinv"]
    corrupted_source = copy(toricbuilder.inputs.source_matrix)
    corrupted_source[1, 1] += one(base_ring(corrupted_source))
    bad_toricbuilder = merge(toricbuilder, (; inputs = merge(toricbuilder.inputs, (; source_matrix = corrupted_source))))
    @test_throws ArgumentError validate_laurent_fixture(bad_toricbuilder)
end
```

- [ ] **Step 2: Run the focused validator and verify it fails for the missing catalog**

Run:

```bash
julia --project=. -e 'include("test/internal/laurent_fixtures.jl")'
```

Expected: FAIL because `test/fixtures/laurent_cases.jl` does not exist yet.

- [ ] **Step 3: Commit the failing test**

Run:

```bash
git add test/internal/laurent_fixtures.jl
git commit -m "test: add laurent fixture validator"
```

---

### Task 2: Add the Catalog and Wire the Validator Into the Suite

**Files:**
- Create: `test/fixtures/laurent_cases.jl`
- Modify: `test/runtests.jl`

**Interfaces:**
- Consumes: `ToricBuilderFactorToricBlock3Fixture.fixture()` from `test/fixtures/toricbuilder_factor_toric_block_3.jl`.
- Produces: `LaurentFixtureCatalog.catalog()` returning `(; ring, cases)`.

- [ ] **Step 1: Create the shared fixture catalog**

Create `test/fixtures/laurent_cases.jl` with this content:

```julia
module LaurentFixtureCatalog

using Oscar
using Suslin

include("toricbuilder_factor_toric_block_3.jl")

function _column_matrix(R, values)
    M = zero_matrix(R, length(values), 1)
    for (i, value) in enumerate(values)
        M[i, 1] = value
    end
    return M
end

function _ring_metadata(R, x, y)
    return (;
        description = "GF(2)[x^+/-1, y^+/-1]",
        object = R,
        generators = (x, y),
    )
end

function _ring_constructor_metadata()
    return (;
        function_name = :suslin_laurent_polynomial_ring,
        coefficient = "GF(2)",
        variables = ("x", "y"),
    )
end

function _synthetic_provenance(description)
    return (;
        source = :synthetic,
        issue = "#6",
        description,
    )
end

function _pinv_toricbuilder_case(ring, ring_constructor)
    toricbuilder_fixture = ToricBuilderFactorToricBlock3Fixture.fixture()
    pinv = only(filter(entry -> entry.toricbuilder_role == "Pinv", toricbuilder_fixture.cases))

    return (;
        id = "toricbuilder-factor-toric-block-3-pinv",
        kind = :toricbuilder_relation,
        ring_constructor,
        ring,
        dimensions = (;
            matrix = pinv.size,
            source_matrix = pinv.size,
        ),
        inputs = (;
            matrix = pinv.matrix,
            source_matrix = pinv.source_matrix,
        ),
        expected_relation = (;
            kind = :toricbuilder_left_inverse,
            description = pinv.relation_description,
        ),
        provenance = (;
            source = :toricbuilder_contract_fixture,
            issue = "#19",
            fixture_id = pinv.name,
            toricbuilder_role = pinv.toricbuilder_role,
            toricbuilder_commit = pinv.provenance.toricbuilder_commit,
            generation_command = pinv.provenance.generation_command,
        ),
        determinant_profile = (;
            relevant = true,
            expected_class = pinv.determinant_classification,
        ),
        consumer_test_ids = ("issue-6-laurent-fixtures", "issue-19-toricbuilder-contract"),
    )
end

function catalog()
    R, (x, y) = Suslin.suslin_laurent_polynomial_ring(GF(2), ["x", "y"])
    ring = _ring_metadata(R, x, y)
    ring_constructor = _ring_constructor_metadata()

    solvable_matrix = matrix(R, [
        one(R) x;
        zero(R) one(R)
    ])
    solvable_solution = _column_matrix(R, [one(R), y^-1])
    solvable_rhs = solvable_matrix * solvable_solution

    unsolvable_matrix = zero_matrix(R, 1, 1)
    unsolvable_rhs = _column_matrix(R, [one(R)])

    negative_vector = _column_matrix(R, [x^-1 * y, x^-2])
    normalization_unit = x^2
    normalized_vector = _column_matrix(R, [x * y, one(R)])

    return (;
        ring,
        cases = [
            (;
                id = "laurent-linear-system-solvable",
                kind = :solvable_linear_system,
                ring_constructor,
                ring,
                dimensions = (;
                    matrix = (2, 2),
                    rhs = (2, 1),
                    solution = (2, 1),
                ),
                inputs = (;
                    matrix = solvable_matrix,
                    rhs = solvable_rhs,
                    expected_solution = solvable_solution,
                ),
                expected_relation = (;
                    kind = :linear_system_solution,
                    description = "matrix * expected_solution == rhs",
                ),
                provenance = _synthetic_provenance("triangular 2x2 Laurent linear system with exact solution"),
                determinant_profile = (;
                    relevant = true,
                    expected_class = "one",
                ),
                consumer_test_ids = ("issue-6-laurent-fixtures", "issue-8-linear-laurent-tests"),
            ),
            (;
                id = "laurent-linear-system-unsolvable",
                kind = :unsolvable_linear_system,
                ring_constructor,
                ring,
                dimensions = (;
                    matrix = (1, 1),
                    rhs = (1, 1),
                ),
                inputs = (;
                    matrix = unsolvable_matrix,
                    rhs = unsolvable_rhs,
                    unsolvability_certificate = (;
                        kind = :zero_matrix_nonzero_rhs,
                        rhs_index = (1, 1),
                    ),
                ),
                expected_relation = (;
                    kind = :linear_system_no_solution,
                    description = "zero matrix cannot produce a nonzero right-hand side",
                ),
                provenance = _synthetic_provenance("zero 1x1 Laurent linear system with nonzero rhs"),
                determinant_profile = (;
                    relevant = true,
                    expected_class = "non-unit",
                ),
                consumer_test_ids = ("issue-6-laurent-fixtures", "issue-9-linear-laurent-rejections"),
            ),
            (;
                id = "laurent-negative-exponent-normalization",
                kind = :negative_exponent_normalization,
                ring_constructor,
                ring,
                dimensions = (;
                    vector = (2, 1),
                    normalized_vector = (2, 1),
                ),
                inputs = (;
                    vector = negative_vector,
                    normalization_unit,
                    normalized_vector,
                ),
                expected_relation = (;
                    kind = :negative_exponent_normalization,
                    description = "x^2 * [x^-1*y, x^-2] == [x*y, 1]",
                ),
                provenance = _synthetic_provenance("minimal vector normalization with negative x exponents"),
                determinant_profile = (;
                    relevant = false,
                    expected_class = nothing,
                ),
                consumer_test_ids = ("issue-6-laurent-fixtures", "issue-12-negative-exponent-normalization"),
            ),
            _pinv_toricbuilder_case(ring, ring_constructor),
        ],
    )
end

end
```

- [ ] **Step 2: Wire the validator into default internal tests**

In `test/runtests.jl`, change the internal group list from:

```julia
    "internal" => [
        "internal/rings.jl",
        "internal/laurent_rings.jl",
        "internal/toricbuilder_contract.jl",
    ],
```

to:

```julia
    "internal" => [
        "internal/rings.jl",
        "internal/laurent_rings.jl",
        "internal/laurent_fixtures.jl",
        "internal/toricbuilder_contract.jl",
    ],
```

- [ ] **Step 3: Run the focused validator and verify it passes**

Run:

```bash
julia --project=. -e 'include("test/internal/laurent_fixtures.jl")'
```

Expected: PASS, including the four required fixtures and all negative controls.

- [ ] **Step 4: Run the package test entry point**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: PASS for the default public and internal groups.

- [ ] **Step 5: Run the documented full suite**

Run:

```bash
julia --project=. test/runtests.jl all
```

Expected: PASS for public, internal, and expert groups.

- [ ] **Step 6: Commit the catalog implementation**

Run:

```bash
git add test/fixtures/laurent_cases.jl test/internal/laurent_fixtures.jl test/runtests.jl
git commit -m "test: add shared laurent fixture catalog"
```

---

## Plan Self-Review

- Spec coverage: Task 1 creates the validator and negative controls; Task 2 creates the loadable catalog, wires it into default internal tests, and verifies focused, package, and full-suite commands.
- Placeholder scan: no placeholder markers are intentionally present.
- Type consistency: the plan consistently uses `LaurentFixtureCatalog.catalog()`, `validate_laurent_fixture(entry)`, and `validate_laurent_fixture_catalog(catalog)`.
