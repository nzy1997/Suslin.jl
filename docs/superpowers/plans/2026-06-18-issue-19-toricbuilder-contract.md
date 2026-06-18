# Issue 19 ToricBuilder Contract Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a small documented ToricBuilder integration contract with exact GF(2) Laurent fixtures and verification tests.

**Architecture:** Keep the contract as documentation and test fixture data. The fixture is a Julia file under `test/fixtures/` that builds exact Oscar Laurent matrices from sparse nonzero entries; internal tests validate metadata, determinant classification, exact inverse relations, and current staged unsupported Suslin behavior.

**Tech Stack:** Julia, Oscar Laurent polynomial rings, Documenter Markdown, Test stdlib.

## Global Constraints

- Do not add ToricBuilder as a Suslin dependency.
- Do not implement full factorization for ToricBuilder matrices in this issue.
- Fixture ring is `GF(2)[x^+/-1, y^+/-1]`.
- Fixture provenance commit is `fa7f82252d42fdc0b2726bc48af24ac4c70a8d73`.
- Fixture generation entry point is `factor_toric_block(3, x, y, R)`.
- Required fixture entries are `Qinv` (`16 x 16`) and `Pinv` (`8 x 8`).
- Required determinant classification for both fixture entries is `one`.
- Required exact relations are `column_transformation * Qinv == I_16` and `row_transformation * Pinv == I_8`.
- Current Suslin behavior for both entries is staged unsupported input through `ArgumentError`; no factorization output is expected yet.
- Final package verification command is `julia --project=. -e 'using Pkg; Pkg.test()'`.
- Broader full-suite verification command is `julia --project=. test/runtests.jl all`.

---

## File Structure

- Create `docs/src/toricbuilder_contract.md`: documents the consumer boundary, provenance, fixture fields, determinant classes, expected current behavior, and eventual ToricBuilder output need.
- Modify `docs/src/index.md`: add a short link to the contract.
- Modify `docs/make.jl`: include the new page in Documenter pages.
- Create `test/internal/toricbuilder_contract.jl`: focused fixture contract test with red/green TDD evidence and negative controls.
- Create `test/fixtures/toricbuilder_factor_toric_block_3.jl`: exact checked-in fixture data.
- Modify `test/runtests.jl`: include the new internal test file.

---

### Task 1: Add the Documented Contract

**Files:**
- Create: `docs/src/toricbuilder_contract.md`
- Modify: `docs/src/index.md`
- Modify: `docs/make.jl`

**Interfaces:**
- Consumes: design spec `docs/superpowers/specs/2026-06-18-issue-19-toricbuilder-contract-design.md`.
- Produces: Documenter-visible contract page with no code dependency on ToricBuilder.

- [ ] **Step 1: Add the contract page**

Create `docs/src/toricbuilder_contract.md` with these sections and exact facts:

```markdown
# ToricBuilder Integration Contract

## Consumer Boundary

The first ToricBuilder boundary Suslin records is the small output of
`factor_toric_block(3, x, y, R)` over `GF(2)[x^+/-1, y^+/-1]`.
ToricBuilder produces exact Laurent transformation matrices; Suslin records
them as contract fixtures before expanding factorization support.

## Fixture Entries

| Entry | Role | Ring | Size | Determinant classification | Current Suslin status |
| --- | --- | --- | --- | --- | --- |
| `factor_toric_block_3_qinv` | `Qinv` | `GF(2)[x^+/-1, y^+/-1]` | `16 x 16` | `one` | staged unsupported input |
| `factor_toric_block_3_pinv` | `Pinv` | `GF(2)[x^+/-1, y^+/-1]` | `8 x 8` | `one` | staged unsupported input |

## Provenance

- ToricBuilder path: `/Users/nzy/jcode/topological-code-decoupling/julia_code/ToricBuilder`
- ToricBuilder commit: `fa7f82252d42fdc0b2726bc48af24ac4c70a8d73`
- Source function: `src/toric_form/toric_factorization.jl:factor_toric_block`
- Generation command: `factor_toric_block(3, x, y, R)` with
  `R, (x, y) = laurent_polynomial_ring(GF(2), ["x", "y"])`
- Returned block coarse-graining size: `(2, 2)`

## Exact Relations

The checked-in fixture stores the source transformation matrix for each
inverse-style entry. Tests verify these equalities by exact multiplication:

- `column_transformation * Qinv == I_16`
- `row_transformation * Pinv == I_8`

## Determinant Classes

Fixture metadata records one of:

- `one`: determinant is exactly `1`.
- `laurent_monomial_unit`: determinant is an invertible Laurent monomial other
  than `1`.
- `other_unit`: determinant is a non-monomial unit if a future coefficient
  ring admits one.
- `non-unit`: determinant is not a unit.

The current fixture entries both use `one`.

## Suslin Output Contract

ToricBuilder ultimately needs a verified transformation certificate. For this
first contract, that certificate is the exact inverse relation plus determinant
classification. Raw elementary factors and normalized factor metadata are not
required until later implementation issues expand Laurent support.

Current Suslin behavior is deliberately unsupported for these matrices:
`elementary_factorization` accepts only a narrow `3 x 3` univariate polynomial
`SL_3` slice, while these fixtures are larger two-variable Laurent matrices.
```

- [ ] **Step 2: Link the page from docs**

In `docs/src/index.md`, add this under `## Scope` after the existing bullets:

```markdown
See [ToricBuilder Integration Contract](@ref) for the first recorded
consumer-boundary fixture contract.
```

In `docs/make.jl`, change the `pages` list to:

```julia
pages=[
    "Home" => "index.md",
    "ToricBuilder Integration Contract" => "toricbuilder_contract.md",
],
```

- [ ] **Step 3: Verify the docs text**

Run:

```bash
sed -n '1,220p' docs/src/toricbuilder_contract.md
sed -n '1,120p' docs/src/index.md
sed -n '1,80p' docs/make.jl
```

Expected: contract page contains the provenance, determinant classes, exact
relations, and current unsupported behavior; index links the page; `make.jl`
lists both pages.

- [ ] **Step 4: Commit**

```bash
git add docs/src/toricbuilder_contract.md docs/src/index.md docs/make.jl
git commit -m "docs: add toricbuilder integration contract"
```

---

### Task 2: Add Exact Fixture and Contract Tests

**Files:**
- Create: `test/internal/toricbuilder_contract.jl`
- Create: `test/fixtures/toricbuilder_factor_toric_block_3.jl`
- Modify: `test/runtests.jl`

**Interfaces:**
- Consumes: Oscar's `laurent_polynomial_ring`, `zero_matrix`, `identity_matrix`, `det`, and exact matrix multiplication.
- Produces: loadable fixture function `ToricBuilderFactorToricBlock3Fixture.fixture()` returning `(; ring, block_size, provenance, cases)`.

- [ ] **Step 1: Write the failing contract test first**

Create `test/internal/toricbuilder_contract.jl` with validators that assert:

```julia
using Test
using Suslin
using Oscar

const TORICBUILDER_FIXTURE_PATH = joinpath(@__DIR__, "..", "fixtures", "toricbuilder_factor_toric_block_3.jl")

function determinant_classification(A)
    R = base_ring(A)
    d = det(A)
    d == one(R) && return "one"
    is_unit(d) && return "laurent_monomial_unit"
    return "non-unit"
end

function assert_required_field(entry, field::Symbol)
    hasproperty(entry, field) || throw(ArgumentError("fixture entry missing field $(field)"))
    return true
end

function assert_toricbuilder_contract_entry_valid(entry)
    for field in (
        :name,
        :toricbuilder_role,
        :matrix,
        :source_matrix,
        :relation_description,
        :ring,
        :size,
        :determinant_classification,
        :expected_suslin_status,
        :expected_output,
        :provenance,
    )
        assert_required_field(entry, field)
    end

    nrows(entry.matrix) == entry.size[1] || throw(ArgumentError("fixture $(entry.name) row count does not match metadata"))
    ncols(entry.matrix) == entry.size[2] || throw(ArgumentError("fixture $(entry.name) column count does not match metadata"))
    entry.size[1] == entry.size[2] || throw(ArgumentError("fixture $(entry.name) must be square"))
    base_ring(entry.matrix) == base_ring(entry.source_matrix) || throw(ArgumentError("fixture $(entry.name) relation matrices use different rings"))
    entry.ring == "GF(2)[x^+/-1, y^+/-1]" || throw(ArgumentError("fixture $(entry.name) has unexpected ring metadata"))

    actual_class = determinant_classification(entry.matrix)
    actual_class == entry.determinant_classification || throw(ArgumentError("fixture $(entry.name) determinant classification $(entry.determinant_classification) does not match $(actual_class)"))

    R = base_ring(entry.matrix)
    expected_identity = identity_matrix(R, entry.size[1])
    entry.source_matrix * entry.matrix == expected_identity || throw(ArgumentError("fixture $(entry.name) exact ToricBuilder relation failed"))
    return true
end

@testset "toricbuilder contract fixture" begin
    @test isfile(TORICBUILDER_FIXTURE_PATH)

    if isfile(TORICBUILDER_FIXTURE_PATH)
        include(TORICBUILDER_FIXTURE_PATH)
        fixture = ToricBuilderFactorToricBlock3Fixture.fixture()

        @test fixture.ring.description == "GF(2)[x^+/-1, y^+/-1]"
        @test fixture.block_size == (2, 2)
        @test fixture.provenance.toricbuilder_commit == "fa7f82252d42fdc0b2726bc48af24ac4c70a8d73"

        roles = Set(entry.toricbuilder_role for entry in fixture.cases)
        @test "Qinv" in roles
        @test "Pinv" in roles

        for entry in fixture.cases
            @test assert_toricbuilder_contract_entry_valid(entry)
            @test entry.expected_suslin_status == :unsupported_now
            @test entry.expected_output == :verified_transformation_certificate

            err = try
                elementary_factorization(entry.matrix)
                nothing
            catch caught
                caught
            end
            @test err isa ArgumentError
            @test occursin("currently supports only 3x3 matrices", sprint(showerror, err))
        end

        qinv = only(filter(entry -> entry.toricbuilder_role == "Qinv", fixture.cases))
        bad_class = merge(qinv, (; determinant_classification = "non-unit"))
        @test_throws ArgumentError assert_toricbuilder_contract_entry_valid(bad_class)

        corrupted_source = copy(qinv.source_matrix)
        corrupted_source[1, 1] += one(base_ring(corrupted_source))
        bad_relation = merge(qinv, (; source_matrix = corrupted_source))
        @test_throws ArgumentError assert_toricbuilder_contract_entry_valid(bad_relation)
    end
end
```

Modify `test/runtests.jl` so the `internal` group includes:

```julia
"internal/toricbuilder_contract.jl",
```

- [ ] **Step 2: Run the focused RED test**

Run:

```bash
julia --project=. test/runtests.jl internal
```

Expected: FAIL because `test/fixtures/toricbuilder_factor_toric_block_3.jl`
does not exist yet, and the `@test isfile(TORICBUILDER_FIXTURE_PATH)` assertion
fails.

- [ ] **Step 3: Create the exact fixture**

Create `test/fixtures/toricbuilder_factor_toric_block_3.jl` with:

```julia
module ToricBuilderFactorToricBlock3Fixture

using Oscar

function _sparse_laurent_matrix(R, rows::Int, cols::Int, entries)
    M = zero_matrix(R, rows, cols)
    for (i, j, value) in entries
        M[i, j] = value
    end
    return M
end

function fixture()
    R, (x, y) = Oscar.laurent_polynomial_ring(Oscar.GF(2), ["x", "y"])

    qinv = _sparse_laurent_matrix(R, 16, 16, [
        (1, 1, one(R)), (1, 3, x), (1, 5, one(R)), (1, 6, y),
        (2, 2, one(R)), (2, 4, x), (2, 5, one(R)), (2, 6, one(R)),
        (3, 3, x + 1), (3, 5, one(R)), (3, 6, y), (3, 7, one(R)), (3, 8, y),
        (4, 4, one(R)), (5, 6, one(R)), (6, 4, one(R)), (6, 7, one(R)),
        (7, 3, one(R)), (7, 4, one(R)), (8, 6, one(R)), (8, 8, one(R)),
        (9, 9, one(R)), (10, 10, one(R)), (11, 9, one(R)), (11, 10, one(R)), (11, 13, one(R)),
        (12, 11, one(R)), (12, 12, one(R)), (12, 13, x^-1), (12, 15, one(R)),
        (13, 9, y^-1), (13, 10, one(R)), (13, 14, one(R)), (13, 16, one(R)),
        (14, 9, one(R)), (14, 10, one(R)), (14, 13, one(R)), (14, 15, one(R)),
        (15, 9, one(R)), (15, 10, 1 + x^-1), (15, 11, one(R)), (15, 13, 1 + x^-1),
        (16, 9, y^-1), (16, 10, y^-1), (16, 13, y^-1), (16, 16, one(R)),
    ])

    column_transformation = _sparse_laurent_matrix(R, 16, 16, [
        (1, 1, one(R)), (1, 3, one(R)), (1, 5, y), (1, 6, one(R)), (1, 7, one(R)), (1, 8, y),
        (2, 2, one(R)), (2, 3, one(R)), (2, 5, one(R)), (2, 6, one(R)), (2, 7, x + 1), (2, 8, y),
        (3, 4, one(R)), (3, 7, one(R)), (4, 4, one(R)),
        (5, 3, one(R)), (5, 4, x), (5, 6, one(R)), (5, 7, x + 1), (5, 8, y),
        (6, 5, one(R)), (7, 4, one(R)), (7, 6, one(R)), (8, 5, one(R)), (8, 8, one(R)),
        (9, 9, one(R)), (10, 10, one(R)), (11, 9, x^-1), (11, 11, 1 + x^-1), (11, 15, one(R)),
        (12, 10, x^-1), (12, 12, one(R)), (12, 14, one(R)), (12, 15, one(R)),
        (13, 9, one(R)), (13, 10, one(R)), (13, 11, one(R)),
        (14, 9, y^-1), (14, 10, one(R)), (14, 11, y^-1), (14, 13, one(R)), (14, 16, one(R)),
        (15, 11, one(R)), (15, 14, one(R)), (16, 11, y^-1), (16, 16, one(R)),
    ])

    pinv = _sparse_laurent_matrix(R, 8, 8, [
        (1, 1, one(R)), (2, 2, one(R)), (3, 1, one(R)), (3, 3, one(R)),
        (4, 2, one(R)), (4, 3, one(R)), (4, 4, one(R)),
        (5, 7, one(R)), (6, 6, one(R)), (7, 5, one(R)),
        (8, 5, one(R)), (8, 6, x^-1), (8, 7, one(R)), (8, 8, one(R)),
    ])

    row_transformation = _sparse_laurent_matrix(R, 8, 8, [
        (1, 1, one(R)), (2, 2, one(R)), (3, 1, one(R)), (3, 3, one(R)),
        (4, 1, one(R)), (4, 2, one(R)), (4, 3, one(R)), (4, 4, one(R)),
        (5, 7, one(R)), (6, 6, one(R)), (7, 5, one(R)),
        (8, 5, one(R)), (8, 6, x^-1), (8, 7, one(R)), (8, 8, one(R)),
    ])

    ring = (;
        description = "GF(2)[x^+/-1, y^+/-1]",
        base_field = "GF(2)",
        variables = ("x", "y"),
        kind = "multivariate Laurent polynomial ring",
    )
    provenance = (;
        toricbuilder_path = "/Users/nzy/jcode/topological-code-decoupling/julia_code/ToricBuilder",
        toricbuilder_commit = "fa7f82252d42fdc0b2726bc48af24ac4c70a8d73",
        source_function = "src/toric_form/toric_factorization.jl:factor_toric_block",
        generation_command = "factor_toric_block(3, x, y, R) with R, (x, y) = laurent_polynomial_ring(GF(2), [\"x\", \"y\"])",
    )

    return (;
        ring,
        block_size = (2, 2),
        provenance,
        cases = [
            (;
                name = "factor_toric_block_3_qinv",
                toricbuilder_role = "Qinv",
                matrix = qinv,
                source_matrix = column_transformation,
                relation_description = "column_transformation * Qinv == I_16",
                ring = ring.description,
                size = (16, 16),
                determinant_classification = "one",
                expected_suslin_status = :unsupported_now,
                expected_output = :verified_transformation_certificate,
                provenance,
            ),
            (;
                name = "factor_toric_block_3_pinv",
                toricbuilder_role = "Pinv",
                matrix = pinv,
                source_matrix = row_transformation,
                relation_description = "row_transformation * Pinv == I_8",
                ring = ring.description,
                size = (8, 8),
                determinant_classification = "one",
                expected_suslin_status = :unsupported_now,
                expected_output = :verified_transformation_certificate,
                provenance,
            ),
        ],
    )
end

end
```

- [ ] **Step 4: Run the focused GREEN test**

Run:

```bash
julia --project=. test/runtests.jl internal
```

Expected: PASS. The output should show the internal group with the existing
ring tests plus the new ToricBuilder contract fixture tests.

- [ ] **Step 5: Commit**

```bash
git add test/internal/toricbuilder_contract.jl test/fixtures/toricbuilder_factor_toric_block_3.jl test/runtests.jl
git commit -m "test: add toricbuilder contract fixture"
```

---

### Task 3: Verify the Contract End to End

**Files:**
- Read: `docs/src/toricbuilder_contract.md`
- Read: `test/internal/toricbuilder_contract.jl`
- Read: `test/fixtures/toricbuilder_factor_toric_block_3.jl`

**Interfaces:**
- Consumes: committed docs and fixture tests from Tasks 1 and 2.
- Produces: final verification evidence for the PR.

- [ ] **Step 1: Run the focused fixture test**

Run:

```bash
julia --project=. test/runtests.jl internal
```

Expected: exits 0.

- [ ] **Step 2: Run the package test entry point**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: exits 0 and runs the default fast tests.

- [ ] **Step 3: Run the documented full suite**

Run:

```bash
julia --project=. test/runtests.jl all
```

Expected: exits 0 and includes public, internal, and expert groups.

- [ ] **Step 4: Check no local-only artifacts are tracked**

Run:

```bash
git status --short
```

Expected: no untracked ToricBuilder temp copy and no `Manifest.toml` staged or
tracked in the branch.
