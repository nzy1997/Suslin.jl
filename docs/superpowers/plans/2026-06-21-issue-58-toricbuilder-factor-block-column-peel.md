# Issue 58 ToricBuilder Factor Block Column-Peel Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Promote the checked-in ToricBuilder Pinv and Qinv determinant-one fixtures through the recursive Laurent column-peel factorization route.

**Architecture:** Add a public acceptance test that exercises the checked-in fixtures through `elementary_factorization` and the existing column-peel replay certificate. Then update the shared ToricBuilder Laurent problem catalog to record both contract fixtures as supported by the `:laurent_column_peel` path and register the new public test.

**Tech Stack:** Julia, Oscar Laurent polynomial rings and matrices, Suslin column-peel internals from Issue #57, Test stdlib.

## Global Constraints

- Do not add a ToricBuilder runtime dependency.
- Use only `test/fixtures/toricbuilder_factor_toric_block_3.jl`.
- `Pinv` must reach final peel block `I_2`.
- `Qinv` must reach final peel block `[y 0; x*y^-1 + y^-1 y^-1]`.
- Both `elementary_factorization(Pinv)` and `elementary_factorization(Qinv)` must return nonempty factor sequences.
- Both returned factor sequences must satisfy `verify_factorization(matrix, factors) == true`.
- Negative controls must fail exact verification rather than treating the fixture as supported.
- Do not support arbitrary determinant-one Laurent matrices beyond the checked-in fixture family.
- Do not commit `Manifest.toml`.
- Focused verification command is `julia --project=. -e 'include("test/public/toricbuilder_factor_toric_block_acceptance.jl")'`.
- Package verification command is `julia --project=. -e 'using Pkg; Pkg.test()'`.

---

## File Structure

- Create `test/public/toricbuilder_factor_toric_block_acceptance.jl`: public acceptance test for Pinv/Qinv exact factorization, peel replay final blocks, and negative controls.
- Modify `test/fixtures/toricbuilder_laurent_problem_catalog.jl`: promote the Pinv/Qinv catalog entries from contract-only status to `:supported_column_peel`.
- Modify `test/internal/toricbuilder_problem_catalog.jl`: allow and assert the promoted catalog status and factorization path metadata.
- Modify `test/runtests.jl`: register the new public acceptance test.

---

### Task 1: Add the Failing Public Acceptance Test

**Files:**
- Create: `test/public/toricbuilder_factor_toric_block_acceptance.jl`

**Interfaces:**
- Consumes: `ToricBuilderFactorToricBlock3Fixture.fixture()`, `elementary_factorization`, `verify_factorization`, `Suslin._factor_laurent_sl_column_peel`, and `Suslin._verify_laurent_column_peel_replay`.
- Produces: focused public coverage for Issue #58 Pinv/Qinv column-peel support.

- [ ] **Step 1: Write the failing test**

Create `test/public/toricbuilder_factor_toric_block_acceptance.jl`:

```julia
using Suslin
using Test
using Oscar

include("../fixtures/toricbuilder_factor_toric_block_3.jl")

function _issue58_product(factors, R, n::Int)
    product = identity_matrix(R, n)
    for factor in factors
        product *= factor
    end
    return product
end

function _issue58_entry_by_role(fixture, role::AbstractString)
    return only(filter(entry -> entry.toricbuilder_role == role, fixture.cases))
end

function _issue58_assert_step_replay(step)
    R = base_ring(step.input_matrix)
    left_product = _issue58_product(step.left_factors, R, step.dimension)
    right_product = _issue58_product(step.right_factors, R, step.dimension)
    @test left_product * step.input_matrix * right_product == step.peeled_matrix
    @test step.peeled_matrix[step.dimension, step.dimension] == one(R)
    @test all(step.peeled_matrix[row, step.dimension] == zero(R) for row in 1:(step.dimension - 1))
    @test all(step.peeled_matrix[step.dimension, col] == zero(R) for col in 1:(step.dimension - 1))
    @test step.next_block == matrix(R, [
        step.peeled_matrix[row, col]
        for row in 1:(step.dimension - 1), col in 1:(step.dimension - 1)
    ])
end

function _issue58_corrupt_source_relation(entry)
    corrupted = copy(entry.source_matrix)
    corrupted[1, 1] += one(base_ring(corrupted))
    return corrupted
end

function _issue58_replace_first_factor_with_identity(factors, R, n::Int)
    corrupted = copy(factors)
    corrupted[1] = identity_matrix(R, n)
    return corrupted
end

function _issue58_assert_column_peel_entry(entry, expected_dimensions, expected_final_block)
    R = base_ring(entry.matrix)
    n = entry.size[1]

    @test entry.expected_suslin_status == :supported_column_peel
    @test entry.expected_suslin_path == :laurent_column_peel
    @test entry.source_matrix * entry.matrix == identity_matrix(R, n)

    certificate = Suslin._factor_laurent_sl_column_peel(entry.matrix)
    @test Suslin._verify_laurent_column_peel_replay(certificate)
    @test [step.dimension for step in certificate.peel_steps] == expected_dimensions
    @test certificate.final_block == expected_final_block
    for step in certificate.peel_steps
        _issue58_assert_step_replay(step)
    end

    factors = elementary_factorization(entry.matrix)
    @test !isempty(factors)
    @test verify_factorization(entry.matrix, factors)

    corrupted_factors = _issue58_replace_first_factor_with_identity(factors, R, n)
    @test !verify_factorization(entry.matrix, corrupted_factors)

    corrupted_source = _issue58_corrupt_source_relation(entry)
    @test corrupted_source * entry.matrix != identity_matrix(R, n)
end

@testset "ToricBuilder factor_toric_block column peel acceptance" begin
    fixture = ToricBuilderFactorToricBlock3Fixture.fixture()
    pinv = _issue58_entry_by_role(fixture, "Pinv")
    qinv = _issue58_entry_by_role(fixture, "Qinv")
    R = base_ring(pinv.matrix)
    x, y = gens(R)

    pinv_expected_final = identity_matrix(R, 2)
    qinv_expected_final = matrix(R, [
        y                  zero(R);
        x * y^-1 + y^-1   y^-1
    ])

    _issue58_assert_column_peel_entry(pinv, collect(8:-1:3), pinv_expected_final)
    _issue58_assert_column_peel_entry(qinv, collect(16:-1:3), qinv_expected_final)
end
```

- [ ] **Step 2: Run the focused test and confirm RED**

Run:

```bash
julia --project=. -e 'include("test/public/toricbuilder_factor_toric_block_acceptance.jl")'
```

Expected: FAIL at the catalog-status assertion in the next task's companion checks until the shared catalog is promoted, or fail before this task because the file does not exist.

- [ ] **Step 3: Commit the failing test and workflow docs**

Run:

```bash
git add docs/superpowers/specs/2026-06-21-issue-58-toricbuilder-factor-block-column-peel-design.md docs/superpowers/plans/2026-06-21-issue-58-toricbuilder-factor-block-column-peel.md test/public/toricbuilder_factor_toric_block_acceptance.jl
git commit -m "test: add ToricBuilder column peel acceptance"
```

Expected: commit succeeds with a focused test in place.

---

### Task 2: Promote Catalog Metadata and Register the Acceptance Test

**Files:**
- Modify: `test/fixtures/toricbuilder_laurent_problem_catalog.jl`
- Modify: `test/internal/toricbuilder_problem_catalog.jl`
- Modify: `test/runtests.jl`

**Interfaces:**
- Consumes: the new public acceptance file path and the fixture metadata fields `expected_suslin_status` and `expected_suslin_path`.
- Produces: catalog entries that record Pinv/Qinv as supported by the Laurent column-peel route.

- [ ] **Step 1: Update the shared catalog entries**

In `test/fixtures/toricbuilder_laurent_problem_catalog.jl`, change `_contract_problem(entry)` so it records:

```julia
expected_current_status = :supported_column_peel,
expected_suslin_path = entry.expected_suslin_path,
verifier = _verifier("public/toricbuilder_factor_toric_block_acceptance.jl", :toricbuilder_factor_toric_block_column_peel_acceptance),
```

Also change the consumer issues to include `"#58"` and the tests tuple to include
`"test/public/toricbuilder_factor_toric_block_acceptance.jl"`.

- [ ] **Step 2: Update catalog validation**

In `test/internal/toricbuilder_problem_catalog.jl`:

Add `:supported_column_peel` to `ALLOWED_TORICBUILDER_PROBLEM_STATUSES`.

Add `:expected_suslin_path` to `REQUIRED_TORICBUILDER_PROBLEM_FIELDS`, and in
`validate_toricbuilder_problem_entry(entry)` require this field only for
`:supported_column_peel` entries:

```julia
if entry.expected_current_status == :supported_column_peel
    hasproperty(entry, :expected_suslin_path) ||
        throw(ArgumentError("problem $(entry.id) missing expected Suslin path"))
    entry.expected_suslin_path == :laurent_column_peel ||
        throw(ArgumentError("problem $(entry.id) expected Laurent column-peel path"))
end
```

Update the public assertions:

```julia
@test by_id["toricbuilder-factor-toric-block-3-qinv"].expected_current_status == :supported_column_peel
@test by_id["toricbuilder-factor-toric-block-3-qinv"].expected_suslin_path == :laurent_column_peel
@test by_id["toricbuilder-factor-toric-block-3-pinv"].expected_current_status == :supported_column_peel
@test by_id["toricbuilder-factor-toric-block-3-pinv"].expected_suslin_path == :laurent_column_peel
```

- [ ] **Step 3: Register the public test**

In `test/runtests.jl`, add the new test to the `"public"` group after
`"public/laurent_large_acceptance.jl"`:

```julia
"public/toricbuilder_factor_toric_block_acceptance.jl",
```

- [ ] **Step 4: Run focused and package verification**

Run:

```bash
julia --project=. -e 'include("test/public/toricbuilder_factor_toric_block_acceptance.jl")'
julia --project=. -e 'include("test/internal/toricbuilder_problem_catalog.jl")'
julia --project=. test/runtests.jl public
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: all commands pass.

- [ ] **Step 5: Commit catalog and registration changes**

Run:

```bash
git add test/fixtures/toricbuilder_laurent_problem_catalog.jl test/internal/toricbuilder_problem_catalog.jl test/runtests.jl
git commit -m "test: promote ToricBuilder factor block column peel"
```

Expected: commit succeeds.

## Self-Review

- Spec coverage: the plan covers public acceptance, exact factorization checks,
  replay final blocks, negative controls, catalog promotion, test registration,
  and package verification.
- Placeholder scan: no unresolved placeholders remain.
- Type consistency: planned names match the existing fixture and Issue #57
  column-peel APIs.
