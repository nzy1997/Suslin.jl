# Issue 17 Laurent Large Acceptance Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an exact 40x40+ Laurent acceptance harness that keeps the ToricBuilder contract fixture in the ladder and verifies supported large Laurent factorizations end to end.

**Architecture:** Put reusable large acceptance cases in `test/fixtures/laurent_large_acceptance_cases.jl`, keep the acceptance gate in `test/public/laurent_large_acceptance.jl`, and register it in the public test group. Enable only the already-designed narrow Laurent block-local `SL_n` reduction path so determinant-one Laurent block-local matrices factor through the public driver, while non-elementary determinant corrections remain staged.

**Tech Stack:** Julia, Oscar, Suslin's existing Laurent normalization, block embedding, `SL_n` reduction, local `SL_3` solver, and `Test`.

## Global Constraints

- The ToricBuilder-motivated fixture from issue #19 must remain in the acceptance catalog and must not be replaced by synthetic large matrices.
- The acceptance catalog must include at least one 40x40 Laurent matrix and one larger enabled case; this plan uses 48x48.
- Large matrices must live in a reusable fixture module, not as opaque inline test data.
- Factorization cases must run `elementary_factorization(A)` and verify exact reconstruction with `verify_factorization(A, factors)`.
- The ToricBuilder case may follow the normalized/contract path: `normalize_laurent_gl_matrix(A)`, `verify_laurent_gl_normalization(A, normalization)`, and exact inverse-relation verification.
- Negative controls must fail by exact verification.
- Do not claim elementary factors for Laurent `GL_n` inputs whose determinant correction is not represented in the returned factor sequence.
- The documented full-suite command is `julia --project=. test/runtests.jl all`.

---

## File Structure

- Create `test/fixtures/laurent_large_acceptance_cases.jl`: reusable issue #17 acceptance catalog and deterministic large Laurent matrix constructors.
- Create `test/public/laurent_large_acceptance.jl`: public acceptance gate for ToricBuilder normalized contract and large exact factorization cases.
- Modify `test/runtests.jl`: register the acceptance gate in the `public` group.
- Modify `src/algorithm/sln_to_sl3_reduction.jl`: allow normalized Laurent block-local reductions with `check_monic=false`.
- Modify `src/algorithm/factorization.jl`: route Laurent candidates through the reduction layer and return factors only after exact verification against the original input.
- Modify `test/internal/gl_laurent_normalization.jl`: update old staged-boundary assertions that become supported by the new Laurent block-local path.

---

### Task 1: Add the Failing Acceptance Harness

**Files:**
- Create: `test/fixtures/laurent_large_acceptance_cases.jl`
- Create: `test/public/laurent_large_acceptance.jl`
- Modify: `test/runtests.jl`

**Interfaces:**
- Consumes: `ToricBuilderFactorToricBlock3Fixture.fixture()`, `block_embedding`, `elementary_factorization`, `reduce_sln_to_sl3`, `verify_factorization`, `normalize_laurent_gl_matrix`, `verify_laurent_gl_normalization`.
- Produces: `LaurentLargeAcceptanceCases.acceptance_catalog()` returning `(; cases = [...])`, with cases carrying `id`, `kind`, `size`, `matrix`, `provenance`, `expected_path`, and negative-control metadata.

- [ ] **Step 1: Write the failing fixture module**

Create `test/fixtures/laurent_large_acceptance_cases.jl` with:

```julia
module LaurentLargeAcceptanceCases

using Oscar
using Suslin

include("toricbuilder_factor_toric_block_3.jl")

const ACCEPTANCE_40_SIZE = 40
const ACCEPTANCE_LARGER_SIZE = 48

function _block_locations(n::Int)
    return [Int[first_idx, first_idx + 1, first_idx + 2] for first_idx in 1:3:(n - 2)]
end

function _local_laurent_block(R, x, y, block_index::Int)
    q = isodd(block_index) ? x * y : x + y^-1
    r = one(R)
    p = one(R) + q * r
    s = one(R)
    return matrix(R, [
        p q zero(R);
        r s zero(R);
        zero(R) zero(R) one(R)
    ])
end

function _large_block_local_case(R, x, y, n::Int, id::String, description::String)
    locations = _block_locations(n)
    A = identity_matrix(R, n)
    for (block_index, indices) in enumerate(locations)
        A *= block_embedding(_local_laurent_block(R, x, y, block_index), n, indices)
    end

    return (;
        id,
        kind = :large_laurent_factorization,
        ring = (;
            description = "GF(2)[x^+/-1, y^+/-1]",
            object = R,
            generators = (x, y),
        ),
        size = (n, n),
        matrix = A,
        block_locations = locations,
        expected_path = :elementary_factorization,
        provenance = (;
            source = :synthetic_supported_block_local,
            issue = "#17",
            description,
            construction = "product of disjoint embedded Laurent local SL3 blocks",
        ),
        negative_control = (;
            kind = :replace_first_factor_with_identity,
            description = "Replacing the first returned factor by identity must break exact reconstruction.",
        ),
    )
end

function _toricbuilder_pinv_case()
    fixture = ToricBuilderFactorToricBlock3Fixture.fixture()
    pinv = only(filter(entry -> entry.toricbuilder_role == "Pinv", fixture.cases))
    return (;
        id = "toricbuilder-factor-toric-block-3-pinv-normalized-contract",
        kind = :toricbuilder_normalized_contract,
        ring = pinv.ring,
        size = pinv.size,
        matrix = pinv.matrix,
        source_matrix = pinv.source_matrix,
        expected_path = :normalized_contract,
        provenance = (;
            source = :toricbuilder_contract_fixture,
            issue = "#19",
            fixture_id = pinv.name,
            toricbuilder_role = pinv.toricbuilder_role,
            toricbuilder_commit = pinv.provenance.toricbuilder_commit,
            generation_command = pinv.provenance.generation_command,
        ),
        negative_control = (;
            kind = :corrupt_inverse_relation_entry,
            row = 1,
            col = 1,
            description = "Toggling one Pinv entry must break the exact ToricBuilder inverse relation.",
        ),
    )
end

function acceptance_catalog()
    R, (x, y) = Suslin.suslin_laurent_polynomial_ring(GF(2), ["x", "y"])
    return (;
        cases = [
            _toricbuilder_pinv_case(),
            _large_block_local_case(
                R,
                x,
                y,
                ACCEPTANCE_40_SIZE,
                "laurent-block-local-40x40",
                "40x40 determinant-one Laurent matrix with disjoint supported local SL3 blocks",
            ),
            _large_block_local_case(
                R,
                x,
                y,
                ACCEPTANCE_LARGER_SIZE,
                "laurent-block-local-48x48",
                "48x48 determinant-one Laurent matrix with disjoint supported local SL3 blocks",
            ),
        ],
    )
end

end
```

- [ ] **Step 2: Write the failing public acceptance test**

Create `test/public/laurent_large_acceptance.jl` with:

```julia
using Suslin
using Test
using Oscar

include("../fixtures/laurent_large_acceptance_cases.jl")
using .LaurentLargeAcceptanceCases

function _acceptance_product(factors, R, n::Int)
    product = identity_matrix(R, n)
    for factor in factors
        product *= factor
    end
    return product
end

function _assert_toricbuilder_contract_case(case)
    R = base_ring(case.matrix)
    normalization = normalize_laurent_gl_matrix(case.matrix)
    @test normalization.determinant_classification == :one
    @test normalization.normalized_matrix == case.matrix
    @test verify_laurent_gl_normalization(case.matrix, normalization)
    @test case.source_matrix * case.matrix == identity_matrix(R, case.size[1])

    corrupted = copy(case.matrix)
    corrupted[case.negative_control.row, case.negative_control.col] += one(R)
    @test case.source_matrix * corrupted != identity_matrix(R, case.size[1])
end

function _assert_large_factorization_case(case)
    R = base_ring(case.matrix)
    n = case.size[1]
    normalization = normalize_laurent_gl_matrix(case.matrix)
    @test normalization.determinant_classification == :one
    @test normalization.normalized_matrix == case.matrix
    @test verify_laurent_gl_normalization(case.matrix, normalization)

    reduction = reduce_sln_to_sl3(case.matrix; block_locations = case.block_locations)
    @test verify_sln_to_sl3_reduction(reduction)
    @test length(reduction.obligations) == length(case.block_locations)
    @test reduction.product == case.matrix

    factors = elementary_factorization(case.matrix)
    @test !isempty(factors)
    @test _acceptance_product(factors, R, n) == case.matrix
    @test verify_factorization(case.matrix, factors)

    corrupted_factors = copy(factors)
    corrupted_factors[1] = identity_matrix(R, n)
    @test !verify_factorization(case.matrix, corrupted_factors)
end

@testset "large Laurent acceptance catalog" begin
    catalog = LaurentLargeAcceptanceCases.acceptance_catalog()
    cases = catalog.cases

    @test any(case -> case.kind == :toricbuilder_normalized_contract, cases)
    @test any(case -> case.size == (40, 40), cases)
    @test any(case -> case.size == (48, 48), cases)

    for case in cases
        if case.expected_path == :normalized_contract
            _assert_toricbuilder_contract_case(case)
        elseif case.expected_path == :elementary_factorization
            _assert_large_factorization_case(case)
        else
            error("unknown acceptance path $(case.expected_path) for $(case.id)")
        end
    end
end
```

- [ ] **Step 3: Register the public test**

In `test/runtests.jl`, add the file to the `public` list immediately after `public/factorization_driver_shell.jl`:

```julia
    "public" => [
        "public/api_surface.jl",
        "public/factorization_driver_shell.jl",
        "public/laurent_large_acceptance.jl",
    ],
```

- [ ] **Step 4: Run the acceptance test to verify RED**

Run:

```bash
julia --project=. -e 'include("test/public/laurent_large_acceptance.jl")'
```

Expected: FAIL because `reduce_sln_to_sl3` or `elementary_factorization`
still reports a staged `SL_n` to local `SL_3` reduction failure for the
40x40 Laurent case.

- [ ] **Step 5: Commit the failing harness**

Commit only the fixture, public test, runner registration, and plan:

```bash
git add test/fixtures/laurent_large_acceptance_cases.jl test/public/laurent_large_acceptance.jl test/runtests.jl docs/superpowers/plans/2026-06-19-issue-17-laurent-large-acceptance.md
git commit -m "test: add laurent large acceptance harness"
```

---

### Task 2: Enable Narrow Laurent Block-Local Reduction

**Files:**
- Modify: `src/algorithm/sln_to_sl3_reduction.jl`
- Modify: `src/algorithm/factorization.jl`
- Modify: `test/internal/gl_laurent_normalization.jl`

**Interfaces:**
- Consumes: failing acceptance test from Task 1.
- Produces: Laurent block-local determinant-one matrices factor through `reduce_sln_to_sl3` and `elementary_factorization`; non-identity Laurent determinant corrections still stage out of the public driver.

- [ ] **Step 1: Update Laurent reduction support**

In `src/algorithm/sln_to_sl3_reduction.jl`, replace `_reduction_generator` and the local-obligation call/constructor with:

```julia
function _reduction_generator(R, ring_profile::Symbol)
    ring_gens = collect(gens(R))
    isempty(ring_gens) && throw(ArgumentError("reduction requires a polynomial or Laurent generator"))

    if ring_profile == :polynomial
        length(ring_gens) == 1 || throw(ArgumentError("ordinary polynomial reduction currently requires a univariate base ring"))
    end

    return ring_gens[1]
end
```

Change the obligation build loop to pass the ring profile:

```julia
        push!(obligations, _build_sl3_local_obligation(normalized_A, normalized_R, indices, X, ring_profile))
```

Replace `_build_sl3_local_obligation` with:

```julia
function _local_obligation_assumptions(ring_profile::Symbol)
    if ring_profile == :laurent
        return Symbol[:disjoint_block_support, :determinant_one_core, :laurent_normalized_check_monic_false]
    end

    return Symbol[:univariate_base_ring, :determinant_one]
end

function _build_sl3_local_obligation(A, R, indices::Vector{Int}, X, ring_profile::Symbol)
    local_target = _principal_submatrix(A, indices)
    local_factors = try
        realize_sl3_local(local_target, X; check_monic = ring_profile != :laurent)
    catch err
        err isa InterruptException && rethrow()
        _throw_staged_sln_to_sl3_failure("failed to solve local SL_3 obligation on block $(indices)")
    end
    embedded_target = block_embedding(local_target, nrows(A), indices)
    embedded_factors = embed_factor_sequence(local_factors, nrows(A), indices)

    reassembly_data = _sl3_obligation_reassembly_data(
        local_factors,
        embedded_factors,
        R,
        nrows(A),
        local_target,
        embedded_target,
    )

    return SL3LocalObligation(
        copy(indices),
        R,
        local_target,
        _local_obligation_assumptions(ring_profile),
        embedded_target,
        local_factors,
        embedded_factors,
        reassembly_data,
    )
end
```

Add a verifier helper and replace the hard-coded assumption check inside `_verify_sl3_local_obligation`:

```julia
function _expected_obligation_assumptions(R)
    return _is_laurent_polynomial_ring(R) ?
        Symbol[:disjoint_block_support, :determinant_one_core, :laurent_normalized_check_monic_false] :
        Symbol[:univariate_base_ring, :determinant_one]
end
```

Use:

```julia
        obligation.required_assumptions == _expected_obligation_assumptions(R) || return false
```

- [ ] **Step 2: Route Laurent candidates through reduction in the public driver**

In `src/algorithm/factorization.jl`, replace:

```julia
    if ring_profile == :polynomial && nrows(normalized_A) > 3
        reduction = reduce_sln_to_sl3(A)
        verify_factorization(A, reduction.factors) && return reduction.factors
    end
```

with:

```julia
    if nrows(normalized_A) > 3 || ring_profile == :laurent
        reduction = reduce_sln_to_sl3(A)
        verify_factorization(A, reduction.factors) && return reduction.factors
    end
```

Keep the existing staged failure after that block so Laurent inputs with determinant corrections still report the Laurent normalization boundary instead of returning factors for only the normalized core.

- [ ] **Step 3: Update old Laurent boundary tests that are now supported**

In `test/internal/gl_laurent_normalization.jl`, replace the `normalized_laurent_sl3` staged-error assertion with exact factorization:

```julia
    normalized_laurent_sl3 = matrix(R, [
        one(R) x zero(R);
        zero(R) one(R) zero(R);
        zero(R) zero(R) one(R)
    ])
    sl3_factors = elementary_factorization(normalized_laurent_sl3)
    @test verify_factorization(normalized_laurent_sl3, sl3_factors)
```

Leave `normalized_then_rejected` as a staged failure because its determinant is `x` and the returned reduction factors would reconstruct only the determinant-one core, not the original `GL_n` input.

- [ ] **Step 4: Run the focused acceptance test to verify GREEN**

Run:

```bash
julia --project=. -e 'include("test/public/laurent_large_acceptance.jl")'
```

Expected: PASS.

- [ ] **Step 5: Run updated internal normalization tests**

Run:

```bash
julia --project=. -e 'include("test/internal/gl_laurent_normalization.jl")'
```

Expected: PASS.

- [ ] **Step 6: Commit the implementation**

```bash
git add src/algorithm/sln_to_sl3_reduction.jl src/algorithm/factorization.jl test/internal/gl_laurent_normalization.jl
git commit -m "feat: enable block-local laurent acceptance factorization"
```

---

### Task 3: Verify Full Acceptance Gate

**Files:**
- No new files; verification-only task.

**Interfaces:**
- Consumes: committed harness and implementation.
- Produces: fresh verification evidence for focused, package, and full-suite commands.

- [ ] **Step 1: Run issue-specific acceptance**

```bash
julia --project=. -e 'include("test/public/laurent_large_acceptance.jl")'
```

Expected: PASS.

- [ ] **Step 2: Run Agent Desk package command**

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: PASS, including the new public acceptance file.

- [ ] **Step 3: Run documented full suite**

```bash
julia --project=. test/runtests.jl all
```

Expected: PASS.

- [ ] **Step 4: Confirm no verification-only commit is needed**

Run:

```bash
git status -sb
```

Expected: no unstaged changes beyond the committed issue #17 work. If this
shows unexpected files, inspect them before finishing; do not create a
verification-only commit unless a concrete issue #17 file needs a correction.

---

## Plan Self-Review

- Spec coverage: the plan covers the real ToricBuilder fixture, reusable large fixture module, 40x40 case, 48x48 case, exact factorization verification, normalized contract verification, negative control, focused command, package command, and full-suite command.
- Placeholder scan: no incomplete implementation markers remain.
- Type consistency: the catalog function is consistently named `acceptance_catalog`, and factorization/reduction verifier names match the existing public API.
