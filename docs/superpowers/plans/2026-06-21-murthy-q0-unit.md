# Murthy q(0)-Unit Recursive Branch Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement the supported univariate Murthy q(0)-unit recursive branch for local `SL_3` special forms.

**Architecture:** Keep the public factor-only API unchanged. Extend `src/algorithm/sl3_local.jl` with an internal q-unit base case, a q(0)-unit reduction replay record, and a `:murthy_q0_unit` certificate branch that composes q-degree normalization, constant-term elimination, split-lemma replay, and recursive child certificates.

**Tech Stack:** Julia, Oscar/AbstractAlgebra polynomial rings, existing Suslin elementary matrix and certificate helpers.

## Global Constraints

- Preserve existing open-slice and unit diagonal pivot branches as fast paths.
- Supported Murthy q(0)-unit recursion is ordinary univariate exact polynomial input with monic `p`.
- Use existing `sl3_local_q_degree_normalization` and `sl3_local_split_lemma_replay` helpers.
- Add a clear recursion guard based on strictly decreasing degree in `p`.
- Do not implement the q(0)-nonunit Bezout/resultant branch.
- Preserve `realize_sl3_local(...)` returning a factor sequence.
- Verification command required by issue: `julia --project=. -e 'include("test/expert/sl3_local_murthy_q_unit.jl")'`.
- Repository verification command required by Agent Desk: `julia --project=. -e 'using Pkg; Pkg.test()'`.

---

### Task 1: Expert Coverage for q(0)-Unit Recursion

**Files:**
- Create: `test/expert/sl3_local_murthy_q_unit.jl`
- Modify: `test/runtests.jl`

**Interfaces:**
- Consumes: `Suslin.realize_sl3_local_certificate`, `Suslin.realize_sl3_local`, `Suslin.verify_factorization`, `Suslin.verify_sl3_local_realization`.
- Produces: failing coverage for `:murthy_q0_unit` certificates, nested q-degree normalization, split replay, recursive child certificates, and q(0)-nonunit staging.

- [ ] **Step 1: Add the failing expert test**

Create `test/expert/sl3_local_murthy_q_unit.jl` with helpers that multiply factors, inspect nested `:murthy_q0_unit` certificates, and assert exact replay:

```julia
using Test
using Suslin
using Oscar

const SL3_Q0_UNIT_FIXTURE_PATH =
    joinpath(@__DIR__, "..", "fixtures", "sl3_murthy_gupta_cases.jl")

function _q0_target(R, p, q, r, s)
    return matrix(R, [
        p q zero(R);
        r s zero(R);
        zero(R) zero(R) one(R)
    ])
end

function _q0_product(factors, R)
    product = identity_matrix(R, 3)
    for factor in factors
        product *= factor
    end
    return product
end

function _q0_degree(value, X)
    var_idx = findfirst(isequal(X), collect(gens(parent(value))))
    return degree(value, var_idx)
end

function _assert_elementary_sequence(target, factors)
    R = base_ring(target)
    @test !isempty(factors)
    for factor in factors
        @test nrows(factor) == 3
        @test ncols(factor) == 3
        @test Suslin._same_base_ring(base_ring(factor), R)
        nonzero_offdiagonal = 0
        for i in 1:3, j in 1:3
            if i == j
                @test factor[i, j] == one(R)
            elseif factor[i, j] != zero(R)
                nonzero_offdiagonal += 1
            end
        end
        @test nonzero_offdiagonal <= 1
    end
    @test _q0_product(factors, R) == target
    @test Suslin.verify_factorization(target, factors)
end

function _assert_reduction_replay(reduction)
    @test Suslin.verify_sl3_local_murthy_q_unit_reduction(reduction)
    split_cert = reduction.split_certificate
    split = split_cert.witness.split
    @test split.original_target == reduction.eliminated_target
    @test split.witness.a == reduction.selected_variable
    @test split.witness.a_prime == reduction.p_prime
    @test split.witness.b == reduction.eliminated_target[1, 2]
    @test _q0_degree(reduction.p_prime, reduction.selected_variable) < reduction.degree_p
    @test Suslin.verify_sl3_local_realization(split_cert)
    @test Suslin.verify_sl3_local_split_lemma_replay(split)
    @test Suslin.verify_sl3_local_realization(split_cert.witness.first_child_certificate)
    @test Suslin.verify_sl3_local_realization(split_cert.witness.second_child_certificate)
end

function _assert_q0_certificate(cert; normalization_expected::Bool)
    @test cert.branch == :murthy_q0_unit
    @test Suslin.verify_sl3_local_realization(cert)
    _assert_elementary_sequence(cert.target, cert.factors)
    if normalization_expected
        @test cert.witness.normalization !== nothing
        @test Suslin.verify_sl3_local_q_degree_normalization(cert.witness.normalization)
        @test cert.witness.normalized_certificate !== nothing
        @test cert.witness.reduction === nothing
        @test cert.witness.normalized_certificate.branch == :murthy_q0_unit
        _assert_q0_certificate(cert.witness.normalized_certificate; normalization_expected = false)
    else
        @test cert.witness.normalization === nothing
        @test cert.witness.normalized_certificate === nothing
        @test cert.witness.reduction !== nothing
        _assert_reduction_replay(cert.witness.reduction)
    end
end

@testset "Murthy q(0)-unit recursive branch for local SL3" begin
    include(SL3_Q0_UNIT_FIXTURE_PATH)
    catalog = SL3MurthyGuptaFixtureCatalog.catalog()
    by_id = Dict(entry.id => entry for entry in catalog.cases)

    fixture = by_id["mg-q0-unit-recursion"]
    fixture_cert = Suslin.realize_sl3_local_certificate(
        fixture.entries.p,
        fixture.entries.q,
        fixture.entries.r,
        fixture.entries.s,
        fixture.variable,
    )
    _assert_q0_certificate(fixture_cert; normalization_expected = false)
    @test fixture_cert.target == fixture.target
    @test Suslin.realize_sl3_local(
        fixture.entries.p,
        fixture.entries.q,
        fixture.entries.r,
        fixture.entries.s,
        fixture.variable,
    ) == fixture_cert.factors

    R, (X,) = Oscar.polynomial_ring(QQ, ["X"])

    p = X^3 + X + 1
    q = X^2 + X + 1
    _, s, minus_r = gcdx(p, q)
    r = -minus_r
    normalized_target = _q0_target(R, p, q, r, s)
    @test det(normalized_target) == one(R)
    @test _q0_degree(q, X) < _q0_degree(p, X)
    @test !is_unit(p)
    @test !is_unit(s)
    normalized_cert = Suslin.realize_sl3_local_certificate(normalized_target, X)
    _assert_q0_certificate(normalized_cert; normalization_expected = false)

    normalizing_fixture = by_id["mg-q-degree-normalization"]
    normalizing_cert = Suslin.realize_sl3_local_certificate(
        normalizing_fixture.entries.p,
        normalizing_fixture.entries.q,
        normalizing_fixture.entries.r,
        normalizing_fixture.entries.s,
        normalizing_fixture.variable,
    )
    _assert_q0_certificate(normalizing_cert; normalization_expected = true)
    @test normalizing_cert.target == normalizing_fixture.target

    nonunit_q0_p = X^2 + one(R)
    nonunit_q0_q = X
    nonunit_q0_s = X + one(R)
    nonunit_q0_r = div(nonunit_q0_p * nonunit_q0_s - one(R), X)
    nonunit_target = _q0_target(R, nonunit_q0_p, nonunit_q0_q, nonunit_q0_r, nonunit_q0_s)
    @test det(nonunit_target) == one(R)
    @test _q0_degree(nonunit_q0_q, X) < _q0_degree(nonunit_q0_p, X)
    @test_throws ArgumentError Suslin.realize_sl3_local_certificate(nonunit_target, X)
end
```

- [ ] **Step 2: Register the expert test**

Add `"expert/sl3_local_murthy_q_unit.jl"` after the q-degree normalization and split-lemma files in `test/runtests.jl`.

- [ ] **Step 3: Verify RED**

Run:

```bash
julia --project=. -e 'include("test/expert/sl3_local_murthy_q_unit.jl")'
```

Expected: FAIL because `verify_sl3_local_murthy_q_unit_reduction` and/or `:murthy_q0_unit` are not implemented yet.

- [ ] **Step 4: Commit Task 1**

```bash
git add test/expert/sl3_local_murthy_q_unit.jl test/runtests.jl
git commit -m "test: cover Murthy q0-unit recursion"
```

### Task 2: Solver and Certificate Replay

**Files:**
- Modify: `src/algorithm/sl3_local.jl`

**Interfaces:**
- Consumes: Task 1 tests, q-degree normalization helper, split-lemma replay helper.
- Produces: `SL3LocalMurthyQUnitReduction`, `verify_sl3_local_murthy_q_unit_reduction`, `:murthy_q0_unit` replay, q-unit base case, and recursive branch dispatch.

- [ ] **Step 1: Add reduction storage and public verifier helper**

Add this struct near the existing local certificate structs:

```julia
struct SL3LocalMurthyQUnitReduction
    target
    q0
    q0_inverse
    p0
    right_e21_coefficient
    eliminated_target
    elimination_factor
    inverse_elimination_factor
    p_prime
    split_certificate
    selected_variable
    degree_p::Int
    degree_p_prime::Int
    locality_witness
end
```

Add `verify_sl3_local_murthy_q_unit_reduction(reduction)::Bool` beside the existing q-degree and split replay verifiers.

- [ ] **Step 2: Add q-unit fast path**

In `_recognize_sl3_local_parameters`, after the existing `p_unit` check, recognize unit top-right entries:

```julia
q_inverse = _unit_inverse_or_nothing(q)
if q_inverse !== nothing
    return (; family = :q_unit, R, p, q, r, s, X, target, pivot_inverse = q_inverse)
end
```

Add `:q_unit` handling in `_sl3_local_form_factors`, `_sl3_local_form_witness`, `_sl3_local_certificate_expected_factors`, `_sl3_local_branch_witness_ok`, and `_sl3_local_witness_keys_ok`.

- [ ] **Step 3: Add Murthy candidate recognition**

Replace the final staged failure in `_recognize_sl3_local_parameters` with an ordinary-polynomial and monic-p check. Return `(; family = :murthy_q0_unit, R, p, q, r, s, X, target, var_idx)` when supported. Throw the staged q(0)-nonunit failure only after q-degree normalization has reduced `q` and the constant coefficient is not a unit.

- [ ] **Step 4: Implement q(0)-unit construction**

Add `_realize_sl3_local_murthy_q0_unit_certificate(form)` and helpers that:

```julia
degree_q = degree(form.q, form.var_idx)
degree_p = degree(form.p, form.var_idx)
if degree_q >= degree_p
    normalization = sl3_local_q_degree_normalization(form.p, form.q, form.r, form.s, form.X)
    normalized_certificate = realize_sl3_local_certificate(normalization.normalized_target, form.X)
    return SL3LocalRealizationCertificate(
        form.target,
        :murthy_q0_unit,
        vcat(normalized_certificate.factors, [normalization.elementary_correction]),
        form.X,
        (; normalization, normalized_certificate, reduction = nothing),
    )
end
```

For normalized input, compute `q0`, `p0`, `lambda`, `eliminated_target`, `p_prime`, split witnesses from `gcdx(form.X, form.q)` and `gcdx(p_prime, form.q)`, child certificates, split certificate, and final factors `vcat(split_certificate.factors, [inverse_elimination_factor])`.

- [ ] **Step 5: Add replay checks**

Extend `_sl3_local_certificate_expected_factors`, `_sl3_local_branch_witness_ok`, and `_sl3_local_witness_keys_ok` for `:murthy_q0_unit`. Add `_sl3_local_murthy_q0_unit_reduction_verification(reduction)` that checks q(0), p(0), elimination factors, `target * elimination_factor == eliminated_target`, `p + lambda*q == X*p_prime`, strict degree decrease, split certificate replay, and exact final factor product.

- [ ] **Step 6: Verify GREEN**

Run:

```bash
julia --project=. -e 'include("test/expert/sl3_local_murthy_q_unit.jl")'
```

Expected: PASS.

- [ ] **Step 7: Commit Task 2**

```bash
git add src/algorithm/sl3_local.jl
git commit -m "feat: implement Murthy q0-unit recursion"
```

### Task 3: Full Verification

**Files:**
- Modify only if verification exposes a focused issue in Task 1 or Task 2 files.

**Interfaces:**
- Consumes: completed Task 1 and Task 2.
- Produces: evidence for PR creation.

- [ ] **Step 1: Run focused expert verification**

```bash
julia --project=. -e 'include("test/expert/sl3_local_murthy_q_unit.jl")'
```

Expected: PASS.

- [ ] **Step 2: Run adjacent expert verification**

```bash
julia --project=. -e 'include("test/expert/sl3_local_q_degree_normalization.jl")'
julia --project=. -e 'include("test/expert/sl3_local_split_lemma.jl")'
julia --project=. -e 'include("test/expert/sl3_local_certificate.jl")'
```

Expected: PASS for all three commands.

- [ ] **Step 3: Run default package verification**

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: PASS.

- [ ] **Step 4: Run diff hygiene**

```bash
git diff --check
```

Expected: no output and exit code 0.
