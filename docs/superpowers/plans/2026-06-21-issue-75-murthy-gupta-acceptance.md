# Issue 75 Murthy-Gupta Acceptance Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add the final #61 acceptance test proving supported Murthy-Gupta local `SL_3` branches factor through `realize_sl3_local`.

**Architecture:** Keep the solver path centered in `src/algorithm/sl3_local.jl`. Add one expert acceptance file that exercises the public factor-returning path and certificate path over the existing fixture catalog, then make the narrow staged-error adjustment exposed by the negative control.

**Tech Stack:** Julia, Oscar/AbstractAlgebra polynomial rings, existing Suslin local `SL_3` certificates, Murthy-Gupta fixture catalog, and expert test runner.

## Global Constraints

- Preserve `realize_sl3_local(...)` returning a factor sequence.
- Preserve existing open-slice, unit-pivot, q-unit, q-degree normalization, split-lemma, q(0)-unit, and q(0)-nonunit branches.
- The focused acceptance test must include at least three hand-checkable special-form matrices over exact polynomial rings where `p` is monic and neither diagonal entry is a unit.
- For each supported case, `realize_sl3_local` or the expert certificate path must return a nonempty elementary factor sequence and `verify_factorization(target, factors) == true`.
- At least one focused case must be shown to fail under the pre-#61 open-slice/unit-pivot solver but succeed through a Murthy branch.
- Replayable certificates must name the Murthy branch and verify the witness relations used along the way.
- A determinant-one special-form matrix with non-monic `p` and no supplied local witness must throw a staged local `SL_3` `ArgumentError`.
- Do not implement Quillen induction, the Elementary Column Property pipeline, arbitrary local rings, factor-count optimization, or the final public Park-Woodburn driver.
- Verification commands required by issue: `julia --project=. -e 'include("test/expert/sl3_local_murthy_gupta.jl")'` and `julia --project=. test/runtests.jl all`.
- Verification command required by Agent Desk: `julia --project=. -e 'using Pkg; Pkg.test()'`.

---

### Task 1: Final Murthy-Gupta Acceptance Test

**Files:**
- Create: `test/expert/sl3_local_murthy_gupta.jl`
- Modify: `test/runtests.jl`

**Interfaces:**
- Consumes: `SL3MurthyGuptaFixtureCatalog.catalog()`, `realize_sl3_local`, `realize_sl3_local_certificate`, `verify_factorization`, and the Murthy certificate verifier helpers.
- Produces: focused #61 acceptance coverage that initially fails because the negative control requires the staged local `SL_3` error shape for non-monic `p`.

- [ ] **Step 1: Create the focused acceptance test**

Create `test/expert/sl3_local_murthy_gupta.jl` with this content:

```julia
using Test
using Suslin
using Oscar

const SL3_MURTHY_GUPTA_ACCEPTANCE_FIXTURE_PATH =
    joinpath(@__DIR__, "..", "fixtures", "sl3_murthy_gupta_cases.jl")

function _mg_acceptance_product(factors, R)
    product = identity_matrix(R, 3)
    for factor in factors
        product *= factor
    end
    return product
end

function _mg_acceptance_degree(value, X)
    var_idx = findfirst(isequal(X), collect(gens(parent(value))))
    return degree(value, var_idx)
end

function _mg_acceptance_special_form_target(R, p, q, r, s)
    return matrix(R, [
        p q zero(R);
        r s zero(R);
        zero(R) zero(R) one(R)
    ])
end

function _mg_acceptance_assert_elementary_sequence(target, factors)
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
    @test _mg_acceptance_product(factors, R) == target
    @test Suslin.verify_factorization(target, factors)
end

function _mg_acceptance_pre_murthy_open_or_unit_pivot(entry)
    R = entry.ring.object
    p, q, r, s = entry.entries.p, entry.entries.q, entry.entries.r, entry.entries.s
    if s == one(R) && p == one(R) + q * r
        return :open_s_one
    elseif p == one(R) && s == one(R) + q * r
        return :open_p_one
    elseif is_unit(s)
        return :s_unit
    elseif is_unit(p)
        return :p_unit
    end
    throw(ArgumentError("pre-#61 open-slice/unit-pivot solver does not support this target"))
end

function _mg_acceptance_assert_supported_case(entry; kwargs...)
    p, q, r, s = entry.entries.p, entry.entries.q, entry.entries.r, entry.entries.s
    X = entry.variable
    target = entry.target
    R = entry.ring.object

    @test target == _mg_acceptance_special_form_target(R, p, q, r, s)
    @test det(target) == one(R)
    @test Suslin._is_monic_in_variable(p, findfirst(isequal(X), collect(gens(R))), R)
    @test !is_unit(p)
    @test !is_unit(s)

    certificate_from_matrix = Suslin.realize_sl3_local_certificate(target, X; kwargs...)
    certificate_from_entries = Suslin.realize_sl3_local_certificate(p, q, r, s, X; kwargs...)
    @test certificate_from_matrix.target == target
    @test certificate_from_entries.target == target
    @test certificate_from_matrix.branch == certificate_from_entries.branch
    @test Suslin.verify_sl3_local_realization(certificate_from_matrix)
    @test Suslin.verify_sl3_local_realization(certificate_from_entries)
    _mg_acceptance_assert_elementary_sequence(target, certificate_from_matrix.factors)

    factors_from_matrix = Suslin.realize_sl3_local(target, X; kwargs...)
    factors_from_entries = Suslin.realize_sl3_local(p, q, r, s, X; kwargs...)
    @test factors_from_matrix == certificate_from_matrix.factors
    @test factors_from_entries == certificate_from_entries.factors
    _mg_acceptance_assert_elementary_sequence(target, factors_from_matrix)
    return certificate_from_matrix
end

function _mg_acceptance_assert_q0_unit_certificate(certificate)
    @test certificate.branch == :murthy_q0_unit
    @test Suslin.verify_sl3_local_realization(certificate)

    if certificate.witness.normalization !== nothing
        normalization = certificate.witness.normalization
        @test Suslin.verify_sl3_local_q_degree_normalization(normalization)
        @test normalization.target == certificate.target
        @test normalization.selected_variable == certificate.selected_variable
        @test certificate.witness.normalized_certificate !== nothing
        @test certificate.witness.reduction === nothing
        _mg_acceptance_assert_q0_unit_certificate(certificate.witness.normalized_certificate)
        return nothing
    end

    reduction = certificate.witness.reduction
    @test reduction !== nothing
    @test Suslin.verify_sl3_local_murthy_q_unit_reduction(reduction)
    @test reduction.target == certificate.target
    @test reduction.selected_variable == certificate.selected_variable
    @test reduction.q0 * reduction.q0_inverse == one(base_ring(certificate.target))
    @test _mg_acceptance_degree(reduction.p_prime, reduction.selected_variable) < reduction.degree_p
    @test reduction.split_certificate.branch == :murthy_split_lemma
    @test Suslin.verify_sl3_local_realization(reduction.split_certificate)
    @test Suslin.verify_sl3_local_split_lemma_replay(reduction.split_certificate.witness.split)
    @test reduction.split_certificate.witness.split.split_id == :murthy_q0_unit_split
    return nothing
end

function _mg_acceptance_assert_resultant_certificate(certificate; expected_source::Symbol)
    @test certificate.branch == :murthy_q0_nonunit_bezout_resultant
    @test Suslin.verify_sl3_local_realization(certificate)
    reduction = certificate.witness.reduction
    @test reduction.witness_source == expected_source
    @test Suslin.verify_sl3_local_murthy_q0_nonunit_reduction(reduction)
    @test reduction.resultant == one(base_ring(certificate.target))
    @test reduction.p_prime * certificate.target[1, 1] - reduction.q_prime * certificate.target[1, 2] ==
        one(base_ring(certificate.target))
    @test reduction.target == reduction.left_factor * reduction.bezout_target
    @test reduction.bezout_target == reduction.first_elementary_factor * reduction.child_link_target
    @test reduction.branch_unit * reduction.branch_unit_inverse == one(base_ring(certificate.target))
    @test _mg_acceptance_degree(reduction.p_prime, reduction.selected_variable) < reduction.degree_q
    @test _mg_acceptance_degree(reduction.q_prime, reduction.selected_variable) < reduction.degree_p
    _mg_acceptance_assert_q0_unit_certificate(reduction.child_certificate)
    return nothing
end

@testset "Issue 61 Murthy-Gupta local SL3 acceptance" begin
    include(SL3_MURTHY_GUPTA_ACCEPTANCE_FIXTURE_PATH)
    catalog = SL3MurthyGuptaFixtureCatalog.catalog()
    by_id = Dict(entry.id => entry for entry in catalog.cases)

    normalization_entry = by_id["mg-q-degree-normalization"]
    q0_unit_entry = by_id["mg-q0-unit-recursion"]
    supplied_entry = by_id["mg-q0-nonunit-normalized-bezout-resultant"]
    extracted_entry = by_id["mg-q0-nonunit-extracted-bezout-resultant"]

    acceptance_entries = (normalization_entry, q0_unit_entry, supplied_entry, extracted_entry)
    @test length(acceptance_entries) >= 3
    @test_throws ArgumentError _mg_acceptance_pre_murthy_open_or_unit_pivot(normalization_entry)
    @test_throws ArgumentError _mg_acceptance_pre_murthy_open_or_unit_pivot(supplied_entry)

    normalization_certificate = _mg_acceptance_assert_supported_case(normalization_entry)
    @test normalization_certificate.branch == :murthy_q0_unit
    @test normalization_certificate.witness.normalization !== nothing
    _mg_acceptance_assert_q0_unit_certificate(normalization_certificate)

    q0_unit_certificate = _mg_acceptance_assert_supported_case(q0_unit_entry)
    @test q0_unit_certificate.witness.normalization === nothing
    _mg_acceptance_assert_q0_unit_certificate(q0_unit_certificate)

    supplied_certificate = _mg_acceptance_assert_supported_case(
        supplied_entry;
        murthy_q0_nonunit_witness = first(supplied_entry.witnesses),
    )
    _mg_acceptance_assert_resultant_certificate(
        supplied_certificate;
        expected_source = :supplied_bezout_witness,
    )

    extracted_certificate = _mg_acceptance_assert_supported_case(extracted_entry)
    _mg_acceptance_assert_resultant_certificate(
        extracted_certificate;
        expected_source = :extracted_bezout_witness,
    )
end

@testset "Issue 61 staged local SL3 unsupported boundary" begin
    R, (X,) = Oscar.polynomial_ring(QQ, ["X"])
    p = 2 * X + one(R)
    q = X
    r = R(2)
    s = one(R)
    target = _mg_acceptance_special_form_target(R, p, q, r, s)
    @test det(target) == one(R)
    @test !Suslin._is_monic_in_variable(p, 1, R)

    err = try
        Suslin.realize_sl3_local(target, X)
        nothing
    catch caught
        caught
    end
    @test err isa ArgumentError
    @test occursin("staged local SL_3 solver failure", sprint(showerror, err))
    @test occursin("p must be monic in X", sprint(showerror, err))
end
```

- [ ] **Step 2: Register the focused test in the expert group**

In `test/runtests.jl`, add the new expert file immediately after the existing Murthy resultants test:

```julia
        "expert/sl3_local_murthy_q_unit.jl",
        "expert/sl3_local_murthy_resultant.jl",
        "expert/sl3_local_murthy_gupta.jl",
        "expert/sln_to_sl3_reduction.jl",
```

- [ ] **Step 3: Verify RED**

Run:

```bash
julia --project=. -e 'include("test/expert/sl3_local_murthy_gupta.jl")'
```

Expected: FAIL in the negative-control test because `realize_sl3_local(target, X)` currently throws `ArgumentError("p must be monic in X")` without the staged local `SL_3` prefix.

- [ ] **Step 4: Commit Task 1**

```bash
git add test/expert/sl3_local_murthy_gupta.jl test/runtests.jl
git commit -m "test: add Murthy-Gupta acceptance coverage"
```

### Task 2: Staged Monicity Boundary

**Files:**
- Modify: `src/algorithm/sl3_local.jl`

**Interfaces:**
- Consumes: Task 1 negative control.
- Produces: a staged local `SL_3` `ArgumentError` for non-monic `p` encountered through `realize_sl3_local`.

- [ ] **Step 1: Change the local-solver recognition monicity failure**

In `_recognize_sl3_local_parameters`, replace:

```julia
        _is_monic_in_variable(p, var_idx, R) || throw(ArgumentError("p must be monic in X"))
```

with:

```julia
        _is_monic_in_variable(p, var_idx, R) || _throw_staged_sl3_local_failure("p must be monic in X")
```

Do not change the q-degree normalization helper or low-level division helper error messages.

- [ ] **Step 2: Verify GREEN**

Run:

```bash
julia --project=. -e 'include("test/expert/sl3_local_murthy_gupta.jl")'
```

Expected: PASS.

- [ ] **Step 3: Run surrounding local solver tests**

Run:

```bash
julia --project=. -e 'include("test/expert/sl3_local_certificate.jl"); include("test/expert/sl3_local_murthy_q_unit.jl"); include("test/expert/sl3_local_murthy_resultant.jl")'
```

Expected: PASS.

- [ ] **Step 4: Commit Task 2**

```bash
git add src/algorithm/sl3_local.jl
git commit -m "fix: stage nonmonic local SL3 boundary"
```

### Task 3: Acceptance Verification

**Files:**
- No code files. This task verifies the branch after Tasks 1 and 2.

**Interfaces:**
- Consumes: Tasks 1 and 2.
- Produces: final evidence for issue #75 and Agent Desk.

- [ ] **Step 1: Run issue verification**

Run:

```bash
julia --project=. -e 'include("test/expert/sl3_local_murthy_gupta.jl")'
julia --project=. test/runtests.jl all
```

Expected: both commands exit 0.

- [ ] **Step 2: Run Agent Desk verification**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: command exits 0.

- [ ] **Step 3: Check formatting and git state**

Run:

```bash
git diff --check
git status -sb
```

Expected: `git diff --check` exits 0 and `git status -sb` shows only intentional committed branch differences from `origin/main`.
