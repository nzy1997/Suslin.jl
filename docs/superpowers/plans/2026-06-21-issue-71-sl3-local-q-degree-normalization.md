# Issue 71 SL3 Local Q-Degree Normalization Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add the Murthy-Gupta q-degree normalization step for local `SL_3` special-form polynomial inputs and make it replayable through the existing certificate verifier.

**Architecture:** Keep the feature inside `src/algorithm/sl3_local.jl` beside the existing local certificate path. Add an internal helper that computes `q = f*p + g` with `divrem`, records the normalized target and elementary correction, and extend `verify_sl3_local_realization` with a partial-normalization branch that checks the exact reconstruction identity instead of full factorization.

**Tech Stack:** Julia, Oscar polynomial rings, AbstractAlgebra/Oscar polynomial division, Suslin elementary matrices, Test stdlib.

## Global Constraints

- Do not implement q(0)-unit recursion, split-lemma recursion, or the resultant/Bezout branch.
- Do not change the public `elementary_factorization` driver.
- Keep the helper polynomial-only unless later work explicitly broadens the Murthy path.
- Reuse at least one fixture from issue 69 where `deg(q) >= deg(p)`.
- The q-degree branch is a partial local realization certificate; it proves normalization replay, not a complete factorization of the original target.
- Focused verification command is `julia --project=. -e 'include("test/expert/sl3_local_q_degree_normalization.jl")'`.
- Required package verification command is `julia --project=. -e 'using Pkg; Pkg.test()'`.

---

## File Structure

- Modify `src/algorithm/sl3_local.jl`: factor out general special-form recognition, add q-degree normalization helper methods, witness construction, polynomial-only guard, degree helper, and branch-specific verifier support.
- Create `test/expert/sl3_local_q_degree_normalization.jl`: focused tests for two exact polynomial examples, fixture reuse, certificate replay, and negative controls.
- Modify `test/runtests.jl`: register the expert test after `expert/sl3_local_certificate.jl`.

---

### Task 1: Focused Q-Degree Normalization Test

**Files:**
- Create: `test/expert/sl3_local_q_degree_normalization.jl`

**Interfaces:**
- Consumes: `SL3MurthyGuptaFixtureCatalog.catalog()` from `test/fixtures/sl3_murthy_gupta_cases.jl`.
- Consumes: existing `Suslin.elementary_matrix`, `Suslin.SL3LocalRealizationCertificate`, and `Suslin.verify_sl3_local_realization`.
- Produces failing expectations for `Suslin.normalize_sl3_local_q_degree(p, q, r, s, X; check_monic=true)` and `Suslin.normalize_sl3_local_q_degree(A, X; check_monic=true)`.

- [ ] **Step 1: Create the failing expert test**

Create `test/expert/sl3_local_q_degree_normalization.jl` with:

```julia
using Test
using Suslin
using Oscar

const SL3_Q_DEGREE_FIXTURE_PATH =
    joinpath(@__DIR__, "..", "fixtures", "sl3_murthy_gupta_cases.jl")

if !isdefined(@__MODULE__, :SL3MurthyGuptaFixtureCatalog)
    include(SL3_Q_DEGREE_FIXTURE_PATH)
end

function _q_degree_target(R, p, q, r, s)
    return matrix(R, [
        p q zero(R);
        r s zero(R);
        zero(R) zero(R) one(R)
    ])
end

function _q_degree_case_data_from_fixture()
    catalog = SL3MurthyGuptaFixtureCatalog.catalog()
    entry = only(filter(case -> case.id == "mg-q-degree-normalization", catalog.cases))
    witness = first(entry.witnesses)
    return (;
        id = entry.id,
        R = entry.ring.object,
        X = entry.variable,
        p = entry.entries.p,
        q = entry.entries.q,
        r = entry.entries.r,
        s = entry.entries.s,
        quotient = witness.quotient,
        remainder = witness.remainder,
        normalized_s = witness.normalized_s,
        target = entry.target,
    )
end

function _q_degree_hand_checked_case()
    R, (X,) = Oscar.polynomial_ring(QQ, ["X"])
    p = X^2 + 1
    f = X^2
    g = X
    r = X
    normalized_s = one(R)
    s = normalized_s + f * r
    q = f * p + g
    target = _q_degree_target(R, p, q, r, s)
    @assert det(target) == one(R)
    return (;
        id = "hand-checked-q-degree",
        R,
        X,
        p,
        q,
        r,
        s,
        quotient = f,
        remainder = g,
        normalized_s,
        target,
    )
end

function _assert_q_degree_certificate(case)
    cert = Suslin.normalize_sl3_local_q_degree(case.p, case.q, case.r, case.s, case.X)
    correction = Suslin.elementary_matrix(3, 1, 2, case.quotient, case.R)
    normalized_target = _q_degree_target(
        case.R,
        case.p,
        case.remainder,
        case.r,
        case.normalized_s,
    )

    @test cert isa Suslin.SL3LocalRealizationCertificate
    @test cert.branch == :q_degree_normalization
    @test cert.target == case.target
    @test cert.selected_variable == case.X
    @test cert.factors == [correction]
    @test cert.witness.p == case.p
    @test cert.witness.q == case.q
    @test cert.witness.r == case.r
    @test cert.witness.s == case.s
    @test cert.witness.quotient == case.quotient
    @test cert.witness.remainder == case.remainder
    @test cert.witness.normalized_target == normalized_target
    @test cert.witness.correction == correction
    @test cert.witness.p_degree == degree(case.p, 1)
    @test cert.witness.remainder_degree == degree(case.remainder, 1)
    @test case.q == case.quotient * case.p + case.remainder
    @test degree(case.remainder, 1) < degree(case.p, 1)
    @test normalized_target * correction == case.target
    @test Suslin.verify_sl3_local_realization(cert)

    matrix_cert = Suslin.normalize_sl3_local_q_degree(case.target, case.X)
    @test matrix_cert.target == cert.target
    @test matrix_cert.witness == cert.witness
    @test Suslin.verify_sl3_local_realization(matrix_cert)
end

@testset "local SL3 q-degree normalization certificate" begin
    cases = (_q_degree_case_data_from_fixture(), _q_degree_hand_checked_case())

    for case in cases
        @test degree(case.q, 1) >= degree(case.p, 1)
        _assert_q_degree_certificate(case)
    end

    case = first(cases)
    non_monic = 2 * case.X^2 + one(case.R)
    non_monic_q = one(case.R)
    non_monic_r = -one(case.R)
    non_monic_s = zero(case.R)
    non_monic_target = _q_degree_target(case.R, non_monic, non_monic_q, non_monic_r, non_monic_s)
    @test det(non_monic_target) == one(case.R)
    @test_throws ArgumentError Suslin.normalize_sl3_local_q_degree(
        non_monic,
        non_monic_q,
        non_monic_r,
        non_monic_s,
        case.X,
    )

    cert = Suslin.normalize_sl3_local_q_degree(case.p, case.q, case.r, case.s, case.X)
    corrupted = Suslin.SL3LocalRealizationCertificate(
        cert.target,
        cert.branch,
        cert.factors,
        cert.selected_variable,
        merge(cert.witness, (; remainder = cert.witness.remainder + one(case.R))),
    )
    @test det(corrupted.target) == one(case.R)
    @test !Suslin.verify_sl3_local_realization(corrupted)

    extra_field = Suslin.SL3LocalRealizationCertificate(
        cert.target,
        cert.branch,
        cert.factors,
        cert.selected_variable,
        merge(cert.witness, (; extra = one(case.R))),
    )
    @test !Suslin.verify_sl3_local_realization(extra_field)
end
```

- [ ] **Step 2: Run focused test to verify RED**

Run:

```bash
julia --project=. -e 'include("test/expert/sl3_local_q_degree_normalization.jl")'
```

Expected: FAIL with `UndefVarError` or equivalent because `normalize_sl3_local_q_degree` does not exist yet.

---

### Task 2: Q-Degree Helper And Replay Branch

**Files:**
- Modify: `src/algorithm/sl3_local.jl`
- Modify: `test/expert/sl3_local_q_degree_normalization.jl`

**Interfaces:**
- Consumes: `_sl3_local_target_entries(target)`, `_sl3_local_factor_sequences_equal(left, right)`, and `elementary_matrix`.
- Produces: `_recognize_sl3_local_special_form_matrix(A, X; check_monic=true)`.
- Produces: `_recognize_sl3_local_special_form_parameters(p, q, r, s, X; check_monic=true)`.
- Produces: `normalize_sl3_local_q_degree(A, X; check_monic::Bool=true)::SL3LocalRealizationCertificate`.
- Produces: `normalize_sl3_local_q_degree(p, q, r, s, X; check_monic::Bool=true)::SL3LocalRealizationCertificate`.
- Extends: `verify_sl3_local_realization(cert)::Bool` for `cert.branch == :q_degree_normalization`.

- [ ] **Step 1: Factor out general special-form recognition**

In `src/algorithm/sl3_local.jl`, replace the bodies of
`_recognize_sl3_local_matrix` and `_recognize_sl3_local_parameters` with a
general special-form validator plus the existing supported-family selection:

```julia
function _recognize_sl3_local_matrix(A, X; check_monic::Bool=true)
    form = _recognize_sl3_local_special_form_matrix(A, X; check_monic)
    return _recognize_sl3_local_supported_family(form)
end

function _recognize_sl3_local_parameters(p, q, r, s, X; check_monic::Bool=true)
    form = _recognize_sl3_local_special_form_parameters(p, q, r, s, X; check_monic)
    return _recognize_sl3_local_supported_family(form)
end

function _recognize_sl3_local_special_form_matrix(A, X; check_monic::Bool=true)
    nrows(A) == 3 || throw(ArgumentError("local SL_3 special-form recognition failed: A must be 3x3"))
    ncols(A) == 3 || throw(ArgumentError("local SL_3 special-form recognition failed: A must be 3x3"))

    R = base_ring(A)
    parent(X) === R || throw(ArgumentError("A and X must lie in the same polynomial or Laurent polynomial ring"))
    if A[1, 3] != zero(R) || A[2, 3] != zero(R) ||
            A[3, 1] != zero(R) || A[3, 2] != zero(R) ||
            A[3, 3] != one(R)
        throw(ArgumentError("local SL_3 special-form recognition failed: A must be an embedded 2x2 block with trailing identity"))
    end

    return _recognize_sl3_local_special_form_parameters(A[1, 1], A[1, 2], A[2, 1], A[2, 2], X; check_monic)
end

function _recognize_sl3_local_special_form_parameters(p, q, r, s, X; check_monic::Bool=true)
    R = parent(p)
    parent(q) === R || throw(ArgumentError("p, q, r, s, and X must lie in the same polynomial or Laurent polynomial ring"))
    parent(r) === R || throw(ArgumentError("p, q, r, s, and X must lie in the same polynomial or Laurent polynomial ring"))
    parent(s) === R || throw(ArgumentError("p, q, r, s, and X must lie in the same polynomial or Laurent polynomial ring"))
    parent(X) === R || throw(ArgumentError("p, q, r, s, and X must lie in the same polynomial or Laurent polynomial ring"))

    ring_gens = collect(gens(R))
    var_idx = findfirst(isequal(X), ring_gens)
    var_idx === nothing && throw(ArgumentError("X must be one of the polynomial or Laurent polynomial ring generators"))

    if check_monic
        if _is_laurent_polynomial_ring(R)
            throw(ArgumentError("p monicity check is only supported for ordinary polynomial local SL_3 inputs; pass check_monic=false only when the caller has discharged the local monicity assumption"))
        end
        _is_monic_in_variable(p, var_idx, R) || throw(ArgumentError("p must be monic in X"))
    end

    target = matrix(R, [
        p q zero(R);
        r s zero(R);
        zero(R) zero(R) one(R)
    ])
    det(target) == one(R) || throw(ArgumentError("constructed matrix must have determinant 1"))

    return (; R, p, q, r, s, X, target)
end

function _recognize_sl3_local_supported_family(form)
    if form.s == one(form.R) && form.p == one(form.R) + form.q * form.r
        return merge(form, (; family = :open_s_one))
    end

    if form.p == one(form.R) && form.s == one(form.R) + form.q * form.r
        return merge(form, (; family = :open_p_one))
    end

    s_inverse = _unit_inverse_or_nothing(form.s)
    if s_inverse !== nothing
        return merge(form, (; family = :s_unit, pivot_inverse = s_inverse))
    end

    p_inverse = _unit_inverse_or_nothing(form.p)
    if p_inverse !== nothing
        return merge(form, (; family = :p_unit, pivot_inverse = p_inverse))
    end

    _throw_staged_sl3_local_failure("supported families require one open unipotent slice or a unit diagonal pivot p or s")
end
```

- [ ] **Step 2: Add q-degree helper methods**

In `src/algorithm/sl3_local.jl`, insert these methods after the existing
`realize_sl3_local_certificate` methods:

```julia
function normalize_sl3_local_q_degree(A, X; check_monic::Bool=true)
    form = _recognize_sl3_local_special_form_matrix(A, X; check_monic)
    return _normalize_sl3_local_q_degree_form(form)
end

function normalize_sl3_local_q_degree(p, q, r, s, X; check_monic::Bool=true)
    form = _recognize_sl3_local_special_form_parameters(p, q, r, s, X; check_monic)
    return _normalize_sl3_local_q_degree_form(form)
end
```

- [ ] **Step 3: Add normalization construction helpers**

In `src/algorithm/sl3_local.jl`, insert these helpers after
`_realize_sl3_local_form(form)`:

```julia
function _normalize_sl3_local_q_degree_form(form)
    _sl3_local_require_ordinary_polynomial_ring(form.R, "q-degree normalization")
    quotient, remainder = divrem(form.q, form.p)
    p_degree = degree(form.p, _sl3_local_variable_index(form.R, form.X))
    remainder_degree = _sl3_local_degree_in_variable(remainder, form.R, form.X)
    remainder_degree < p_degree ||
        throw(ArgumentError("q-degree normalization failed to reduce q below p in X"))

    normalized_s = form.s - quotient * form.r
    normalized_target = matrix(form.R, [
        form.p remainder zero(form.R);
        form.r normalized_s zero(form.R);
        zero(form.R) zero(form.R) one(form.R)
    ])
    correction = elementary_matrix(3, 1, 2, quotient, form.R)
    normalized_target * correction == form.target ||
        error("internal q-degree normalization reconstruction failed")

    certificate = SL3LocalRealizationCertificate(
        form.target,
        :q_degree_normalization,
        [correction],
        form.X,
        (;
            p = form.p,
            q = form.q,
            r = form.r,
            s = form.s,
            quotient,
            remainder,
            normalized_target,
            correction,
            p_degree,
            remainder_degree,
        ),
    )
    verify_sl3_local_realization(certificate) ||
        error("internal q-degree normalization certificate verification failed")
    return certificate
end

function _sl3_local_require_ordinary_polynomial_ring(R, label::AbstractString)
    _is_laurent_polynomial_ring(R) &&
        throw(ArgumentError("local SL_3 $(label) is only supported for ordinary polynomial rings"))
    return nothing
end

function _sl3_local_variable_index(R, X)
    var_idx = findfirst(isequal(X), collect(gens(R)))
    var_idx === nothing &&
        throw(ArgumentError("X must be one of the polynomial ring generators"))
    return var_idx
end

function _sl3_local_degree_in_variable(value, R, X)
    return degree(value, _sl3_local_variable_index(R, X))
end
```

- [ ] **Step 4: Extend verifier for partial q-degree certificates**

In `_sl3_local_realization_verification(certificate)`, replace:

```julia
factors_ok = size_ok && verify_factorization(target, certificate.factors)
```

with:

```julia
factors_ok = size_ok && (
    certificate.branch == :q_degree_normalization ?
    _sl3_local_q_degree_reconstruction_ok(certificate) :
    verify_factorization(target, certificate.factors)
)
```

In `_sl3_local_certificate_expected_factors(certificate)`, add this branch before
the final `end`:

```julia
elseif certificate.branch == :q_degree_normalization
    return [witness.correction]
```

In `_sl3_local_branch_witness_ok(certificate)`, add this branch before the final
`end`:

```julia
elseif certificate.branch == :q_degree_normalization
    return _sl3_local_q_degree_witness_ok(certificate)
```

In `_sl3_local_witness_keys_ok(branch, witness)`, add this branch before the
final `return false`:

```julia
elseif branch == :q_degree_normalization
    return keys(witness) == (
        :p,
        :q,
        :r,
        :s,
        :quotient,
        :remainder,
        :normalized_target,
        :correction,
        :p_degree,
        :remainder_degree,
    )
```

Then add these helpers after `_sl3_local_witness_keys_ok`:

```julia
function _sl3_local_q_degree_witness_ok(certificate)
    entries = _sl3_local_target_entries(certificate.target)
    entries === nothing && return false

    R = base_ring(certificate.target)
    witness = certificate.witness
    witness.p == entries.p || return false
    witness.q == entries.q || return false
    witness.r == entries.r || return false
    witness.s == entries.s || return false
    parent(witness.quotient) === R || return false
    parent(witness.remainder) === R || return false
    witness.q == witness.quotient * witness.p + witness.remainder || return false
    witness.p_degree == _sl3_local_degree_in_variable(witness.p, R, certificate.selected_variable) ||
        return false
    witness.remainder_degree ==
        _sl3_local_degree_in_variable(witness.remainder, R, certificate.selected_variable) ||
        return false
    witness.remainder_degree < witness.p_degree || return false

    expected_normalized = matrix(R, [
        witness.p witness.remainder zero(R);
        witness.r witness.s - witness.quotient * witness.r zero(R);
        zero(R) zero(R) one(R)
    ])
    witness.normalized_target == expected_normalized || return false
    witness.correction == elementary_matrix(3, 1, 2, witness.quotient, R) || return false
    return _sl3_local_q_degree_reconstruction_ok(certificate)
end

function _sl3_local_q_degree_reconstruction_ok(certificate)
    witness = certificate.witness
    return witness.normalized_target * witness.correction == certificate.target &&
        _sl3_local_factor_sequences_equal(certificate.factors, [witness.correction])
end
```

- [ ] **Step 5: Run focused test to verify GREEN**

Run:

```bash
julia --project=. -e 'include("test/expert/sl3_local_q_degree_normalization.jl")'
```

Expected: PASS with the `local SL3 q-degree normalization certificate` testset.

---

### Task 3: Expert Test Registration And Verification

**Files:**
- Modify: `test/runtests.jl`
- Modify: `docs/superpowers/plans/2026-06-21-issue-71-sl3-local-q-degree-normalization.md`

**Interfaces:**
- Consumes: focused test from Task 1 and helper/verifier implementation from Task 2.
- Produces: registered expert test and committed implementation.

- [ ] **Step 1: Register the expert test**

In `test/runtests.jl`, add this entry immediately after
`"expert/sl3_local_certificate.jl",`:

```julia
"expert/sl3_local_q_degree_normalization.jl",
```

- [ ] **Step 2: Run expert and package verification**

Run:

```bash
julia --project=. -e 'include("test/expert/sl3_local_q_degree_normalization.jl")'
julia --project=. test/runtests.jl expert
julia --project=. -e 'using Pkg; Pkg.test()'
git diff --check
```

Expected:

- Focused q-degree test exits 0.
- Expert suite exits 0.
- Package test exits 0.
- `git diff --check` exits 0.

- [ ] **Step 3: Commit implementation**

Stage and commit:

```bash
git add \
    docs/superpowers/plans/2026-06-21-issue-71-sl3-local-q-degree-normalization.md \
    src/algorithm/sl3_local.jl \
    test/expert/sl3_local_q_degree_normalization.jl \
    test/runtests.jl
git commit -m "Add SL3 local q-degree normalization"
```

Expected: commit succeeds with only the plan, source, focused expert test, and test registration staged.

---

## Self-Review

- Spec coverage: the tasks cover a polynomial-only helper, fixture reuse, two exact q-degree examples, quotient/remainder degree checks, normalized target and correction recording, exact replay through `verify_sl3_local_realization`, non-monic rejection, corrupted witness rejection, focused verification, package verification, and no public driver change.
- Placeholder scan: no deferred implementation markers remain.
- Type consistency: `normalize_sl3_local_q_degree`, `SL3LocalRealizationCertificate`, `verify_sl3_local_realization`, `:q_degree_normalization`, and witness field names are used consistently across test and implementation steps.
