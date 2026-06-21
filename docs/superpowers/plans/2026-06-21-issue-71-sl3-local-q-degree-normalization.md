# Issue 71 SL3 Local q-Degree Normalization Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an internal, replayable Murthy q-degree normalization step for local polynomial `SL_3` special-form inputs.

**Architecture:** Add a small normalization record in `src/algorithm/sl3_local.jl`, a selected-variable monic division helper, and a certificate branch that replays `[normalized_target, E_12(f)]`. Keep `realize_sl3_local` and `elementary_factorization` behavior unchanged.

**Tech Stack:** Julia, Oscar/AbstractAlgebra polynomial rings and matrices, existing Suslin elementary matrix and local certificate helpers.

## Global Constraints

- Do not implement q(0)-unit recursion, split-lemma recursion, or the resultant/Bezout branch in this issue.
- Do not change the public `elementary_factorization` driver.
- Keep the helper polynomial-only; Laurent inputs must be rejected.
- Do not export the new expert helpers from `src/Suslin.jl`.
- Every stored quotient, remainder, normalized target, elementary correction, and selected variable must participate in replay verification.
- Use exact polynomial operations and exponent-vector coefficient extraction; do not parse or manipulate polynomial strings.

---

### Task 1: Focused Expert Test

**Files:**
- Create: `test/expert/sl3_local_q_degree_normalization.jl`

**Interfaces:**
- Consumes: `Suslin.sl3_local_q_degree_normalization`, `Suslin.sl3_local_q_degree_normalization_certificate`, `Suslin.verify_sl3_local_q_degree_normalization`, `Suslin.verify_sl3_local_realization`, `Suslin.SL3LocalQDegreeNormalization`.
- Produces: Failing coverage for the exact API and replay behavior.

- [ ] **Step 1: Write the failing test**

```julia
using Test
using Suslin
using Oscar

const SL3_Q_DEGREE_FIXTURE_PATH =
    joinpath(@__DIR__, "..", "fixtures", "sl3_murthy_gupta_cases.jl")

function _qdegree_target(R, p, q, r, s)
    return matrix(R, [
        p q zero(R);
        r s zero(R);
        zero(R) zero(R) one(R)
    ])
end

function _qdegree_product(factors, R)
    product = identity_matrix(R, 3)
    for factor in factors
        product *= factor
    end
    return product
end

function _degree_in_variable(value, X)
    var_idx = findfirst(isequal(X), collect(gens(parent(value))))
    return degree(value, var_idx)
end

function _assert_qdegree_record(record, p, q, r, s, X, expected_f, expected_g)
    R = parent(p)
    @test record.target == _qdegree_target(R, p, q, r, s)
    @test record.quotient == expected_f
    @test record.remainder == expected_g
    @test q == record.quotient * p + record.remainder
    @test _degree_in_variable(record.remainder, X) < _degree_in_variable(p, X)
    @test record.normalized_target ==
        _qdegree_target(R, p, record.remainder, r, s - record.quotient * r)
    @test record.elementary_correction ==
        Suslin.elementary_matrix(3, 1, 2, record.quotient, R)
    @test record.normalized_target * record.elementary_correction == record.target
    @test Suslin.verify_sl3_local_q_degree_normalization(record)
end

@testset "Murthy q-degree normalization replay for local SL3" begin
    include(SL3_Q_DEGREE_FIXTURE_PATH)
    catalog = SL3MurthyGuptaFixtureCatalog.catalog()
    by_id = Dict(entry.id => entry for entry in catalog.cases)

    fixture = by_id["mg-q-degree-normalization"]
    witness = first(fixture.witnesses)
    record = Suslin.sl3_local_q_degree_normalization(
        fixture.entries.p,
        fixture.entries.q,
        fixture.entries.r,
        fixture.entries.s,
        fixture.variable,
    )
    _assert_qdegree_record(
        record,
        fixture.entries.p,
        fixture.entries.q,
        fixture.entries.r,
        fixture.entries.s,
        fixture.variable,
        witness.quotient,
        witness.remainder,
    )
    @test record.normalized_target[2, 2] == witness.normalized_s

    cert = Suslin.sl3_local_q_degree_normalization_certificate(record)
    @test cert.branch == :murthy_q_degree_normalization
    @test cert.witness.normalization == record
    @test _qdegree_product(cert.factors, fixture.ring.object) == record.target
    @test Suslin.verify_sl3_local_realization(cert)

    matrix_record = Suslin.sl3_local_q_degree_normalization(fixture.target, fixture.variable)
    @test matrix_record == record

    R, (X,) = Oscar.polynomial_ring(QQ, ["X"])
    p = X^2 + X + 1
    r = one(R)
    s = X^2
    q = p * s - one(R)
    second = Suslin.sl3_local_q_degree_normalization(p, q, r, s, X)
    _assert_qdegree_record(second, p, q, r, s, X, X^2, -one(R))
    second_cert = Suslin.sl3_local_q_degree_normalization_certificate(p, q, r, s, X)
    @test Suslin.verify_sl3_local_realization(second_cert)

    nonmonic_p = 2 * X + one(R)
    nonmonic_q = 2 * X
    @test det(_qdegree_target(R, nonmonic_p, nonmonic_q, one(R), one(R))) == one(R)
    @test_throws ArgumentError Suslin.sl3_local_q_degree_normalization(
        nonmonic_p,
        nonmonic_q,
        one(R),
        one(R),
        X,
    )

    bad_record = Suslin.SL3LocalQDegreeNormalization(
        second.target,
        second.quotient + one(R),
        second.remainder,
        second.normalized_target,
        second.elementary_correction,
        second.selected_variable,
    )
    @test det(bad_record.target) == one(R)
    @test !Suslin.verify_sl3_local_q_degree_normalization(bad_record)
    bad_cert = Suslin.SL3LocalRealizationCertificate(
        bad_record.target,
        :murthy_q_degree_normalization,
        [bad_record.normalized_target, bad_record.elementary_correction],
        bad_record.selected_variable,
        (; normalization = bad_record),
    )
    @test !Suslin.verify_sl3_local_realization(bad_cert)
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project=. -e 'include("test/expert/sl3_local_q_degree_normalization.jl")'`

Expected: FAIL with `UndefVarError` for `sl3_local_q_degree_normalization` or `SL3LocalQDegreeNormalization`.

- [ ] **Step 3: Commit**

```bash
git add test/expert/sl3_local_q_degree_normalization.jl
git commit -m "test: cover SL3 q-degree normalization replay"
```

### Task 2: Normalization Record and Certificate Replay

**Files:**
- Modify: `src/algorithm/sl3_local.jl`

**Interfaces:**
- Consumes: failing test from Task 1.
- Produces:
  - `SL3LocalQDegreeNormalization`
  - `sl3_local_q_degree_normalization`
  - `sl3_local_q_degree_normalization_certificate`
  - `verify_sl3_local_q_degree_normalization`
  - `:murthy_q_degree_normalization` branch support in `verify_sl3_local_realization`

- [ ] **Step 1: Add the record type and construction entry points**

Add the struct near `SL3LocalRealizationCertificate`:

```julia
struct SL3LocalQDegreeNormalization
    target
    quotient
    remainder
    normalized_target
    elementary_correction
    selected_variable
end
```

Add matrix and parameter helpers:

```julia
function sl3_local_q_degree_normalization(A, X; check_monic::Bool=true)
    entries = _sl3_local_target_entries(A)
    entries === nothing && throw(ArgumentError("q-degree normalization requires a local SL_3 special-form matrix"))
    return sl3_local_q_degree_normalization(entries.p, entries.q, entries.r, entries.s, X; check_monic)
end

function sl3_local_q_degree_normalization(p, q, r, s, X; check_monic::Bool=true)
    form = _recognize_sl3_local_q_degree_normalization_parameters(p, q, r, s, X; check_monic)
    quotient, remainder = _sl3_local_divrem_monic_in_variable(form.q, form.p, form.var_idx, form.R)
    normalized_target = _sl3_local_special_form_target(
        form.R,
        form.p,
        remainder,
        form.r,
        form.s - quotient * form.r,
    )
    record = SL3LocalQDegreeNormalization(
        form.target,
        quotient,
        remainder,
        normalized_target,
        elementary_matrix(3, 1, 2, quotient, form.R),
        form.X,
    )
    verify_sl3_local_q_degree_normalization(record) ||
        error("internal Murthy q-degree normalization verification failed")
    return record
end
```

- [ ] **Step 2: Add selected-variable division and recognition helpers**

Implement helpers in `src/algorithm/sl3_local.jl` using existing monicity and
target helpers:

```julia
function _recognize_sl3_local_q_degree_normalization_parameters(p, q, r, s, X; check_monic::Bool=true)
    R = parent(p)
    parent(q) === R || throw(ArgumentError("p, q, r, s, and X must lie in the same polynomial ring"))
    parent(r) === R || throw(ArgumentError("p, q, r, s, and X must lie in the same polynomial ring"))
    parent(s) === R || throw(ArgumentError("p, q, r, s, and X must lie in the same polynomial ring"))
    parent(X) === R || throw(ArgumentError("p, q, r, s, and X must lie in the same polynomial ring"))
    _is_laurent_polynomial_ring(R) &&
        throw(ArgumentError("q-degree normalization is only supported for ordinary polynomial local SL_3 inputs"))

    ring_gens = collect(gens(R))
    var_idx = findfirst(isequal(X), ring_gens)
    var_idx === nothing && throw(ArgumentError("X must be one of the polynomial ring generators"))
    if check_monic
        _is_monic_in_variable(p, var_idx, R) || throw(ArgumentError("p must be monic in X"))
    end

    target = _sl3_local_special_form_target(R, p, q, r, s)
    det(target) == one(R) || throw(ArgumentError("constructed matrix must have determinant 1"))
    return (; R, p, q, r, s, X, var_idx, target)
end

function _sl3_local_divrem_monic_in_variable(q, p, var_idx::Int, R)
    degree_p = degree(p, var_idx)
    degree_p >= 0 || throw(ArgumentError("p must have nonnegative degree in X"))
    quotient = zero(R)
    remainder = q
    X = collect(gens(R))[var_idx]
    while !iszero(remainder) && degree(remainder, var_idx) >= degree_p
        degree_gap = degree(remainder, var_idx) - degree_p
        leading = _sl3_local_coefficient_in_variable_degree(remainder, var_idx, degree(remainder, var_idx), R)
        term = leading * X^degree_gap
        quotient += term
        remainder -= term * p
    end
    return quotient, remainder
end

function _sl3_local_coefficient_in_variable_degree(value, var_idx::Int, target_degree::Int, R)
    ring_gens = collect(gens(R))
    total = zero(R)
    for (coeff, exponents) in zip(AbstractAlgebra.coefficients(value), AbstractAlgebra.exponent_vectors(value))
        exponents[var_idx] == target_degree || continue
        term = R(coeff)
        for idx in eachindex(ring_gens)
            idx == var_idx && continue
            exponent = exponents[idx]
            exponent == 0 && continue
            term *= ring_gens[idx]^exponent
        end
        total += term
    end
    return total
end
```

Then update `_is_monic_in_variable` to call
`_sl3_local_coefficient_in_variable_degree(p, var_idx, target_degree, R) == one(R)`.

- [ ] **Step 3: Add record verification and certificate helpers**

Add:

```julia
function verify_sl3_local_q_degree_normalization(record)::Bool
    try
        return _sl3_local_q_degree_normalization_verification(record).overall_ok
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end

function sl3_local_q_degree_normalization_certificate(record::SL3LocalQDegreeNormalization)
    verify_sl3_local_q_degree_normalization(record) ||
        throw(ArgumentError("Murthy q-degree normalization record must verify before certificate construction"))
    certificate = SL3LocalRealizationCertificate(
        record.target,
        :murthy_q_degree_normalization,
        [record.normalized_target, record.elementary_correction],
        record.selected_variable,
        (; normalization = record),
    )
    verify_sl3_local_realization(certificate) ||
        error("internal Murthy q-degree normalization certificate verification failed")
    return certificate
end

function sl3_local_q_degree_normalization_certificate(args...; check_monic::Bool=true)
    return sl3_local_q_degree_normalization_certificate(
        sl3_local_q_degree_normalization(args...; check_monic),
    )
end
```

Extend `_sl3_local_witness_keys_ok`, `_sl3_local_branch_witness_ok`, and
`_sl3_local_certificate_expected_factors` for `:murthy_q_degree_normalization`.

- [ ] **Step 4: Run focused test to verify it passes**

Run: `julia --project=. -e 'include("test/expert/sl3_local_q_degree_normalization.jl")'`

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/algorithm/sl3_local.jl
git commit -m "feat: add SL3 q-degree normalization replay"
```

### Task 3: Test Registration and Verification

**Files:**
- Modify: `test/runtests.jl`

**Interfaces:**
- Consumes: passing focused expert test.
- Produces: expert suite coverage for the new normalization test.

- [ ] **Step 1: Register the expert test**

Add `"expert/sl3_local_q_degree_normalization.jl"` after
`"expert/sl3_local_certificate.jl"` in `TEST_GROUP_FILES["expert"]`.

- [ ] **Step 2: Run focused and package verification**

Run:

```bash
julia --project=. -e 'include("test/expert/sl3_local_q_degree_normalization.jl")'
julia --project=. -e 'using Pkg; Pkg.test()'
git diff --check
```

Expected: all commands exit 0.

- [ ] **Step 3: Commit**

```bash
git add test/runtests.jl
git commit -m "test: register SL3 q-degree normalization expert test"
```

## Plan Self-Review

- The plan covers the issue objective, the #69 fixture reuse requirement, the
  #70 certificate replay requirement, negative controls, and public API
  guardrails.
- No task solves recursive Murthy branches or changes `elementary_factorization`.
- All new replay data has an explicit verification check.
