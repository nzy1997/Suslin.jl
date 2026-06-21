# Issue 70 SL3 Local Certificates Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add replayable local `SL_3` realization certificates for the existing local solver families while keeping `realize_sl3_local(...)` factor-returning.

**Architecture:** Reuse the existing local `SL_3` recognition path and add a thin non-exported certificate/result path in `src/algorithm/sl3_local.jl`. The legacy API delegates to the certificate path and returns `certificate.factors`; expert tests exercise the certificate path and replay verifier directly through `Suslin.<name>`.

**Tech Stack:** Julia, Oscar polynomial rings, Suslin elementary matrices, Test stdlib.

## Global Constraints

- Do not implement Murthy-Gupta recursive branches.
- Do not change `elementary_factorization` routing.
- Do not optimize factor counts.
- Preserve existing public behavior: `realize_sl3_local(...)` returns only the factor sequence.
- Keep the certificate layer thin and internal; do not export the new names from `src/Suslin.jl`.
- Every certificate field must participate in `verify_sl3_local_realization`.
- Supported certificate branches for this issue are exactly `:open_s_one`, `:open_p_one`, `:s_unit`, and `:p_unit`.
- Focused verification command is `julia --project=. -e 'include("test/expert/sl3_local_certificate.jl")'`.
- Required package verification command is `julia --project=. -e 'using Pkg; Pkg.test()'`.

---

## File Structure

- Modify `src/algorithm/sl3_local.jl`: add `SL3LocalRealizationCertificate`, expert certificate constructors, branch witnesses, and verifier; refactor legacy factor path to delegate to certificate path.
- Create `test/expert/sl3_local_certificate.jl`: focused expert tests for open-slice and unit-pivot certificates, legacy API preservation, and tamper rejection.
- Modify `test/runtests.jl`: register `expert/sl3_local_certificate.jl` in the expert group.

---

### Task 1: Replayable Local SL3 Certificate Path

**Files:**
- Modify: `src/algorithm/sl3_local.jl`
- Create: `test/expert/sl3_local_certificate.jl`
- Modify: `test/runtests.jl`

**Interfaces:**
- Consumes: existing `_recognize_sl3_local_matrix`, `_recognize_sl3_local_parameters`, `_sl3_diagonal_unit_factors`, `_verify_sl3_local_factorization`, `verify_factorization`, and `elementary_matrix`.
- Produces: non-exported `SL3LocalRealizationCertificate`.
- Produces: non-exported `realize_sl3_local_certificate(A, X; check_monic=true)`.
- Produces: non-exported `realize_sl3_local_certificate(p, q, r, s, X; check_monic=true)`.
- Produces: non-exported `verify_sl3_local_realization(cert)::Bool`.

- [ ] **Step 1: Write the failing expert test**

Create `test/expert/sl3_local_certificate.jl` with this structure:

```julia
using Test
using Suslin
using Oscar

function _sl3_certificate_target(p, q, r, s, R)
    return matrix(R, [
        p q zero(R);
        r s zero(R);
        zero(R) zero(R) one(R)
    ])
end

function _sl3_certificate_product(factors, R)
    product = identity_matrix(R, 3)
    for factor in factors
        product *= factor
    end
    return product
end

function _assert_sl3_certificate_replays(cert)
    R = base_ring(cert.target)
    @test _sl3_certificate_product(cert.factors, R) == cert.target
    @test Suslin.verify_factorization(cert.target, cert.factors)
    @test Suslin.verify_sl3_local_realization(cert)
end

function _tamper_first_factor(cert)
    R = base_ring(cert.target)
    factors = copy(cert.factors)
    factors[1] = identity_matrix(R, 3)
    return Suslin.SL3LocalRealizationCertificate(
        cert.target,
        cert.branch,
        factors,
        cert.selected_variable,
        cert.witness,
    )
end

function _tamper_witness_q(cert)
    R = base_ring(cert.target)
    return Suslin.SL3LocalRealizationCertificate(
        cert.target,
        cert.branch,
        cert.factors,
        cert.selected_variable,
        merge(cert.witness, (; q = cert.witness.q + one(R))),
    )
end

@testset "local SL3 realization certificates" begin
    R, (X,) = Oscar.polynomial_ring(QQ, ["X"])

    q_open = one(R)
    r_open = X
    p_open = one(R) + q_open * r_open
    s_open = one(R)
    open_target = _sl3_certificate_target(p_open, q_open, r_open, s_open, R)
    open_cert = Suslin.realize_sl3_local_certificate(p_open, q_open, r_open, s_open, X)
    @test open_cert.target == open_target
    @test open_cert.branch == :open_s_one
    @test open_cert.selected_variable == X
    @test open_cert.witness.q == q_open
    @test open_cert.witness.r == r_open
    _assert_sl3_certificate_replays(open_cert)

    legacy_factors = Suslin.realize_sl3_local(p_open, q_open, r_open, s_open, X)
    @test legacy_factors isa Vector
    @test _sl3_certificate_product(legacy_factors, R) == open_target
    @test Suslin.verify_factorization(open_target, legacy_factors)

    q_dual = X + one(R)
    r_dual = X
    p_dual = one(R)
    s_dual = one(R) + q_dual * r_dual
    dual_cert = Suslin.realize_sl3_local_certificate(p_dual, q_dual, r_dual, s_dual, X)
    @test dual_cert.branch == :open_p_one
    _assert_sl3_certificate_replays(dual_cert)

    s_unit = R(2)
    q_s_unit = one(R)
    r_s_unit = 2 * X
    p_s_unit = X + R(1 // 2)
    s_unit_target = _sl3_certificate_target(p_s_unit, q_s_unit, r_s_unit, s_unit, R)
    s_unit_cert = Suslin.realize_sl3_local_certificate(s_unit_target, X)
    @test s_unit_cert.branch == :s_unit
    @test s_unit_cert.witness.pivot == s_unit
    @test s_unit_cert.witness.pivot_inverse == inv(s_unit)
    _assert_sl3_certificate_replays(s_unit_cert)

    p_unit = R(2)
    q_p_unit = X
    r_p_unit = one(R)
    s_p_unit = (X + one(R)) * R(1 // 2)
    p_unit_target = _sl3_certificate_target(p_unit, q_p_unit, r_p_unit, s_p_unit, R)
    p_unit_cert = Suslin.realize_sl3_local_certificate(
        p_unit,
        q_p_unit,
        r_p_unit,
        s_p_unit,
        X;
        check_monic = false,
    )
    @test p_unit_cert.target == p_unit_target
    @test p_unit_cert.branch == :p_unit
    @test p_unit_cert.witness.pivot == p_unit
    @test p_unit_cert.witness.pivot_inverse == inv(p_unit)
    _assert_sl3_certificate_replays(p_unit_cert)

    @test !Suslin.verify_sl3_local_realization(_tamper_first_factor(open_cert))
    @test !Suslin.verify_sl3_local_realization(_tamper_witness_q(open_cert))
end
```

- [ ] **Step 2: Run focused test to verify RED**

Run:

```bash
julia --project=. -e 'include("test/expert/sl3_local_certificate.jl")'
```

Expected: FAIL with `UndefVarError` or equivalent because `realize_sl3_local_certificate` and `SL3LocalRealizationCertificate` do not exist yet.

- [ ] **Step 3: Implement the certificate path**

In `src/algorithm/sl3_local.jl`, add this struct near the top of the file:

```julia
struct SL3LocalRealizationCertificate
    target
    branch::Symbol
    factors::Vector
    selected_variable
    witness
end
```

Add expert constructors next to the existing `realize_sl3_local` methods:

```julia
function realize_sl3_local_certificate(A, X; check_monic::Bool=true)
    form = _recognize_sl3_local_matrix(A, X; check_monic)
    return _realize_sl3_local_certificate_form(form)
end

function realize_sl3_local_certificate(p, q, r, s, X; check_monic::Bool=true)
    form = _recognize_sl3_local_parameters(p, q, r, s, X; check_monic)
    return _realize_sl3_local_certificate_form(form)
end
```

Refactor the legacy methods to delegate without changing their return type:

```julia
function realize_sl3_local(A, X; check_monic::Bool=true)
    return realize_sl3_local_certificate(A, X; check_monic).factors
end

function realize_sl3_local(p, q, r, s, X; check_monic::Bool=true)
    return realize_sl3_local_certificate(p, q, r, s, X; check_monic).factors
end
```

Move the branch factor logic from `_realize_sl3_local_form` into a helper:

```julia
function _sl3_local_form_factors(form)
    R = form.R
    if form.family == :open_s_one
        return [
            elementary_matrix(3, 1, 2, form.q, R),
            elementary_matrix(3, 2, 1, form.r, R),
        ]
    elseif form.family == :open_p_one
        return [
            elementary_matrix(3, 2, 1, form.r, R),
            elementary_matrix(3, 1, 2, form.q, R),
        ]
    elseif form.family == :s_unit
        s_inverse = form.pivot_inverse
        return vcat(
            [elementary_matrix(3, 1, 2, form.q * s_inverse, R)],
            _sl3_diagonal_unit_factors(s_inverse, R),
            [elementary_matrix(3, 2, 1, form.r * s_inverse, R)],
        )
    elseif form.family == :p_unit
        p_inverse = form.pivot_inverse
        return vcat(
            [elementary_matrix(3, 2, 1, form.r * p_inverse, R)],
            _sl3_diagonal_unit_factors(form.p, R),
            [elementary_matrix(3, 1, 2, form.q * p_inverse, R)],
        )
    end
    _throw_staged_sl3_local_failure("unrecognized local solver family $(form.family)")
end
```

Build branch witnesses from the recognized form:

```julia
function _sl3_local_form_witness(form)
    if form.family in (:open_s_one, :open_p_one)
        return (; q = form.q, r = form.r)
    elseif form.family == :s_unit
        return (; pivot = form.s, pivot_inverse = form.pivot_inverse)
    elseif form.family == :p_unit
        return (; pivot = form.p, pivot_inverse = form.pivot_inverse)
    end
    _throw_staged_sl3_local_failure("unrecognized local solver family $(form.family)")
end
```

Add certificate construction and legacy wrapper:

```julia
function _realize_sl3_local_certificate_form(form)
    factors = _sl3_local_form_factors(form)
    _verify_sl3_local_factorization(form.target, factors)
    certificate = SL3LocalRealizationCertificate(
        form.target,
        form.family,
        factors,
        form.X,
        _sl3_local_form_witness(form),
    )
    verify_sl3_local_realization(certificate) ||
        error("internal local SL_3 realization certificate verification failed")
    return certificate
end

function _realize_sl3_local_form(form)
    return _realize_sl3_local_certificate_form(form).factors
end
```

The existing recognition named tuple must include the selected variable. In both
recognized branches in `_recognize_sl3_local_parameters`, add `X` to every
returned named tuple:

```julia
return (; family = :open_s_one, R, p, q, r, s, X, target)
return (; family = :open_p_one, R, p, q, r, s, X, target)
return (; family = :s_unit, R, p, q, r, s, X, target, pivot_inverse = s_inverse)
return (; family = :p_unit, R, p, q, r, s, X, target, pivot_inverse = p_inverse)
```

Add replay helpers below `_verify_sl3_local_factorization`:

```julia
function verify_sl3_local_realization(certificate)::Bool
    try
        return _sl3_local_realization_verification(certificate).overall_ok
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end

function _sl3_local_realization_verification(certificate)
    target = certificate.target
    size_ok = nrows(target) == 3 && ncols(target) == 3
    R = size_ok ? base_ring(target) : nothing
    variable_ok =
        size_ok &&
        parent(certificate.selected_variable) === R &&
        certificate.selected_variable in collect(gens(R))
    shape_ok = size_ok && _sl3_local_target_entries(target) !== nothing
    determinant_ok = size_ok && det(target) == one(R)
    expected_factors =
        size_ok && variable_ok && shape_ok && determinant_ok ?
        _sl3_local_certificate_expected_factors(certificate) :
        nothing
    witness_ok = expected_factors !== nothing
    factors_match_ok = witness_ok && _sl3_local_factor_sequences_equal(certificate.factors, expected_factors)
    factors_ok = size_ok && verify_factorization(target, certificate.factors)
    overall_ok = size_ok && variable_ok && shape_ok && determinant_ok &&
        witness_ok && factors_match_ok && factors_ok
    return (;
        overall_ok,
        size_ok,
        variable_ok,
        shape_ok,
        determinant_ok,
        witness_ok,
        factors_match_ok,
        factors_ok,
    )
end
```

Implement `_sl3_local_realization_verification` with these exact checks:

- `target` is `3 x 3`.
- `selected_variable` lies in the same base ring as `target` and is one of
  `gens(R)`.
- `target` has local special-form shape and determinant one.
- Branch-specific witness fields match target entries and branch equations.
- Expected factors recomputed from `branch`, `target`, and `witness` are equal
  to stored `factors` element-by-element.
- Stored factors multiply exactly to `target`.

Use helper functions with concrete names:

```julia
_sl3_local_target_entries(target)
_sl3_local_factor_sequences_equal(left, right)
_sl3_local_certificate_expected_factors(certificate)
_sl3_local_branch_witness_ok(certificate)
```

`_sl3_local_certificate_expected_factors` must use the same formulas as
`_sl3_local_form_factors`, but must read only from the certificate target and
witness:

```julia
if certificate.branch == :open_s_one
    witness.q == q && witness.r == r && s == one(R) && p == one(R) + q * r || return nothing
    return [elementary_matrix(3, 1, 2, witness.q, R), elementary_matrix(3, 2, 1, witness.r, R)]
elseif certificate.branch == :open_p_one
    witness.q == q && witness.r == r && p == one(R) && s == one(R) + q * r || return nothing
    return [elementary_matrix(3, 2, 1, witness.r, R), elementary_matrix(3, 1, 2, witness.q, R)]
elseif certificate.branch == :s_unit
    witness.pivot == s && witness.pivot * witness.pivot_inverse == one(R) || return nothing
    return vcat([elementary_matrix(3, 1, 2, q * witness.pivot_inverse, R)], _sl3_diagonal_unit_factors(witness.pivot_inverse, R), [elementary_matrix(3, 2, 1, r * witness.pivot_inverse, R)])
elseif certificate.branch == :p_unit
    witness.pivot == p && witness.pivot * witness.pivot_inverse == one(R) || return nothing
    return vcat([elementary_matrix(3, 2, 1, r * witness.pivot_inverse, R)], _sl3_diagonal_unit_factors(witness.pivot, R), [elementary_matrix(3, 1, 2, q * witness.pivot_inverse, R)])
end
return nothing
```

- [ ] **Step 4: Run focused test to verify GREEN**

Run:

```bash
julia --project=. -e 'include("test/expert/sl3_local_certificate.jl")'
```

Expected: PASS with the `local SL3 realization certificates` testset.

- [ ] **Step 5: Register the expert test**

Modify `test/runtests.jl` and add:

```julia
"expert/sl3_local_certificate.jl",
```

after:

```julia
"expert/sl3_local_extended.jl",
```

- [ ] **Step 6: Run expert group and package verification**

Run:

```bash
julia --project=. test/runtests.jl expert
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: both commands exit 0. The package command should still run the default public and internal groups.

- [ ] **Step 7: Commit implementation**

Stage exactly the design/plan/test/source changes that belong to this issue:

```bash
git add \
    docs/superpowers/plans/2026-06-21-issue-70-sl3-local-certificates.md \
    src/algorithm/sl3_local.jl \
    test/expert/sl3_local_certificate.jl \
    test/runtests.jl
git commit -m "Add replayable SL3 local certificates"
```

Expected: commit succeeds with the plan, certificate implementation, focused test, and test registration.

---

## Self-Review

- Spec coverage: the task preserves legacy factors, adds an expert certificate result and verifier, records target/branch/factors/selected variable/witnesses, checks branch-specific relations, rejects tampered factors and witnesses, and registers focused expert tests.
- Placeholder scan: no incomplete markers or deferred implementation steps remain.
- Type consistency: the plan consistently names `SL3LocalRealizationCertificate`, `realize_sl3_local_certificate`, and `verify_sl3_local_realization`.
