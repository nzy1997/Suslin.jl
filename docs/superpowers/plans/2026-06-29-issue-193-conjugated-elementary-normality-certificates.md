# Issue 193 Conjugated Elementary Normality Certificates Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add replayable ordinary-polynomial normality certificates for conjugates of elementary matrices.

**Architecture:** `src/algorithm/normality.jl` owns validation, convention metadata, vector extraction, and replay verification for `A * E_ij(a) * inv(A)`. It delegates the rank-one factor proof to `realize_rank_one_normality_certificate`, so this layer does not duplicate Cohn-type or rank-one factor formulas.

**Tech Stack:** Julia, Oscar/AbstractAlgebra ordinary polynomial matrices, existing Suslin elementary matrix helpers, #190 polynomial normality fixtures, #192 rank-one normality certificates, `Test`.

## Global Constraints

- Certificate inputs are ordinary-polynomial square matrices `A` with `nrows(A) == ncols(A) >= 3`, `det(A) == one(base_ring(A))`, distinct indices `i` and `j`, and a coefficient `a` coercible into `base_ring(A)`.
- The stored convention must be explicit as `conjugation_convention = :A_E_invA`.
- The stored target must be named `conjugation_target` and must equal `A * elementary_matrix(n, i, j, a, R) * inv(A)`.
- The certificate must store `A`, `i`, `j`, `a`, `ring`, `determinant`, `inverse_A`, `elementary_matrix`, `conjugation_convention`, `conjugation_target`, `v`, `w`, `g`, a child rank-one certificate, final `factors`, reconstructed `product`, and verification metadata.
- The derived vectors must be `v = A[:, i]`, `w = a * inv(A)[j, :]`, and `g = inv(A)[i, :]`.
- The child rank-one certificate must be produced by `realize_rank_one_normality_certificate(v, w, g, R)`.
- `verify_conjugate_elementary_certificate(cert)` must accept only exact replay, including stored verification metadata equality.
- Keep `realize_conjugate_elementary(B, i, j, a)` working for current callers, including existing Laurent behavior.
- Tests must include non-fixture ordinary-polynomial `SL_3` and `SL_4` examples where certificate factors multiply exactly to the stored convention target.
- Negative controls must reject determinant-not-one matrices, equal indices, coefficients from the wrong ring, and a certificate with a tampered factor.
- Do not implement Quillen patching, Murthy local solving, ECP reducers, recursive `SL_n` factorization, Laurent certificates, ToricBuilder support, or Steinberg factor-count optimization.

---

## File Structure

- Modify `test/expert/normality.jl`: add failing certificate-focused tests around the existing raw helper tests.
- Modify `test/public/api_surface.jl`: add failing public export checks for the new type and functions.
- Modify `src/algorithm/normality.jl`: add `ConjugatedElementaryNormalityCertificate`, checked-data helper, constructor, core verifier, and public verifier before the existing raw factor helper.
- Modify `src/Suslin.jl`: export `ConjugatedElementaryNormalityCertificate`, `realize_conjugate_elementary_certificate`, and `verify_conjugate_elementary_certificate`.

### Task 1: Certificate Tests

**Files:**
- Modify: `test/expert/normality.jl`
- Modify: `test/public/api_surface.jl`

**Interfaces:**
- Consumes: `Suslin.realize_conjugate_elementary_certificate(A, i, j, a)`, `Suslin.verify_conjugate_elementary_certificate(cert)`, `Suslin.ConjugatedElementaryNormalityCertificate`, `Suslin.verify_rank_one_normality_certificate(cert)`, and the existing `realize_conjugate_elementary` raw helper.
- Produces: Focused tests that fail before the certificate API exists.

- [ ] **Step 1: Add test helpers in `test/expert/normality.jl`**

After the existing `product_of_factors(factors)` helper, add:

```julia
function product_of_factors(factors, R, n::Int)
    product = identity_matrix(R, n)
    for factor in factors
        product *= factor
    end
    return product
end

function replace_conjugate_certificate_field(cert, field::Symbol, value)
    names = fieldnames(typeof(cert))
    idx = findfirst(==(field), names)
    idx === nothing && throw(ArgumentError("certificate field $(field) not found"))
    values = ntuple(k -> k == idx ? value : getfield(cert, names[k]), length(names))
    return typeof(cert)(values...)
end

function assert_conjugate_certificate(cert, A, i::Int, j::Int, a)
    R = base_ring(A)
    n = nrows(A)
    inverse_A = inv(A)
    elementary = elementary_matrix(n, i, j, a, R)
    target = A * elementary * inverse_A
    coerced_a = R(a)
    expected_v = [A[row, i] for row in 1:n]
    expected_w = [coerced_a * inverse_A[j, col] for col in 1:n]
    expected_g = [inverse_A[i, col] for col in 1:n]

    @test cert isa Suslin.ConjugatedElementaryNormalityCertificate
    @test cert.n == n
    @test cert.A == A
    @test cert.i == i
    @test cert.j == j
    @test cert.a == coerced_a
    @test cert.ring == R
    @test cert.determinant == one(R)
    @test cert.inverse_A == inverse_A
    @test cert.elementary_matrix == elementary
    @test cert.conjugation_convention == :A_E_invA
    @test cert.conjugation_target == target
    @test cert.v == expected_v
    @test cert.w == expected_w
    @test cert.g == expected_g
    @test Suslin.verify_rank_one_normality_certificate(cert.rank_one_certificate)
    @test cert.rank_one_certificate.v == expected_v
    @test cert.rank_one_certificate.w == expected_w
    @test cert.rank_one_certificate.g == expected_g
    @test cert.rank_one_certificate.target == target
    @test cert.factors == cert.rank_one_certificate.factors
    @test cert.product == product_of_factors(cert.factors, R, n)
    @test cert.product == cert.conjugation_target
    @test cert.verification.target_matches_product_ok
    @test Suslin.verify_conjugate_elementary_certificate(cert)
end
```

- [ ] **Step 2: Add the failing certificate tests in `test/expert/normality.jl`**

Inside the existing `"constructive normality"` testset, after the current raw-helper checks, add:

```julia
    cert = Suslin.realize_conjugate_elementary_certificate(B, 1, 3, x + 1)
    assert_conjugate_certificate(cert, B, 1, 3, x + 1)
    @test Suslin.realize_conjugate_elementary(B, 1, 3, x + 1) == cert.factors

    cert4 = Suslin.realize_conjugate_elementary_certificate(B4, 1, 4, x + 1)
    assert_conjugate_certificate(cert4, B4, 1, 4, x + 1)
    @test Suslin.realize_conjugate_elementary(B4, 1, 4, x + 1) == cert4.factors

    fixture_path = joinpath(@__DIR__, "..", "fixtures", "polynomial_normality_cases.jl")
    if !isdefined(Main, :ParkWoodburnPolynomialNormalityFixtureCatalog)
        include(fixture_path)
    end
    fixture = Main.ParkWoodburnPolynomialNormalityFixtureCatalog.cases_by_id()["pw-section2-conjugated-elementary-qq"]
    fixture_cert = Suslin.realize_conjugate_elementary_certificate(
        fixture.inputs.B,
        fixture.inputs.i,
        fixture.inputs.j,
        fixture.inputs.a,
    )
    assert_conjugate_certificate(
        fixture_cert,
        fixture.inputs.B,
        fixture.inputs.i,
        fixture.inputs.j,
        fixture.inputs.a,
    )
    @test fixture_cert.conjugation_target == fixture.target_matrix

    non_sl = matrix(R, [
        2 0 0;
        0 1 0;
        0 0 1
    ])
    @test_throws ArgumentError Suslin.realize_conjugate_elementary_certificate(non_sl, 1, 2, one(R))
    @test_throws ArgumentError Suslin.realize_conjugate_elementary_certificate(B, 1, 1, x + 1)

    L, (u,) = Oscar.laurent_polynomial_ring(QQ, ["u"])
    @test_throws ArgumentError Suslin.realize_conjugate_elementary_certificate(B, 1, 3, u)

    tampered_factor_cert = Suslin.realize_conjugate_elementary_certificate(B, 1, 3, x + 1)
    tampered_factor_cert.factors[1][1, 1] += one(R)
    @test !Suslin.verify_conjugate_elementary_certificate(tampered_factor_cert)

    tampered_convention_cert = replace_conjugate_certificate_field(cert, :conjugation_convention, :invA_E_A)
    @test !Suslin.verify_conjugate_elementary_certificate(tampered_convention_cert)
```

- [ ] **Step 3: Add failing public API checks**

In `test/public/api_surface.jl`, add these `isdefined` checks near the existing normality exports:

```julia
    @test isdefined(Suslin, :ConjugatedElementaryNormalityCertificate)
    @test isdefined(Suslin, :realize_conjugate_elementary_certificate)
    @test isdefined(Suslin, :verify_conjugate_elementary_certificate)
```

Add these equality checks near the existing `realize_conjugate_elementary` equality check:

```julia
    @test Suslin.ConjugatedElementaryNormalityCertificate === ConjugatedElementaryNormalityCertificate
    @test Suslin.realize_conjugate_elementary_certificate === realize_conjugate_elementary_certificate
    @test Suslin.verify_conjugate_elementary_certificate === verify_conjugate_elementary_certificate
```

- [ ] **Step 4: Run the focused test to verify RED**

Run:

```bash
julia --project=. -e 'include("test/expert/normality.jl")'
```

Expected: FAIL with an undefined `realize_conjugate_elementary_certificate` or `ConjugatedElementaryNormalityCertificate` error.

- [ ] **Step 5: Commit the red tests**

Run:

```bash
git add test/expert/normality.jl test/public/api_surface.jl
git commit -m "test: cover conjugated elementary certificates"
```

Expected: commit succeeds with only test files staged.

### Task 2: Certificate API

**Files:**
- Modify: `src/algorithm/normality.jl`
- Modify: `src/Suslin.jl`
- Test: `test/expert/normality.jl`
- Test: `test/public/api_surface.jl`

**Interfaces:**
- Consumes: `elementary_matrix(n, i, j, a, R)`, `_require_ordinary_polynomial_certificate_ring(R)`, `_coerce_into_ring(R, value, label)`, `_same_base_ring(left, right)`, `_cohn_type_factor_product(factors, R, n)`, `realize_rank_one_normality_certificate(v, w, g, R)`, and `verify_rank_one_normality_certificate(cert)`.
- Produces: `ConjugatedElementaryNormalityCertificate`, `realize_conjugate_elementary_certificate(A, i, j, a)`, and `verify_conjugate_elementary_certificate(cert)::Bool`.

- [ ] **Step 1: Add exports in `src/Suslin.jl`**

Add these exports near the existing normality certificate exports:

```julia
export ConjugatedElementaryNormalityCertificate
export realize_conjugate_elementary_certificate
export verify_conjugate_elementary_certificate
```

- [ ] **Step 2: Add the certificate type in `src/algorithm/normality.jl`**

Add this code after `verify_rank_one_normality_certificate` and before `realize_conjugate_elementary`:

```julia
struct ConjugatedElementaryNormalityCertificate
    n::Int
    A
    i::Int
    j::Int
    a
    ring
    determinant
    inverse_A
    elementary_matrix
    conjugation_convention::Symbol
    conjugation_target
    v::Vector
    w::Vector
    g::Vector
    rank_one_certificate::RankOneNormalityCertificate
    factors::Vector
    product
    verification
end

function _conjugate_elementary_checked_data(A, i::Int, j::Int, a)
    nrows(A) == ncols(A) || throw(ArgumentError("A must be square"))

    n = nrows(A)
    n >= 3 || throw(ArgumentError("A must have size at least 3"))
    1 <= i <= n || throw(ArgumentError("i must be between 1 and n"))
    1 <= j <= n || throw(ArgumentError("j must be between 1 and n"))
    i == j && throw(ArgumentError("i and j must differ"))

    R = base_ring(A)
    _require_ordinary_polynomial_certificate_ring(R)
    coerced_a = _coerce_into_ring(R, a, "a")
    determinant = det(A)
    determinant == one(R) || throw(ArgumentError("A must have determinant one"))
    inverse_A = inv(A)
    identity = identity_matrix(R, n)
    A * inverse_A == identity || throw(ArgumentError("A inverse replay failed"))
    inverse_A * A == identity || throw(ArgumentError("A inverse replay failed"))

    elementary = elementary_matrix(n, i, j, coerced_a, R)
    conjugation_target = A * elementary * inverse_A
    v = [A[row, i] for row in 1:n]
    w = [coerced_a * inverse_A[j, col] for col in 1:n]
    g = [inverse_A[i, col] for col in 1:n]
    rank_one_certificate = realize_rank_one_normality_certificate(v, w, g, R)
    rank_one_certificate.target == conjugation_target ||
        error("internal conjugated elementary target did not match rank-one target")

    return (;
        n,
        A,
        i,
        j,
        a = coerced_a,
        ring = R,
        determinant,
        inverse_A,
        elementary_matrix = elementary,
        conjugation_convention = :A_E_invA,
        conjugation_target,
        v,
        w,
        g,
        rank_one_certificate,
    )
end
```

- [ ] **Step 3: Add constructor and verifier in `src/algorithm/normality.jl`**

Add this code after `_conjugate_elementary_checked_data`:

```julia
function realize_conjugate_elementary_certificate(A, i::Int, j::Int, a)
    data = _conjugate_elementary_checked_data(A, i, j, a)
    factors = data.rank_one_certificate.factors
    product = _cohn_type_factor_product(factors, data.ring, data.n)
    provisional = ConjugatedElementaryNormalityCertificate(
        data.n,
        data.A,
        data.i,
        data.j,
        data.a,
        data.ring,
        data.determinant,
        data.inverse_A,
        data.elementary_matrix,
        data.conjugation_convention,
        data.conjugation_target,
        data.v,
        data.w,
        data.g,
        data.rank_one_certificate,
        factors,
        product,
        nothing,
    )
    verification = _conjugate_elementary_certificate_core_verification(provisional)
    verification.overall_core_ok ||
        error("internal conjugated elementary certificate verification failed")
    return ConjugatedElementaryNormalityCertificate(
        provisional.n,
        provisional.A,
        provisional.i,
        provisional.j,
        provisional.a,
        provisional.ring,
        provisional.determinant,
        provisional.inverse_A,
        provisional.elementary_matrix,
        provisional.conjugation_convention,
        provisional.conjugation_target,
        provisional.v,
        provisional.w,
        provisional.g,
        provisional.rank_one_certificate,
        provisional.factors,
        provisional.product,
        verification,
    )
end

function _conjugate_elementary_certificate_core_verification(cert)
    ring_ok = false
    checked_inputs_ok = false
    determinant_ok = false
    inverse_replay_ok = false
    elementary_matrix_ok = false
    convention_ok = false
    target_replay_ok = false
    vector_replay_ok = false
    rank_one_certificate_ok = false
    factor_sequence_ok = false
    factor_count = 0
    product_replay_ok = false
    product_matches_stored_ok = false
    target_matches_product_ok = false

    data = _conjugate_elementary_checked_data(cert.A, cert.i, cert.j, cert.a)
    checked_inputs_ok = true
    ring_ok = _same_base_ring(cert.ring, data.ring)
    determinant_ok = cert.determinant == one(data.ring) && cert.determinant == data.determinant
    identity = identity_matrix(data.ring, data.n)
    inverse_replay_ok =
        cert.inverse_A == data.inverse_A &&
        cert.A * cert.inverse_A == identity &&
        cert.inverse_A * cert.A == identity
    elementary_matrix_ok = cert.elementary_matrix == data.elementary_matrix
    convention_ok =
        cert.conjugation_convention == :A_E_invA &&
        cert.conjugation_convention == data.conjugation_convention
    target_replay_ok = cert.conjugation_target == data.conjugation_target
    vector_replay_ok = cert.v == data.v && cert.w == data.w && cert.g == data.g
    rank_one_certificate_ok =
        verify_rank_one_normality_certificate(cert.rank_one_certificate) &&
        cert.rank_one_certificate.n == data.n &&
        cert.rank_one_certificate.v == data.v &&
        cert.rank_one_certificate.w == data.w &&
        cert.rank_one_certificate.g == data.g &&
        _same_base_ring(cert.rank_one_certificate.ring, data.ring) &&
        cert.rank_one_certificate.target == data.conjugation_target &&
        cert.rank_one_certificate.factors == data.rank_one_certificate.factors &&
        cert.rank_one_certificate.product == data.rank_one_certificate.product &&
        cert.rank_one_certificate.verification == data.rank_one_certificate.verification
    factor_sequence_ok = cert.factors == cert.rank_one_certificate.factors
    factor_count = length(cert.factors)
    replayed_product = _cohn_type_factor_product(cert.factors, data.ring, data.n)
    product_replay_ok = true
    product_matches_stored_ok = replayed_product == cert.product
    target_matches_product_ok = target_replay_ok && replayed_product == cert.conjugation_target

    overall_core_ok =
        ring_ok &&
        checked_inputs_ok &&
        determinant_ok &&
        inverse_replay_ok &&
        elementary_matrix_ok &&
        convention_ok &&
        target_replay_ok &&
        vector_replay_ok &&
        rank_one_certificate_ok &&
        factor_sequence_ok &&
        product_replay_ok &&
        product_matches_stored_ok &&
        target_matches_product_ok

    return (;
        ring_ok,
        checked_inputs_ok,
        determinant_ok,
        inverse_replay_ok,
        elementary_matrix_ok,
        convention_ok,
        target_replay_ok,
        vector_replay_ok,
        rank_one_certificate_ok,
        factor_sequence_ok,
        factor_count,
        product_replay_ok,
        product_matches_stored_ok,
        target_matches_product_ok,
        overall_core_ok,
    )
end

function _conjugate_elementary_certificate_verification(cert)
    core = _conjugate_elementary_certificate_core_verification(cert)
    stored_verification_ok = cert.verification == core
    return merge(core, (;
        stored_verification_ok,
        overall_ok = core.overall_core_ok && stored_verification_ok,
    ))
end

function verify_conjugate_elementary_certificate(cert)::Bool
    try
        return _conjugate_elementary_certificate_verification(cert).overall_ok
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end
```

- [ ] **Step 4: Run focused tests to verify GREEN**

Run:

```bash
julia --project=. -e 'include("test/expert/normality.jl")'
```

Expected: PASS.

- [ ] **Step 5: Run public API tests**

Run:

```bash
julia --project=. -e 'include("test/public/api_surface.jl")'
```

Expected: PASS.

- [ ] **Step 6: Commit implementation**

Run:

```bash
git add src/algorithm/normality.jl src/Suslin.jl
git commit -m "feat: add conjugated elementary certificates"
```

Expected: commit succeeds with only implementation files staged.

## Final Verification

After both tasks are complete, run:

```bash
julia --project=. -e 'include("test/expert/normality.jl")'
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: both commands exit 0. If the full package command fails on a test group unrelated to this issue, capture the failure output before deciding whether it blocks the PR.

## Plan Self-Review

- The plan covers the issue's positive `SL_3` and `SL_4` certificate examples, explicit convention metadata, rank-one child certificate reuse, raw helper compatibility, and requested negative controls.
- The plan keeps the certificate route ordinary-polynomial only and does not alter the existing Laurent raw helper route.
- The plan uses TDD: Task 1 adds failing tests before Task 2 adds production code.
- No task implements out-of-scope Quillen, Murthy, ECP, recursive `SL_n`, Laurent certificate, ToricBuilder, or Steinberg optimization work.
