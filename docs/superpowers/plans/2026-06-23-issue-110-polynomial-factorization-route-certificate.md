# Issue 110 Polynomial Factorization Route Certificate Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an internal replayable route certificate for currently supported ordinary-polynomial factorization routes.

**Architecture:** Keep the certificate in `src/algorithm/factorization.jl` beside `elementary_factorization` and `verify_factorization`. The constructor builds route-specific evidence using existing local `SL_3` and `SL_n` reduction APIs, while the verifier recomputes product and evidence checks instead of trusting stored metadata.

**Tech Stack:** Julia, Oscar polynomial matrices, Suslin internal factorization helpers, Test stdlib.

## Global Constraints

- Preserve `elementary_factorization(A)` as factor-returning.
- Keep route certificate APIs non-exported.
- Support only ordinary-polynomial univariate local `SL_3`, block-local `reduce_sln_to_sl3`, and staged failure records.
- Do not add recursive column-peel support.
- Do not consume Quillen patches.
- For successful certificates, replay must check route tag, determinant-one precondition, route-specific evidence, stored factor product, and exact equality to `A`.
- Required focused verification command: `julia --project=. -e 'include("test/expert/park_woodburn_route_certificate.jl")'`.
- Required package verification command: `julia --project=. -e 'using Pkg; Pkg.test()'`.

---

## File Structure

- Modify `src/algorithm/factorization.jl`: define `PolynomialFactorizationRouteCertificate`, `_polynomial_factorization_route_certificate`, `_verify_polynomial_factorization_route_certificate`, and internal replay helpers.
- Create `test/expert/park_woodburn_route_certificate.jl`: test supported fixture certificates and tamper controls.
- Modify `test/runtests.jl`: register the expert test.

---

### Task 1: Add Expert Route Certificate Tests First

**Files:**
- Create: `test/expert/park_woodburn_route_certificate.jl`
- Modify: `test/runtests.jl`

**Interfaces:**
- Consumes: `ParkWoodburnPolynomialFixtureCatalog.cases_by_id()`.
- Consumes planned internal helpers:
  - `Suslin._polynomial_factorization_route_certificate(A; route=nothing)`
  - `Suslin._verify_polynomial_factorization_route_certificate(cert)::Bool`
  - `Suslin.PolynomialFactorizationRouteCertificate`
- Produces: failing expert coverage for issue 110 before implementation.

- [ ] **Step 1: Write the failing test file**

Create `test/expert/park_woodburn_route_certificate.jl` with helper functions that:

```julia
using Test
using Suslin
using Oscar

const PARK_WOODBURN_ROUTE_CATALOG_PATH =
    joinpath(@__DIR__, "..", "fixtures", "park_woodburn_polynomial_cases.jl")

function _pw_route_product(factors, R, n::Int)
    product = identity_matrix(R, n)
    for factor in factors
        product *= factor
    end
    return product
end

function _pw_replace_certificate(
        cert;
        matrix = cert.matrix,
        route = cert.route,
        factors = cert.factors,
        product = cert.product,
        evidence = cert.evidence,
        status = cert.status,
        verification = cert.verification)
    return Suslin.PolynomialFactorizationRouteCertificate(
        matrix,
        route,
        factors,
        product,
        evidence,
        status,
        verification,
    )
end
```

The positive assertions must:

```julia
include(PARK_WOODBURN_ROUTE_CATALOG_PATH)
entries = ParkWoodburnPolynomialFixtureCatalog.cases_by_id()
fast_entry = entries["pw-poly-univariate-sl3-fast-local-qq"]
block_entry = entries["pw-poly-univariate-sln-disjoint-blocks-qq"]

fast_cert = Suslin._polynomial_factorization_route_certificate(
    fast_entry.matrix;
    route = fast_entry.route,
)
@test fast_cert.route == :fast_local_sl3
@test Suslin._verify_polynomial_factorization_route_certificate(fast_cert)
@test Suslin.verify_sl3_local_realization(fast_cert.evidence)
@test _pw_route_product(fast_cert.factors, base_ring(fast_entry.matrix), nrows(fast_entry.matrix)) == fast_entry.matrix

block_cert = Suslin._polynomial_factorization_route_certificate(
    block_entry.matrix;
    route = block_entry.route,
)
@test block_cert.route == :disjoint_local_blocks
@test Suslin._verify_polynomial_factorization_route_certificate(block_cert)
@test Suslin.verify_sln_to_sl3_reduction(block_cert.evidence)
@test _pw_route_product(block_cert.factors, base_ring(block_entry.matrix), nrows(block_entry.matrix)) == block_entry.matrix
```

The negative controls must construct separate tampered certificates for a wrong route tag, a wrong factor, a wrong stored product, and wrong evidence:

```julia
R = base_ring(fast_cert.matrix)
n = nrows(fast_cert.matrix)
bad_route = _pw_replace_certificate(fast_cert; route = :quillen_patched_substitution)
@test !Suslin._verify_polynomial_factorization_route_certificate(bad_route)

bad_factors = copy(fast_cert.factors)
bad_factors[1] = identity_matrix(R, n)
bad_factor_cert = _pw_replace_certificate(fast_cert; factors = bad_factors)
@test !Suslin._verify_polynomial_factorization_route_certificate(bad_factor_cert)

bad_product = identity_matrix(R, n)
bad_product_cert = _pw_replace_certificate(fast_cert; product = bad_product)
@test !Suslin._verify_polynomial_factorization_route_certificate(bad_product_cert)

bad_evidence = Suslin.SL3LocalRealizationCertificate(
    fast_cert.evidence.target,
    fast_cert.evidence.branch,
    fast_cert.evidence.factors,
    fast_cert.evidence.selected_variable,
    merge(fast_cert.evidence.witness, (; q = fast_cert.evidence.witness.q + one(R))),
)
bad_evidence_cert = _pw_replace_certificate(fast_cert; evidence = bad_evidence)
@test Suslin.verify_factorization(bad_evidence_cert.matrix, bad_evidence_cert.factors)
@test !Suslin._verify_polynomial_factorization_route_certificate(bad_evidence_cert)
```

- [ ] **Step 2: Register the expert test**

Add `"expert/park_woodburn_route_certificate.jl"` to `TEST_GROUP_FILES["expert"]` in `test/runtests.jl` immediately after `"expert/sln_to_sl3_reduction.jl"`.

- [ ] **Step 3: Run focused test and confirm RED**

Run:

```bash
julia --project=. -e 'include("test/expert/park_woodburn_route_certificate.jl")'
```

Expected: FAIL because `Suslin._polynomial_factorization_route_certificate` and `Suslin.PolynomialFactorizationRouteCertificate` do not exist yet.

---

### Task 2: Implement Internal Route Certificate Replay

**Files:**
- Modify: `src/algorithm/factorization.jl`

**Interfaces:**
- Produces:
  - `struct PolynomialFactorizationRouteCertificate`
  - `_polynomial_factorization_route_certificate(A; route=nothing)`
  - `_verify_polynomial_factorization_route_certificate(cert)::Bool`

- [ ] **Step 1: Add the certificate type and constructor dispatch**

Add the struct and route constants near the top of `src/algorithm/factorization.jl`:

```julia
struct PolynomialFactorizationRouteCertificate
    matrix
    route::Symbol
    factors::Vector
    product
    evidence
    status::Symbol
    verification
end

const _POLYNOMIAL_FACTORIZATION_ROUTE_TAGS = Set([
    :fast_local_sl3,
    :disjoint_local_blocks,
    :staged_failure,
])
```

Add `_polynomial_factorization_route_certificate(A; route=nothing)` that validates square ordinary-polynomial determinant-one input, then builds `:fast_local_sl3`, `:disjoint_local_blocks`, or `:staged_failure` certificates.

- [ ] **Step 2: Add product and evidence replay helpers**

Add helpers that recompute products, compare factor sequences, and verify evidence:

```julia
function _polynomial_route_factor_product(factors, R, n::Int)
    product = identity_matrix(R, n)
    for factor in factors
        nrows(factor) == n || throw(ArgumentError("route certificate factor has wrong row count"))
        ncols(factor) == n || throw(ArgumentError("route certificate factor has wrong column count"))
        _same_base_ring(base_ring(factor), R) ||
            throw(ArgumentError("route certificate factor has wrong base ring"))
        product *= factor
    end
    return product
end
```

Successful route evidence must require:

```julia
cert.route == :fast_local_sl3 &&
    cert.evidence isa SL3LocalRealizationCertificate &&
    cert.evidence.target == cert.matrix &&
    verify_sl3_local_realization(cert.evidence) &&
    _polynomial_route_factor_sequences_equal(cert.factors, cert.evidence.factors)
```

and:

```julia
cert.route == :disjoint_local_blocks &&
    cert.evidence isa SLNToSL3Reduction &&
    cert.evidence.original_matrix == cert.matrix &&
    verify_sln_to_sl3_reduction(cert.evidence) &&
    _polynomial_route_factor_sequences_equal(cert.factors, cert.evidence.factors)
```

- [ ] **Step 3: Add fail-closed verifier**

Add `_verify_polynomial_factorization_route_certificate(cert)::Bool` that catches non-interrupt exceptions and returns false. It must recompute a core verification record and require the stored `cert.verification` to match the recomputed core record. Successful route verification requires determinant one, product equals stored product, product equals matrix, route-specific evidence verifies, and `verify_factorization(cert.matrix, cert.factors)` returns true.

- [ ] **Step 4: Run focused test and confirm GREEN**

Run:

```bash
julia --project=. -e 'include("test/expert/park_woodburn_route_certificate.jl")'
```

Expected: PASS.

---

### Task 3: Final Verification and Commit

**Files:**
- Modify: implementation and test files from Tasks 1-2.

**Interfaces:**
- Consumes completed route certificate implementation.
- Produces a committed branch ready for PR.

- [ ] **Step 1: Run focused issue verification**

Run:

```bash
julia --project=. -e 'include("test/expert/park_woodburn_route_certificate.jl")'
```

Expected: PASS.

- [ ] **Step 2: Run full package verification**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: PASS.

- [ ] **Step 3: Commit implementation**

Run:

```bash
git add src/algorithm/factorization.jl test/expert/park_woodburn_route_certificate.jl test/runtests.jl docs/superpowers/plans/2026-06-23-issue-110-polynomial-factorization-route-certificate.md
git commit -m "feat: add polynomial route certificates"
```

Expected: commit succeeds.

---

## Plan Self-Review

- The plan covers the design spec and issue 110 verification command.
- The test is written and run before production implementation.
- The route certificate stays non-exported and does not alter public API return types.
- The negative controls cover route tag, factor, product, and evidence tampering.
- Recursive column peel and Quillen support remain out of scope.
