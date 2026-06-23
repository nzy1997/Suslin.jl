# Issue 112 Polynomial Column-Peel Certificate Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a replayable internal ordinary-polynomial recursive column-peel certificate through ECP for fixture-backed Park-Woodburn cases.

**Architecture:** Add a focused internal certificate implementation in `src/algorithm/polynomial_column_peel.jl`, mirroring the Laurent column-peel replay while enforcing ordinary-polynomial preconditions. Extend the existing issue 110 route certificate with a `:recursive_column_peel` evidence variant, but do not change public `elementary_factorization(A)`.

**Tech Stack:** Julia, Oscar exact polynomial matrices, Suslin ECP column reducer, existing polynomial route certificates, Test stdlib.

## Global Constraints

- Do not wire the route into public `elementary_factorization`.
- Do not implement Quillen patching.
- Do not handle Laurent matrices in the polynomial column-peel certificate.
- Keep new certificate APIs non-exported.
- Use `reduce_unimodular_column(v, R)` for the selected last column.
- Successful certificates must store peel steps, ECP factors, right-clearing factors, recursive next blocks, final local data, final factor sequence, product, and verification metadata.
- Required focused verification command: `julia --project=. -e 'include("test/expert/park_woodburn_polynomial_column_peel.jl")'`.
- Required package verification command: `julia --project=. -e 'using Pkg; Pkg.test()'`.

---

## File Structure

- Create `src/algorithm/polynomial_column_peel.jl`: internal peel step/certificate structs, constructor, recursion, replay, and verifier.
- Modify `src/Suslin.jl`: include `algorithm/polynomial_column_peel.jl` after the route certificate and `SL_n` reduction helpers.
- Modify `src/algorithm/factorization.jl`: add the `:recursive_column_peel` route tag, route constructor branch, and evidence verifier.
- Create `test/expert/park_woodburn_polynomial_column_peel.jl`: focused acceptance and tamper coverage.
- Modify `test/runtests.jl`: register the expert test after the route-certificate test.

---

### Task 1: Add Expert Polynomial Column-Peel Tests First

**Files:**
- Create: `test/expert/park_woodburn_polynomial_column_peel.jl`
- Modify: `test/runtests.jl`

**Interfaces:**
- Consumes planned internals:
  - `Suslin.PolynomialColumnPeelStep`
  - `Suslin.PolynomialColumnPeelCertificate`
  - `Suslin._polynomial_column_peel_certificate(A; final_route=nothing)`
  - `Suslin._verify_polynomial_column_peel_certificate(cert)::Bool`
  - `Suslin._polynomial_factorization_route_certificate(A; route=:recursive_column_peel)`
- Produces failing expert coverage for issue 112 before implementation.

- [ ] **Step 1: Create the failing expert test**

Add `test/expert/park_woodburn_polynomial_column_peel.jl` with:

```julia
using Test
using Suslin
using Oscar

const PARK_WOODBURN_POLY_PEEL_CATALOG_PATH =
    joinpath(@__DIR__, "..", "fixtures", "park_woodburn_polynomial_cases.jl")

function _pw_poly_peel_product(factors, R, n::Int)
    product = identity_matrix(R, n)
    for factor in factors
        product *= factor
    end
    return product
end

function _pw_poly_peel_target_column(R, n::Int)
    target = zero_matrix(R, n, 1)
    target[n, 1] = one(R)
    return target
end

function _pw_poly_peel_wrap_final_block(final_block, tail_coeffs)
    R = base_ring(final_block)
    n = nrows(final_block) + 1
    A = block_embedding(final_block, n, collect(1:(n - 1)))
    for row in 1:(n - 1)
        A[row, n] = tail_coeffs[row]
    end
    return A
end

function _pw_poly_replace_step(
        step;
        dimension = step.dimension,
        input_matrix = step.input_matrix,
        last_column = step.last_column,
        left_factors = step.left_factors,
        after_left_matrix = step.after_left_matrix,
        right_factors = step.right_factors,
        peeled_matrix = step.peeled_matrix,
        next_block = step.next_block)
    return Suslin.PolynomialColumnPeelStep(
        dimension,
        input_matrix,
        last_column,
        left_factors,
        after_left_matrix,
        right_factors,
        peeled_matrix,
        next_block,
    )
end

function _pw_poly_replace_certificate(
        cert;
        original_matrix = cert.original_matrix,
        peel_steps = cert.peel_steps,
        final_block = cert.final_block,
        final_certificate = cert.final_certificate,
        final_factors = cert.final_factors,
        factors = cert.factors,
        product = cert.product,
        verification = cert.verification)
    return Suslin.PolynomialColumnPeelCertificate(
        original_matrix,
        peel_steps,
        final_block,
        final_certificate,
        final_factors,
        factors,
        product,
        verification,
    )
end

function _pw_assert_real_peel_certificate(cert, A)
    R = base_ring(A)
    @test !isempty(cert.peel_steps)
    @test Suslin._verify_polynomial_column_peel_certificate(cert)
    @test cert.product == A
    @test _pw_poly_peel_product(cert.factors, R, nrows(A)) == A
    @test verify_factorization(A, cert.factors)
    for step in cert.peel_steps
        left_product = _pw_poly_peel_product(step.left_factors, R, step.dimension)
        recorded_column = matrix(R, step.dimension, 1, step.last_column)
        @test left_product * recorded_column == _pw_poly_peel_target_column(R, step.dimension)
        right_product = _pw_poly_peel_product(step.right_factors, R, step.dimension)
        @test step.after_left_matrix * right_product == step.peeled_matrix
        @test step.peeled_matrix == block_embedding(step.next_block, step.dimension, collect(1:(step.dimension - 1)))
    end
end

@testset "Park-Woodburn ordinary polynomial column-peel certificates" begin
    if !isdefined(Main, :ParkWoodburnPolynomialFixtureCatalog)
        include(PARK_WOODBURN_POLY_PEEL_CATALOG_PATH)
    end
    entries = ParkWoodburnPolynomialFixtureCatalog.cases_by_id()

    recursive_entry = entries["pw-poly-recursive-column-peel-gf2"]
    recursive_cert = Suslin._polynomial_column_peel_certificate(recursive_entry.matrix)
    _pw_assert_real_peel_certificate(recursive_cert, recursive_entry.matrix)
    @test recursive_cert.final_certificate.route == :fast_local_sl3

    fast_entry = entries["pw-poly-univariate-sl3-fast-local-qq"]
    R = base_ring(fast_entry.matrix)
    X = only(gens(R))
    wrapped = _pw_poly_peel_wrap_final_block(fast_entry.matrix, [X, X + one(R), X^2 + X])
    wrapped_cert = Suslin._polynomial_column_peel_certificate(wrapped)
    _pw_assert_real_peel_certificate(wrapped_cert, wrapped)
    @test wrapped_cert.final_block == fast_entry.matrix
    @test wrapped_cert.final_certificate.route == :fast_local_sl3

    route_cert = Suslin._polynomial_factorization_route_certificate(
        wrapped;
        route = :recursive_column_peel,
    )
    @test route_cert.route == :recursive_column_peel
    @test route_cert.evidence isa Suslin.PolynomialColumnPeelCertificate
    @test Suslin._verify_polynomial_factorization_route_certificate(route_cert)

    first_step = recursive_cert.peel_steps[1]
    bad_last = copy(first_step.last_column)
    bad_last[1] += one(base_ring(recursive_entry.matrix))
    bad_last_step = _pw_poly_replace_step(first_step; last_column = bad_last)
    bad_last_cert = _pw_poly_replace_certificate(
        recursive_cert;
        peel_steps = Suslin.PolynomialColumnPeelStep[bad_last_step],
        product = recursive_cert.product,
    )
    @test !Suslin._verify_polynomial_column_peel_certificate(bad_last_cert)

    bad_left = copy(first_step.left_factors)
    bad_left[1] = identity_matrix(base_ring(recursive_entry.matrix), first_step.dimension)
    bad_left_step = _pw_poly_replace_step(first_step; left_factors = bad_left)
    bad_left_cert = _pw_poly_replace_certificate(
        recursive_cert;
        peel_steps = Suslin.PolynomialColumnPeelStep[bad_left_step],
        product = recursive_cert.product,
    )
    @test !Suslin._verify_polynomial_column_peel_certificate(bad_left_cert)

    bad_right = copy(first_step.right_factors)
    push!(bad_right, elementary_matrix(first_step.dimension, first_step.dimension, 1, one(base_ring(recursive_entry.matrix)), base_ring(recursive_entry.matrix)))
    bad_right_step = _pw_poly_replace_step(first_step; right_factors = bad_right)
    bad_right_cert = _pw_poly_replace_certificate(
        recursive_cert;
        peel_steps = Suslin.PolynomialColumnPeelStep[bad_right_step],
        product = recursive_cert.product,
    )
    @test !Suslin._verify_polynomial_column_peel_certificate(bad_right_cert)

    bad_block = copy(first_step.next_block)
    bad_block[1, 1] += one(base_ring(recursive_entry.matrix))
    bad_block_step = _pw_poly_replace_step(first_step; next_block = bad_block)
    bad_block_cert = _pw_poly_replace_certificate(
        recursive_cert;
        peel_steps = Suslin.PolynomialColumnPeelStep[bad_block_step],
        product = recursive_cert.product,
    )
    @test !Suslin._verify_polynomial_column_peel_certificate(bad_block_cert)

    @test_throws ArgumentError Suslin._polynomial_column_peel_certificate(identity_matrix(base_ring(recursive_entry.matrix), 2))
    @test_throws ArgumentError Suslin._polynomial_column_peel_certificate(first_step.peeled_matrix)
end
```

- [ ] **Step 2: Register the expert test**

Add `"expert/park_woodburn_polynomial_column_peel.jl"` to `TEST_GROUP_FILES["expert"]` in `test/runtests.jl` immediately after `"expert/park_woodburn_route_certificate.jl"`.

- [ ] **Step 3: Run focused test and confirm RED**

Run:

```bash
julia --project=. -e 'include("test/expert/park_woodburn_polynomial_column_peel.jl")'
```

Expected: FAIL because `Suslin.PolynomialColumnPeelStep` and `Suslin._polynomial_column_peel_certificate` are not defined yet.

---

### Task 2: Implement Internal Polynomial Column-Peel Replay

**Files:**
- Create: `src/algorithm/polynomial_column_peel.jl`
- Modify: `src/Suslin.jl`
- Modify: `src/algorithm/factorization.jl`

**Interfaces:**
- Consumes: `reduce_unimodular_column`, `_polynomial_factorization_route_certificate`, `verify_factorization`, `block_embedding`, and elementary factor helpers.
- Produces:
  - `PolynomialColumnPeelStep`
  - `PolynomialColumnPeelCertificate`
  - `_polynomial_column_peel_certificate(A; final_route=nothing)`
  - `_verify_polynomial_column_peel_certificate(cert)::Bool`

- [ ] **Step 1: Add the focused implementation file**

Create `src/algorithm/polynomial_column_peel.jl` with structs and methods matching the test interfaces. Core functions must:

```julia
function _validate_polynomial_column_peel_input(A)
    n = _validate_factorization_matrix(A)
    R = base_ring(A)
    _factorization_ring_profile(R) == :polynomial ||
        throw(ArgumentError("polynomial column-peel certificates require an ordinary polynomial ring"))
    _require_polynomial_sl_determinant(A)
    return n
end

function _polynomial_column_peel_step(current)
    R = base_ring(current)
    d = nrows(current)
    last_column = [current[row, d] for row in 1:d]
    left_factors = reduce_unimodular_column(last_column, R)
    left_product = _factor_product(left_factors, R, d)
    after_left = left_product * current
    right_factors = _expected_column_peel_right_factors(after_left, d, R)
    right_product = _factor_product(right_factors, R, d)
    peeled = after_left * right_product
    next_block = matrix(R, [peeled[row, col] for row in 1:(d - 1), col in 1:(d - 1)])
    return PolynomialColumnPeelStep(d, current, last_column, left_factors, after_left, right_factors, peeled, next_block)
end
```

Use `_inverse_elementary_sequence`, `_embed_upper_left_factors`,
`_expected_column_peel_right_factors`, `_column_peel_target_column`, and
`_factor_product` consistently with the Laurent implementation.

- [ ] **Step 2: Add recursion and final-route stopping**

Implement `_polynomial_column_peel_recursive(current, final_route)` so it first
tries `_polynomial_column_peel_final_certificate(current, final_route)`. That
helper must call `_polynomial_factorization_route_certificate(current; route =
final_route, allow_recursive_column_peel = false)` and accept only
`:fast_local_sl3` and `:disjoint_local_blocks` supported certificates. If no
final route is available and the dimension is `3`, throw an `ArgumentError`.
Otherwise record one peel step and recurse on `step.next_block`.

- [ ] **Step 3: Add replay verifier**

Add `_is_valid_polynomial_column_peel_step_data`, `_replay_polynomial_column_peel_factors`, and `_polynomial_column_peel_verification`. The verifier must fail closed and return a NamedTuple containing at least `overall_ok`, `preconditions_ok`, `step_chain_ok`, `steps_ok`, `final_certificate_ok`, `factor_sequence_ok`, `product_ok`, and `factors_ok`.

- [ ] **Step 4: Include the file**

In `src/Suslin.jl`, add:

```julia
include("algorithm/polynomial_column_peel.jl")
```

after `include("algorithm/sln_to_sl3_reduction.jl")`.

- [ ] **Step 5: Extend route certificates**

In `src/algorithm/factorization.jl`:

```julia
const _POLYNOMIAL_FACTORIZATION_ROUTE_TAGS = Set([
    :fast_local_sl3,
    :disjoint_local_blocks,
    :recursive_column_peel,
    :staged_failure,
])
```

Change `_polynomial_factorization_route_certificate(A; route=nothing)` to
`_polynomial_factorization_route_certificate(A; route=nothing,
allow_recursive_column_peel::Bool=true)`. After the disjoint-local attempt
fails in auto mode, try `_polynomial_recursive_column_peel_route_certificate(A)`
when `allow_recursive_column_peel` is true. Add an explicit
`:recursive_column_peel` branch that also requires the flag.

Add:

```julia
function _polynomial_recursive_column_peel_route_certificate(A)
    evidence = _polynomial_column_peel_certificate(A)
    factors = copy(evidence.factors)
    product = _polynomial_route_factor_product(factors, base_ring(A), nrows(A))
    return _polynomial_route_certificate(A, :recursive_column_peel, factors, product, evidence, :supported)
end
```

and update `_polynomial_route_evidence_ok(cert)` to verify recursive evidence
with `_verify_polynomial_column_peel_certificate(cert.evidence)` and exact
factor sequence equality.

- [ ] **Step 6: Run focused test and confirm GREEN**

Run:

```bash
julia --project=. -e 'include("test/expert/park_woodburn_polynomial_column_peel.jl")'
```

Expected: PASS.

---

### Task 3: Final Verification And Commit

**Files:**
- Modify: files from Tasks 1-2.

**Interfaces:**
- Consumes completed certificate implementation.
- Produces a committed branch ready for PR.

- [ ] **Step 1: Run the focused issue verification**

Run:

```bash
julia --project=. -e 'include("test/expert/park_woodburn_polynomial_column_peel.jl")'
```

Expected: PASS.

- [ ] **Step 2: Run package verification**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: PASS.

- [ ] **Step 3: Commit implementation**

Run:

```bash
git add src/Suslin.jl src/algorithm/factorization.jl src/algorithm/polynomial_column_peel.jl test/expert/park_woodburn_polynomial_column_peel.jl test/runtests.jl docs/superpowers/plans/2026-06-23-issue-112-polynomial-column-peel.md
git commit -m "feat: add polynomial column peel certificates"
```

Expected: commit succeeds.

---

## Plan Self-Review

- The plan covers every issue 112 output field and negative control.
- The test file is created and run before implementation.
- The public `elementary_factorization(A)` route is not changed.
- The route-certificate extension is internal and reuses the new peel certificate as evidence.
- Laurent and Quillen support remain out of scope.
