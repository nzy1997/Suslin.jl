# Issue 191 Cohn-Type Replay Certificates Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a replayable Cohn-type certificate API for ordinary polynomial rings while preserving the existing factor-only helper.

**Architecture:** Keep Cohn-type factor generation in `src/algorithm/cohn_type.jl` as the primitive. Add a concrete certificate struct plus constructor and verifier that recompute the target and factor product from stored data and compare them with stored replay metadata.

**Tech Stack:** Julia, Oscar/AbstractAlgebra polynomial matrices, existing Suslin elementary matrix helpers, `Test`.

## Global Constraints

- Input is `n >= 3`, distinct one-based `i` and `j`, coefficient `a`, one-based vector `v` of length `n`, and an ordinary polynomial ring `R`.
- Output certificate stores coerced input data, target matrix `I + a*v*(v_j*e_i - v_i*e_j)`, elementary factor sequence, reconstructed product, and verification metadata.
- `verify_cohn_type_certificate(cert)` returns `true` only when stored factors exactly replay the stored and recomputed target.
- Reject Laurent rings in the certificate route.
- Keep `realize_cohn_type` available for existing callers.
- Use at least one fixture-backed Cohn-type case from `test/fixtures/polynomial_normality_cases.jl`.
- Negative controls must reject a tampered factor, tampered target entry, or changed stored `a`.
- Do not implement rank-one normality, conjugated elementary normality, Murthy, Quillen, ECP, Laurent, ToricBuilder, or Steinberg factor-count optimization.

---

## File Structure

- Modify `src/algorithm/cohn_type.jl`: certificate struct, input normalization helpers, target/product replay helpers, constructor, verifier, and a factor-only wrapper path.
- Modify `src/Suslin.jl`: export `CohnTypeRealizationCertificate`, `realize_cohn_type_certificate`, and `verify_cohn_type_certificate`.
- Modify `test/expert/cohn_type.jl`: write red-green certificate tests and negative controls.
- Modify `test/expert/polynomial_normality_fixtures.jl`: guard fixture module inclusion so expert tests can include the same fixture file once per Julia session without module replacement warnings.

### Task 1: Cohn-Type Certificate API

**Files:**
- Modify: `test/expert/cohn_type.jl`
- Modify: `test/expert/polynomial_normality_fixtures.jl`
- Modify: `src/algorithm/cohn_type.jl`
- Modify: `src/Suslin.jl`

**Interfaces:**
- Consumes: `elementary_matrix(n, i, j, a, R)`, `_coerce_into_ring(R, value, label)`, `_is_laurent_polynomial_ring(R)`, `_same_base_ring(left, right)`.
- Produces: `CohnTypeRealizationCertificate`, `realize_cohn_type_certificate(n::Int, i::Int, j::Int, a, v::AbstractVector, R)`, `verify_cohn_type_certificate(cert)::Bool`.

- [ ] **Step 1: Write the failing certificate tests**

Add fixture loading and helpers to `test/expert/cohn_type.jl`:

```julia
const PARK_WOODBURN_POLYNOMIAL_NORMALITY_CATALOG_PATH =
    joinpath(@__DIR__, "..", "fixtures", "polynomial_normality_cases.jl")

if !isdefined(Main, :ParkWoodburnPolynomialNormalityFixtureCatalog)
    include(PARK_WOODBURN_POLYNOMIAL_NORMALITY_CATALOG_PATH)
end

function replace_cohn_certificate_field(cert, field::Symbol, value)
    names = fieldnames(typeof(cert))
    idx = findfirst(==(field), names)
    idx === nothing && throw(ArgumentError("certificate field $(field) not found"))
    values = ntuple(k -> k == idx ? value : getfield(cert, names[k]), length(names))
    return typeof(cert)(values...)
end
```

Add a certificate testset that exercises the new public API:

```julia
@testset "cohn-type replay certificate" begin
    cases = Main.ParkWoodburnPolynomialNormalityFixtureCatalog.cases_by_id()
    fixture = cases["pw-section2-cohn-type-qq"]
    R = fixture.ring.object
    inputs = fixture.inputs
    n = nrows(fixture.target_matrix)

    cert = Suslin.realize_cohn_type_certificate(n, inputs.i, inputs.j, inputs.a, inputs.v, R)

    @test cert isa Suslin.CohnTypeRealizationCertificate
    @test cert.n == n
    @test cert.i == inputs.i
    @test cert.j == inputs.j
    @test cert.a == R(inputs.a)
    @test cert.v == [R(value) for value in inputs.v]
    @test cert.target == fixture.target_matrix
    @test cert.product == cert.target
    @test cert.verification.target_matches_product_ok
    @test Suslin.verify_cohn_type_certificate(cert)
    @test Suslin.realize_cohn_type(n, inputs.i, inputs.j, inputs.a, inputs.v, R) == cert.factors

    tampered_factor_cert = Suslin.realize_cohn_type_certificate(n, inputs.i, inputs.j, inputs.a, inputs.v, R)
    tampered_factor_cert.factors[1][1, 1] += one(R)
    @test !Suslin.verify_cohn_type_certificate(tampered_factor_cert)

    tampered_target_cert = Suslin.realize_cohn_type_certificate(n, inputs.i, inputs.j, inputs.a, inputs.v, R)
    tampered_target_cert.target[1, 1] += one(R)
    @test !Suslin.verify_cohn_type_certificate(tampered_target_cert)

    changed_a_cert = replace_cohn_certificate_field(cert, :a, cert.a + one(R))
    @test !Suslin.verify_cohn_type_certificate(changed_a_cert)

    L, (lx, ly) = Oscar.laurent_polynomial_ring(QQ, ["lx", "ly"])
    @test_throws ArgumentError Suslin.realize_cohn_type_certificate(
        3,
        1,
        2,
        lx,
        [one(L), lx, ly],
        L,
    )
end
```

Guard fixture inclusion in `test/expert/polynomial_normality_fixtures.jl`:

```julia
if !isdefined(Main, :ParkWoodburnPolynomialNormalityFixtureCatalog)
    include(PARK_WOODBURN_POLYNOMIAL_NORMALITY_CATALOG_PATH)
end
catalog = Main.ParkWoodburnPolynomialNormalityFixtureCatalog.catalog()
```

- [ ] **Step 2: Run the focused test to verify RED**

Run:

```bash
julia --project=. -e 'include("test/expert/cohn_type.jl")'
```

Expected: FAIL because `Suslin.realize_cohn_type_certificate` and `Suslin.CohnTypeRealizationCertificate` are not defined yet.

- [ ] **Step 3: Implement the certificate API**

In `src/algorithm/cohn_type.jl`, add the concrete certificate and helpers before `realize_cohn_type`:

```julia
struct CohnTypeRealizationCertificate
    n::Int
    i::Int
    j::Int
    a
    v::Vector
    ring
    auxiliary_index::Int
    target
    factors::Vector
    product
    verification
end
```

Add helper functions with these exact responsibilities:

```julia
function _require_ordinary_polynomial_certificate_ring(R)
    _is_laurent_polynomial_ring(R) &&
        throw(ArgumentError("Cohn-type certificates require an ordinary polynomial ring"))
    try
        collect(gens(R))
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("Cohn-type certificates require an ordinary polynomial ring"))
    end
    return R
end

function _cohn_type_checked_data(n::Int, i::Int, j::Int, a, v::AbstractVector, R)
    n >= 3 || throw(ArgumentError("n must be at least 3"))
    1 <= i <= n || throw(ArgumentError("i must be between 1 and n"))
    1 <= j <= n || throw(ArgumentError("j must be between 1 and n"))
    i == j && throw(ArgumentError("i and j must differ"))
    Base.require_one_based_indexing(v)
    length(v) == n || throw(ArgumentError("v must contain exactly n entries"))
    t = findfirst(k -> k != i && k != j, 1:n)
    t === nothing && throw(ArgumentError("could not choose an auxiliary index distinct from i and j"))
    coerced_v = [_coerce_into_ring(R, v[idx], "v[$idx]") for idx in 1:n]
    coerced_a = _coerce_into_ring(R, a, "a")
    return (; n, i, j, a = coerced_a, v = coerced_v, ring = R, auxiliary_index = t)
end
```

Move the existing factor construction into `_cohn_type_factors_from_checked_data(data)` and make `realize_cohn_type` return that helper's factors. Add:

```julia
function _cohn_type_target_from_checked_data(data)
    target = identity_matrix(data.ring, data.n)
    vi = data.v[data.i]
    vj = data.v[data.j]
    for row in 1:data.n
        target[row, data.i] += data.a * data.v[row] * vj
        target[row, data.j] -= data.a * data.v[row] * vi
    end
    return target
end

function _cohn_type_factor_product(factors, R, n::Int)
    product = identity_matrix(R, n)
    for factor in factors
        nrows(factor) == n || throw(ArgumentError("Cohn-type factor has wrong row count"))
        ncols(factor) == n || throw(ArgumentError("Cohn-type factor has wrong column count"))
        _same_base_ring(base_ring(factor), R) ||
            throw(ArgumentError("Cohn-type factor has wrong base ring"))
        product *= factor
    end
    return product
end
```

Implement `realize_cohn_type_certificate` by requiring the ordinary polynomial
ring, building checked data, factors, target, product, and a verification pass.
Implement `verify_cohn_type_certificate(cert)::Bool` by recomputing fresh
verification and returning `verification.overall_ok`, catching non-interrupt
errors as `false`.

- [ ] **Step 4: Export the public API**

In `src/Suslin.jl`, add:

```julia
export CohnTypeRealizationCertificate
export realize_cohn_type_certificate
export verify_cohn_type_certificate
```

- [ ] **Step 5: Run the focused test to verify GREEN**

Run:

```bash
julia --project=. -e 'include("test/expert/cohn_type.jl")'
```

Expected: PASS. The output should include both `cohn-type realization` and `cohn-type replay certificate` testsets.

- [ ] **Step 6: Run required package verification**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: PASS for the default package test entry point.

- [ ] **Step 7: Commit the implementation**

Run:

```bash
git add src/Suslin.jl src/algorithm/cohn_type.jl test/expert/cohn_type.jl test/expert/polynomial_normality_fixtures.jl docs/superpowers/plans/2026-06-29-issue-191-cohn-type-replay-certificates.md
git commit -m "feat: add cohn-type replay certificates"
```

Expected: a commit containing the plan, source changes, and tests.

## Plan Self-Review

- Spec coverage: the task covers the constructor, verifier, exported API,
  fixture-backed test, negative controls, ordinary-polynomial boundary, and
  preservation of `realize_cohn_type`.
- Placeholder scan: no deferred markers or unspecified "add tests" steps remain.
- Type consistency: produced names match the design and exports exactly:
  `CohnTypeRealizationCertificate`, `realize_cohn_type_certificate`, and
  `verify_cohn_type_certificate`.
