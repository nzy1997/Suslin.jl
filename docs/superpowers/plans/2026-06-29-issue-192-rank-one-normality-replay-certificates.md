# Issue 192 Rank-One Normality Replay Certificates Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add replayable ordinary-polynomial certificates for Park-Woodburn orthogonal rank-one normality decompositions.

**Architecture:** `src/algorithm/normality.jl` owns rank-one input validation, coefficient-table replay, child certificate assembly, and full product verification. It delegates every Cohn-type child to `realize_cohn_type_certificate` from `src/algorithm/cohn_type.jl`, so this issue does not duplicate the Cohn elementary factor formula.

**Tech Stack:** Julia, Oscar/AbstractAlgebra ordinary polynomial matrices, existing Suslin elementary matrix helpers, #190 polynomial normality fixtures, #191 Cohn-type certificate API, `Test`.

## Global Constraints

- Input vectors `v`, `w`, and `g` must be one-based vectors of the same length `n >= 3` over an ordinary polynomial ring `R`.
- The constructor must reject inputs unless `w*v == 0` and `g*v == 1` after coercion into `R`.
- The coefficient table must be deterministic lexicographic `i < j` with `a_ij = w_i*g_j - w_j*g_i`, including zero coefficients.
- Child Cohn-type certificates must be produced by `realize_cohn_type_certificate`; do not duplicate the Cohn-type elementary factor formula.
- Child certificates and concatenated elementary factors skip zero `a_ij` entries, preserving nonzero table order.
- The certificate must store input vectors, orthogonality and Bezout replay values, coefficient table, child certificates, concatenated factor sequence, rank-one target `I + v*w`, reconstructed product, and verification metadata.
- `verify_rank_one_normality_certificate(cert)` must accept only exact replay, including stored verification metadata equality.
- Use at least one multivariate ordinary-polynomial fixture from `test/fixtures/polynomial_normality_cases.jl`.
- Negative controls must reject bad orthogonality, bad Bezout data, and a tampered child Cohn-type factor.
- Do not implement conjugated elementary normality certificates, Murthy, Quillen, ECP, Laurent, ToricBuilder, or Steinberg factor-count optimization.

---

## File Structure

- Modify `src/algorithm/normality.jl`: add `RankOneNormalityCertificate`, checked-data helpers, coefficient-table/target builders, constructor, core verifier, public verifier, and small factor-product replay wrapper.
- Modify `src/Suslin.jl`: export `RankOneNormalityCertificate`, `realize_rank_one_normality_certificate`, and `verify_rank_one_normality_certificate`.
- Create `test/expert/normality_rank_one.jl`: fixture-backed red-green tests and negative controls.
- Modify `test/runtests.jl`: add `expert/normality_rank_one.jl` to the expert group near `expert/normality.jl`.

### Task 1: Rank-One Certificate Tests

**Files:**
- Create: `test/expert/normality_rank_one.jl`

**Interfaces:**
- Consumes: `Main.ParkWoodburnPolynomialNormalityFixtureCatalog.cases_by_id()`, `Suslin.realize_rank_one_normality_certificate(v, w, g, R)`, `Suslin.verify_rank_one_normality_certificate(cert)`, and `Suslin.verify_cohn_type_certificate(child)`.
- Produces: A focused expert test file that fails before the rank-one certificate API exists.

- [ ] **Step 1: Write the failing focused test**

Create `test/expert/normality_rank_one.jl` with this content:

```julia
using Test
using Suslin
using Oscar

const PARK_WOODBURN_POLYNOMIAL_NORMALITY_CATALOG_PATH =
    joinpath(@__DIR__, "..", "fixtures", "polynomial_normality_cases.jl")

if !isdefined(Main, :ParkWoodburnPolynomialNormalityFixtureCatalog)
    include(PARK_WOODBURN_POLYNOMIAL_NORMALITY_CATALOG_PATH)
end

function rank_one_product_of_factors(factors, R, n::Int)
    product = identity_matrix(R, n)
    for factor in factors
        product *= factor
    end
    return product
end

function expected_rank_one_cohn_coefficients(v, w, g)
    n = length(v)
    return [
        (; i, j, a = w[i] * g[j] - w[j] * g[i])
        for i in 1:(n - 1) for j in (i + 1):n
    ]
end

function replace_rank_one_certificate_field(cert, field::Symbol, value)
    names = fieldnames(typeof(cert))
    idx = findfirst(==(field), names)
    idx === nothing && throw(ArgumentError("certificate field $(field) not found"))
    values = ntuple(k -> k == idx ? value : getfield(cert, names[k]), length(names))
    return typeof(cert)(values...)
end

@testset "orthogonal rank-one normality replay certificate" begin
    cases = Main.ParkWoodburnPolynomialNormalityFixtureCatalog.cases_by_id()
    fixture = cases["pw-section2-orthogonal-rank-one-qq"]
    R = fixture.ring.object
    inputs = fixture.inputs
    n = nrows(fixture.target_matrix)

    cert = Suslin.realize_rank_one_normality_certificate(inputs.v, inputs.w, inputs.g, R)

    @test cert isa Suslin.RankOneNormalityCertificate
    @test cert.n == n
    @test cert.v == [R(value) for value in inputs.v]
    @test cert.w == [R(value) for value in inputs.w]
    @test cert.g == [R(value) for value in inputs.g]
    @test cert.orthogonality == zero(R)
    @test cert.bezout == one(R)
    @test cert.cohn_coefficients == expected_rank_one_cohn_coefficients(cert.v, cert.w, cert.g)
    @test length(cert.child_certificates) == count(entry -> entry.a != zero(R), cert.cohn_coefficients)
    @test all(Suslin.verify_cohn_type_certificate, cert.child_certificates)
    @test cert.factors == reduce(vcat, [child.factors for child in cert.child_certificates]; init = Any[])
    @test cert.product == rank_one_product_of_factors(cert.factors, R, n)
    @test cert.target == fixture.target_matrix
    @test cert.product == cert.target
    @test cert.verification.target_matches_product_ok
    @test Suslin.verify_rank_one_normality_certificate(cert)

    bad_w = copy(inputs.w)
    bad_w[end] += one(R)
    @test_throws ArgumentError Suslin.realize_rank_one_normality_certificate(inputs.v, bad_w, inputs.g, R)

    bad_g = copy(inputs.g)
    bad_g[1] += one(R)
    @test_throws ArgumentError Suslin.realize_rank_one_normality_certificate(inputs.v, inputs.w, bad_g, R)

    tampered_child_cert = Suslin.realize_rank_one_normality_certificate(inputs.v, inputs.w, inputs.g, R)
    tampered_child_cert.child_certificates[1].factors[1][1, 1] += one(R)
    @test !Suslin.verify_rank_one_normality_certificate(tampered_child_cert)

    tampered_factor_cert = Suslin.realize_rank_one_normality_certificate(inputs.v, inputs.w, inputs.g, R)
    tampered_factor_cert.factors[1][1, 1] += one(R)
    @test !Suslin.verify_rank_one_normality_certificate(tampered_factor_cert)

    changed_table = copy(cert.cohn_coefficients)
    changed_table[1] = (; changed_table[1].i, changed_table[1].j, a = changed_table[1].a + one(R))
    @test !Suslin.verify_rank_one_normality_certificate(
        replace_rank_one_certificate_field(cert, :cohn_coefficients, changed_table),
    )
end
```

- [ ] **Step 2: Run the focused test to verify RED**

Run:

```bash
julia --project=. -e 'include("test/expert/normality_rank_one.jl")'
```

Expected: FAIL because `Suslin.realize_rank_one_normality_certificate` and `Suslin.RankOneNormalityCertificate` are not defined.

- [ ] **Step 3: Commit the red test**

Run:

```bash
git add test/expert/normality_rank_one.jl
git commit -m "test: cover rank-one normality certificates"
```

Expected: commit succeeds with only the new focused test file staged.

### Task 2: Rank-One Certificate API

**Files:**
- Modify: `src/algorithm/normality.jl`
- Modify: `src/Suslin.jl`
- Modify: `test/runtests.jl`
- Test: `test/expert/normality_rank_one.jl`

**Interfaces:**
- Consumes: `realize_cohn_type_certificate(n, i, j, a, v, R)`, `verify_cohn_type_certificate(cert)`, `_require_ordinary_polynomial_certificate_ring(R)`, `_cohn_type_factor_product(factors, R, n)`, `_coerce_into_ring(R, value, label)`, `_same_base_ring(left, right)`, `_dot(values_a, values_b, R)`.
- Produces: `RankOneNormalityCertificate`, `realize_rank_one_normality_certificate(v::AbstractVector, w::AbstractVector, g::AbstractVector, R)`, `verify_rank_one_normality_certificate(cert)::Bool`.

- [ ] **Step 1: Implement rank-one certificate replay in `src/algorithm/normality.jl`**

Add this code before `realize_conjugate_elementary`:

```julia
struct RankOneNormalityCertificate
    n::Int
    v::Vector
    w::Vector
    g::Vector
    ring
    orthogonality
    bezout
    cohn_coefficients::Vector
    child_certificates::Vector{CohnTypeRealizationCertificate}
    factors::Vector
    target
    product
    verification
end

function _rank_one_checked_data(v::AbstractVector, w::AbstractVector, g::AbstractVector, R)
    Base.require_one_based_indexing(v)
    Base.require_one_based_indexing(w)
    Base.require_one_based_indexing(g)
    n = length(v)
    n >= 3 || throw(ArgumentError("v, w, and g must have length at least 3"))
    length(w) == n || throw(ArgumentError("w must have the same length as v"))
    length(g) == n || throw(ArgumentError("g must have the same length as v"))
    _require_ordinary_polynomial_certificate_ring(R)
    coerced_v = [_coerce_into_ring(R, v[idx], "v[$idx]") for idx in 1:n]
    coerced_w = [_coerce_into_ring(R, w[idx], "w[$idx]") for idx in 1:n]
    coerced_g = [_coerce_into_ring(R, g[idx], "g[$idx]") for idx in 1:n]
    orthogonality = _dot(coerced_w, coerced_v, R)
    bezout = _dot(coerced_g, coerced_v, R)
    return (;
        n,
        v = coerced_v,
        w = coerced_w,
        g = coerced_g,
        ring = R,
        orthogonality,
        bezout,
    )
end

function _rank_one_target_from_checked_data(data)
    target = identity_matrix(data.ring, data.n)
    for row in 1:data.n, col in 1:data.n
        target[row, col] += data.v[row] * data.w[col]
    end
    return target
end

function _rank_one_cohn_coefficients_from_checked_data(data)
    return [
        (; i, j, a = data.w[i] * data.g[j] - data.w[j] * data.g[i])
        for i in 1:(data.n - 1) for j in (i + 1):data.n
    ]
end

function _rank_one_child_certificates_from_checked_data(data, coefficients)
    children = CohnTypeRealizationCertificate[]
    for entry in coefficients
        entry.a == zero(data.ring) && continue
        push!(children, realize_cohn_type_certificate(data.n, entry.i, entry.j, entry.a, data.v, data.ring))
    end
    return children
end

function _rank_one_factors_from_children(children)
    factors = Any[]
    for child in children
        append!(factors, child.factors)
    end
    return factors
end
```

Then add the constructor and verifier:

```julia
function realize_rank_one_normality_certificate(v::AbstractVector, w::AbstractVector, g::AbstractVector, R)
    data = _rank_one_checked_data(v, w, g, R)
    data.orthogonality == zero(R) || throw(ArgumentError("rank-one inputs must satisfy w*v == 0"))
    data.bezout == one(R) || throw(ArgumentError("rank-one inputs must satisfy g*v == 1"))
    coefficients = _rank_one_cohn_coefficients_from_checked_data(data)
    children = _rank_one_child_certificates_from_checked_data(data, coefficients)
    factors = _rank_one_factors_from_children(children)
    target = _rank_one_target_from_checked_data(data)
    product = _cohn_type_factor_product(factors, R, data.n)
    provisional = RankOneNormalityCertificate(
        data.n,
        data.v,
        data.w,
        data.g,
        data.ring,
        data.orthogonality,
        data.bezout,
        coefficients,
        children,
        factors,
        target,
        product,
        nothing,
    )
    verification = _rank_one_certificate_core_verification(provisional)
    verification.overall_core_ok ||
        error("internal rank-one normality certificate verification failed")
    return RankOneNormalityCertificate(
        provisional.n,
        provisional.v,
        provisional.w,
        provisional.g,
        provisional.ring,
        provisional.orthogonality,
        provisional.bezout,
        provisional.cohn_coefficients,
        provisional.child_certificates,
        provisional.factors,
        provisional.target,
        provisional.product,
        verification,
    )
end

function verify_rank_one_normality_certificate(cert)::Bool
    try
        return _rank_one_certificate_verification(cert).overall_ok
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end
```

Implement `_rank_one_certificate_core_verification(cert)` and
`_rank_one_certificate_verification(cert)` so fresh replay checks all Global
Constraints:

```julia
function _rank_one_certificate_core_verification(cert)
    data = _rank_one_checked_data(cert.v, cert.w, cert.g, cert.ring)
    orthogonality_ok = cert.orthogonality == zero(cert.ring) && cert.orthogonality == data.orthogonality
    bezout_ok = cert.bezout == one(cert.ring) && cert.bezout == data.bezout
    expected_coefficients = _rank_one_cohn_coefficients_from_checked_data(data)
    coefficient_table_ok = cert.cohn_coefficients == expected_coefficients
    expected_children = _rank_one_child_certificates_from_checked_data(data, expected_coefficients)
    child_count = length(cert.child_certificates)
    child_certificates_ok =
        child_count == length(expected_children) &&
        all(verify_cohn_type_certificate, cert.child_certificates) &&
        all(zip(cert.child_certificates, expected_children)) do (actual, expected)
            actual.n == expected.n &&
                actual.i == expected.i &&
                actual.j == expected.j &&
                actual.a == expected.a &&
                actual.v == expected.v &&
                _same_base_ring(actual.ring, expected.ring) &&
                actual.auxiliary_index == expected.auxiliary_index &&
                actual.target == expected.target &&
                actual.factors == expected.factors &&
                actual.product == expected.product &&
                actual.verification == expected.verification
        end
    expected_factors = _rank_one_factors_from_children(cert.child_certificates)
    factor_sequence_ok = cert.factors == expected_factors
    replayed_product = _cohn_type_factor_product(cert.factors, cert.ring, cert.n)
    product_replay_ok = true
    product_matches_stored_ok = replayed_product == cert.product
    expected_target = _rank_one_target_from_checked_data(data)
    target_replay_ok = cert.target == expected_target
    target_matches_product_ok = target_replay_ok && replayed_product == cert.target
    overall_core_ok =
        orthogonality_ok &&
        bezout_ok &&
        coefficient_table_ok &&
        child_certificates_ok &&
        factor_sequence_ok &&
        product_replay_ok &&
        product_matches_stored_ok &&
        target_matches_product_ok
    return (;
        ring_ok = true,
        checked_inputs_ok = true,
        orthogonality_ok,
        bezout_ok,
        coefficient_table_ok,
        child_count,
        child_certificates_ok,
        factor_sequence_ok,
        factor_count = length(cert.factors),
        product_replay_ok,
        product_matches_stored_ok,
        target_replay_ok,
        target_matches_product_ok,
        overall_core_ok,
    )
end

function _rank_one_certificate_verification(cert)
    core = _rank_one_certificate_core_verification(cert)
    stored_verification_ok = cert.verification == core
    return merge(core, (;
        stored_verification_ok,
        overall_ok = core.overall_core_ok && stored_verification_ok,
    ))
end
```

- [ ] **Step 2: Export the new public API**

In `src/Suslin.jl`, add:

```julia
export RankOneNormalityCertificate
export realize_rank_one_normality_certificate
export verify_rank_one_normality_certificate
```

- [ ] **Step 3: Register the expert test**

In `test/runtests.jl`, add this entry immediately after `"expert/normality.jl"`:

```julia
"expert/normality_rank_one.jl",
```

- [ ] **Step 4: Run the focused test to verify GREEN**

Run:

```bash
julia --project=. -e 'include("test/expert/normality_rank_one.jl")'
```

Expected: PASS with the focused rank-one certificate testset.

- [ ] **Step 5: Commit the implementation**

Run:

```bash
git add src/algorithm/normality.jl src/Suslin.jl test/runtests.jl
git commit -m "feat: add rank-one normality certificates"
```

Expected: commit succeeds with the implementation and test registration.

## Self-Review

- Spec coverage: the two tasks cover test-first focused coverage, constructor rejection of bad orthogonality/Bezout data, exact replay verification, child Cohn certificate reuse, deterministic coefficient metadata, exports, and expert test registration.
- Placeholder scan: no task contains placeholder markers or an unspecified implementation step.
- Type consistency: the public type and function names match between tests, implementation, exports, and verification commands.
