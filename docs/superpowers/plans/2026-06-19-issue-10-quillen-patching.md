# Constructive Quillen Patching Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build an inspectable constructive Quillen patching layer with exact denominator coverage tracking and exact product verification.

**Architecture:** Extend the existing Quillen scaffolding in `src/core/groebner_tools.jl` with a small data model for denominator metadata, elementary local corrections, local contributions, global patches, and verification records. Patch construction remains explicit: callers supply coverage coefficients, the constructor verifies that the denominator-weighted coverage sum is one, assembles elementary factors, multiplies them exactly, and checks the product against the requested target matrix.

**Tech Stack:** Julia, Oscar.jl exact polynomial and Laurent polynomial rings, existing Suslin elementary matrix helpers, `Test`.

## Global Constraints

- No `AGENTS.md`, `CLAUDE.md`, or `GEMINI.md` file is present in this checkout.
- Issue #10 dependencies #8, #9, #12, and #21 are closed as completed on GitHub.
- Use the explicit coverage identity `sum(coverage_multiplier_i * denominator_i) == 1`; do not add internal Bezout solving.
- Keep patch assembly exact-ring generic for ordinary polynomial and Laurent polynomial parents.
- Keep polynomial-only substeps out of this patch layer; Laurent inputs must remain valid exact-ring inputs.
- The patch object must record denominators, local certificates, substitutions, factors, product, target, and exact verification data.
- A certificate set whose denominators do not provide coverage must throw `ArgumentError("denominator coverage must sum to one")`.
- A tampered denominator in an otherwise valid patch must make `verify_quillen_patch(patch)` return `false`.
- Required focused verification command: `julia --project=. -e 'include("test/expert/quillen_patching_exact.jl")'`.
- Required package verification command: `julia --project=. -e 'using Pkg; Pkg.test()'`.
- Documented full-suite command from #21: `julia --project=. test/runtests.jl all`.

---

## File Structure

- Modify `src/core/groebner_tools.jl`: promote `LocalCertificate` and add the Quillen patch data model, constructor helpers, exact factor assembly, and patch verification.
- Modify `src/Suslin.jl`: export the new patching records and helper functions, plus the existing Quillen scaffolding names that are now part of the public patching interface.
- Create `test/expert/quillen_patching_exact.jl`: focused red/green coverage for exact multi-certificate patching, exact multiplication, denominator stability, missing coverage, tampered denominator verification, certificate mismatch, substitution validation, and Laurent exact-ring support.
- Modify `test/runtests.jl`: include the new expert test in the expert group.
- Modify `test/public/api_surface.jl`: assert that the exported patching names are present and identical to the `Suslin.` bindings.

### Task 1: Constructive Quillen Patch Data Model

**Files:**
- Modify: `src/core/groebner_tools.jl`
- Modify: `src/Suslin.jl`
- Create: `test/expert/quillen_patching_exact.jl`
- Modify: `test/runtests.jl`
- Modify: `test/public/api_surface.jl`

**Interfaces:**
- Consumes:
  - `LocalCertificate(indices, denominators)`
  - `elementary_matrix(n, i, j, a, R)`
  - `_coerce_into_ring(R, value, label)`
  - `_is_laurent_polynomial_ring(R)`
  - `_require_square_matrix(M, label)`
- Produces:
  - `QuillenDenominatorData(denominator, coverage_multiplier)`
  - `QuillenElementaryCorrection(row::Int, col::Int, entry)`
  - `QuillenLocalContribution(certificate::LocalCertificate, denominator, coverage_multiplier, correction::QuillenElementaryCorrection)`
  - `QuillenPatchVerification(denominator_data_ok::Bool, coverage_sum, coverage_ok::Bool, product, target, product_ok::Bool)`
  - `QuillenPatch(ring, size::Int, substitution_variable, denominator_data::Vector{QuillenDenominatorData}, local_contributions::Vector{QuillenLocalContribution}, factors::Vector, product, target, verification::QuillenPatchVerification)`
  - `construct_quillen_patch(n::Int, X, contributions; target)::QuillenPatch`
  - `verify_quillen_patch(patch::QuillenPatch)::Bool`

- [ ] **Step 1: Write the failing expert test**

Create `test/expert/quillen_patching_exact.jl` with:

```julia
using Test
using Suslin
using Oscar

function quillen_patch_product(factors, R, n)
    product = identity_matrix(R, n)
    for factor in factors
        product *= factor
    end
    return product
end

@testset "constructive quillen patching exact" begin
    R, (X, r, g) = Oscar.polynomial_ring(QQ, ["X", "r", "g"])
    n = 3
    target_entry = X + g + 1
    target = elementary_matrix(n, 1, 2, target_entry, R)
    contributions = [
        QuillenLocalContribution(
            LocalCertificate([1, 2], [r, r]),
            r,
            one(R),
            QuillenElementaryCorrection(1, 2, target_entry),
        ),
        QuillenLocalContribution(
            LocalCertificate([1, 2], [one(R) - r, one(R) - r]),
            one(R) - r,
            one(R),
            QuillenElementaryCorrection(1, 2, target_entry),
        ),
    ]

    patch = construct_quillen_patch(n, X, contributions; target)

    @test verify_quillen_patch(patch)
    @test patch.ring == R
    @test patch.size == n
    @test patch.substitution_variable == X
    @test patch.product == target
    @test quillen_patch_product(patch.factors, R, n) == target
    @test patch.verification.product == target
    @test patch.verification.target == target
    @test patch.verification.coverage_sum == one(R)
    @test patch.verification.coverage_ok
    @test patch.verification.product_ok
    @test patch.verification.denominator_data_ok

    base_matrix = matrix(R, [
        one(R)  X       g;
        zero(R) one(R)  r;
        zero(R) zero(R) one(R)
    ])
    @test base_matrix * patch.product == base_matrix * target

    @test length(patch.denominator_data) == 2
    @test [data.denominator for data in patch.denominator_data] == [r, one(R) - r]
    @test [data.coverage_multiplier for data in patch.denominator_data] == [one(R), one(R)]
    @test [contribution.certificate.indices for contribution in patch.local_contributions] == [[1, 2], [1, 2]]
    @test [contribution.certificate.denominators for contribution in patch.local_contributions] == [
        [r, r],
        [one(R) - r, one(R) - r],
    ]

    uncovered_contributions = [
        contributions[1],
        QuillenLocalContribution(
            LocalCertificate([1, 2], [one(R) - r, one(R) - r]),
            one(R) - r,
            r,
            QuillenElementaryCorrection(1, 2, target_entry),
        ),
    ]
    @test_throws ArgumentError construct_quillen_patch(n, X, uncovered_contributions; target)

    tampered_denominator_data = copy(patch.denominator_data)
    tampered_denominator_data[2] = QuillenDenominatorData(one(R) - r + X, one(R))
    tampered_patch = QuillenPatch(
        patch.ring,
        patch.size,
        patch.substitution_variable,
        tampered_denominator_data,
        patch.local_contributions,
        patch.factors,
        patch.product,
        patch.target,
        patch.verification,
    )
    @test !verify_quillen_patch(tampered_patch)

    certificate_mismatch = [
        QuillenLocalContribution(
            LocalCertificate([1, 2], [r + 1, r + 1]),
            r,
            one(R),
            QuillenElementaryCorrection(1, 2, target_entry),
        ),
        contributions[2],
    ]
    @test_throws ArgumentError construct_quillen_patch(n, X, certificate_mismatch; target)
    @test_throws ArgumentError construct_quillen_patch(n, one(R), contributions; target)

    L, (x, u) = suslin_laurent_polynomial_ring(QQ, ["x", "u"])
    laurent_target_entry = x^-1 + u
    laurent_target = elementary_matrix(n, 2, 3, laurent_target_entry, L)
    laurent_contributions = [
        QuillenLocalContribution(
            LocalCertificate([2, 3], [u, u]),
            u,
            one(L),
            QuillenElementaryCorrection(2, 3, laurent_target_entry),
        ),
        QuillenLocalContribution(
            LocalCertificate([2, 3], [one(L) - u, one(L) - u]),
            one(L) - u,
            one(L),
            QuillenElementaryCorrection(2, 3, laurent_target_entry),
        ),
    ]
    laurent_patch = construct_quillen_patch(n, x, laurent_contributions; target = laurent_target)
    @test verify_quillen_patch(laurent_patch)
    @test laurent_patch.product == laurent_target
    @test quillen_patch_product(laurent_patch.factors, L, n) == laurent_target
    @test laurent_patch.verification.coverage_sum == one(L)
end
```

- [ ] **Step 2: Register public API and expert runner expectations**

In `test/public/api_surface.jl`, add these `isdefined` checks inside the existing testset, next to the other public API checks:

```julia
    @test isdefined(Suslin, :LocalCertificate)
    @test isdefined(Suslin, :common_denominator_factor)
    @test isdefined(Suslin, :patched_substitution)
    @test isdefined(Suslin, :QuillenDenominatorData)
    @test isdefined(Suslin, :QuillenElementaryCorrection)
    @test isdefined(Suslin, :QuillenLocalContribution)
    @test isdefined(Suslin, :QuillenPatchVerification)
    @test isdefined(Suslin, :QuillenPatch)
    @test isdefined(Suslin, :construct_quillen_patch)
    @test isdefined(Suslin, :verify_quillen_patch)
```

Add these identity checks after the existing identity checks:

```julia
    @test Suslin.LocalCertificate === LocalCertificate
    @test Suslin.common_denominator_factor === common_denominator_factor
    @test Suslin.patched_substitution === patched_substitution
    @test Suslin.QuillenDenominatorData === QuillenDenominatorData
    @test Suslin.QuillenElementaryCorrection === QuillenElementaryCorrection
    @test Suslin.QuillenLocalContribution === QuillenLocalContribution
    @test Suslin.QuillenPatchVerification === QuillenPatchVerification
    @test Suslin.QuillenPatch === QuillenPatch
    @test Suslin.construct_quillen_patch === construct_quillen_patch
    @test Suslin.verify_quillen_patch === verify_quillen_patch
```

In `test/runtests.jl`, add the new file to the expert group immediately after `expert/quillen_induction.jl`:

```julia
        "expert/quillen_patching_exact.jl",
```

- [ ] **Step 3: Run the focused expert test to verify RED**

Run:

```bash
julia --project=. -e 'include("test/expert/quillen_patching_exact.jl")'
```

Expected: FAIL before production changes with `UndefVarError` for `QuillenLocalContribution` or another new patching API name. That proves the test is exercising missing behavior.

- [ ] **Step 4: Implement the patching data model and constructor**

In `src/core/groebner_tools.jl`, replace the `LocalCertificate` definition with an abstract-vector friendly constructor and add the new structs and helpers after `_denominator_factor`:

```julia
# Backend-sensitive scaffolding for the later Quillen patching layer.
struct LocalCertificate
    indices::Vector{Int}
    denominators::Vector

    function LocalCertificate(indices::AbstractVector{<:Integer}, denominators::AbstractVector)
        Base.require_one_based_indexing(indices)
        Base.require_one_based_indexing(denominators)
        length(indices) == length(denominators) || throw(ArgumentError("indices and denominators must have the same length"))
        return new(Int.(collect(indices)), collect(denominators))
    end
end
```

Append:

```julia
struct QuillenDenominatorData
    denominator
    coverage_multiplier
end

struct QuillenElementaryCorrection
    row::Int
    col::Int
    entry
end

struct QuillenLocalContribution
    certificate::LocalCertificate
    denominator
    coverage_multiplier
    correction::QuillenElementaryCorrection
end

struct QuillenPatchVerification
    denominator_data_ok::Bool
    coverage_sum
    coverage_ok::Bool
    product
    target
    product_ok::Bool
end

struct QuillenPatch
    ring
    size::Int
    substitution_variable
    denominator_data::Vector{QuillenDenominatorData}
    local_contributions::Vector{QuillenLocalContribution}
    factors::Vector
    product
    target
    verification::QuillenPatchVerification
end

function _require_supported_quillen_ring(R)
    _is_laurent_polynomial_ring(R) && return R
    try
        gens(R)
        coefficient_ring(R)
        return R
    catch err
        err isa MethodError || err isa ErrorException || rethrow()
        throw(ArgumentError("target base ring must be a supported exact polynomial or Laurent polynomial ring"))
    end
end

function _require_quillen_target(target, n::Int)
    n >= 2 || throw(ArgumentError("patch size must be at least 2"))
    size = _require_square_matrix(target, "target correction")
    size == n || throw(DimensionMismatch("target correction size must match the requested patch size"))
    return _require_supported_quillen_ring(base_ring(target))
end

function _require_substitution_generator(R, X)
    coerced = _coerce_into_ring(R, X, "substitution variable")
    any(gen -> gen == coerced, collect(gens(R))) ||
        throw(ArgumentError("substitution variable must be a generator of the target ring"))
    return coerced
end

function _require_elementary_indices(n::Int, row::Int, col::Int)
    1 <= row <= n || throw(ArgumentError("elementary correction row must be between 1 and the patch size"))
    1 <= col <= n || throw(ArgumentError("elementary correction column must be between 1 and the patch size"))
    row != col || throw(ArgumentError("elementary correction row and column must differ"))
    return row, col
end

function _coerced_certificate_denominators(certificate::LocalCertificate, R)
    return [_coerce_into_ring(R, denominator, "certificate denominator") for denominator in certificate.denominators]
end

function _normalize_quillen_contribution(contribution::QuillenLocalContribution, R, n::Int)
    correction = contribution.correction
    row, col = _require_elementary_indices(n, correction.row, correction.col)
    row in contribution.certificate.indices && col in contribution.certificate.indices ||
        throw(ArgumentError("local certificate indices must include the correction row and column"))

    denominator = _coerce_into_ring(R, contribution.denominator, "denominator")
    certificate_denominators = _coerced_certificate_denominators(contribution.certificate, R)
    any(certificate_denominator -> certificate_denominator == denominator, certificate_denominators) ||
        throw(ArgumentError("local contribution denominator must be present in its certificate denominators"))

    coverage_multiplier = _coerce_into_ring(R, contribution.coverage_multiplier, "coverage multiplier")
    entry = _coerce_into_ring(R, correction.entry, "correction entry")
    normalized_correction = QuillenElementaryCorrection(row, col, entry)
    normalized_certificate = LocalCertificate(contribution.certificate.indices, certificate_denominators)
    return QuillenLocalContribution(normalized_certificate, denominator, coverage_multiplier, normalized_correction)
end

function _quillen_denominator_data(local_contributions)
    return [
        QuillenDenominatorData(contribution.denominator, contribution.coverage_multiplier)
        for contribution in local_contributions
    ]
end

function _quillen_coverage_sum(R, denominator_data)
    total = zero(R)
    for data in denominator_data
        total += data.coverage_multiplier * data.denominator
    end
    return total
end

function _quillen_factors(R, n::Int, local_contributions)
    factor_type = typeof(identity_matrix(R, n))
    factors = factor_type[]
    for contribution in local_contributions
        correction = contribution.correction
        weighted_entry = contribution.coverage_multiplier * contribution.denominator * correction.entry
        push!(factors, elementary_matrix(n, correction.row, correction.col, weighted_entry, R))
    end
    return factors
end

function _quillen_product(R, n::Int, factors)
    product = identity_matrix(R, n)
    for factor in factors
        product *= factor
    end
    return product
end

function _same_quillen_denominator_data(left, right)::Bool
    length(left) == length(right) || return false
    for idx in eachindex(left)
        left[idx] == right[idx] || return false
    end
    return true
end

function _quillen_patch_verification(R, n::Int, denominator_data, local_contributions, factors, product, target)
    expected_denominator_data = _quillen_denominator_data(local_contributions)
    denominator_data_ok = _same_quillen_denominator_data(denominator_data, expected_denominator_data)
    coverage_sum = _quillen_coverage_sum(R, denominator_data)
    coverage_ok = coverage_sum == one(R)
    actual_product = _quillen_product(R, n, factors)
    product_ok = actual_product == target && product == actual_product
    return QuillenPatchVerification(denominator_data_ok, coverage_sum, coverage_ok, actual_product, target, product_ok)
end

function construct_quillen_patch(n::Int, X, contributions; target)
    collected = collect(contributions)
    isempty(collected) && throw(ArgumentError("local contributions must be nonempty"))

    R = _require_quillen_target(target, n)
    substitution_variable = _require_substitution_generator(R, X)
    local_contributions = [
        _normalize_quillen_contribution(contribution, R, n)
        for contribution in collected
    ]
    denominator_data = _quillen_denominator_data(local_contributions)
    coverage_sum = _quillen_coverage_sum(R, denominator_data)
    coverage_sum == one(R) || throw(ArgumentError("denominator coverage must sum to one"))

    factors = _quillen_factors(R, n, local_contributions)
    product = _quillen_product(R, n, factors)
    product == target || throw(ArgumentError("constructed Quillen patch does not multiply to the target correction"))
    verification = _quillen_patch_verification(R, n, denominator_data, local_contributions, factors, product, target)
    return QuillenPatch(R, n, substitution_variable, denominator_data, local_contributions, factors, product, target, verification)
end

function verify_quillen_patch(patch::QuillenPatch)::Bool
    try
        verification = _quillen_patch_verification(
            patch.ring,
            patch.size,
            patch.denominator_data,
            patch.local_contributions,
            patch.factors,
            patch.product,
            patch.target,
        )
        return verification.denominator_data_ok && verification.coverage_ok && verification.product_ok
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end
```

- [ ] **Step 5: Export the patching API**

In `src/Suslin.jl`, add exports near the existing core/factor exports:

```julia
export LocalCertificate
export common_denominator_factor
export patched_substitution
export QuillenDenominatorData
export QuillenElementaryCorrection
export QuillenLocalContribution
export QuillenPatchVerification
export QuillenPatch
export construct_quillen_patch
export verify_quillen_patch
```

- [ ] **Step 6: Run focused tests to verify GREEN**

Run:

```bash
julia --project=. -e 'include("test/expert/quillen_patching_exact.jl")'
```

Expected: PASS, including the negative controls for missing denominator coverage and tampered denominator metadata.

- [ ] **Step 7: Run integration checks**

Run:

```bash
julia --project=. -e 'include("test/expert/quillen_induction.jl")'
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected:
- Existing Quillen induction scaffolding remains green.
- Package entry point passes public and internal groups.

- [ ] **Step 8: Commit**

Run:

```bash
git add src/core/groebner_tools.jl src/Suslin.jl test/expert/quillen_patching_exact.jl test/runtests.jl test/public/api_surface.jl
git commit -m "feat: add constructive quillen patching"
```

Expected: one feature commit containing the implementation and tests.

## Plan Self-Review

- Spec coverage: the task covers the data model, exact denominator coverage, patch assembly, verification data, exact product equality, negative controls, Laurent support, exports, focused tests, package test wiring, and full-suite compatibility.
- Placeholder scan: the plan contains no incomplete-work markers.
- Type consistency: all produced names match the design spec and the public API assertions.
- Execution choice: under the standing answer policy, choose **Subagent-Driven (recommended)** when prompted because it is marked recommended by the writing-plans workflow.
