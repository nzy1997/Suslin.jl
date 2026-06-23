# Issue 114 Quillen Patch Route Adapter Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an internal route adapter that consumes verified Quillen patch output and embeds it as replayable ordinary-polynomial factorization route evidence.

**Architecture:** Extend the existing `PolynomialFactorizationRouteCertificate` model with an explicit `:quillen_patch` route. The new adapter verifies the supplied Quillen patch, copies its global elementary factors, recomputes the product against the requested target matrix, and is replayed by the route certificate verifier.

**Tech Stack:** Julia, Oscar exact polynomial matrices, existing Suslin Quillen patch verifiers, Test stdlib.

## Global Constraints

- No `AGENTS.md`, `CLAUDE.md`, `GEMINI.md`, or `CODEX.md` file is present in this checkout.
- The worker branch is `agent/issue-114-add-a-quillen-patch-route-adapter-for-polynomial-run-1`.
- Reuse #99 fixture ids and #105 constructive helpers; do not invent local-cover data.
- The adapter must call `verify_quillen_patch(patch)` before accepting factors.
- Keep the adapter internal/expert; do not export new names from `src/Suslin.jl`.
- Do not implement Quillen patch internals.
- Do not solve arbitrary local realizability.
- Do not wire public `elementary_factorization(A)` to multivariate inputs yet.
- Automatic polynomial route selection must not choose Quillen without an explicit patch.
- Focused verification command: `julia --project=. -e 'include("test/expert/park_woodburn_quillen_route_adapter.jl")'`.
- Required package verification command: `julia --project=. -e 'using Pkg; Pkg.test()'`.

---

## File Structure

- Modify `src/algorithm/factorization.jl`: define `PolynomialQuillenPatchRouteAdapter`, `_polynomial_quillen_patch_route_adapter`, `_verify_polynomial_quillen_patch_route_adapter`, and the `:quillen_patch` route-certificate integration.
- Create `test/expert/park_woodburn_quillen_route_adapter.jl`: focused expert test using the #109 Park-Woodburn Quillen fixture and the #105 constructive patch helper path.
- Modify `test/runtests.jl`: register the new expert test immediately after `expert/quillen_induction_constructive.jl`.
- Keep `src/Suslin.jl` and `test/public/api_surface.jl` unchanged.

---

### Task 1: Add RED Coverage for Quillen Patch Route Adaptation

**Files:**
- Create: `test/expert/park_woodburn_quillen_route_adapter.jl`
- Modify: `test/runtests.jl`

**Interfaces:**
- Consumes existing #105 helpers from `test/expert/quillen_induction_constructive.jl`, especially `constructive_patch`.
- Consumes planned internal helpers:
  - `Suslin.PolynomialQuillenPatchRouteAdapter`
  - `Suslin._polynomial_quillen_patch_route_adapter(target, patch)`
  - `Suslin._verify_polynomial_quillen_patch_route_adapter(adapter)::Bool`
  - `Suslin._polynomial_factorization_route_certificate(A; route=:quillen_patch, quillen_patch=patch)`
- Produces failing focused coverage before implementation.

- [ ] **Step 1: Add the focused expert test**

Create `test/expert/park_woodburn_quillen_route_adapter.jl`:

```julia
using Test
using Suslin
using Oscar

const PW_QUILLEN_ROUTE_CATALOG_PATH =
    joinpath(@__DIR__, "..", "fixtures", "park_woodburn_polynomial_cases.jl")

if !isdefined(Main, :ParkWoodburnPolynomialFixtureCatalog)
    include(PW_QUILLEN_ROUTE_CATALOG_PATH)
end
if !isdefined(Main, :constructive_patch)
    include(joinpath(@__DIR__, "quillen_induction_constructive.jl"))
end

function _pwq_product(factors, R, n::Int)
    product = identity_matrix(R, n)
    for factor in factors
        product *= factor
    end
    return product
end

function _pwq_rebuild(record; kwargs...)
    overrides = Dict{Symbol,Any}()
    for pair in kwargs
        overrides[pair.first] = pair.second
    end
    values = map(fieldnames(typeof(record))) do name
        get(overrides, name, getproperty(record, name))
    end
    return typeof(record)(values...)
end

function _pwq_replace_adapter(
        adapter;
        target = adapter.target,
        route = adapter.route,
        quillen_patch = adapter.quillen_patch,
        global_elementary_factors = adapter.global_elementary_factors,
        product = adapter.product,
        target_matrix = adapter.target_matrix,
        replay_metadata = adapter.replay_metadata,
        verification = adapter.verification)
    return Suslin.PolynomialQuillenPatchRouteAdapter(
        target,
        route,
        quillen_patch,
        global_elementary_factors,
        product,
        target_matrix,
        replay_metadata,
        verification,
    )
end

function _pwq_replace_route_certificate(
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

function _pwq_adapter_accepts(target, patch)::Bool
    try
        adapter = Suslin._polynomial_quillen_patch_route_adapter(target, patch)
        return Suslin._verify_polynomial_quillen_patch_route_adapter(adapter)
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end

function _pwq_route_accepts(target, patch)::Bool
    try
        cert = Suslin._polynomial_factorization_route_certificate(
            target;
            route = :quillen_patch,
            quillen_patch = patch,
        )
        return Suslin._verify_polynomial_factorization_route_certificate(cert)
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end

@testset "Park-Woodburn Quillen patch route adapter" begin
    route_entries = ParkWoodburnPolynomialFixtureCatalog.cases_by_id()
    route_entry = route_entries["quillen-patched-substitution-witness-qq"]
    quillen_entries = Main.QuillenPatchFixtureCatalog.cases_by_id()
    quillen_entry = quillen_entries[route_entry.provenance.quillen_fixture_id]

    _, _, _, patch = constructive_patch(quillen_entry)
    @test Suslin.verify_quillen_patch(patch)
    @test patch.target == route_entry.matrix

    adapter = Suslin._polynomial_quillen_patch_route_adapter(route_entry.matrix, patch)
    @test adapter isa Suslin.PolynomialQuillenPatchRouteAdapter
    @test adapter.route == :quillen_patch
    @test adapter.target_matrix == route_entry.matrix
    @test adapter.product == route_entry.matrix
    @test adapter.global_elementary_factors == patch.global_elementary_factors
    @test _pwq_product(adapter.global_elementary_factors, base_ring(route_entry.matrix), nrows(route_entry.matrix)) == route_entry.matrix
    @test Suslin._verify_polynomial_quillen_patch_route_adapter(adapter)

    cert = Suslin._polynomial_factorization_route_certificate(
        route_entry.matrix;
        route = :quillen_patch,
        quillen_patch = patch,
    )
    @test cert.route == :quillen_patch
    @test cert.status == :supported
    @test cert.evidence isa Suslin.PolynomialQuillenPatchRouteAdapter
    @test cert.factors == adapter.global_elementary_factors
    @test cert.product == route_entry.matrix
    @test verify_factorization(route_entry.matrix, cert.factors)
    @test Suslin._verify_polynomial_factorization_route_certificate(cert)

    bad_cover = Suslin.QuillenDenominatorCoverCertificate(
        patch.cover_certificate.ring,
        patch.cover_certificate.denominators,
        [
            patch.cover_certificate.coverage_multipliers[1] + one(patch.ring),
            patch.cover_certificate.coverage_multipliers[2],
        ],
        patch.cover_certificate.coverage_sum + patch.cover_certificate.denominators[1],
        patch.cover_certificate.verification,
    )
    tampered_cover_patch = _pwq_rebuild(patch; cover_certificate = bad_cover)
    @test !Suslin.verify_quillen_patch(tampered_cover_patch)
    @test !_pwq_adapter_accepts(route_entry.matrix, tampered_cover_patch)
    @test !_pwq_route_accepts(route_entry.matrix, tampered_cover_patch)

    tampered_local_certificates = copy(patch.local_certificates)
    tampered_factors = copy(tampered_local_certificates[1].factors)
    tampered_factors[1] =
        tampered_factors[1] *
        elementary_matrix(patch.size, 1, 3, one(patch.ring), patch.ring)
    tampered_local_certificates[1] = _pwq_rebuild(
        tampered_local_certificates[1];
        factors = tampered_factors,
    )
    tampered_local_patch = _pwq_rebuild(
        patch;
        local_certificates = tampered_local_certificates,
    )
    @test !Suslin.verify_quillen_patch(tampered_local_patch)
    @test !_pwq_adapter_accepts(route_entry.matrix, tampered_local_patch)
    @test !_pwq_route_accepts(route_entry.matrix, tampered_local_patch)

    overwritten_patch = _pwq_rebuild(
        tampered_local_patch;
        patched_product = route_entry.matrix,
        target = route_entry.matrix,
    )
    @test overwritten_patch.patched_product == route_entry.matrix
    @test overwritten_patch.target == route_entry.matrix
    @test !Suslin.verify_quillen_patch(overwritten_patch)
    @test !_pwq_adapter_accepts(route_entry.matrix, overwritten_patch)

    bad_adapter_factors = copy(adapter.global_elementary_factors)
    bad_adapter_factors[1] =
        bad_adapter_factors[1] *
        elementary_matrix(patch.size, 1, 3, one(patch.ring), patch.ring)
    bad_adapter = _pwq_replace_adapter(
        adapter;
        global_elementary_factors = bad_adapter_factors,
    )
    @test !Suslin._verify_polynomial_quillen_patch_route_adapter(bad_adapter)

    bad_route_cert = _pwq_replace_route_certificate(
        cert;
        factors = bad_adapter_factors,
        evidence = bad_adapter,
    )
    @test !Suslin._verify_polynomial_factorization_route_certificate(bad_route_cert)

    @test_throws ArgumentError Suslin._polynomial_factorization_route_certificate(
        route_entry.matrix;
        route = :quillen_patch,
    )
end
```

- [ ] **Step 2: Register the expert test**

In `test/runtests.jl`, add this entry immediately after
`"expert/quillen_induction_constructive.jl",`:

```julia
        "expert/park_woodburn_quillen_route_adapter.jl",
```

- [ ] **Step 3: Run focused RED verification**

Run:

```bash
julia --project=. -e 'include("test/expert/park_woodburn_quillen_route_adapter.jl")'
```

Expected: FAIL because `Suslin._polynomial_quillen_patch_route_adapter` and
`Suslin.PolynomialQuillenPatchRouteAdapter` do not exist yet.

---

### Task 2: Implement the Internal Quillen Patch Route Adapter

**Files:**
- Modify: `src/algorithm/factorization.jl`

**Interfaces:**
- Produces `PolynomialQuillenPatchRouteAdapter`.
- Produces `_polynomial_quillen_patch_route_adapter(target, patch)`.
- Produces `_verify_polynomial_quillen_patch_route_adapter(adapter)::Bool`.
- Extends `_polynomial_factorization_route_certificate(A; route=nothing, quillen_patch=nothing)`.

- [ ] **Step 1: Add the adapter type and route tag**

Near `PolynomialFactorizationRouteCertificate`, add:

```julia
struct PolynomialQuillenPatchRouteAdapter
    target
    route::Symbol
    quillen_patch
    global_elementary_factors::Vector
    product
    target_matrix
    replay_metadata
    verification
end
```

Add `:quillen_patch` to `_POLYNOMIAL_FACTORIZATION_ROUTE_TAGS`.

- [ ] **Step 2: Extend the route certificate constructor**

Change the constructor signature to:

```julia
function _polynomial_factorization_route_certificate(A; route=nothing, quillen_patch=nothing)
```

Inside the automatic `route === nothing` branch, reject accidental patch
arguments without an explicit route:

```julia
        quillen_patch === nothing ||
            throw(ArgumentError("Quillen patch route certificates require route = :quillen_patch"))
```

Add an explicit route branch before the final unsupported-route error:

```julia
    elseif route == :quillen_patch
        quillen_patch !== nothing ||
            throw(ArgumentError("Quillen patch route requires a supplied verified patch"))
        return _polynomial_quillen_patch_route_certificate(A, quillen_patch)
    end
```

- [ ] **Step 3: Add adapter construction and replay helpers**

Add these helpers before `_polynomial_route_evidence_ok`:

```julia
function _polynomial_quillen_patch_route_certificate(A, patch)
    adapter = _polynomial_quillen_patch_route_adapter(A, patch)
    factors = copy(adapter.global_elementary_factors)
    return _polynomial_route_certificate(
        adapter.target_matrix,
        :quillen_patch,
        factors,
        adapter.product,
        adapter,
        :supported,
    )
end

function _polynomial_quillen_patch_route_adapter(target, patch)
    verify_quillen_patch(patch) ||
        throw(ArgumentError("Quillen patch must verify before route adaptation"))
    target_matrix = _polynomial_quillen_route_target_matrix(target, patch)
    patch_target = _polynomial_quillen_patch_target_matrix(patch)
    patch_target == target_matrix ||
        throw(ArgumentError("Quillen patch target does not match route target"))
    factors = copy(collect(_polynomial_quillen_patch_factors(patch)))
    product = _polynomial_route_factor_product(
        factors,
        base_ring(target_matrix),
        nrows(target_matrix),
    )
    product == target_matrix ||
        throw(ArgumentError("Quillen patch route factors do not multiply to the target"))
    _polynomial_quillen_patch_product(patch) == product ||
        throw(ArgumentError("Quillen patch stored product does not match the adapted product"))
    replay_metadata = _polynomial_quillen_patch_route_metadata(patch)
    raw = PolynomialQuillenPatchRouteAdapter(
        target,
        :quillen_patch,
        patch,
        factors,
        product,
        target_matrix,
        replay_metadata,
        nothing,
    )
    verification = _polynomial_quillen_patch_route_core_verification(raw)
    adapter = PolynomialQuillenPatchRouteAdapter(
        target,
        :quillen_patch,
        patch,
        factors,
        product,
        target_matrix,
        replay_metadata,
        verification,
    )
    _verify_polynomial_quillen_patch_route_adapter(adapter) ||
        error("internal Quillen patch route adapter verification failed")
    return adapter
end
```

Add field accessors that support both `QuillenGlobalPatchAssembly` and older
`QuillenPatch` objects:

```julia
_polynomial_quillen_patch_factors(patch) =
    hasproperty(patch, :global_elementary_factors) ? patch.global_elementary_factors :
    hasproperty(patch, :factors) ? patch.factors :
    throw(ArgumentError("Quillen patch has no global elementary factors"))

_polynomial_quillen_patch_product(patch) =
    hasproperty(patch, :patched_product) ? patch.patched_product :
    hasproperty(patch, :product) ? patch.product :
    throw(ArgumentError("Quillen patch has no patched product"))

function _polynomial_quillen_patch_target_matrix(patch)
    R = _require_supported_quillen_ring(patch.ring)
    return _quillen_global_target_matrix(
        patch.target;
        ring = R,
        size = patch.size,
        label = "Quillen patch target",
    )
end
```

Add `_polynomial_quillen_route_target_matrix(target, patch)` that accepts a
matrix target or `QuillenElementaryCorrection`, checks ordinary-polynomial base
ring support, exact size equality with `patch.size`, and determinant one.

Add `_polynomial_quillen_patch_route_metadata(patch)` returning:

```julia
(;
    patch_size = patch.size,
    substitution_variable = patch.substitution_variable,
    denominator_data = copy(collect(patch.denominator_data)),
    local_certificate_count = hasproperty(patch, :local_certificates) ?
        length(patch.local_certificates) :
        length(patch.local_contributions),
    normalized_contribution_count = hasproperty(patch, :normalized_local_contributions) ?
        length(patch.normalized_local_contributions) :
        0,
    patch_replay_metadata = hasproperty(patch, :replay_metadata) ?
        patch.replay_metadata :
        nothing,
)
```

- [ ] **Step 4: Add fail-closed adapter verification**

Add `_polynomial_quillen_patch_route_core_verification(adapter)` that recomputes
the adapter fields from `adapter.quillen_patch` and returns a named tuple with:

```julia
(;
    route_tag_ok,
    patch_verified_ok,
    target_matrix_ok,
    patch_target_ok,
    factors_ok,
    product_ok,
    replay_metadata_ok,
    overall_core_ok,
)
```

`overall_core_ok` must require:

- `adapter.route == :quillen_patch`,
- `verify_quillen_patch(adapter.quillen_patch)`,
- the recomputed target matrix equals `adapter.target_matrix`,
- the patch target equals `adapter.target_matrix`,
- factors copied from the patch equal `adapter.global_elementary_factors`,
- the recomputed product equals both `adapter.product` and `adapter.target_matrix`,
- the patch stored product equals the recomputed product, and
- replay metadata matches.

Add:

```julia
function _verify_polynomial_quillen_patch_route_adapter(adapter)::Bool
    try
        verification = _polynomial_quillen_patch_route_verification(adapter)
        return verification.overall_ok
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end
```

where `_polynomial_quillen_patch_route_verification(adapter)` merges the core
record with `stored_verification_ok = adapter.verification == core` and
`overall_ok = core.overall_core_ok && stored_verification_ok`.

- [ ] **Step 5: Extend route evidence replay**

In `_polynomial_route_evidence_ok(cert)`, add:

```julia
        elseif cert.route == :quillen_patch
            return cert.evidence isa PolynomialQuillenPatchRouteAdapter &&
                cert.evidence.target_matrix == cert.matrix &&
                _verify_polynomial_quillen_patch_route_adapter(cert.evidence) &&
                _polynomial_route_factor_sequences_equal(
                    cert.factors,
                    cert.evidence.global_elementary_factors,
                )
```

- [ ] **Step 6: Run focused GREEN verification**

Run:

```bash
julia --project=. -e 'include("test/expert/park_woodburn_quillen_route_adapter.jl")'
```

Expected: PASS.

---

### Task 3: Final Verification, Commit, and PR Readiness

**Files:**
- Modify: implementation and tests from Tasks 1-2.

**Interfaces:**
- Produces a committed branch ready to push and open as a PR.

- [ ] **Step 1: Run focused issue verification**

Run:

```bash
julia --project=. -e 'include("test/expert/park_woodburn_quillen_route_adapter.jl")'
```

Expected: PASS.

- [ ] **Step 2: Run required package verification**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: PASS.

- [ ] **Step 3: Run whitespace check**

Run:

```bash
git diff --check origin/main..HEAD
```

Expected: no output, exit 0.

- [ ] **Step 4: Commit implementation**

Run:

```bash
git add src/algorithm/factorization.jl test/expert/park_woodburn_quillen_route_adapter.jl test/runtests.jl docs/superpowers/plans/2026-06-23-issue-114-quillen-patch-route-adapter.md
git commit -m "feat: add quillen patch route adapter"
```

Expected: commit succeeds.

---

## Plan Self-Review

- Spec coverage: the plan covers adapter construction, explicit route certificate embedding, #99/#105 fixture reuse, no public export, focused verification, and package verification.
- Placeholder scan: no incomplete markers remain.
- Type consistency: adapter, verifier, and route constructor names match the design document and planned tests.
- Scope check: public `elementary_factorization(A)` remains unchanged for multivariate inputs.

## Execution Choice

Under the Standing Answer Policy, choose **Subagent-Driven (recommended)** from
the writing-plans execution options. Use `superpowers:subagent-driven-development`
to execute this plan task-by-task.
