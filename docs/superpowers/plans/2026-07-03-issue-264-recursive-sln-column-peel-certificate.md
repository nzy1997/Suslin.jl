# Issue 264 Recursive SLn Column-Peel Certificate Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Assemble and verify the #186 recursive ordinary-polynomial `SL_n` column-peel certificate.

**Architecture:** Upgrade `PolynomialColumnPeelCertificate` with recomputed descent and #186 mainline-support metadata, while preserving existing constructors and route behavior. Add a narrow supplied-Quillen `SL_3` route evidence wrapper so column peel can accept #260 final blocks without accepting raw adapter-only finals.

**Tech Stack:** Julia, Oscar, existing Suslin internal certificate records, `Test`.

## Global Constraints

- Do not change public `elementary_factorization` dispatch.
- Keep legacy `:disjoint_local_blocks` and fast-local recursive routes verifiable for regression coverage, but do not mark them as #186 mainline support.
- #186 mainline support is true only when every peel step verifies through #262 ECP evidence, strict dimension descent reaches a `3 x 3` final block, and the final route verifies through #263/#184 evidence-backed `SL_3` route evidence.
- Preserve #263 adapter-only rejection: raw `PolynomialQuillenPatchRouteAdapter` evidence is not sufficient for a column-peel final route.
- Required focused verification commands:
  `julia --project=. -e 'include("test/expert/park_woodburn_sln_recursive_driver.jl")'`
  `julia --project=. -e 'include("test/expert/park_woodburn_polynomial_column_peel.jl")'`
- Required package verification command:
  `julia --project=. -e 'using Pkg; Pkg.test()'`

---

### Task 1: Add Recursive Driver Red Tests

**Files:**
- Create: `test/expert/park_woodburn_sln_recursive_driver.jl`
- Modify: `test/runtests.jl`

**Interfaces:**
- Consumes: `Suslin._polynomial_column_peel_certificate`, `Suslin._verify_polynomial_column_peel_certificate`, `Suslin._polynomial_column_peel_core_verification`, `verify_factorization`.
- Produces: Failing coverage for the new recursive metadata and supplied final-route evidence requirements.

- [ ] **Step 1: Write the failing test file**

Create `test/expert/park_woodburn_sln_recursive_driver.jl` with helper functions matching the local style in `park_woodburn_polynomial_column_peel.jl`:

```julia
using Test
using Oscar
using Suslin

const PARK_WOODBURN_SLN_RECURSIVE_DRIVER_CATALOG_PATH =
    joinpath(@__DIR__, "..", "fixtures", "park_woodburn_sln_driver_cases.jl")

function _sln_recursive_product(factors, R, n::Int)
    product = identity_matrix(R, n)
    for factor in factors
        product *= factor
    end
    return product
end

function _sln_recursive_replace_certificate(
        cert;
        original_matrix = cert.original_matrix,
        peel_steps = cert.peel_steps,
        final_block = cert.final_block,
        final_certificate = cert.final_certificate,
        final_factors = cert.final_factors,
        factors = cert.factors,
        product = cert.product,
        verification = cert.verification,
        descent_metadata = cert.descent_metadata,
        mainline_support_metadata = cert.mainline_support_metadata,
        final_route_provenance = cert.final_route_provenance)
    return Suslin.PolynomialColumnPeelCertificate(
        original_matrix,
        peel_steps,
        final_block,
        final_certificate,
        final_factors,
        factors,
        product,
        verification,
        descent_metadata,
        mainline_support_metadata,
        final_route_provenance,
    )
end

function _sln_recursive_replace_route_certificate(
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

function _sln_recursive_rebuild(record; kwargs...)
    overrides = Dict{Symbol,Any}(pair.first => pair.second for pair in kwargs)
    values = map(fieldnames(typeof(record))) do name
        get(overrides, name, getproperty(record, name))
    end
    return typeof(record)(values...)
end

function _sln_recursive_assert_mainline_certificate(cert, entry)
    @test Suslin._verify_polynomial_column_peel_certificate(cert)
    @test cert.original_matrix == entry.matrix
    @test length(cert.peel_steps) == entry.expected_peel_count
    @test tuple((step.dimension for step in cert.peel_steps)..., nrows(cert.final_block)) ==
          entry.descent_dimensions
    @test cert.final_block == entry.final_route.matrix
    @test nrows(cert.final_block) == 3
    @test ncols(cert.final_block) == 3
    @test cert.final_certificate.route == :quillen_patch
    @test cert.final_route_provenance == :issue184_evidence_backed_sl3
    @test cert.descent_metadata.strict_dimension_descent
    @test cert.descent_metadata.final_block_is_sl3
    @test cert.mainline_support_metadata.supported
    @test cert.mainline_support_metadata.marker == :issue186_mainline
    @test cert.mainline_support_metadata.final_route_issue184_ok
    @test cert.mainline_support_metadata.peel_steps_ecp_verified
    @test cert.mainline_support_metadata.factor_replay_ok
    @test cert.mainline_support_metadata.reconstruction_ok
    @test cert.product == entry.matrix
    @test _sln_recursive_product(cert.factors, base_ring(entry.matrix), nrows(entry.matrix)) == entry.matrix
    @test verify_factorization(entry.matrix, cert.factors)
    for step in cert.peel_steps
        @test Suslin._polynomial_column_peel_step_verification(step).overall_ok
        @test Suslin.verify_ecp_column_reduction(step.ecp_evidence)
    end
end

@testset "Park-Woodburn recursive SLn column-peel certificate" begin
    if !isdefined(Main, :ParkWoodburnSLnDriverFixtureCatalog)
        include(PARK_WOODBURN_SLN_RECURSIVE_DRIVER_CATALOG_PATH)
    end
    entries = ParkWoodburnSLnDriverFixtureCatalog.cases_by_id()

    sl4 = entries["sln-driver-sl4-gf2-ecp-mainline"]
    sl4_cert = Suslin._polynomial_column_peel_certificate(sl4.matrix)
    _sln_recursive_assert_mainline_certificate(sl4_cert, sl4)

    sl5 = entries["sln-driver-sl5-gf2-two-step"]
    sl5_cert = Suslin._polynomial_column_peel_certificate(sl5.matrix)
    _sln_recursive_assert_mainline_certificate(sl5_cert, sl5)
    @test length(sl5_cert.peel_steps) >= 2

    legacy = entries["sln-driver-legacy-recursive-column-peel-qq"]
    legacy_cert = Suslin._polynomial_column_peel_certificate(legacy.matrix)
    @test Suslin._verify_polynomial_column_peel_certificate(legacy_cert)
    @test !legacy_cert.mainline_support_metadata.supported
    @test legacy_cert.mainline_support_metadata.marker == :not_issue186_mainline
    @test :missing_issue184_final_sl3_route in legacy_cert.mainline_support_metadata.reason_codes

    reordered_steps = reverse(copy(sl5_cert.peel_steps))
    reordered = _sln_recursive_replace_certificate(sl5_cert; peel_steps = reordered_steps)
    @test !Suslin._verify_polynomial_column_peel_certificate(reordered)
    @test !Suslin._polynomial_column_peel_core_verification(reordered).descent_metadata_ok

    duplicated_steps = vcat(copy(sl5_cert.peel_steps), [last(sl5_cert.peel_steps)])
    duplicated = _sln_recursive_replace_certificate(sl5_cert; peel_steps = duplicated_steps)
    @test !Suslin._verify_polynomial_column_peel_certificate(duplicated)

    bad_final_evidence = _sln_recursive_rebuild(
        sl5_cert.final_certificate.evidence;
        replay_metadata = merge(
            sl5_cert.final_certificate.evidence.replay_metadata,
            (; source = :tampered_recursive_final_route),
        ),
    )
    bad_final_certificate = _sln_recursive_replace_route_certificate(
        sl5_cert.final_certificate;
        evidence = bad_final_evidence,
    )
    bad_final = _sln_recursive_replace_certificate(sl5_cert; final_certificate = bad_final_certificate)
    @test !Suslin._verify_polynomial_column_peel_certificate(bad_final)

    bad_mainline_metadata = merge(
        sl5_cert.mainline_support_metadata,
        (; marker = :not_issue186_mainline),
    )
    bad_metadata = _sln_recursive_replace_certificate(
        sl5_cert;
        mainline_support_metadata = bad_mainline_metadata,
    )
    @test verify_factorization(bad_metadata.original_matrix, bad_metadata.factors)
    @test !Suslin._verify_polynomial_column_peel_certificate(bad_metadata)

    bad_factors = copy(sl5_cert.factors)
    bad_factors[1] =
        bad_factors[1] *
        elementary_matrix(nrows(sl5.matrix), 1, 2, one(base_ring(sl5.matrix)), base_ring(sl5.matrix))
    bad_factor_cert = _sln_recursive_replace_certificate(sl5_cert; factors = bad_factors)
    @test !Suslin._verify_polynomial_column_peel_certificate(bad_factor_cert)

    bad_product = _sln_recursive_replace_certificate(
        sl5_cert;
        product = identity_matrix(base_ring(sl5.matrix), nrows(sl5.matrix)),
    )
    @test !Suslin._verify_polynomial_column_peel_certificate(bad_product)
end
```

- [ ] **Step 2: Register the expert test**

Add `"expert/park_woodburn_sln_recursive_driver.jl"` immediately after
`"expert/park_woodburn_polynomial_column_peel.jl"` in `test/runtests.jl`.

- [ ] **Step 3: Run red test**

Run:

```bash
julia --project=. -e 'include("test/expert/park_woodburn_sln_recursive_driver.jl")'
```

Expected: FAIL before implementation, because `_polynomial_column_peel_certificate` cannot yet build the #260 supplied-Quillen final route and the certificate lacks `descent_metadata` and `mainline_support_metadata`.

- [ ] **Step 4: Commit red tests**

```bash
git add test/expert/park_woodburn_sln_recursive_driver.jl test/runtests.jl
git commit -m "test: cover recursive sln column peel certificate"
```

### Task 2: Add Supplied-Quillen SL3 Route Evidence Wrapper

**Files:**
- Modify: `src/algorithm/factorization.jl`
- Modify: `src/algorithm/polynomial_column_peel.jl`
- Test: `test/expert/park_woodburn_sln_recursive_driver.jl`
- Test: `test/expert/park_woodburn_polynomial_column_peel.jl`

**Interfaces:**
- Produces: `PolynomialSL3SuppliedQuillenRouteEvidence`, `_verify_polynomial_sl3_supplied_quillen_route_evidence`, and a route certificate accepted by column peel as #184/#263 evidence.
- Consumes: `_polynomial_quillen_supplied_evidence_data`, `_polynomial_quillen_supplied_evidence_patch`, `_polynomial_quillen_patch_route_adapter`.

- [ ] **Step 1: Add the new evidence type near existing route evidence structs**

In `src/algorithm/factorization.jl`, place this after `PolynomialSL3QuillenMurthyRouteEvidence`:

```julia
struct PolynomialSL3SuppliedQuillenRouteEvidence
    target
    route::Symbol
    supplied_evidence
    quillen_route_adapter::PolynomialQuillenPatchRouteAdapter
    replay_metadata::NamedTuple
    verification
end
```

- [ ] **Step 2: Add constructor and verifier helpers**

Add functions that build the wrapper from supplied Quillen evidence:

```julia
function _polynomial_sl3_supplied_quillen_route_fields(A; metadata = (;))
    nrows(A) == 3 && ncols(A) == 3 ||
        throw(ArgumentError("SL_3 supplied-Quillen route requires a 3 x 3 input"))
    data = _polynomial_quillen_supplied_evidence_data(A)
    data === nothing &&
        throw(ArgumentError("SL_3 supplied-Quillen route requires supplied Quillen evidence"))
    patch = _polynomial_quillen_supplied_evidence_patch(A)
    adapter = _polynomial_quillen_patch_route_adapter(A, patch)
    replay_metadata = (;
        source = :sl3_supplied_quillen_polynomial_route,
        route_issue_id = "#263",
        driver_issue_id = "#184",
        patch_replay_metadata = adapter.quillen_patch.replay_metadata,
        adapter_replay_metadata = adapter.replay_metadata,
        metadata,
    )
    return (;
        target = adapter.target_matrix,
        route = :quillen_patch,
        supplied_evidence = data,
        quillen_route_adapter = adapter,
        replay_metadata,
    )
end

function _polynomial_sl3_supplied_quillen_route_core_verification(evidence)
    route_ok = evidence.route == :quillen_patch
    fields = route_ok ? _polynomial_sl3_supplied_quillen_route_fields(evidence.target) : nothing
    supplied_evidence_ok = fields !== nothing && evidence.supplied_evidence == fields.supplied_evidence
    adapter_ok =
        fields !== nothing &&
        evidence.quillen_route_adapter == fields.quillen_route_adapter &&
        _verify_polynomial_quillen_patch_route_adapter(evidence.quillen_route_adapter)
    replay_metadata_ok = fields !== nothing && evidence.replay_metadata == fields.replay_metadata
    overall_core_ok = route_ok && supplied_evidence_ok && adapter_ok && replay_metadata_ok
    return (;
        route_ok,
        supplied_evidence_ok,
        adapter_ok,
        replay_metadata_ok,
        overall_core_ok,
    )
end

function _polynomial_sl3_supplied_quillen_route_verification(evidence)
    core = _polynomial_sl3_supplied_quillen_route_core_verification(evidence)
    stored_verification_ok = evidence.verification == core
    return merge(core, (; stored_verification_ok, overall_ok = core.overall_core_ok && stored_verification_ok))
end

function _verify_polynomial_sl3_supplied_quillen_route_evidence(evidence)::Bool
    try
        return _polynomial_sl3_supplied_quillen_route_verification(evidence).overall_ok
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end

function _polynomial_sl3_supplied_quillen_route_evidence(A; metadata = (;))
    fields = _polynomial_sl3_supplied_quillen_route_fields(A; metadata)
    raw = PolynomialSL3SuppliedQuillenRouteEvidence(values(merge(fields, (; verification = nothing,)))...)
    verification = _polynomial_sl3_supplied_quillen_route_core_verification(raw)
    evidence = PolynomialSL3SuppliedQuillenRouteEvidence(values(merge(fields, (; verification,)))...)
    _verify_polynomial_sl3_supplied_quillen_route_evidence(evidence) ||
        error("internal SL_3 supplied-Quillen route evidence verification failed")
    return evidence
end

function _polynomial_sl3_supplied_quillen_route_certificate(A)
    evidence = _polynomial_sl3_supplied_quillen_route_evidence(A)
    adapter = evidence.quillen_route_adapter
    factors = copy(adapter.global_elementary_factors)
    return _polynomial_route_certificate(
        adapter.target_matrix,
        :quillen_patch,
        factors,
        adapter.product,
        evidence,
        :supported,
    )
end
```

Wrap `_polynomial_sl3_supplied_quillen_route_core_verification` in a `try/catch`
that returns all four booleans as `false` when non-interrupt errors are raised.

- [ ] **Step 3: Teach route verification about the wrapper**

In `_polynomial_route_evidence_ok`, add a `cert.evidence isa PolynomialSL3SuppliedQuillenRouteEvidence` branch under `cert.route == :quillen_patch` that verifies:

```julia
cert.evidence.target == cert.matrix &&
_verify_polynomial_sl3_supplied_quillen_route_evidence(cert.evidence) &&
_polynomial_route_factor_sequences_equal(
    cert.factors,
    cert.evidence.quillen_route_adapter.global_elementary_factors,
)
```

- [ ] **Step 4: Keep adapter-only rejection but accept the wrapper**

In `src/algorithm/polynomial_column_peel.jl`, update
`_polynomial_column_peel_quillen_issue184_final_route_ok` to accept either
`PolynomialSL3QuillenMurthyRouteEvidence` or
`PolynomialSL3SuppliedQuillenRouteEvidence`.

In `_polynomial_column_peel_try_final_route`, when `route == :quillen_patch`,
try `_polynomial_sl3_supplied_quillen_route_certificate(current)` before the
generic `_polynomial_factorization_route_certificate` fallback.

- [ ] **Step 5: Run focused tests**

Run:

```bash
julia --project=. -e 'include("test/expert/park_woodburn_sln_recursive_driver.jl")'
```

Expected after Task 2: it may still FAIL on missing certificate metadata from Task 3, but the failure should no longer be "requires a supported final route at size 3".

- [ ] **Step 6: Commit route wrapper**

```bash
git add src/algorithm/factorization.jl src/algorithm/polynomial_column_peel.jl
git commit -m "feat: add supplied quillen sl3 route evidence"
```

### Task 3: Add Recursive Descent and Mainline Metadata

**Files:**
- Modify: `src/algorithm/polynomial_column_peel.jl`
- Test: `test/expert/park_woodburn_sln_recursive_driver.jl`
- Test: `test/expert/park_woodburn_polynomial_column_peel.jl`

**Interfaces:**
- Produces: `descent_metadata` and `mainline_support_metadata` on `PolynomialColumnPeelCertificate`.
- Extends: `_polynomial_column_peel_core_verification` with `descent_metadata_ok` and `mainline_support_metadata_ok`.

- [ ] **Step 1: Extend the certificate struct and compatibility constructors**

Add `descent_metadata::NamedTuple` and `mainline_support_metadata::NamedTuple`
before `final_route_provenance::Symbol`.

Keep old constructor arities by defining:

```julia
function PolynomialColumnPeelCertificate(
    original_matrix,
    peel_steps::Vector{PolynomialColumnPeelStep},
    final_block,
    final_certificate,
    final_factors::Vector,
    factors::Vector,
    product,
    verification,
)
    final_route_provenance = _polynomial_column_peel_final_route_provenance(final_certificate)
    descent_metadata = _polynomial_column_peel_certificate_descent_metadata(
        peel_steps,
        original_matrix,
        final_block,
        final_route_provenance,
    )
    mainline_support_metadata = _polynomial_column_peel_mainline_support_metadata(
        original_matrix,
        peel_steps,
        final_block,
        final_certificate,
        final_factors,
        factors,
        product,
        final_route_provenance,
        descent_metadata,
    )
    return PolynomialColumnPeelCertificate(
        original_matrix,
        peel_steps,
        final_block,
        final_certificate,
        final_factors,
        factors,
        product,
        verification,
        descent_metadata,
        mainline_support_metadata,
        final_route_provenance,
    )
end
```

Add a matching nine-argument constructor that accepts explicit
`final_route_provenance` and computes both metadata records.

- [ ] **Step 2: Add descent metadata helpers**

Implement `_polynomial_column_peel_certificate_descent_metadata`:

```julia
function _polynomial_column_peel_certificate_descent_metadata(
    peel_steps,
    original_matrix,
    final_block,
    final_route_provenance::Symbol,
)
    step_dimensions = tuple((step.dimension for step in peel_steps)...)
    next_dimensions = tuple((nrows(step.next_block) for step in peel_steps)...)
    descent_dimensions = tuple(nrows(original_matrix), next_dimensions...)
    expected_step_count = nrows(original_matrix) - nrows(final_block)
    strict_dimension_descent =
        length(peel_steps) == expected_step_count &&
        all(idx -> next_dimensions[idx] == step_dimensions[idx] - 1, eachindex(step_dimensions)) &&
        all(idx -> idx == 1 || step_dimensions[idx] == next_dimensions[idx - 1], eachindex(step_dimensions))
    final_block_is_sl3 = nrows(final_block) == 3 && ncols(final_block) == 3
    return (;
        route = :park_woodburn_recursive_column_peel,
        input_dimension = nrows(original_matrix),
        final_dimension = nrows(final_block),
        step_count = length(peel_steps),
        expected_step_count,
        step_dimensions,
        next_dimensions,
        descent_dimensions,
        strict_dimension_descent,
        final_block_is_sl3,
        final_route_provenance,
    )
end
```

Implement `_polynomial_column_peel_descent_metadata_ok(cert)::Bool` to compare
stored metadata with a recomputation and require `strict_dimension_descent`.

- [ ] **Step 3: Add mainline-support metadata helpers**

Implement `_polynomial_column_peel_mainline_support_metadata` to recompute:

```julia
peel_steps_ecp_verified = all(_polynomial_column_peel_left_certificate_ok, peel_steps)
final_route_issue184_ok = _polynomial_column_peel_quillen_issue184_final_route_ok(final_certificate)
factor_replay_ok = _factor_sequences_equal(
    factors,
    _replay_polynomial_column_peel_factors(peel_steps, final_factors, base_ring(original_matrix)),
)
product_replay_ok = product == _factor_product(factors, base_ring(original_matrix), nrows(original_matrix))
reconstruction_ok = verify_factorization(original_matrix, factors)
supported = peel_steps_ecp_verified &&
    descent_metadata.strict_dimension_descent &&
    descent_metadata.final_block_is_sl3 &&
    final_route_issue184_ok &&
    final_route_provenance == _POLYNOMIAL_COLUMN_PEEL_ISSUE184_FINAL_ROUTE_PROVENANCE &&
    factor_replay_ok &&
    product_replay_ok &&
    reconstruction_ok
```

Build `reason_codes` from failed conditions using symbols:
`:missing_ecp_peel_evidence`, `:non_strict_dimension_descent`,
`:final_block_not_sl3`, `:missing_issue184_final_sl3_route`,
`:factor_replay_mismatch`, `:product_replay_mismatch`,
`:factorization_reconstruction_mismatch`.

Return:

```julia
(;
    issue_id = "#186",
    marker = supported ? :issue186_mainline : :not_issue186_mainline,
    supported,
    reason_codes = tuple(reason_codes...),
    peel_steps_ecp_verified,
    strict_dimension_descent = descent_metadata.strict_dimension_descent,
    final_block_is_sl3 = descent_metadata.final_block_is_sl3,
    final_route_issue184_ok,
    final_route_provenance,
    factor_replay_ok,
    product_replay_ok,
    reconstruction_ok,
)
```

Implement `_polynomial_column_peel_mainline_support_metadata_ok(cert)::Bool`
to compare stored metadata with recomputation and to require
`marker == :issue186_mainline` exactly when `supported == true`.

- [ ] **Step 4: Extend core verification**

In `_polynomial_column_peel_core_verification`, compute:

```julia
descent_metadata_ok = _polynomial_column_peel_descent_metadata_ok(cert)
mainline_support_metadata_ok = _polynomial_column_peel_mainline_support_metadata_ok(cert)
```

Include both booleans in `overall_core_ok` and in the returned named tuple.

- [ ] **Step 5: Run focused tests**

Run:

```bash
julia --project=. -e 'include("test/expert/park_woodburn_sln_recursive_driver.jl")'
julia --project=. -e 'include("test/expert/park_woodburn_polynomial_column_peel.jl")'
```

Expected: both exit 0.

- [ ] **Step 6: Commit metadata implementation**

```bash
git add src/algorithm/polynomial_column_peel.jl test/expert/park_woodburn_sln_recursive_driver.jl test/runtests.jl
git commit -m "feat: add recursive sln column peel metadata"
```

### Task 4: Final Verification and Cleanup

**Files:**
- Modify only files required by failures found in this task.

**Interfaces:**
- Consumes: all previous tasks.
- Produces: verified branch ready for PR.

- [ ] **Step 1: Run required focused verification**

```bash
julia --project=. -e 'include("test/expert/park_woodburn_sln_recursive_driver.jl")'
julia --project=. -e 'include("test/expert/park_woodburn_polynomial_column_peel.jl")'
```

Expected: both exit 0.

- [ ] **Step 2: Run required package verification**

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: exits 0.

- [ ] **Step 3: Run diff hygiene**

```bash
git diff --check
```

Expected: exits 0.

- [ ] **Step 4: Commit any final fixes**

If Step 1, Step 2, or Step 3 requires a fix, commit only the touched files:

```bash
git add src/algorithm/factorization.jl src/algorithm/polynomial_column_peel.jl test/expert/park_woodburn_sln_recursive_driver.jl test/expert/park_woodburn_polynomial_column_peel.jl test/runtests.jl
git commit -m "fix: stabilize recursive sln column peel certificate"
```

No commit is needed if the working tree is clean after verification.
