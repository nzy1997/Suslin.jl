# Issue 263 Polynomial Column-Peel SL3 Final Route Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Allow ordinary-polynomial column peel to finish at a size-3 final block through the verified #184 evidence-backed `SL_3` `:quillen_patch` route.

**Architecture:** Keep `PolynomialFactorizationRouteCertificate` unchanged. Extend `PolynomialColumnPeelCertificate` with a stable final-route provenance marker and accept `:quillen_patch` as a peel final route only when the normal route certificate verifier succeeds and the evidence is `PolynomialSL3QuillenMurthyRouteEvidence`.

**Tech Stack:** Julia, Oscar, existing Suslin route certificates, existing Superpowers docs/test workflow.

## Global Constraints

- Preserve existing `:fast_local_sl3` and `:disjoint_local_blocks` final-route behavior.
- Do not accept supplied Quillen patch adapters, raw factor lists, or unverified patches as column-peel final routes.
- `:quillen_patch` final route provenance must be exactly `:issue184_evidence_backed_sl3` when accepted.
- The peel verifier must call `_verify_polynomial_factorization_route_certificate` for final route certificates.
- Focus changes on `src/algorithm/polynomial_column_peel.jl` and `test/expert/park_woodburn_polynomial_column_peel.jl`.
- Do not implement public recursive `SL_n` dispatch, Laurent/ToricBuilder support, #187 final acceptance, ECP-backed peel public routing, or Steinberg optimization.

---

### Task 1: Add Failing Column-Peel Tests For #184 Final Route

**Files:**
- Modify: `test/expert/park_woodburn_polynomial_column_peel.jl`

**Interfaces:**
- Consumes: existing `Suslin._polynomial_column_peel_certificate`, `Suslin._verify_polynomial_column_peel_certificate`, `Suslin._polynomial_factorization_route_certificate`, and existing test helper constructors.
- Produces: failing tests that require `cert.final_route_provenance`, `:quillen_patch` final-route acceptance, and rejection of tampered provenance/evidence/final factors.

- [ ] **Step 1: Add helper rebuild and #184 wrapper fixtures**

Add these helpers near the existing `_pw_poly_replace_route_certificate` helper:

```julia
function _pw_poly_rebuild(record; kwargs...)
    overrides = Dict{Symbol,Any}()
    for pair in kwargs
        overrides[pair.first] = pair.second
    end
    values = map(fieldnames(typeof(record))) do name
        get(overrides, name, getproperty(record, name))
    end
    return typeof(record)(values...)
end

function _pw_poly_issue184_sl3_route_case()
    R, (X, r, g) = Oscar.polynomial_ring(QQ, ["X", "r", "g"])
    p = X + r * g + one(R)
    q = one(R)
    s = one(R)
    lower = X + r * g
    A = matrix(R, [
        p q zero(R);
        lower s zero(R);
        zero(R) zero(R) one(R)
    ])
    @assert det(A) == one(R)
    return (; R, X, r, g, p, q, s, lower, A)
end

function _pw_poly_wrap_sl4_final_block(final_block, tail_entries)
    R = base_ring(final_block)
    length(tail_entries) == 3 || throw(ArgumentError("SL4 wrapper needs three tail entries"))
    wrapped = block_embedding(final_block, 4, [1, 2, 3])
    for row in 1:3
        wrapped[row, 4] = tail_entries[row]
    end
    return wrapped
end

function _pw_poly_certificate_with_provenance(
        cert,
        final_route_provenance)
    return Suslin.PolynomialColumnPeelCertificate(
        cert.original_matrix,
        cert.peel_steps,
        cert.final_block,
        cert.final_certificate,
        cert.final_factors,
        cert.factors,
        cert.product,
        cert.verification,
        final_route_provenance,
    )
end
```

- [ ] **Step 2: Add the RED test block**

Add this test block inside the existing `"Park-Woodburn ordinary polynomial column-peel certificates"` testset after the disjoint-block regression:

```julia
    issue184 = _pw_poly_issue184_sl3_route_case()
    issue184_wrapped = _pw_poly_wrap_sl4_final_block(
        issue184.A,
        [issue184.X + issue184.r, issue184.g + one(issue184.R), issue184.r * issue184.g],
    )
    issue184_cert = Suslin._polynomial_column_peel_certificate(issue184_wrapped)
    explicit_issue184_cert = Suslin._polynomial_column_peel_certificate(
        issue184_wrapped;
        final_route = :quillen_patch,
    )
    for cert in (issue184_cert, explicit_issue184_cert)
        @test cert.final_block == issue184.A
        @test cert.final_certificate.route == :quillen_patch
        @test cert.final_certificate.evidence isa Suslin.PolynomialSL3QuillenMurthyRouteEvidence
        @test cert.final_route_provenance == :issue184_evidence_backed_sl3
        @test cert.verification.final_route_provenance_ok
        _pw_poly_assert_real_peel_certificate(cert, issue184_wrapped)
    end

    adapter_only_final = _pw_poly_replace_route_certificate(
        issue184_cert.final_certificate;
        evidence = issue184_cert.final_certificate.evidence.quillen_route_adapter,
    )
    @test Suslin._verify_polynomial_factorization_route_certificate(adapter_only_final)
    adapter_only_peel = _pw_poly_replace_certificate(
        issue184_cert;
        final_certificate = adapter_only_final,
    )
    @test verify_factorization(adapter_only_peel.original_matrix, adapter_only_peel.factors)
    @test !Suslin._verify_polynomial_column_peel_certificate(adapter_only_peel)

    tampered_provenance =
        _pw_poly_certificate_with_provenance(issue184_cert, :tampered_quillen_patch)
    @test verify_factorization(tampered_provenance.original_matrix, tampered_provenance.factors)
    @test !Suslin._verify_polynomial_column_peel_certificate(tampered_provenance)

    tampered_witness = merge(
        issue184_cert.final_certificate.evidence.context.local_form_witness,
        (; monic_entry_position = (1, 2)),
    )
    tampered_context = _pw_poly_rebuild(
        issue184_cert.final_certificate.evidence.context;
        local_form_witness = tampered_witness,
    )
    tampered_evidence = _pw_poly_rebuild(
        issue184_cert.final_certificate.evidence;
        context = tampered_context,
    )
    tampered_final_certificate = _pw_poly_replace_route_certificate(
        issue184_cert.final_certificate;
        evidence = tampered_evidence,
    )
    tampered_final_evidence_peel = _pw_poly_replace_certificate(
        issue184_cert;
        final_certificate = tampered_final_certificate,
    )
    @test verify_factorization(
        tampered_final_evidence_peel.original_matrix,
        tampered_final_evidence_peel.factors,
    )
    @test !Suslin._verify_polynomial_column_peel_certificate(tampered_final_evidence_peel)

    tampered_final_factors = copy(issue184_cert.final_factors)
    tampered_final_factors[1] =
        tampered_final_factors[1] *
        elementary_matrix(3, 1, 2, one(issue184.R), issue184.R)
    tampered_final_factor_peel = _pw_poly_replace_certificate(
        issue184_cert;
        final_factors = tampered_final_factors,
    )
    @test verify_factorization(tampered_final_factor_peel.original_matrix, tampered_final_factor_peel.factors)
    @test !Suslin._verify_polynomial_column_peel_certificate(tampered_final_factor_peel)
```

- [ ] **Step 3: Run RED verification**

Run:

```bash
julia --project=. -e 'include("test/expert/park_woodburn_polynomial_column_peel.jl")'
```

Expected before implementation: FAIL because `:quillen_patch` is not an accepted column-peel final route and/or `final_route_provenance` is missing.

- [ ] **Step 4: Commit the failing tests**

```bash
git add test/expert/park_woodburn_polynomial_column_peel.jl
git commit -m "test: cover issue 263 column peel sl3 final route"
```

### Task 2: Implement Provenance-Gated Quillen Final Routes

**Files:**
- Modify: `src/algorithm/polynomial_column_peel.jl`

**Interfaces:**
- Consumes: `PolynomialSL3QuillenMurthyRouteEvidence`, `_verify_polynomial_sl3_quillen_murthy_route_evidence`, `_polynomial_factorization_route_certificate`, and `_verify_polynomial_factorization_route_certificate`.
- Produces: `PolynomialColumnPeelCertificate.final_route_provenance` and final-route verification that accepts `:quillen_patch` only for #184 evidence-backed `SL_3`.

- [ ] **Step 1: Add provenance constants and helpers**

Add near the `PolynomialColumnPeelCertificate` definition:

```julia
const _POLYNOMIAL_COLUMN_PEEL_FAST_FINAL_ROUTE_PROVENANCE = :fast_local_sl3
const _POLYNOMIAL_COLUMN_PEEL_BLOCK_FINAL_ROUTE_PROVENANCE = :disjoint_local_blocks
const _POLYNOMIAL_COLUMN_PEEL_ISSUE184_FINAL_ROUTE_PROVENANCE = :issue184_evidence_backed_sl3
const _POLYNOMIAL_COLUMN_PEEL_UNSUPPORTED_FINAL_ROUTE_PROVENANCE = :unsupported_final_route

function _polynomial_column_peel_quillen_issue184_final_route_ok(final_certificate)::Bool
    try
        final_certificate.route == :quillen_patch || return false
        final_certificate.evidence isa PolynomialSL3QuillenMurthyRouteEvidence || return false
        _verify_polynomial_sl3_quillen_murthy_route_evidence(final_certificate.evidence) || return false
        return true
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end

function _polynomial_column_peel_supported_final_route_ok(final_certificate)::Bool
    final_certificate.route in (:fast_local_sl3, :disjoint_local_blocks) && return true
    return _polynomial_column_peel_quillen_issue184_final_route_ok(final_certificate)
end

function _polynomial_column_peel_final_route_provenance(final_certificate)::Symbol
    try
        final_certificate.route == :fast_local_sl3 &&
            return _POLYNOMIAL_COLUMN_PEEL_FAST_FINAL_ROUTE_PROVENANCE
        final_certificate.route == :disjoint_local_blocks &&
            return _POLYNOMIAL_COLUMN_PEEL_BLOCK_FINAL_ROUTE_PROVENANCE
        _polynomial_column_peel_quillen_issue184_final_route_ok(final_certificate) &&
            return _POLYNOMIAL_COLUMN_PEEL_ISSUE184_FINAL_ROUTE_PROVENANCE
    catch err
        err isa InterruptException && rethrow()
    end
    return _POLYNOMIAL_COLUMN_PEEL_UNSUPPORTED_FINAL_ROUTE_PROVENANCE
end
```

- [ ] **Step 2: Extend the certificate struct with compatibility constructor**

Change the struct and add an outer constructor:

```julia
struct PolynomialColumnPeelCertificate
    original_matrix
    peel_steps::Vector{PolynomialColumnPeelStep}
    final_block
    final_certificate
    final_factors::Vector
    factors::Vector
    product
    verification
    final_route_provenance::Symbol
end

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
    return PolynomialColumnPeelCertificate(
        original_matrix,
        peel_steps,
        final_block,
        final_certificate,
        final_factors,
        factors,
        product,
        verification,
        _polynomial_column_peel_final_route_provenance(final_certificate),
    )
end
```

- [ ] **Step 3: Update final-route validation and candidate selection**

Change `_validate_polynomial_column_peel_final_route` to allow
`:quillen_patch`, and update the candidate tuple in
`_polynomial_column_peel_try_final_route`:

```julia
final_route in (:fast_local_sl3, :disjoint_local_blocks, :quillen_patch) ||
    throw(ArgumentError("polynomial column-peel final route must be :fast_local_sl3, :disjoint_local_blocks, or :quillen_patch"))
```

```julia
candidate_routes =
    final_route === nothing ?
    (:fast_local_sl3, :disjoint_local_blocks, :quillen_patch) :
    (final_route,)
```

- [ ] **Step 4: Build route certificates for `:quillen_patch` through normal route selection but provenance-gate them**

Replace the route-certificate construction in `_polynomial_column_peel_try_final_route` with:

```julia
certificate = try
    if route == :quillen_patch
        _polynomial_factorization_route_certificate(
            current;
            allow_recursive_column_peel=false,
        )
    else
        _polynomial_factorization_route_certificate(
            current;
            route=route,
            allow_recursive_column_peel=false,
        )
    end
catch err
    err isa InterruptException && rethrow()
    err isa ArgumentError || rethrow()
    nothing
end
```

Then accept only:

```julia
if certificate !== nothing &&
        certificate.status == :supported &&
        certificate.route == route &&
        _polynomial_column_peel_supported_final_route_ok(certificate) &&
        certificate.matrix != identity_matrix(base_ring(current), nrows(current))
    return certificate
end
```

- [ ] **Step 5: Add provenance verification to core verification**

Add helper:

```julia
function _polynomial_column_peel_final_route_provenance_ok(cert)::Bool
    try
        hasproperty(cert, :final_route_provenance) || return false
        expected = _polynomial_column_peel_final_route_provenance(cert.final_certificate)
        expected == _POLYNOMIAL_COLUMN_PEEL_UNSUPPORTED_FINAL_ROUTE_PROVENANCE && return false
        return cert.final_route_provenance == expected
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end
```

In `_polynomial_column_peel_core_verification`, compute
`final_route_provenance_ok` and include it in `overall_core_ok` and the returned
named tuple.

- [ ] **Step 6: Tighten `_polynomial_column_peel_final_certificate_ok`**

Replace the hard-coded route tuple check with:

```julia
_polynomial_column_peel_supported_final_route_ok(final_certificate) || return false
```

Keep the existing normal route verifier call:

```julia
_verify_polynomial_factorization_route_certificate(final_certificate) || return false
```

- [ ] **Step 7: Run GREEN verification**

Run:

```bash
julia --project=. -e 'include("test/expert/park_woodburn_polynomial_column_peel.jl")'
```

Expected after implementation: PASS.

- [ ] **Step 8: Commit implementation**

```bash
git add src/algorithm/polynomial_column_peel.jl
git commit -m "feat: allow evidence-backed sl3 peel finals"
```

### Task 3: Route Certificate Regression And Full Verification

**Files:**
- Test: `test/expert/park_woodburn_route_certificate.jl`
- Test: `test/runtests.jl`

**Interfaces:**
- Consumes: implementation from Task 2.
- Produces: verified package state and no accidental route-certificate regression.

- [ ] **Step 1: Run focused route certificate regression**

Run:

```bash
julia --project=. -e 'include("test/expert/park_woodburn_route_certificate.jl")'
```

Expected: PASS.

- [ ] **Step 2: Run package test entry point**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: PASS.

- [ ] **Step 3: Run whitespace check**

Run:

```bash
git diff --check
```

Expected: exit 0.

- [ ] **Step 4: Commit any final test/cleanup changes**

If Task 3 required edits, commit them:

```bash
git add test/runtests.jl test/expert/park_woodburn_route_certificate.jl test/expert/park_woodburn_polynomial_column_peel.jl src/algorithm/polynomial_column_peel.jl
git commit -m "test: verify issue 263 route regressions"
```

If no files changed, do not create an empty commit.

## Plan Self-Review

Spec coverage is complete: the plan adds the positive #184 final route, the
stable provenance marker, normal final route certificate verification, and
negative controls for adapter-only evidence, provenance tampering, witness
tampering, and final factor tampering. The route tags and field names are
consistent across tasks. No marker text remains.
