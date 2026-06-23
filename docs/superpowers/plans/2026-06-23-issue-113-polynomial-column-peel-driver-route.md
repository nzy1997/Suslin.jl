# Issue 113 Polynomial Column-Peel Driver Route Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Integrate the ordinary-polynomial recursive column-peel certificate into automatic verified public factorization routing.

**Architecture:** Keep the existing route shell in `src/algorithm/factorization.jl`. Add `:polynomial_column_peel` as the automatic route tag, keep `:recursive_column_peel` as an explicit alias, and enable recursive column peel in public polynomial route selection after narrower routes fail.

**Tech Stack:** Julia, Oscar exact polynomial matrices, existing Suslin route certificates, Test stdlib.

## Global Constraints

- Preserve public `elementary_factorization(A)` as a factor-returning API.
- Keep route order: `:fast_local_sl3`, then `:disjoint_local_blocks`, then `:polynomial_column_peel`, then `:staged_failure`.
- Preserve explicit `route=:recursive_column_peel` compatibility for issue 112 tests and fixtures.
- Route certificate verification must replay nested `PolynomialColumnPeelCertificate` data, not only check final multiplication.
- Preserve staged errors for determinant-one matrices whose selected column or recursive final block is outside implemented witness families.
- Do not add multivariate Quillen routing.
- Do not optimize factor count.
- Do not change Laurent behavior.
- Required focused commands: `julia --project=. -e 'include("test/expert/park_woodburn_route_certificate.jl")'` and `julia --project=. -e 'include("test/public/factorization_driver_shell.jl")'`.
- Required package command: `julia --project=. -e 'using Pkg; Pkg.test()'`.

---

## File Structure

- Modify `src/algorithm/factorization.jl`: add the `:polynomial_column_peel` route tag, route-tag-aware recursive constructor, automatic route selection, public shell enabling, and route verifier branch.
- Modify `test/expert/park_woodburn_route_certificate.jl`: add route-tag and nested tamper coverage.
- Modify `test/public/factorization_driver_shell.jl`: add public success and staged-failure checks for catalog recursive column-peel cases.

---

### Task 1: Add RED Public And Expert Route Tests

**Files:**
- Modify: `test/expert/park_woodburn_route_certificate.jl`
- Modify: `test/public/factorization_driver_shell.jl`

**Interfaces:**
- Consumes existing `ParkWoodburnPolynomialFixtureCatalog`.
- Consumes planned automatic route tag `:polynomial_column_peel`.
- Produces failing tests before implementation.

- [ ] **Step 1: Add route-certificate helper and expert assertions**

In `test/expert/park_woodburn_route_certificate.jl`, add a helper that corrupts
the nested peel evidence:

```julia
function _pw_corrupt_route_peel_evidence(cert)
    evidence = cert.evidence
    first_step = first(evidence.peel_steps)
    bad_last_column = copy(first_step.last_column)
    bad_last_column[1] += one(base_ring(evidence.original_matrix))
    bad_step = Suslin.PolynomialColumnPeelStep(
        first_step.dimension,
        first_step.input_matrix,
        bad_last_column,
        first_step.left_factors,
        first_step.after_left_matrix,
        first_step.right_factors,
        first_step.peeled_matrix,
        first_step.next_block,
    )
    bad_evidence = Suslin.PolynomialColumnPeelCertificate(
        evidence.original_matrix,
        Suslin.PolynomialColumnPeelStep[bad_step; evidence.peel_steps[2:end]],
        evidence.final_block,
        evidence.final_certificate,
        evidence.final_factors,
        evidence.factors,
        evidence.product,
        evidence.verification,
    )
    return _pw_replace_certificate(cert; evidence = bad_evidence)
end
```

Inside the route testset after `auto_staged_cert` assertions, add:

```julia
    recursive_supported_entry = entries["pw-poly-recursive-column-peel-sl3-qq"]
    auto_peel_cert = Suslin._polynomial_factorization_route_certificate(
        recursive_supported_entry.matrix,
    )
    @test auto_peel_cert.route == :polynomial_column_peel
    @test auto_peel_cert.status == :supported
    @test auto_peel_cert.evidence isa Suslin.PolynomialColumnPeelCertificate
    @test Suslin._verify_polynomial_factorization_route_certificate(auto_peel_cert)
    @test verify_factorization(auto_peel_cert.matrix, auto_peel_cert.factors)

    alias_peel_cert = Suslin._polynomial_factorization_route_certificate(
        recursive_supported_entry.matrix;
        route = :recursive_column_peel,
    )
    @test alias_peel_cert.route == :recursive_column_peel
    @test Suslin._verify_polynomial_factorization_route_certificate(alias_peel_cert)

    bad_peel_route_cert = _pw_corrupt_route_peel_evidence(auto_peel_cert)
    @test verify_factorization(bad_peel_route_cert.matrix, bad_peel_route_cert.factors)
    @test !Suslin._verify_polynomial_factorization_route_certificate(bad_peel_route_cert)
```

- [ ] **Step 2: Add public driver catalog checks**

In `test/public/factorization_driver_shell.jl`, define:

```julia
const PARK_WOODBURN_DRIVER_CATALOG_PATH =
    joinpath(@__DIR__, "..", "fixtures", "park_woodburn_polynomial_cases.jl")
```

At the end of the public driver testset, add:

```julia
    if !isdefined(Main, :ParkWoodburnPolynomialFixtureCatalog)
        include(PARK_WOODBURN_DRIVER_CATALOG_PATH)
    end
    entries = ParkWoodburnPolynomialFixtureCatalog.cases_by_id()
    recursive_supported = entries["pw-poly-recursive-column-peel-sl3-qq"].matrix
    recursive_factors = elementary_factorization(recursive_supported)
    @test verify_factorization(recursive_supported, recursive_factors)
    recursive_cert = Suslin._polynomial_factorization_route_certificate(recursive_supported)
    @test recursive_cert.route == :polynomial_column_peel
    @test recursive_factors == recursive_cert.factors
    @test recursive_cert.evidence isa Suslin.PolynomialColumnPeelCertificate
    @test Suslin._verify_polynomial_column_peel_certificate(recursive_cert.evidence)

    recursive_unsupported = entries["pw-poly-recursive-column-peel-gf2"].matrix
    recursive_unsupported_err =
        _captured_error(() -> elementary_factorization(recursive_unsupported))
    @test recursive_unsupported_err isa ArgumentError
    @test occursin("polynomial column-peel", sprint(showerror, recursive_unsupported_err)) ||
        occursin("SL_n reduction layer", sprint(showerror, recursive_unsupported_err))
```

- [ ] **Step 3: Run focused tests and confirm RED**

Run:

```bash
julia --project=. -e 'include("test/expert/park_woodburn_route_certificate.jl")'
julia --project=. -e 'include("test/public/factorization_driver_shell.jl")'
```

Expected: at least one failure because automatic route selection still returns
`:staged_failure` instead of `:polynomial_column_peel`, and the public driver
does not yet enable recursive column peel.

---

### Task 2: Enable Polynomial Column-Peel Auto Routing

**Files:**
- Modify: `src/algorithm/factorization.jl`

**Interfaces:**
- Produces automatic route tag `:polynomial_column_peel`.
- Preserves explicit `route=:recursive_column_peel`.
- Updates `_polynomial_verified_route_factors(A)` to enable recursive column peel.

- [ ] **Step 1: Add the new tag and helper predicates**

Update the route tag set and add helpers:

```julia
const _POLYNOMIAL_FACTORIZATION_ROUTE_TAGS = Set([
    :fast_local_sl3,
    :disjoint_local_blocks,
    :polynomial_column_peel,
    :recursive_column_peel,
    :staged_failure,
])

_is_polynomial_column_peel_route(route::Symbol) =
    route in (:polynomial_column_peel, :recursive_column_peel)
```

- [ ] **Step 2: Enable automatic selection before staged failure**

Change `_polynomial_factorization_route_certificate` default to
`allow_recursive_column_peel::Bool=true`. In automatic mode, call:

```julia
if allow_recursive_column_peel
    try
        return _polynomial_recursive_column_peel_route_certificate(
            A;
            route_tag = :polynomial_column_peel,
        )
    catch err
        err isa InterruptException && rethrow()
        err isa ArgumentError || rethrow()
    end
end
```

For explicit route handling, accept both tags:

```julia
elseif _is_polynomial_column_peel_route(route)
    return _polynomial_recursive_column_peel_route_certificate(A; route_tag = route)
```

- [ ] **Step 3: Store the requested route tag in the recursive certificate**

Change:

```julia
function _polynomial_recursive_column_peel_route_certificate(A)
```

to:

```julia
function _polynomial_recursive_column_peel_route_certificate(
    A;
    route_tag::Symbol = :polynomial_column_peel,
)
    _is_polynomial_column_peel_route(route_tag) ||
        throw(ArgumentError("unsupported polynomial column-peel route tag $(route_tag)"))
    evidence = _polynomial_column_peel_certificate(A)
    factors = copy(evidence.factors)
    product = _polynomial_route_factor_product(factors, base_ring(A), nrows(A))
    return _polynomial_route_certificate(A, route_tag, factors, product, evidence, :supported)
end
```

- [ ] **Step 4: Verify both route tags as successful column-peel routes**

In route verification, change:

```julia
successful_route = route in (:fast_local_sl3, :disjoint_local_blocks, :recursive_column_peel)
```

to include `_is_polynomial_column_peel_route(route)`. In
`_polynomial_route_evidence_ok`, replace the `cert.route == :recursive_column_peel`
branch with `_is_polynomial_column_peel_route(cert.route)`.

- [ ] **Step 5: Keep final recursive blocks non-recursive**

In `src/algorithm/polynomial_column_peel.jl`,
`_polynomial_column_peel_try_final_route` must keep recursive route selection
disabled:

```julia
_polynomial_factorization_route_certificate(
    current;
    route=final_route,
    allow_recursive_column_peel=false,
)
```

- [ ] **Step 6: Run focused tests and confirm GREEN**

Run:

```bash
julia --project=. -e 'include("test/expert/park_woodburn_route_certificate.jl")'
julia --project=. -e 'include("test/public/factorization_driver_shell.jl")'
```

Expected: both exit 0.

---

### Task 3: Final Verification And PR Preparation

**Files:**
- Modify: files from Tasks 1-2.

**Interfaces:**
- Produces a committed branch ready for PR.

- [ ] **Step 1: Run issue-required focused commands**

Run:

```bash
julia --project=. -e 'include("test/expert/park_woodburn_route_certificate.jl")'
julia --project=. -e 'include("test/public/factorization_driver_shell.jl")'
```

Expected: both exit 0.

- [ ] **Step 2: Run package verification**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: exits 0.

- [ ] **Step 3: Commit implementation**

Run:

```bash
git add docs/superpowers/specs/2026-06-23-issue-113-polynomial-column-peel-driver-route-design.md docs/superpowers/plans/2026-06-23-issue-113-polynomial-column-peel-driver-route.md src/algorithm/factorization.jl test/expert/park_woodburn_route_certificate.jl test/public/factorization_driver_shell.jl
git commit -m "feat: route polynomial column peel through driver"
```

Expected: commit succeeds.

---

## Plan Self-Review

- The plan covers automatic route selection, public dispatch, route verifier replay, and staged unsupported behavior.
- Tests are written before production changes and should fail for the expected missing automatic route.
- Existing `:recursive_column_peel` behavior remains covered by explicit alias tests.
- Laurent and Quillen behavior remain out of scope.
