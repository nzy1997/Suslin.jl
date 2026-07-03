# Issue 265 Public SLn Recursive Route Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Route supported ordinary-polynomial `SL_n`, `n > 3`, public factorization through the #264 recursive driver only when its #186 mainline evidence verifies.

**Architecture:** Reuse `:polynomial_column_peel` as the public route tag and make the nested `PolynomialColumnPeelCertificate.mainline_support_metadata` the public support gate. Automatic routing tries the recursive driver before legacy disjoint blocks for real recursive candidates, and staged route evidence reports recursive reason codes when the recursive proof boundary fails.

**Tech Stack:** Julia, Oscar, `Test`, existing Suslin route-certificate structs and Park-Woodburn fixtures.

## Global Constraints

- Keep the public route tag as `:polynomial_column_peel`; do not add a new tag unless the implementation cannot verify provenance with existing metadata.
- Public recursive `SL_n`, `n > 3`, support requires `mainline_support_metadata.issue_id == "#186"`, `marker == :issue186_mainline`, `supported == true`, and `final_route_provenance == :issue184_evidence_backed_sl3`.
- Legacy `:disjoint_local_blocks`, fast-final column peel, and old recursive-column-peel evidence must not verify as public #186 proof.
- Staged public failures must distinguish `:missing_ecp_evidence`, `:missing_final_sl3_route`, determinant failure, and unsupported coefficient rings.
- Do not broaden Laurent, ToricBuilder, or arbitrary coefficient-ring support.

---

### Task 1: Add Failing Public Recursive Route Tests

**Files:**
- Modify: `test/expert/park_woodburn_route_certificate.jl`
- Modify: `test/public/factorization_driver_shell.jl`
- Modify: `test/public/park_woodburn_polynomial_factorization.jl`

**Interfaces:**
- Consumes: `ParkWoodburnSLnDriverFixtureCatalog.cases_by_id()`, `Suslin._polynomial_factorization_route_certificate`, `Suslin._verify_polynomial_factorization_route_certificate`, `elementary_factorization`, `verify_factorization`.
- Produces: failing tests for public #186 mainline acceptance, legacy disjoint/fast-final rejection, and tampered mainline metadata rejection.

- [ ] **Step 1: Add SLn fixture constants and helper assertions**

Add this constant next to the existing fixture path constants in each touched test file that needs #260/#264 fixtures:

```julia
const PARK_WOODBURN_SLN_DRIVER_CATALOG_PATH =
    joinpath(@__DIR__, "..", "fixtures", "park_woodburn_sln_driver_cases.jl")
```

In `test/expert/park_woodburn_route_certificate.jl`, add a helper:

```julia
function _pw_replace_column_peel_certificate(
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
```

- [ ] **Step 2: Update expert route-certificate tests**

In the existing route-certificate test set, include the SLn driver catalog and use `entries["sln-driver-sl4-gf2-ecp-mainline"]` for the public recursive success case:

```julia
if !isdefined(Main, :ParkWoodburnSLnDriverFixtureCatalog)
    include(PARK_WOODBURN_SLN_DRIVER_CATALOG_PATH)
end
sln_entries = ParkWoodburnSLnDriverFixtureCatalog.cases_by_id()
mainline_entry = sln_entries["sln-driver-sl4-gf2-ecp-mainline"]
mainline_cert = Suslin._polynomial_factorization_route_certificate(mainline_entry.matrix)
@test mainline_cert.route == :polynomial_column_peel
@test mainline_cert.status == :supported
@test mainline_cert.evidence.mainline_support_metadata.marker == :issue186_mainline
@test mainline_cert.evidence.mainline_support_metadata.supported
@test mainline_cert.evidence.final_route_provenance == :issue184_evidence_backed_sl3
@test Suslin._verify_polynomial_factorization_route_certificate(mainline_cert)
@test verify_factorization(mainline_entry.matrix, mainline_cert.factors)
```

For the legacy polynomial recursive fixtures already in this file, assert they no longer verify as public column-peel route evidence:

```julia
legacy_evidence = Suslin._polynomial_column_peel_certificate(recursive_supported_entry.matrix)
@test !legacy_evidence.mainline_support_metadata.supported
legacy_route = Suslin.PolynomialFactorizationRouteCertificate(
    recursive_supported_entry.matrix,
    :polynomial_column_peel,
    copy(legacy_evidence.factors),
    legacy_evidence.product,
    legacy_evidence,
    :supported,
    Suslin._polynomial_factorization_route_core_verification(
        Suslin.PolynomialFactorizationRouteCertificate(
            recursive_supported_entry.matrix,
            :polynomial_column_peel,
            copy(legacy_evidence.factors),
            legacy_evidence.product,
            legacy_evidence,
            :supported,
            nothing,
        ),
    ),
)
@test verify_factorization(legacy_route.matrix, legacy_route.factors)
@test !Suslin._verify_polynomial_factorization_route_certificate(legacy_route)
@test_throws ArgumentError Suslin._polynomial_factorization_route_certificate(
    recursive_supported_entry.matrix;
    route = :polynomial_column_peel,
)
```

Add a tampered-mainline negative control using `mainline_cert`:

```julia
bad_mainline_metadata = merge(
    mainline_cert.evidence.mainline_support_metadata,
    (; marker = :not_issue186_mainline),
)
bad_mainline_evidence = _pw_replace_column_peel_certificate(
    mainline_cert.evidence;
    mainline_support_metadata = bad_mainline_metadata,
)
bad_mainline_route = _pw_replace_certificate(mainline_cert; evidence = bad_mainline_evidence)
@test verify_factorization(bad_mainline_route.matrix, bad_mainline_route.factors)
@test !Suslin._verify_polynomial_factorization_route_certificate(bad_mainline_route)
```

- [ ] **Step 3: Update public driver tests**

In `test/public/factorization_driver_shell.jl`, replace the old recursive supported fixture assertion with the SLn driver mainline fixture:

```julia
if !isdefined(Main, :ParkWoodburnSLnDriverFixtureCatalog)
    include(PARK_WOODBURN_SLN_DRIVER_CATALOG_PATH)
end
sln_entries = ParkWoodburnSLnDriverFixtureCatalog.cases_by_id()
recursive_supported = sln_entries["sln-driver-sl4-gf2-ecp-mainline"].matrix
recursive_factors = elementary_factorization(recursive_supported)
@test verify_factorization(recursive_supported, recursive_factors)
recursive_cert = Suslin._polynomial_factorization_route_certificate(recursive_supported)
@test recursive_cert.route == :polynomial_column_peel
@test recursive_factors == recursive_cert.factors
@test recursive_cert.evidence isa Suslin.PolynomialColumnPeelCertificate
@test recursive_cert.evidence.mainline_support_metadata.marker == :issue186_mainline
@test recursive_cert.evidence.mainline_support_metadata.supported
@test recursive_cert.evidence.final_route_provenance == :issue184_evidence_backed_sl3
@test Suslin._verify_polynomial_column_peel_certificate(recursive_cert.evidence)
```

Keep the identity/disjoint local-block test unchanged so legacy explicit support is preserved.

- [ ] **Step 4: Update public Park-Woodburn acceptance tests**

In `test/public/park_woodburn_polynomial_factorization.jl`, include the SLn driver fixture and replace the old `pw-poly-recursive-column-peel-sln-block-qq` public acceptance case with `sln-driver-sl4-gf2-ecp-mainline`. Assert the same #186 marker and final-route provenance as in Step 3.

- [ ] **Step 5: Run red tests**

Run:

```bash
julia --project=. -e 'include("test/expert/park_woodburn_route_certificate.jl")'
julia --project=. -e 'include("test/public/factorization_driver_shell.jl")'
julia --project=. -e 'include("test/public/park_woodburn_polynomial_factorization.jl")'
```

Expected before implementation: at least one command fails because legacy column-peel route evidence still verifies as public route proof and/or explicit legacy public route certificates still return factors.

### Task 2: Gate Public Column-Peel Route Evidence On #186 Mainline Metadata

**Files:**
- Modify: `src/algorithm/polynomial_column_peel.jl`
- Modify: `src/algorithm/factorization.jl`

**Interfaces:**
- Consumes: `PolynomialColumnPeelCertificate.mainline_support_metadata`, `_verify_polynomial_column_peel_certificate`.
- Produces: `_polynomial_column_peel_public_mainline_supported`, staged failure metadata helpers, and public route verification that rejects non-mainline column-peel evidence.

- [ ] **Step 1: Add public mainline helper**

In `src/algorithm/polynomial_column_peel.jl`, after `_polynomial_column_peel_mainline_support_metadata_ok`, add:

```julia
function _polynomial_column_peel_public_mainline_supported(cert)::Bool
    try
        cert isa PolynomialColumnPeelCertificate || return false
        _verify_polynomial_column_peel_certificate(cert) || return false
        hasproperty(cert, :mainline_support_metadata) || return false
        metadata = cert.mainline_support_metadata
        return hasproperty(metadata, :issue_id) &&
            hasproperty(metadata, :marker) &&
            hasproperty(metadata, :supported) &&
            hasproperty(metadata, :final_route_provenance) &&
            metadata.issue_id == "#186" &&
            metadata.marker == :issue186_mainline &&
            metadata.supported == true &&
            metadata.final_route_provenance ==
                _POLYNOMIAL_COLUMN_PEEL_ISSUE184_FINAL_ROUTE_PROVENANCE
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end
```

- [ ] **Step 2: Add staged reason/message helpers**

In `src/algorithm/polynomial_column_peel.jl`, add helpers near the mainline helper:

```julia
const _POLYNOMIAL_COLUMN_PEEL_PUBLIC_NOT_APPLICABLE =
    :not_polynomial_column_peel_candidate

function _polynomial_column_peel_public_reason_code(cert)::Symbol
    try
        cert isa PolynomialColumnPeelCertificate || return :missing_final_sl3_route
        hasproperty(cert, :mainline_support_metadata) || return :missing_final_sl3_route
        reasons = cert.mainline_support_metadata.reason_codes
        :missing_ecp_peel_evidence in reasons && return :missing_ecp_evidence
        :missing_issue184_final_sl3_route in reasons && return :missing_final_sl3_route
        :final_block_not_sl3 in reasons && return :missing_final_sl3_route
    catch err
        err isa InterruptException && rethrow()
    end
    return :missing_final_sl3_route
end

function _polynomial_column_peel_public_reason_code(err::ArgumentError)::Symbol
    message = sprint(showerror, err)
    occursin("ECP", message) && return :missing_ecp_evidence
    occursin("left factors", message) && return :missing_ecp_evidence
    occursin("last column", message) && return :missing_ecp_evidence
    return :missing_final_sl3_route
end

function _polynomial_column_peel_public_staged_message(reason_code::Symbol)
    if reason_code == :missing_ecp_evidence
        return "SL_n recursive column-peel route is staged: missing verified #185/#262 ECP peel evidence"
    elseif reason_code == :missing_final_sl3_route
        return "SL_n recursive column-peel route is staged: missing verified #184/#263 final SL_3 route evidence"
    end
    return "SL_n recursive column-peel route is staged: unsupported recursive route evidence"
end

function _polynomial_column_peel_public_staged_failure_evidence(reason_code::Symbol)
    return (;
        error_type = :ArgumentError,
        reason_code,
        message = _polynomial_column_peel_public_staged_message(reason_code),
    )
end

function _polynomial_column_peel_public_non_candidate_error(err::ArgumentError)::Bool
    return occursin(
        "polynomial column-peel certificate requires at least one real peel step",
        sprint(showerror, err),
    )
end
```

- [ ] **Step 3: Expose recursive staged evidence**

Add:

```julia
function _polynomial_recursive_column_peel_public_staged_failure_evidence(A)
    try
        evidence = _polynomial_column_peel_certificate(A)
        _polynomial_column_peel_public_mainline_supported(evidence) && return nothing
        return _polynomial_column_peel_public_staged_failure_evidence(
            _polynomial_column_peel_public_reason_code(evidence),
        )
    catch err
        err isa InterruptException && rethrow()
        err isa ArgumentError || rethrow()
        _polynomial_column_peel_public_non_candidate_error(err) &&
            return _POLYNOMIAL_COLUMN_PEEL_PUBLIC_NOT_APPLICABLE
        return _polynomial_column_peel_public_staged_failure_evidence(
            _polynomial_column_peel_public_reason_code(err),
        )
    end
end
```

- [ ] **Step 4: Require public mainline evidence in route construction**

In `_polynomial_recursive_column_peel_route_certificate`, after the identity-final guard, add:

```julia
if !_polynomial_column_peel_public_mainline_supported(evidence)
    reason_code = _polynomial_column_peel_public_reason_code(evidence)
    throw(ArgumentError(_polynomial_column_peel_public_staged_message(reason_code)))
end
```

- [ ] **Step 5: Require public mainline evidence in route verification**

In `_polynomial_route_evidence_ok`, in the `_is_polynomial_column_peel_route(cert.route)` branch, add:

```julia
_polynomial_column_peel_public_mainline_supported(cert.evidence) &&
```

between `_verify_polynomial_column_peel_certificate(cert.evidence)` and the factor-sequence equality check.

### Task 3: Route Automatic `n > 3` Polynomial Inputs Through The Recursive Driver

**Files:**
- Modify: `src/algorithm/factorization.jl`

**Interfaces:**
- Consumes: `_polynomial_recursive_column_peel_public_staged_failure_evidence`, `_POLYNOMIAL_COLUMN_PEEL_PUBLIC_NOT_APPLICABLE`.
- Produces: automatic route ordering and staged route evidence that preserves recursive reason codes.

- [ ] **Step 1: Reorder automatic route selection**

In `_polynomial_factorization_route_certificate`, for `route === nothing`, move recursive-column-peel before the `:disjoint_local_blocks` attempt for `n > 3`:

```julia
if n > 3 && allow_recursive_column_peel
    recursive_staged_evidence = _polynomial_recursive_column_peel_public_staged_failure_evidence(A)
    if recursive_staged_evidence === nothing
        return _polynomial_recursive_column_peel_route_certificate(
            A;
            route_tag = :polynomial_column_peel,
        )
    elseif recursive_staged_evidence !== _POLYNOMIAL_COLUMN_PEEL_PUBLIC_NOT_APPLICABLE
        return _polynomial_staged_failure_route_certificate(
            A;
            allow_recursive_column_peel = allow_recursive_column_peel,
        )
    end
end
```

Leave the existing disjoint local-block fallback after this block.

- [ ] **Step 2: Preserve recursive staged evidence**

In `_polynomial_staged_failure_evidence`, inside the `nrows(A) > 3` branch, check recursive staged evidence before `reduce_sln_to_sl3(A)`:

```julia
if allow_recursive_column_peel
    recursive_staged_evidence =
        _polynomial_recursive_column_peel_public_staged_failure_evidence(A)
    if recursive_staged_evidence === nothing
        return (; error_type = :none, message = "")
    elseif recursive_staged_evidence !== _POLYNOMIAL_COLUMN_PEEL_PUBLIC_NOT_APPLICABLE
        return recursive_staged_evidence
    end
end
```

The identity/no-real-peel case returns `_POLYNOMIAL_COLUMN_PEEL_PUBLIC_NOT_APPLICABLE`, so the existing disjoint route and existing identity public tests remain valid.

- [ ] **Step 3: Run green tests**

Run the same three focused commands from Task 1. Expected: all exit 0.

### Task 4: Final Verification, Review, Commit, And PR

**Files:**
- Modify only files changed by Tasks 1-3 unless a test reveals a necessary adjacent update.

**Interfaces:**
- Consumes: all modified tests and implementation.
- Produces: passing focused verification, passing package test if feasible, review notes, commits, pushed branch, and PR URL.

- [ ] **Step 1: Run required verification**

Run:

```bash
julia --project=. -e 'include("test/expert/park_woodburn_route_certificate.jl")'
julia --project=. -e 'include("test/public/factorization_driver_shell.jl")'
julia --project=. -e 'include("test/public/park_woodburn_polynomial_factorization.jl")'
julia --project=. -e 'using Pkg; Pkg.test()'
```

- [ ] **Step 2: Review the branch**

Use the required code-review workflow. Fix Critical and Important findings, then rerun the focused tests affected by any fix.

- [ ] **Step 3: Commit implementation**

Commit all implementation and test changes:

```bash
git add src/algorithm/factorization.jl src/algorithm/polynomial_column_peel.jl \
    test/expert/park_woodburn_route_certificate.jl \
    test/public/factorization_driver_shell.jl \
    test/public/park_woodburn_polynomial_factorization.jl \
    docs/superpowers/plans/2026-07-03-issue-265-public-sln-recursive-route.md
git commit -m "Implement #265 public SLn recursive route"
```

- [ ] **Step 4: Push and create pull request**

Use the finishing-a-development-branch flow and choose option 2, "Push and create a Pull Request". Open the PR against `main` with a summary of the recursive route gate, staged diagnostics, and tests. Stop after PR creation.
