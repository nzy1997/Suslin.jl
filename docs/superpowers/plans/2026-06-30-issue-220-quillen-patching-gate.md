# Issue 220 Quillen Patching Gate Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Close the #183 Quillen patching gate by making ordinary-polynomial Quillen route acceptance evidence-backed and documenting the staged boundary.

**Architecture:** Reuse the existing supplied-evidence and Murthy-adapter Quillen assemblers. Replace the automatic fixture-only polynomial Quillen route with a narrow supplied-sequence builder for supported elementary ordinary-polynomial cases, then adapt verified patches through the existing route-certificate adapter.

**Tech Stack:** Julia, Oscar.jl polynomial rings, existing Suslin Quillen certificate structs and test runner.

## Global Constraints

- Automatic Quillen route acceptance must pass through verified `QuillenLocalFactorSequenceCertificate` evidence and `assemble_quillen_patch_from_local_evidence`.
- The automatic supported boundary is exact ordinary polynomial rings, `3 x 3` elementary matrices, first generator as selected variable, second generator as the two-open cover `s, 1-s`, and trivial or supplied elementary `A(0)` base-term handling.
- Unsupported determinant-one multivariate matrices must fail before returning factors with a clear missing Quillen/local evidence or missing base-term evidence message.
- Do not implement #184 general `SL_3`, #185 ECP, #186 recursive `SL_n`, #187 full public acceptance, Laurent/ToricBuilder mainline support, or factor-count optimization.
- Tests must include positive supplied-evidence route coverage that is not accepted because of the old exact fixture-id shortcut, and negative controls for tampered sequence, substitution chain, base-term evidence, and patch certificates.
- README and docs must state #183 support precisely without claiming #184-#187 or Laurent/ToricBuilder support.

---

## File Structure

- Modify `src/algorithm/factorization.jl`: add the automatic supplied-evidence Quillen patch builder, update automatic route selection, and keep staged failures clear.
- Modify `test/expert/park_woodburn_route_certificate.jl`: add #220 expert route assertions and tamper controls.
- Modify `test/public/factorization_driver_shell.jl`: add a public non-fixture Quillen acceptance and staged missing-evidence negative control.
- Modify `test/public/park_woodburn_polynomial_factorization.jl`: assert public acceptance uses supplied-evidence Quillen route internals.
- Modify `README.md` and `docs/src/index.md`: update the current scope wording from staged #183 to supported #183 boundary while keeping later issues staged.
- `test/runtests.jl` already includes the touched test files; update only if a new test file is created.

### Task 1: Evidence-Backed Automatic Quillen Route

**Files:**
- Modify: `src/algorithm/factorization.jl`
- Modify: `test/expert/park_woodburn_route_certificate.jl`
- Modify: `test/public/factorization_driver_shell.jl`
- Modify: `test/public/park_woodburn_polynomial_factorization.jl`

**Interfaces:**
- Produces: `_polynomial_quillen_supplied_evidence_data(A)` returning either `nothing` or a named tuple with `selected_variable`, `cover_generator`, `row`, `col`, `entry`, `base_entry`, `delta_entry`, `base_term_policy`, `base_term_factors`, and `local_certificates`.
- Produces: `_polynomial_quillen_supplied_evidence_patch(A)` returning a verified `QuillenSuppliedEvidencePatchAssembly` or `nothing`.
- Consumes: `quillen_local_realization_certificate`, `quillen_local_factor_sequence_certificate`, `assemble_quillen_patch_from_local_evidence`, `_polynomial_quillen_patch_route_certificate`.

- [ ] **Step 1: Write failing expert route tests**

Add tests to `test/expert/park_woodburn_route_certificate.jl` in the existing `"Park-Woodburn polynomial route certificates"` testset:

```julia
    S = base_ring(quillen_entry.matrix)
    X, r, g = collect(gens(S))
    nonfixture_quillen = elementary_matrix(
        3,
        1,
        3,
        X * r + g + one(S),
        S,
    )
    nonfixture_quillen_cert =
        Suslin._polynomial_factorization_route_certificate(nonfixture_quillen)
    @test nonfixture_quillen_cert.route == :quillen_patch
    @test nonfixture_quillen_cert.evidence.quillen_patch isa
          Suslin.QuillenSuppliedEvidencePatchAssembly
    @test all(
        Suslin.verify_quillen_local_factor_sequence_certificate,
        nonfixture_quillen_cert.evidence.quillen_patch.local_certificates,
    )
    @test nonfixture_quillen_cert.evidence.quillen_patch.base_term_policy == :supplied
    @test nonfixture_quillen_cert.evidence.quillen_patch.base_term ==
          elementary_matrix(3, 1, 3, g + one(S), S)
    @test nonfixture_quillen_cert.evidence.quillen_patch.substitution_chain.verification.telescope_ok
    @test verify_factorization(nonfixture_quillen, nonfixture_quillen_cert.factors)
```

Also update the existing old fixture Quillen assertions so `quillen_cert.evidence.quillen_patch isa Suslin.QuillenSuppliedEvidencePatchAssembly`, and add tamper checks by rebuilding the patch with a modified first local sequence factor and a modified substitution-chain sign convention. Each tampered patch must make `verify_quillen_patch` false and `_polynomial_quillen_patch_route_adapter` reject it.

- [ ] **Step 2: Run the expert test and capture RED**

Run:

```bash
julia --project=. -e 'include("test/expert/park_woodburn_route_certificate.jl")'
```

Expected before implementation: FAIL because the non-fixture multivariate elementary matrix has no automatic Quillen/local evidence route, or because the old fixture still produces a deterministic `QuillenGlobalPatchAssembly`.

- [ ] **Step 3: Implement the automatic supplied-evidence route**

In `src/algorithm/factorization.jl`, replace the fixture-only automatic path with helpers shaped like:

```julia
function _polynomial_quillen_elementary_entry(A)
    nrows(A) == 3 && ncols(A) == 3 || return nothing
    R = base_ring(A)
    _factorization_ring_profile(R) == :polynomial || return nothing
    ring_gens = collect(gens(R))
    length(ring_gens) >= 2 || return nothing
    row = 0
    col = 0
    entry = zero(R)
    for i in 1:3, j in 1:3
        if i == j
            A[i, j] == one(R) || return nothing
        elseif A[i, j] != zero(R)
            row == 0 || return nothing
            row = i
            col = j
            entry = A[i, j]
        end
    end
    row == 0 && return nothing
    return (; row, col, entry, selected_variable = ring_gens[1], cover_generator = ring_gens[2])
end
```

Use `evaluate(entry, [gen == selected ? zero(R) : gen for gen in ring_gens])` to compute the base entry, use `delta_entry = entry - base_entry`, reject zero `delta_entry` as missing local Quillen evidence, and build two local sequence certificates over denominators `cover_generator` and `one(R) - cover_generator`.

Use `base_term_policy = :trivial` when the base entry is zero. Otherwise use `base_term_policy = :supplied` and `base_term_factors = [elementary_matrix(3, row, col, base_entry, R)]`. Call:

```julia
assemble_quillen_patch_from_local_evidence(
    A,
    selected_variable,
    local_certificates;
    exponent = 1,
    coverage_multipliers = [one(R), one(R)],
    base_term_policy = base_term_policy,
    base_term_factors = base_term_factors,
    metadata = (; source = :automatic_quillen_supplied_evidence, consumer_issue_id = "#220"),
)
```

Verify the patch before returning it. Update automatic route selection and staged-evidence checks to use `_polynomial_quillen_supplied_evidence_route_certificate(A)` instead of the old fixture function.

- [ ] **Step 4: Run the expert test and capture GREEN**

Run:

```bash
julia --project=. -e 'include("test/expert/park_woodburn_route_certificate.jl")'
```

Expected after implementation: all tests in that file pass.

- [ ] **Step 5: Add and run public acceptance coverage**

Update `test/public/factorization_driver_shell.jl` and `test/public/park_woodburn_polynomial_factorization.jl` to assert:

- the existing fixture-shaped Quillen matrix now uses `QuillenSuppliedEvidencePatchAssembly`;
- a non-fixture multivariate elementary matrix over `QQ[X,r,g]` returns factors from `elementary_factorization`;
- `verify_factorization(A, factors) == true`;
- a multivariate determinant-one non-elementary matrix still fails with a missing Quillen/local evidence message.

Run:

```bash
julia --project=. -e 'include("test/public/factorization_driver_shell.jl")'
julia --project=. -e 'include("test/public/park_woodburn_polynomial_factorization.jl")'
```

Expected: both commands pass.

- [ ] **Step 6: Commit**

```bash
git add src/algorithm/factorization.jl test/expert/park_woodburn_route_certificate.jl test/public/factorization_driver_shell.jl test/public/park_woodburn_polynomial_factorization.jl
git commit -m "feat: gate quillen route on supplied evidence"
```

### Task 2: Murthy Route Gate And Staged Base-Term Controls

**Files:**
- Modify: `test/expert/park_woodburn_route_certificate.jl`
- Modify: `test/expert/park_woodburn_quillen_route_adapter.jl`

**Interfaces:**
- Consumes: `_murthy_quillen_local_adapter`, `consume_murthy_quillen_adapters_for_patch`, `_polynomial_factorization_route_certificate(...; route = :quillen_patch, quillen_patch = patch)`.
- Produces: additional #220 acceptance/negative controls without changing public APIs.

- [ ] **Step 1: Add Murthy explicit route acceptance**

In `test/expert/park_woodburn_quillen_route_adapter.jl`, add a focused case using a denominator-one elementary `SL3LocalRealizationCertificate` for `E_12(X)` over `QQ[X]`. Convert it through `_murthy_quillen_local_adapter`, call `consume_murthy_quillen_adapters_for_patch` with `base_term_policy = :trivial`, and pass the resulting patch to `_polynomial_factorization_route_certificate(target; route = :quillen_patch, quillen_patch = result.patch)`. Assert the route verifies, the evidence patch verifies, and `verify_factorization(target, cert.factors)`.

- [ ] **Step 2: Add missing base-term evidence and tampered patch controls**

In `test/expert/park_woodburn_route_certificate.jl`, use `_polynomial_quillen_supplied_evidence_data(nonfixture_quillen)` to get the route-generated local certificates. Call `assemble_quillen_patch_from_local_evidence` with those certificates but without `base_term_policy` or `base_term_factors`; assert it throws an `ArgumentError` containing `A(0)` or `base-term evidence`. Rebuild a valid automatic patch with a modified `product` or modified `replay_metadata`; assert `verify_quillen_patch` is false and explicit route adaptation rejects it.

- [ ] **Step 3: Run focused expert tests**

Run:

```bash
julia --project=. -e 'include("test/expert/park_woodburn_route_certificate.jl")'
julia --project=. -e 'include("test/expert/park_woodburn_quillen_route_adapter.jl")'
```

Expected: both commands pass.

- [ ] **Step 4: Commit**

```bash
git add test/expert/park_woodburn_route_certificate.jl test/expert/park_woodburn_quillen_route_adapter.jl
git commit -m "test: cover quillen patching acceptance gate"
```

### Task 3: Documentation And Suite Integration

**Files:**
- Modify: `README.md`
- Modify: `docs/src/index.md`
- Modify: `test/runtests.jl` only if Task 1 or Task 2 created a new test file.

**Interfaces:**
- Produces: documented #183 boundary matching the implemented automatic route and explicit supplied/Murthy adapter route support.

- [ ] **Step 1: Update README scope**

In `README.md`, change the current-scope bullet that says Quillen automatic patching remains staged. State that #183 automatic Quillen patching is supported only for ordinary-polynomial cases with verified supplied or Murthy-adapter local evidence, exact cover replay, sequence replay, substitution-chain replay, and trivial or supplied `A(0)` base-term handling. Keep #184-#187, coefficient-ring broadening, Laurent/ToricBuilder mainline acceptance, and factor-count optimization staged.

- [ ] **Step 2: Mirror docs index wording**

Apply the same scope boundary to `docs/src/index.md`.

- [ ] **Step 3: Run docs/test focused checks**

Run:

```bash
git diff --check
julia --project=. -e 'include("test/expert/park_woodburn_route_certificate.jl")'
julia --project=. -e 'include("test/public/park_woodburn_polynomial_factorization.jl")'
```

Expected: all commands pass.

- [ ] **Step 4: Commit**

```bash
git add README.md docs/src/index.md test/runtests.jl
git commit -m "docs: state quillen patching support boundary"
```
