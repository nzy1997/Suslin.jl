# Issue 266 Recursive SLn Acceptance Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add parent-level #186 acceptance coverage and documentation for the supported recursive ordinary-polynomial `SL_n`, `n > 3`, public route without claiming #187 final acceptance.

**Architecture:** Treat the existing `:polynomial_column_peel` route plus `PolynomialColumnPeelCertificate.mainline_support_metadata` as the #186 public proof boundary. Add tests that exercise the two-step public `SL_5` route and documentation smoke checks that force README/Documenter/audit wording to distinguish supported, staged, legacy-regression, and out-of-scope boundaries.

**Tech Stack:** Julia, Oscar, `Test`, existing Suslin route-certificate structs and Park-Woodburn fixture catalogs.

## Global Constraints

- Supported #186 public scope is exact field-backed ordinary-polynomial `SL_n`, `n > 3`, whose recursive peel steps verify #185 ECP evidence and whose final `SL_3` block verifies #184 route evidence.
- Public #186 route evidence must carry `mainline_support_metadata.issue_id == "#186"`, `marker == :issue186_mainline`, `supported == true`, and `final_route_provenance == :issue184_evidence_backed_sl3`.
- Staged determinant-one `SL_n` inputs missing ECP or final `SL_3` evidence must report stable reason codes including `:missing_ecp_evidence` and `:missing_final_sl3_route` before factors are returned.
- Legacy fast-local/disjoint-block examples may still verify factor products but must not count as #186 mainline support by themselves.
- Do not close #187, add broad Laurent/ToricBuilder support, add arbitrary coefficient-ring support, or optimize factor counts.

---

### Task 1: Public Two-Step Recursive Acceptance Coverage

**Files:**
- Modify: `test/public/factorization_driver_shell.jl`
- Modify: `test/public/park_woodburn_polynomial_factorization.jl`

**Interfaces:**
- Consumes: `ParkWoodburnSLnDriverFixtureCatalog.cases_by_id()`, `elementary_factorization`, `verify_factorization`, `Suslin._polynomial_factorization_route_certificate`, `Suslin._verify_polynomial_factorization_route_certificate`, `Suslin.verify_ecp_column_reduction`.
- Produces: public acceptance checks proving that the `sln-driver-sl5-gf2-two-step` fixture returns factors through two recursive peel steps with #186 mainline metadata.

- [ ] **Step 1: Add a two-step recursive public assertion helper to `test/public/factorization_driver_shell.jl`**

Add the following helper after `_captured_error`:

```julia
function _assert_public_issue186_recursive_route(A, expected_step_count::Int)
    factors = elementary_factorization(A)
    @test verify_factorization(A, factors)
    cert = Suslin._polynomial_factorization_route_certificate(A)
    @test cert.route == :polynomial_column_peel
    @test cert.status == :supported
    @test factors == cert.factors
    @test Suslin._verify_polynomial_factorization_route_certificate(cert)
    @test cert.evidence isa Suslin.PolynomialColumnPeelCertificate
    @test length(cert.evidence.peel_steps) == expected_step_count
    @test cert.evidence.descent_metadata.descent_dimensions ==
          tuple((nrows(A):-1:3)...)
    @test cert.evidence.mainline_support_metadata.issue_id == "#186"
    @test cert.evidence.mainline_support_metadata.marker == :issue186_mainline
    @test cert.evidence.mainline_support_metadata.supported
    @test cert.evidence.mainline_support_metadata.peel_steps_ecp_verified
    @test cert.evidence.mainline_support_metadata.final_route_issue184_ok
    @test cert.evidence.final_route_provenance ==
          :issue184_evidence_backed_sl3
    @test all(step -> Suslin.verify_ecp_column_reduction(step.ecp_evidence),
        cert.evidence.peel_steps)
    @test cert.evidence.final_certificate.route == :quillen_patch
    @test !(cert.evidence.final_certificate.evidence isa
        Suslin.SL3LocalRealizationCertificate)

    corrupted = copy(factors)
    corrupted[1] =
        corrupted[1] *
        elementary_matrix(nrows(A), 1, 2, one(base_ring(A)), base_ring(A))
    @test !verify_factorization(A, corrupted)

    return cert
end
```

- [ ] **Step 2: Replace the existing single-step recursive assertions in `test/public/factorization_driver_shell.jl` with helper calls**

In the existing section after `sln_entries = ParkWoodburnSLnDriverFixtureCatalog.cases_by_id()`, replace the manual `recursive_supported` assertions with:

```julia
recursive_supported = sln_entries["sln-driver-sl4-gf2-ecp-mainline"].matrix
recursive_cert = _assert_public_issue186_recursive_route(recursive_supported, 1)
@test Suslin._verify_polynomial_column_peel_certificate(recursive_cert.evidence)

recursive_two_step = sln_entries["sln-driver-sl5-gf2-two-step"].matrix
recursive_two_step_cert = _assert_public_issue186_recursive_route(recursive_two_step, 2)
@test Suslin._verify_polynomial_column_peel_certificate(recursive_two_step_cert.evidence)
```

- [ ] **Step 3: Add a matching helper to `test/public/park_woodburn_polynomial_factorization.jl`**

Add the following helper after `_pw_acceptance_result_or_error`:

```julia
function _pw_assert_public_issue186_recursive_acceptance(A, expected_step_count::Int)
    factors, err = _pw_acceptance_result_or_error(A)
    @test err === nothing
    @test factors !== nothing
    @test verify_factorization(A, factors)
    cert = Suslin._polynomial_factorization_route_certificate(A)
    @test cert.route == :polynomial_column_peel
    @test cert.status == :supported
    @test factors == cert.factors
    @test cert.evidence isa Suslin.PolynomialColumnPeelCertificate
    @test length(cert.evidence.peel_steps) == expected_step_count
    @test cert.evidence.descent_metadata.descent_dimensions ==
          tuple((nrows(A):-1:3)...)
    @test cert.evidence.mainline_support_metadata.issue_id == "#186"
    @test cert.evidence.mainline_support_metadata.marker == :issue186_mainline
    @test cert.evidence.mainline_support_metadata.supported
    @test cert.evidence.mainline_support_metadata.peel_steps_ecp_verified
    @test cert.evidence.mainline_support_metadata.final_route_issue184_ok
    @test cert.evidence.final_route_provenance ==
          :issue184_evidence_backed_sl3
    @test all(step -> Suslin.verify_ecp_column_reduction(step.ecp_evidence),
        cert.evidence.peel_steps)
    @test cert.evidence.final_certificate.route == :quillen_patch
    @test Suslin._verify_polynomial_column_peel_certificate(cert.evidence)
    @test Suslin._verify_polynomial_factorization_route_certificate(cert)

    corrupted = copy(factors)
    corrupted[1] =
        corrupted[1] *
        elementary_matrix(nrows(A), 1, 2, one(base_ring(A)), base_ring(A))
    @test !verify_factorization(A, corrupted)

    return cert
end
```

- [ ] **Step 4: Replace the existing recursive acceptance assertions in `test/public/park_woodburn_polynomial_factorization.jl`**

Replace the manual `recursive = sln_entries["sln-driver-sl4-gf2-ecp-mainline"].matrix` block with:

```julia
recursive = sln_entries["sln-driver-sl4-gf2-ecp-mainline"].matrix
_pw_assert_public_issue186_recursive_acceptance(recursive, 1)

recursive_two_step = sln_entries["sln-driver-sl5-gf2-two-step"].matrix
_pw_assert_public_issue186_recursive_acceptance(recursive_two_step, 2)
```

- [ ] **Step 5: Run focused public tests**

Run:

```bash
julia --project=. -e 'include("test/public/factorization_driver_shell.jl")'
julia --project=. -e 'include("test/public/park_woodburn_polynomial_factorization.jl")'
```

Expected: both commands exit 0. These are test-only acceptance additions for already-implemented #265 behavior, so no production-code red phase is required.

- [ ] **Step 6: Commit Task 1**

```bash
git add test/public/factorization_driver_shell.jl test/public/park_woodburn_polynomial_factorization.jl
git commit -m "test: cover issue 186 public recursive SLn acceptance"
```

---

### Task 2: Boundary Documentation Red-Green And Audit Note

**Files:**
- Modify: `test/expert/polynomial_normality_support_boundary.jl`
- Modify: `README.md`
- Modify: `docs/src/index.md`
- Create: `docs/audits/2026-07-03-issue-186-recursive-sln-acceptance.md`

**Interfaces:**
- Consumes: repo docs as plain text, existing expert documentation smoke style.
- Produces: failing-then-passing documentation checks for supported/staged/legacy/out-of-scope #186 wording and a parent audit note.

- [ ] **Step 1: Write failing documentation boundary tests**

In `test/expert/polynomial_normality_support_boundary.jl`, replace `_normality_boundary_required_phrases()` with:

```julia
function _normality_boundary_required_phrases()
    return [
        "ordinary-polynomial normality/conjugation certificates",
        "ordinary-polynomial ECP unimodular-column reducer (#185) is accepted",
        "recursive ordinary-polynomial `SL_n` driver (#186) is supported for exact field-backed ordinary-polynomial `SL_n`, `n > 3`, inputs",
        "whose recursive peel steps verify #185 ECP evidence",
        "whose final `SL_3` block verifies #184 route evidence",
        "public route certificate carries #186 mainline provenance",
        "determinant-one `SL_n` inputs missing ECP, final `SL_3`, variable, local-form, Quillen/Murthy, or unsupported-ring evidence remain staged",
        "`:missing_ecp_evidence`",
        "`:missing_final_sl3_route`",
        "legacy fast-local/disjoint-block examples may still verify factors but do not count as #186 mainline support by themselves",
        "full public Park-Woodburn acceptance (#187), coefficient-ring support beyond exact field-backed ordinary polynomial rings",
        "Laurent/ToricBuilder mainline acceptance, and Steinberg factor-count optimization remain staged boundaries",
    ]
end
```

Add this helper after `_normality_boundary_contains`:

```julia
function _normality_boundary_audit_path()
    return joinpath(
        POLYNOMIAL_NORMALITY_BOUNDARY_REPO_ROOT,
        "docs",
        "audits",
        "2026-07-03-issue-186-recursive-sln-acceptance.md",
    )
end
```

In the negative-phrase loop, add:

```julia
@test !_normality_boundary_contains(text, "#187 final mainline acceptance is supported")
@test !_normality_boundary_contains(text, "arbitrary Laurent `GL_n` support is complete")
@test !_normality_boundary_contains(text, "ToricBuilder support is complete")
@test !_normality_boundary_contains(text, "Steinberg factor-count optimization is supported")
```

Add a new testset before the registration testset:

```julia
@testset "issue 186 recursive SLn acceptance audit note" begin
    audit_path = _normality_boundary_audit_path()
    @test isfile(audit_path)
    audit = read(audit_path, String)
    for phrase in (
            "#260",
            "#261",
            "#262",
            "#263",
            "#264",
            "#265",
            "#266",
            "Park-Woodburn Section 3",
            "Park-Woodburn Section 4",
            "Park-Woodburn Section 5",
            ":missing_ecp_evidence",
            ":missing_final_sl3_route",
            "does not close #187")
        @test _normality_boundary_contains(audit, phrase)
    end
end
```

- [ ] **Step 2: Run the red documentation test**

Run:

```bash
julia --project=. -e 'include("test/expert/polynomial_normality_support_boundary.jl")'
```

Expected before docs are updated: FAIL because the README/docs still describe recursive `SL_n` (#186) as staged and the new audit file does not exist.

- [ ] **Step 3: Update README scope wording**

In `README.md`, replace the ECP bullet and the staged `SL_3` boundary bullet with wording that includes this exact support boundary:

```markdown
- The ordinary-polynomial ECP unimodular-column reducer (#185) is accepted for
  exact field-backed ordinary polynomial rings through
  `reduce_unimodular_column`, `ecp_column_reduction_certificate`, and replaying
  ECP certificate verifiers. The route covers the checked input context (#243),
  monicity normalization (#244), link witness extraction (#245), link-step
  replay (#246), induction/normality composition (#247), and public reducer
  dispatch (#248).
- The recursive ordinary-polynomial `SL_n` driver (#186) is supported for exact
  field-backed ordinary-polynomial `SL_n`, `n > 3`, inputs whose recursive peel
  steps verify #185 ECP evidence, whose final `SL_3` block verifies #184 route
  evidence, and whose public route certificate carries #186 mainline
  provenance. Determinant-one `SL_n` inputs missing ECP, final `SL_3`,
  variable, local-form, Quillen/Murthy, or unsupported-ring evidence remain
  staged with stable reason codes such as `:missing_ecp_evidence` and
  `:missing_final_sl3_route`. Legacy fast-local/disjoint-block examples may
  still verify factors but do not count as #186 mainline support by themselves.
- Staged ordinary-polynomial `SL_3` inputs include determinant-one matrices with
  no supported local-form, variable-change, normality/conjugation, Murthy, or
  Quillen evidence path. Outside the evidence-backed #184 and #186 slices
  above, Quillen automatic patching (#183), general `SL_3` (#184), full public
  Park-Woodburn acceptance (#187), coefficient-ring support beyond exact
  field-backed ordinary polynomial rings, arbitrary Laurent `GL_n` determinant
  correction, Laurent/ToricBuilder mainline acceptance, and Steinberg
  factor-count optimization remain staged boundaries.
```

- [ ] **Step 4: Apply the same scope wording to `docs/src/index.md`**

Mirror the README bullet wording in `docs/src/index.md` so the Documenter page and README share the same support boundary.

- [ ] **Step 5: Add the issue #186 audit note**

Create `docs/audits/2026-07-03-issue-186-recursive-sln-acceptance.md` with:

```markdown
# Issue 186 Recursive SLn Acceptance Coverage

Date: 2026-07-03

Issue #266 is the parent-level closeout gate for #186. It records that the
ordinary-polynomial recursive `SL_n` driver is accepted for the implemented
Park-Woodburn path when #185 ECP peel evidence and #184 final `SL_3` route
evidence both replay, without closing #187.

## Stage Map

| Issue | Stage | Evidence boundary |
| --- | --- | --- |
| #260 | Driver catalog | Representative ordinary-polynomial `SL_n`, `n > 3`, fixtures for Park-Woodburn Section 3 recursion |
| #261 | Driver context | Checked determinant-one, exact field-backed ordinary-polynomial input context and staged reason vocabulary |
| #262 | ECP peel step | Park-Woodburn Section 4 last-column peel records verified #185 ECP certificates |
| #263 | Final route | Park-Woodburn Section 5 final `SL_3` block routes through verified #184 evidence |
| #264 | Recursive certificate | Full recursive column-peel certificate records descent, factors, nested ECP evidence, final route evidence, and #186 metadata |
| #265 | Public route | `elementary_factorization` accepts only #186-mainline recursive route certificates and reports staged boundaries otherwise |
| #266 | Parent acceptance gate | Public/expert tests and docs prove the accepted #186 boundary |

## Acceptance Evidence

- `test/expert/park_woodburn_sln_recursive_driver.jl` checks one-step and
  two-step recursive certificates, final #184 route evidence, verified #185 ECP
  peel evidence, tampered nested evidence rejection, tampered factor rejection,
  and tampered #186 provenance rejection.
- `test/expert/park_woodburn_route_certificate.jl` checks public route
  certificate verification, legacy route rejection, staged reason-code mapping,
  and forged provenance rejection.
- `test/public/factorization_driver_shell.jl` and
  `test/public/park_woodburn_polynomial_factorization.jl` check that
  `elementary_factorization(A)` returns factors for representative `SL_4` and
  two-step `SL_5` #186 mainline fixtures, and that corrupted returned factors
  fail `verify_factorization`.

## Staged Boundaries

- `:missing_ecp_evidence` covers determinant-one recursive candidates whose
  peel step cannot replay verified #185 ECP evidence.
- `:missing_final_sl3_route` covers determinant-one recursive candidates whose
  final `SL_3` block lacks verified #184 route evidence.
- Determinant-not-one and unsupported coefficient-ring inputs fail before
  factors are returned.
- Legacy fast-local or disjoint-block examples can remain regression fixtures,
  but they do not count as #186 mainline support by themselves.

## Non-Claims

- This gate does not close #187 final public Park-Woodburn acceptance.
- It does not add arbitrary Laurent `GL_n`, ToricBuilder, or unsupported
  coefficient-ring support.
- It does not optimize Steinberg factor counts (#188).
```

- [ ] **Step 6: Run the green documentation test**

Run:

```bash
julia --project=. -e 'include("test/expert/polynomial_normality_support_boundary.jl")'
```

Expected after docs/audit updates: exit 0.

- [ ] **Step 7: Commit Task 2**

```bash
git add test/expert/polynomial_normality_support_boundary.jl README.md docs/src/index.md docs/audits/2026-07-03-issue-186-recursive-sln-acceptance.md
git commit -m "docs: record issue 186 recursive SLn acceptance boundary"
```
