# Issue 266 Recursive SLn Acceptance Gate Design

Issue #266 closes parent #186 at the acceptance and documentation layer. The
recursive ordinary-polynomial `SL_n`, `n > 3`, path is already public after
#265 when the nested column-peel certificate verifies #185 ECP evidence at each
peel step and #184 route evidence for the final `SL_3` block. This gate proves
that boundary without claiming #187 final public Park-Woodburn acceptance.

## Context

There is no repository `AGENTS.md`; the README test instructions apply. Live
GitHub issue fetching is unavailable in this Agent Desk sandbox because
`gh issue view` cannot reach the configured proxy, so the supplied #266 issue
body, local merge history for #260-#265, and repository docs/tests are the
source of truth.

The current implementation already exposes:

- `PolynomialColumnPeelCertificate.mainline_support_metadata` with
  `issue_id == "#186"` and `marker == :issue186_mainline` for verified
  recursive support.
- Public route certificates that reject legacy column-peel evidence lacking
  #186 mainline metadata, even when factor multiplication succeeds.
- Stable staged reason codes such as `:missing_ecp_evidence`,
  `:missing_final_sl3_route`, `:determinant_not_one`, and
  `:unsupported_coefficient_ring`.

## Chosen Approach

Use a focused parent gate instead of changing the public API:

1. Strengthen public acceptance coverage so at least one multivariate
   ordinary-polynomial `SL_5` case performs two real recursive peel steps
   through `elementary_factorization`.
2. Keep expert coverage on the certificate verifier, nested #185 ECP evidence,
   final #184 route evidence, tampered factors, tampered nested evidence, and
   tampered #186 provenance.
3. Add documentation smoke coverage for the new support boundary and update
   README/Documenter scope wording.
4. Add an audit note mapping #260-#266 to the Park-Woodburn Section 3/4/5
   dependency boundary.

This is preferable to adding a dedicated route tag because the existing
`:polynomial_column_peel` route plus #186 metadata is already the public proof
boundary. It is also preferable to broadening production code because #266 is a
closeout gate and explicitly leaves #187, Laurent/ToricBuilder, arbitrary
coefficient rings, and Steinberg optimization out of scope.

## Acceptance Surface

The gate covers these calls and evidence objects:

- `elementary_factorization(A)` returns factors for representative
  exact-field-backed ordinary-polynomial `SL_n`, `n > 3`, inputs whose recursive
  proof is #186-mainline supported.
- `verify_factorization(A, factors)` confirms exact multiplication, including a
  negative control with one corrupted returned factor.
- `_polynomial_factorization_route_certificate(A)` exposes
  `:polynomial_column_peel` route evidence with #186 mainline metadata.
- `_verify_polynomial_factorization_route_certificate` and
  `_verify_polynomial_column_peel_certificate` reject tampered factors,
  corrupted peel evidence, corrupted final `SL_3` evidence, and forged #186
  provenance.
- Staged public route certificates return no factors for determinant-not-one,
  unsupported coefficient-ring, missing-ECP, and missing-final-route inputs,
  with stable reason codes.

## Documentation Boundary

README and `docs/src/index.md` should distinguish:

- Supported: exact field-backed ordinary-polynomial `SL_n`, `n > 3`, inputs
  whose recursive peel steps verify #185 ECP evidence and whose final `SL_3`
  block verifies #184 route evidence, with public #186 mainline provenance.
- Staged: determinant-one `SL_n` candidates missing ECP, final `SL_3`,
  variable, local-form, Quillen/Murthy, or unsupported-ring evidence.
- Legacy regression: older fast-local or disjoint-block examples may still
  verify factors, but they do not count as #186 mainline support by themselves.
- Out of scope: #187 final mainline acceptance, Laurent/ToricBuilder support,
  arbitrary coefficient rings, and Steinberg factor-count optimization.

## Tests

Focused verification commands:

```bash
julia --project=. -e 'include("test/expert/park_woodburn_sln_recursive_driver.jl")'
julia --project=. -e 'include("test/expert/park_woodburn_route_certificate.jl")'
julia --project=. -e 'include("test/public/park_woodburn_polynomial_factorization.jl")'
julia --project=. test/runtests.jl all
julia --project=. -e 'using Pkg; Pkg.test()'
```

The implementation plan should add one documentation smoke test that fails
before docs are updated, then update docs/audit text to pass it. Existing #265
tests already exercise much of the core verifier behavior; this issue should
avoid duplicating production code unless the new parent gate exposes a real
gap.

## Self-Review

The design is scoped to acceptance tests and docs unless a failing test reveals
a production defect. It does not add public APIs, change certificate
structures, close #187, or claim broad Laurent, ToricBuilder, arbitrary
coefficient-ring, or factor-count support.
