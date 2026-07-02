# Issue 263 Polynomial Column-Peel SL3 Final Route Design

Issue #263 lets ordinary-polynomial column peel finish a size-3 final block
through the evidence-backed #184 `SL_3` route. The existing peel final-route
gate accepts only `:fast_local_sl3` and `:disjoint_local_blocks`; this change
adds a narrow `:quillen_patch` final route only when the route certificate is
the #184 `PolynomialSL3QuillenMurthyRouteEvidence` chain and its normal
certificate verifier passes.

## Context

There is no repository `AGENTS.md`, `CLAUDE.md`, or `GEMINI.md`; the README and
CI test commands apply. GitHub issue fetching is blocked by the sandbox proxy,
so the issue body supplied by Agent Desk, local merge history, and merged
predecessor specs for #238, #260, and #261 are the source of truth.

Relevant current behavior:

- `src/algorithm/factorization.jl` already returns supported
  `:quillen_patch` certificates for two distinct evidence shapes: supplied
  Quillen patch adapters and the #184 `PolynomialSL3QuillenMurthyRouteEvidence`
  route.
- `src/algorithm/polynomial_column_peel.jl` currently tries only
  `:fast_local_sl3` and `:disjoint_local_blocks` final routes and verifies that
  final route certificate with `_verify_polynomial_factorization_route_certificate`.
- #263 requires column peel to reject supplied patches, raw factor lists, and
  tampered provenance even when the final factors still multiply exactly.

## Chosen Approach

Use the existing `:quillen_patch` route tag, but add a column-peel-specific
provenance gate. A peel certificate will carry
`final_route_provenance = :issue184_evidence_backed_sl3` only when its
`:quillen_patch` final certificate has verified
`PolynomialSL3QuillenMurthyRouteEvidence`. Legacy final routes keep stable
provenance values matching their route names.

This is preferable to adding a new public route tag because the #184 path
already produces a normal `PolynomialFactorizationRouteCertificate` whose
factors, product, Quillen patch adapter, local evidence provider, context, and
witness chain are verified. The new peel layer only needs to distinguish the
trusted #184 driver path from other valid `:quillen_patch` certificates that
are not acceptable as recursive peel finals.

## Implementation

- Extend `PolynomialColumnPeelCertificate` with a trailing
  `final_route_provenance::Symbol` field and keep an outer constructor for the
  previous positional shape so older test helpers remain usable.
- Add helpers in `polynomial_column_peel.jl` to classify and verify supported
  final routes:
  - `:fast_local_sl3`
  - `:disjoint_local_blocks`
  - `:quillen_patch` only with verified #184
    `PolynomialSL3QuillenMurthyRouteEvidence`
- Allow `final_route = :quillen_patch` through input validation and include it
  in automatic final-route candidates.
- For the `:quillen_patch` candidate, obtain a normal polynomial route
  certificate without recursive column peel and accept it only when the
  evidence-backed #184 provenance helper returns
  `:issue184_evidence_backed_sl3`.
- Add a `final_route_provenance_ok` check to column-peel core verification, so
  tampering the marker fails before factors are returned.

## Tests

Extend `test/expert/park_woodburn_polynomial_column_peel.jl` with a one-step
`SL_4` wrapper around a non-fixture multivariate #184 `SL_3` special-form
matrix. The test asserts:

- automatic and explicit `:quillen_patch` final routing finish with a supported
  `:quillen_patch` final certificate;
- the final evidence is `PolynomialSL3QuillenMurthyRouteEvidence`;
- `final_route_provenance == :issue184_evidence_backed_sl3`;
- full peel factors multiply exactly to the original `SL_4` matrix.

Negative controls mutate final-route provenance, replace the #184 evidence with
the verified raw Quillen patch adapter, tamper the #236 local-form witness, and
tamper one final-route factor while leaving unrelated products where possible.
The peel verifier must reject each corrupted certificate.

Required verification commands:

```bash
julia --project=. -e 'include("test/expert/park_woodburn_polynomial_column_peel.jl")'
julia --project=. -e 'include("test/expert/park_woodburn_route_certificate.jl")'
julia --project=. -e 'using Pkg; Pkg.test()'
git diff --check
```

## Non-Claims

This issue does not implement ECP-backed public dispatch, broad `SL_n`
recursive factorization, Laurent or ToricBuilder support, #187 final
acceptance, general multivariate `SL_3` search, or Steinberg optimization.
Supplied Quillen patch adapters remain valid route certificates in
`factorization.jl`, but they are not accepted as column-peel final routes.

## Automatic Approval

This Agent Desk run is non-interactive. Under the Standing Answer Policy, this
design is approved automatically because it is the narrowest implementation
that preserves existing legacy routes, keeps all final factors under the normal
route verifier, and adds the required #184 provenance gate before recursive
peel factors are accepted.

## Self-Review

No unresolved markers remain. The design keeps public route tags unchanged,
records the requested stable provenance marker, states how supplied Quillen
patches are rejected at the peel layer, and covers the positive and negative
controls requested by the issue.
