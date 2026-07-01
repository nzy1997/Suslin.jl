# Issue 238 SL3 Evidence Route Design

## Context

Issue #238 connects the already-merged #235, #236, #237, #219, and #220
evidence layers to the public ordinary-polynomial factorization route. The
current route order accepts the cheap univariate `:fast_local_sl3` path, selected
older `n > 3` paths, and the #220 supplied-evidence `:quillen_patch` route for
elementary matrices. It still does not ask the #184 `SL_3` driver context to
produce a Murthy/Quillen provider and global patch before public
`elementary_factorization` returns factors.

The checkout has no `AGENTS.md`, `CLAUDE.md`, or `GEMINI.md`. GitHub CLI access
is blocked by the sandbox proxy, but the GitHub connector showed no issue
comments on #238 and confirmed that #237 is merged. The local source contains
the checked `SL3RealizationInputContext`, `SL3LocalFormWitnessSelection`,
`SL3MurthyQuillenLocalEvidenceProvider`, Murthy adapter consumption, and
Quillen supplied-evidence patch verification APIs needed for this issue.

## Approaches Considered

1. Reuse the existing `:quillen_patch` route tag and make the automatic route
   build its patch through the SL3 provider path. This is the chosen approach.
   The route already verifies global patch factors, exact patch replay, base
   term evidence, and route-product equality; the missing piece is metadata and
   construction that prove the patch came from #235-#237 rather than fixture
   matching.
2. Add a new `:sl3_quillen_murthy` route tag. This would make the route name
   more explicit, but it would duplicate the existing patch route verification
   unless it still wrapped the same `PolynomialQuillenPatchRouteAdapter`.
3. Return Murthy local factors or provider local sequences directly. This is
   rejected because the issue requires verified global Quillen patch assembly
   before returning factors.

## Design

Add an internal SL3 evidence route helper in `src/algorithm/factorization.jl`.
For a `3 x 3` determinant-one matrix over an exact field-backed ordinary
polynomial ring, the helper will:

1. build a checked `_sl3_realization_input_context` with selected variable
   metadata and a #236 local-form witness payload;
2. select a checked `_select_sl3_local_form_witness`;
3. construct `_sl3_murthy_quillen_local_evidence_provider` and require its
   staged diagnostic to be supported;
4. call `consume_murthy_quillen_adapters_for_patch` with the verified Murthy
   adapter, carrying route metadata that names the #235 context, #236 witness,
   #237 provider, #219 adapter consumption, and #220 patch gate;
5. adapt the resulting verified patch through
   `_polynomial_quillen_patch_route_certificate`.

The automatic route order remains conservative:

- `:fast_local_sl3` stays first for the existing univariate local path.
- `n > 3` disjoint-block and recursive column-peel behavior stays unchanged.
- The existing #220 supplied-evidence elementary route remains available.
- The new SL3 evidence route is tried for `3 x 3` ordinary-polynomial inputs
  before recursive column peel and staged failure.

The first supported slice is deliberately narrow: already-special-form
multivariate `SL_3` contexts whose local-form payload replays exactly against
the original matrix, whose provider emits ordinary Quillen local sequences, and
whose global patch verifies through #219/#220. Transformed local-form replay,
localized denominator-cleared handoff, arbitrary coordinate-change search, ECP,
recursive `SL_n`, Laurent/ToricBuilder support, and broad multivariate `SL_3`
search remain staged.

## Error Boundary

Every stage must fail before factors are returned:

- context construction failures identify determinant, coefficient ring,
  ordinary polynomial, selected-variable, or context-boundary problems;
- witness selection/provider failures identify missing #236 local-form payloads,
  local-form mismatch, non-monic selected entries, transformed local-form replay
  that is not yet globally composed, and localized Murthy handoff;
- provider and adapter verification failures identify Murthy local evidence or
  Murthy-to-Quillen adapter replay problems;
- patch consumption failures identify denominator cover, base-term,
  substitution-chain, local sequence, or global patch verification problems.

The staged failure certificate should report the earliest stable diagnostic from
these attempts. Public `elementary_factorization` continues to throw
`ArgumentError` for staged certificates, so unsupported determinant-one inputs do
not receive factors.

## Verification

Focused tests will extend:

- `test/expert/park_woodburn_route_certificate.jl` with a non-fixture
  evidence-backed SL3 route, provider/consumption metadata checks, and corrupted
  route-certificate rejection.
- `test/public/factorization_driver_shell.jl` with public acceptance for the
  same SL3 evidence route and staged negative controls for unsupported
  determinant-one, unsupported coefficient ring, missing #236 witness, missing
  #237 evidence, missing #220 patch evidence, and tampered patch evidence.
- `test/public/park_woodburn_polynomial_factorization.jl` with acceptance that
  `elementary_factorization(A)` returns factors whose product verifies for a
  #234-style multivariate `SL_3` ordinary-polynomial matrix.

Required commands:

```bash
julia --project=. -e 'include("test/expert/park_woodburn_route_certificate.jl")'
julia --project=. -e 'include("test/public/factorization_driver_shell.jl")'
julia --project=. -e 'include("test/public/park_woodburn_polynomial_factorization.jl")'
julia --project=. -e 'using Pkg; Pkg.test()'
git diff --check
```

## Automatic Approval

This Agent Desk run is non-interactive. Under the Standing Answer Policy, the
design is approved automatically because it reuses verified local and global
evidence APIs, preserves existing cheap supported routes, returns factors only
after global patch verification, and keeps unsupported general `SL_3` inputs
staged with diagnostics.

## Spec Self-Review

- No placeholder markers remain.
- The chosen route tag and replay metadata boundary are explicit.
- The design does not trust fixture ids or raw matrix matching.
- All issue-required positive and negative controls are represented.
- Out-of-scope ECP, recursive `SL_n`, Laurent/ToricBuilder, general coordinate
  search, and direct Murthy-factor return paths remain excluded.
