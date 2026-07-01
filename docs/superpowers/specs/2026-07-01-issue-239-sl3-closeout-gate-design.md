# Issue 239 SL3 Closeout Gate Design

## Context

Issue #239 closes parent #184 after #234-#238 by proving the public ordinary
polynomial `SL_3` driver accepts only the evidence-backed Park-Woodburn route
that replays the implemented upstream layers:

- #235 checked `SL_3` realization input context,
- #236 local-form or selected-variable witness,
- #237 ordinary Quillen local evidence provider,
- #183/#220 global Quillen patch evidence.

The current implementation already routes the narrow multivariate special-form
case through `PolynomialSL3QuillenMurthyRouteEvidence`. The closeout should make
that support boundary visible in public tests and docs without adding broad
local-form search or claiming later issues.

## Non-Interactive Decisions

This Agent Desk run is non-interactive, so the standing policy selects the
conservative option whenever the brainstorming skill would normally ask for
approval. The approved design is the narrow closeout gate below because it
matches the issue body and avoids irreversible or broad public API changes.

## Approaches Considered

1. **Recommended: acceptance/docs closeout only.** Add final public acceptance
   assertions, expert tamper controls, and README/docs wording. This proves #184
   over the already implemented #234-#238 evidence path while preserving staged
   boundaries.
2. **Broaden the driver to discover more local forms.** This would attempt to
   close more of arbitrary determinant-one multivariate `SL_3`, but it conflicts
   with the issue's out-of-scope list and risks claiming #185-#187 early.
3. **Docs-only closeout.** This would be low-risk but would not add the required
   acceptance and negative-control evidence for parent #184.

## Design

### Public Acceptance

Extend the public Park-Woodburn polynomial acceptance tests to name the #239
gate explicitly. The suite must keep the existing univariate fast-local route,
the legacy Quillen fixture regression, and the non-fixture multivariate
ordinary-polynomial `SL_3` case that returns factors through replayed #235-#237
and #220 evidence.

The non-fixture `SL_3` case will assert:

- `elementary_factorization(A)` returns factors,
- `verify_factorization(A, factors) == true`,
- the route certificate is `:quillen_patch`,
- the evidence type is `PolynomialSL3QuillenMurthyRouteEvidence`,
- the replayed context, witness, local evidence provider, adapter consumption,
  and global patch verifier all pass,
- corrupting one returned factor makes `verify_factorization(A, corrupted) ==
  false`.

The public shell tests already cover determinant-not-one and unsupported
coefficient-ring rejection before factors are returned. The final closeout will
keep or sharpen staged diagnostics for determinant-one `SL_3` inputs with no
supported local-form or ordinary Quillen evidence.

### Expert Route Evidence

Extend the route-certificate expert test with an explicit #236 witness tamper
control. Starting from the supported multivariate `SL_3` route evidence, mutate
the replayed local-form witness in the context or witness selection and assert
that `_verify_polynomial_sl3_quillen_murthy_route_evidence` and the enclosing
route certificate verifier reject the forged evidence.

Existing #238 controls already cover provider, adapter-consumption, global
patch, and forged metadata rejection; the #239 addition makes the local-form
witness boundary directly visible.

### Documentation

Update `README.md` and `docs/src/index.md` so the scope says:

- supported: exact ordinary-polynomial `3 x 3` determinant-one inputs whose
  #235 context, #236 witness, #237 local evidence, and #183/#220 global patch
  evidence replay before factors are returned;
- staged: determinant-one `SL_3` inputs without a supported local-form,
  variable-change, normality/conjugation, Murthy, or Quillen evidence path;
- out of scope: ECP, recursive `SL_n`, full Park-Woodburn public acceptance,
  Laurent/ToricBuilder mainline support, and factor-count optimization.

The wording must not claim arbitrary determinant-one multivariate `SL_3`
realization.

## Testing

Run the issue's required commands:

```bash
julia --project=. -e 'include("test/expert/park_woodburn_route_certificate.jl")'
julia --project=. -e 'include("test/public/park_woodburn_polynomial_factorization.jl")'
julia --project=. test/runtests.jl all
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Spec Self-Review

- No incomplete markers remain.
- The design is scoped to #239 acceptance/docs and does not broaden the public
  API.
- The supported/staged/out-of-scope docs boundary matches the issue body.
- The test additions are focused on observable public behavior and replay
  verification, not fixture-id equality.
