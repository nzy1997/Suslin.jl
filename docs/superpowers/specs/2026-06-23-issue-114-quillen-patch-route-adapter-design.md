# Issue 114 Quillen Patch Route Adapter Design

## Context

Issue 114 connects two already-merged expert layers. Issue 105 completed the
constructive Quillen patch path through verified denominator covers, local
certificates, normalized contributions, deterministic global patch assembly,
and `verify_quillen_patch(patch)`. Issue 110 added internal ordinary-polynomial
route certificates, and issue 111 made the public polynomial factorization
driver return factors only through verified route certificates.

The adapter should be a narrow bridge. It must consume an already-built
Quillen patch and replay the Quillen verifier before exposing factors to the
polynomial route layer. It should not rebuild denominator covers, local
certificate replay, normalized contribution replay, or global patch assembly.

The repository has no `AGENTS.md`, `CLAUDE.md`, `GEMINI.md`, or `CODEX.md` file
in this checkout. The Agent Desk branch is already an isolated worktree on
`agent/issue-114-add-a-quillen-patch-route-adapter-for-polynomial-run-1`.

## Approaches Considered

Recommended: add an internal `PolynomialQuillenPatchRouteAdapter` beside the
existing polynomial route certificate code in `src/algorithm/factorization.jl`.
The adapter accepts a target matrix or correction plus a Quillen patch object,
requires `verify_quillen_patch(patch)`, copies the verified global elementary
factors, recomputes their product, stores replay metadata, and can be embedded
as route-specific evidence in `PolynomialFactorizationRouteCertificate` with a
new `:quillen_patch` route tag. This reuses both existing replay systems and
keeps the public factor-returning API unchanged.

Alternative: expose a public Quillen route inspection API now. That would make
the adapter easier to call from outside the package, but it would widen public
API surface before issue 64 decides the final certificate inspection shape.

Alternative: keep all Quillen-to-route adaptation in the new expert test. That
would prove the fixture path can be replayed, but later polynomial route work
would still lack a reusable internal evidence record.

The chosen design is the recommended internal adapter because it is the
narrowest reusable bridge and follows the existing issue 110 certificate model.

## API Surface

No new exported names.

Add internal names accessed as `Suslin.<name>` by expert tests:

- `PolynomialQuillenPatchRouteAdapter`
- `_polynomial_quillen_patch_route_adapter(target, patch)`
- `_verify_polynomial_quillen_patch_route_adapter(adapter)::Bool`

Extend `_polynomial_factorization_route_certificate(A; route=nothing)` with an
optional `quillen_patch = nothing` keyword. Existing automatic route selection
does not choose Quillen. Explicit `route = :quillen_patch` requires a supplied
patch and returns a supported route certificate whose evidence is the adapter.

## Adapter Semantics

The adapter stores:

- `target`: the ordinary-polynomial matrix to be reconstructed,
- `patch`: the supplied Quillen patch object,
- `route`: `:quillen_patch`,
- `global_elementary_factors`: copied from the verified patch,
- `product`: the exact product recomputed from those factors,
- `target_matrix`: the expected target matrix, and
- `replay_metadata`: a small record binding the patch size, selected variable,
  denominator data, local certificate count, normalized contribution count, and
  patch replay metadata.

Construction fails unless `verify_quillen_patch(patch)` is true and the
adapter product equals the target matrix. Verification fails closed and
recomputes all adapter fields from the stored patch. It also calls
`verify_quillen_patch(adapter.patch)` again, so overwriting `patch.patched_product`
or adapter product fields is not enough to pass.

For `QuillenGlobalPatchAssembly`, the adapter reads
`global_elementary_factors`, `patched_product`, `target`, `substitution_variable`,
`denominator_data`, `local_certificates`, `normalized_local_contributions`, and
`replay_metadata`. For the older `QuillenPatch` scaffold, the same adapter
reads `factors`, `product`, `target`, `substitution_variable`,
`denominator_data`, and `local_contributions`. This keeps the bridge compatible
with both #63-style patch objects without adding patch internals.

## Route Certificate Integration

Add `:quillen_patch` to `_POLYNOMIAL_FACTORIZATION_ROUTE_TAGS`. The route
certificate builder handles only explicit Quillen routing:

```julia
Suslin._polynomial_factorization_route_certificate(
    A;
    route = :quillen_patch,
    quillen_patch = patch,
)
```

The resulting certificate has status `:supported`, factors copied from the
adapter, product equal to the target matrix, and adapter evidence. Its replay
requires the adapter verifier, exact factor sequence equality with the adapter,
stored product equality, and `verify_factorization(A, cert.factors)`.

Automatic public `elementary_factorization(A)` remains unchanged and will not
dispatch multivariate inputs through Quillen until issue 64 intentionally wires
that public route.

## Test Design

Add `test/expert/park_woodburn_quillen_route_adapter.jl` and register it in the
expert group after `expert/quillen_induction_constructive.jl`, so the full
expert runner loads the constructive #105 helpers once before the adapter test.

The test reuses the #109 Park-Woodburn multivariate Quillen fixture id
`quillen-patched-substitution-witness-qq` and the #105 constructive helpers by
including `test/expert/quillen_induction_constructive.jl`. It builds a verified
constructive patch from the corresponding #99 fixture, adapts it, embeds the
adapter in a route certificate, and checks:

- `verify_quillen_patch(patch)` is true before adaptation,
- adapter verification is true,
- adapter factors multiply exactly to the Park-Woodburn target matrix,
- the route certificate verifier accepts the untampered route, and
- `verify_factorization(target, cert.factors)` is true.

Negative controls cover separate tampering points:

- one Quillen denominator multiplier,
- one local certificate factor,
- one adapted global factor, and
- a patch whose `patched_product` and `target` are overwritten to the target
  while upstream Quillen evidence remains tampered.

All must be rejected by construction or by the adapter/route verifier.

## Files

- Modify `src/algorithm/factorization.jl`: add the adapter type, adapter
  verifier, Quillen route builder, and route evidence replay branch.
- Create `test/expert/park_woodburn_quillen_route_adapter.jl`: focused expert
  positive and negative coverage.
- Modify `test/runtests.jl`: register the new expert test.
- Do not modify `src/Suslin.jl` or `test/public/api_surface.jl`.

## Verification

Focused command:

```bash
julia --project=. -e 'include("test/expert/park_woodburn_quillen_route_adapter.jl")'
```

Required package command:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Automatic Approval

This Agent Desk run is non-interactive. Under the Standing Answer Policy, the
design is approved automatically because it chooses the conservative
expert/internal adapter, reuses #99/#105/#110 artifacts, avoids public API
changes, and leaves public multivariate dispatch out of scope.

## Spec Self-Review

- No incomplete markers remain.
- The adapter only consumes verified Quillen output and does not reimplement
  Quillen patch internals.
- The route certificate integration is explicit-only and keeps public dispatch
  unchanged.
- Negative controls cover denominator, local certificate, global factor, and
  overwritten final-product false positives.
