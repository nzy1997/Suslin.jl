# Issue 115 Quillen Public Route Design

## Context

Issue 114 added an internal `PolynomialQuillenPatchRouteAdapter` that turns a
verified constructive Quillen patch into replayable ordinary-polynomial route
evidence. The public `elementary_factorization(A)` path still treats
multivariate ordinary-polynomial matrices as staged failures, even when the
matrix is the deterministic Park-Woodburn Quillen fixture from issues 99, 105,
and 109.

This issue wires one narrow public route. It should accept only the
fixture-backed Quillen witness currently present in the catalog and continue to
reject determinant-one multivariate matrices outside that witness family with a
clear staged error naming the missing Quillen/local realizability layer.

The repository has no `AGENTS.md`, `CLAUDE.md`, `GEMINI.md`, or `CODEX.md` file
in this checkout. The Agent Desk branch is already an isolated worktree on
`agent/issue-115-wire-multivariate-quillen-routing-into-elementar-run-1`.

## Approaches Considered

Recommended: add a narrow internal fixture lookup near the polynomial route
certificate code. The lookup recognizes the `quillen-patched-substitution-witness-qq`
matrix shape over `QQ[X, r, g]`, reconstructs the same local certificates and
denominator cover used by the issue 99 and 105 fixture path, assembles a
verified Quillen patch, and then reuses the issue 114 adapter. This keeps the
public route conservative while proving the real driver can consume verified
Quillen evidence.

Alternative: expose a public API that lets callers pass a Quillen patch to
`elementary_factorization`. That would avoid matching, but it changes public API
shape before issue 64 settles certificate inspection.

Alternative: make tests monkey-patch the route certificate builder with a
supplied patch. That would exercise the adapter, but it would not wire the
ordinary public driver requested by this issue.

The chosen design is the recommended fixture lookup because it is deterministic,
fail-closed, and reuses the verified Quillen adapter rather than adding new
local-to-global logic.

## Design

Extend `_polynomial_factorization_route_certificate(A)` automatic selection so
it tries the Quillen fixture route after the univariate fast/disjoint routes and
before recursive column peel or staged failure. The route returns a normal
`PolynomialFactorizationRouteCertificate` with route `:quillen_patch`, status
`:supported`, adapter evidence, and factors copied from the verified patch.

Add private helpers in `src/algorithm/factorization.jl`:

- `_polynomial_quillen_fixture_route_certificate(A)` returns a route certificate
  or `nothing`.
- `_polynomial_quillen_fixture_patch(A)` returns a verified patch or `nothing`.
- `_polynomial_quillen_patched_substitution_witness_fixture(A)` performs the
  exact matrix/ring match for `QQ[X, r, g]`.
- `_polynomial_quillen_fixture_local_certificate(fixture; local_index)` rebuilds
  each local certificate using existing Quillen certificate constructors.

The fixture match is intentionally narrow: square `3 x 3`, ordinary polynomial
ring, coefficient ring `QQ`, generator names `X`, `r`, `g`, and exact equality
with `E_12(X + r^2*g + g + 1)`. Matrices that are determinant one and close to
that shape but not equal are not accepted.

Staged errors remain split by precondition:

- determinant not one still throws the existing determinant/unit precondition
  message;
- determinant-one multivariate ordinary-polynomial inputs outside the fixture
  route throw an `ArgumentError` mentioning the missing Quillen/local
  realizability witness.

## Catalog Metadata

Update the Park-Woodburn polynomial fixture catalog so
`quillen-patched-substitution-witness-qq` now records route `:quillen_patch` and
status `:supported`. The internal catalog validator should keep checking that
the Park-Woodburn matrix exactly matches the referenced issue 99 Quillen target
matrix and that the fixture references issues 99, 105, and 115.

## Tests

Extend `test/expert/park_woodburn_route_certificate.jl` with automatic Quillen
route assertions and a nested adapter tamper check. The untampered route must
verify, its factors must satisfy `verify_factorization`, and mutating nested
adapter evidence while keeping top-level factors valid must make the route
certificate verifier return `false`.

Update `test/expert/park_woodburn_quillen_route_adapter.jl` so automatic route
selection for the Quillen fixture now returns `:quillen_patch` rather than the
old staged failure.

Extend `test/public/factorization_driver_shell.jl` with a public multivariate
Quillen success case and a close determinant-one negative control. The success
case must call `elementary_factorization(A)` directly and verify the returned
factors. The negative control must throw a staged `ArgumentError` that names the
missing Quillen/local realizability witness.

## Verification

Focused commands:

```bash
julia --project=. -e 'include("test/expert/park_woodburn_quillen_route_adapter.jl")'
julia --project=. -e 'include("test/expert/park_woodburn_route_certificate.jl")'
julia --project=. -e 'include("test/public/factorization_driver_shell.jl")'
```

Required package command:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Automatic Approval

This Agent Desk run is non-interactive. Under the Standing Answer Policy, the
design is approved automatically because it chooses the conservative
catalog-backed route, leaves public API shape unchanged, and keeps unsupported
multivariate matrices staged.

## Spec Self-Review

- No incomplete markers remain.
- The only accepted multivariate Quillen route is the deterministic fixture.
- Unsupported determinant-one matrices get a Quillen/local realizability staged
  error rather than an unchecked fallback.
- The design reuses existing Quillen verification and the issue 114 adapter.
