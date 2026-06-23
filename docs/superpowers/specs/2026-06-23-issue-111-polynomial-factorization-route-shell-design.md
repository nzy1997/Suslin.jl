# Issue 111 Polynomial Factorization Route Shell Design

## Context

Issue 110 added an internal replayable
`PolynomialFactorizationRouteCertificate` for ordinary-polynomial factorization
routes. The public `elementary_factorization(A)` driver still contains direct
factor-returning branches for polynomial inputs. Issue 111 moves those ordinary
polynomial branches behind the verified certificate shell, while preserving the
public API: callers still receive only a factor sequence.

The existing supported ordinary-polynomial routes are:

- `:fast_local_sl3`, selected first for the univariate local `SL_3` shape.
- `:disjoint_local_blocks`, selected next through `reduce_sln_to_sl3(A)` for
  currently supported block-local `SL_n` matrices.
- `:staged_failure`, used only to preserve clear staged errors for unsupported
  determinant-one ordinary-polynomial inputs.

Laurent behavior remains outside this change. Laurent determinant-one fallback
and Laurent `GL_n` determinant-correction errors must stay as they are.

## Approaches Considered

Recommended: add a small ordinary-polynomial route-shell helper in
`src/algorithm/factorization.jl`, and have `elementary_factorization(A)` call it
after polynomial validation. The helper builds a route certificate, verifies it,
then returns `cert.factors` only for supported certificates. If the certificate
is staged or fails verification, the helper throws the same staged
`ArgumentError` layer the previous driver would have thrown.

Alternative: inline the certificate construction directly inside
`elementary_factorization(A)`. This would satisfy the issue, but it would keep
the public driver responsible for route details and make later route additions
harder to sequence.

Alternative: return the certificate from the public API and let callers extract
factors. This would violate the issue's requirement to preserve the public
factor-returning API and would widen API surface before it is stable.

## Design

Add `_polynomial_verified_route_factors(A)` as the ordinary-polynomial shell.
It calls `_polynomial_factorization_route_certificate(A)` with automatic route
selection. If the certificate status is `:supported` and
`_verify_polynomial_factorization_route_certificate(cert)` is true, it returns
`cert.factors`. If verification fails for a supported route, it throws an
internal verification error instead of returning unchecked factors. If the
certificate is a staged failure, the helper rethrows the existing staged layer
by calling `_throw_staged_factorization_failure(A, :polynomial, nothing)`.

The public driver keeps its current top-level validation and Laurent
normalization behavior. For ordinary polynomials it still checks determinant
one first, then delegates to `_polynomial_verified_route_factors(normalized_A)`.
For Laurent inputs, the current direct local `SL_3`, determinant-one fallback,
and determinant-correction behavior remain unchanged.

The route ordering remains local `SL_3`, then block-local `SL_n`, then staged
failure, because that ordering is already encoded by
`_polynomial_factorization_route_certificate(A)`.

## Files

- Modify `src/algorithm/factorization.jl`: add the verified polynomial route
  shell and route ordinary-polynomial public dispatch through it.
- Modify `test/public/factorization_driver_shell.jl`: assert public polynomial
  supported cases return verified factors from certificates and unsupported
  cases still throw staged `ArgumentError` messages.
- Modify `test/expert/park_woodburn_route_certificate.jl`: add a focused
  negative control proving the public dispatcher refuses a corrupted supported
  route certificate.

## Verification

Focused commands:

```bash
julia --project=. -e 'include("test/public/factorization_driver_shell.jl")'
julia --project=. -e 'include("test/expert/park_woodburn_route_certificate.jl")'
```

Required public command:

```bash
julia --project=. test/runtests.jl public
```

Required package command:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

The public tests must show that supported polynomial cases still return only
factor lists satisfying `verify_factorization(A, factors) == true`, and that
unsupported polynomial cases still throw staged `ArgumentError` messages. The
expert negative control must temporarily force a supported certificate to carry
a corrupted factor sequence and confirm the public dispatcher refuses to return
it.

## Out of Scope

- No new polynomial routes.
- No Quillen routing.
- No recursive column-peel routing.
- No Laurent behavior changes.
- No public return type changes.

## Spec Self-Review

- The design preserves the public factor-returning API.
- Supported ordinary-polynomial dispatch is gated by certificate replay and
  `verify_factorization`.
- Unsupported ordinary-polynomial inputs keep staged `ArgumentError` messages.
- Laurent dispatch is intentionally left untouched.
- The negative control covers corrupted supported certificate factors at the
  public dispatcher boundary.
