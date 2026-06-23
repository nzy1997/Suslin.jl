# Issue 113 Polynomial Column-Peel Driver Route Design

## Context

Issue 111 routed ordinary-polynomial public factorization through the internal
verified route certificate shell. Issue 112 added replayable ordinary-polynomial
recursive column-peel certificates, but automatic public route selection still
keeps that path disabled. This issue connects those two pieces.

The existing narrow routes stay first: local univariate `SL_3`, then the
disjoint local-block `SL_n` reduction. The new public route is tried only after
those fail and before staged failure evidence is produced.

## Approaches Considered

Recommended: add `:polynomial_column_peel` as the automatic route tag and keep
the older `:recursive_column_peel` tag as an explicit compatibility alias. This
matches the issue's public route wording without breaking the issue 112 expert
tests and fixture metadata that already exercise `:recursive_column_peel`.

Alternative: rename `:recursive_column_peel` everywhere. That would make the
route vocabulary cleaner, but it would churn the catalog and expert tests that
already describe issue 112 behavior.

Alternative: only make `_polynomial_verified_route_factors(A)` pass
`allow_recursive_column_peel=true` while preserving the old route tag. That is a
smaller code change, but it does not satisfy the issue's requirement that the
route shell record `:polynomial_column_peel` evidence.

## Design

Change `_polynomial_factorization_route_certificate(A)` automatic route
selection so it tries:

1. `:fast_local_sl3`;
2. `:disjoint_local_blocks`;
3. `:polynomial_column_peel`;
4. `:staged_failure`.

The recursive route constructor will build a `PolynomialColumnPeelCertificate`
and record it under the requested route tag. Both `:polynomial_column_peel` and
`:recursive_column_peel` are accepted route tags, but automatic selection uses
`:polynomial_column_peel`.

The route verifier treats both route tags as successful polynomial column-peel
routes and verifies the nested peel certificate with
`_verify_polynomial_column_peel_certificate`. It also keeps the exact factor
sequence equality check so tampered nested evidence is rejected even when the
top-level factors still multiply to the original matrix.

The public `elementary_factorization(A)` remains factor-returning. Its
ordinary-polynomial route shell will request automatic route selection with
recursive column peel enabled. Supported catalog inputs return verified factors;
unsupported determinant-one inputs whose peel step or final recursive block is
outside the implemented witness families still throw staged `ArgumentError`s.

## Tests

Add public coverage in `test/public/factorization_driver_shell.jl` using the
Park-Woodburn polynomial catalog:

- a supported recursive column-peel case succeeds through
  `elementary_factorization(A)`;
- the automatic route certificate records `:polynomial_column_peel`;
- a nearby unsupported determinant-one recursive catalog case still throws a
  staged error through the public driver.

Extend `test/expert/park_woodburn_route_certificate.jl`:

- automatic route selection for a supported recursive input records
  `:polynomial_column_peel`;
- explicit `route=:recursive_column_peel` remains accepted for compatibility;
- corrupting the nested polynomial column-peel evidence while keeping the final
  factors unchanged makes the top-level route verifier return `false`.

## Verification

Focused commands:

```bash
julia --project=. -e 'include("test/expert/park_woodburn_route_certificate.jl")'
julia --project=. -e 'include("test/public/factorization_driver_shell.jl")'
```

Required package command:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Out Of Scope

- No multivariate Quillen routing.
- No factor-count optimization.
- No Laurent behavior changes.
- No new exported names.

## Spec Self-Review

- The design preserves existing narrower route ordering before column peel.
- The new automatic route tag matches the issue text.
- Existing explicit issue 112 route tests can keep using the compatibility tag.
- Unsupported determinant-one polynomial inputs still fail through staged errors
  when the recursive peel certificate cannot be built.
