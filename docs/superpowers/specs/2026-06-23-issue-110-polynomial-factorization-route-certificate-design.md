# Issue 110 Polynomial Factorization Route Certificate Design

## Context

Issue 64 needs the ordinary-polynomial factorization driver to compose multiple
algorithm blocks while keeping route evidence that expert tests can replay. The
current public `elementary_factorization(A)` API returns only factor matrices,
and that behavior must not change. Existing supported polynomial paths are the
univariate local `SL_3` slice and the `SL_n` to local `SL_3` block reduction.
Issue 109 added Park-Woodburn fixture entries for these supported routes plus
staged recursive and Quillen entries that must remain out of scope here.

The existing replay building blocks are:

- `realize_sl3_local_certificate` and `verify_sl3_local_realization` for local
  `SL_3` evidence.
- `reduce_sln_to_sl3` and `verify_sln_to_sl3_reduction` for block-local
  `SL_n` evidence.
- `verify_factorization(A, factors)` for exact factor products.

## Approaches Considered

Recommended: add a non-exported certificate type and helper constructors in
`src/algorithm/factorization.jl`. This keeps the public factor-returning API
unchanged, reuses existing route verifiers, and gives expert tests direct access
through `Suslin._...` names without widening `src/Suslin.jl` exports.

Alternative: expose a public route-inspection API. This would require API
surface changes and would make an internal test-support certificate harder to
evolve before the final Park-Woodburn driver lands.

Alternative: keep certificates entirely in expert tests. That would not give
later route orchestration code a reusable internal object and would make route
metadata easier to fake.

## Design

Add an internal `PolynomialFactorizationRouteCertificate` that stores the
original matrix, route tag, factors, stored product, route-specific evidence,
and a verification record. Add
`_polynomial_factorization_route_certificate(A; route=nothing)` to build a
certificate for currently supported ordinary-polynomial routes, and
`_verify_polynomial_factorization_route_certificate(cert)::Bool` to replay it.

Supported certificate routes for this issue:

- `:fast_local_sl3`: requires a univariate ordinary-polynomial `3 x 3` local
  shape. Evidence is the `SL3LocalRealizationCertificate` returned by
  `realize_sl3_local_certificate`.
- `:disjoint_local_blocks`: requires an ordinary-polynomial `n > 3` input that
  `reduce_sln_to_sl3` can reduce. Evidence is the resulting
  `SLNToSL3Reduction`.
- `:staged_failure`: records unsupported ordinary-polynomial determinant-one
  inputs with the diagnostic error message and no factors. This gives later
  route orchestration a replayable failure record without accepting it as a
  successful factorization route.

The verifier must fail closed. It checks that the matrix is square, ordinary
polynomial, determinant one, and that the route tag is one of the supported
certificate tags. For successful routes it recomputes the factor product from
stored factors, requires it to match both `cert.product` and `cert.matrix`, and
then verifies route-specific evidence through the existing local or reduction
verifier. It also requires the stored factors to match the route evidence, so a
tampered product or a manually restored `verify_factorization(A, factors)` is
not enough to pass.

For staged failures, verification requires an empty factor list, an identity
stored product, status metadata marking the record as staged, and evidence that
records a non-empty failure message. Staged records never satisfy the successful
factor-product route checks.

## Files

- Modify `src/algorithm/factorization.jl`: add the certificate type,
  constructor helpers, route-specific builders, and verifier.
- Create `test/expert/park_woodburn_route_certificate.jl`: focused expert tests
  using the issue 109 Park-Woodburn fixture catalog.
- Modify `test/runtests.jl`: register the new expert test.

## Verification

Focused command:

```bash
julia --project=. -e 'include("test/expert/park_woodburn_route_certificate.jl")'
```

Required package command:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

The focused test must build certificates for one univariate `SL_3` fixture and
one univariate `n > 3` block-local fixture from issue 109. Each certificate must
verify, its factors must multiply exactly to the original matrix, and its
route-specific evidence must replay through the existing verifier. Separate
negative controls must tamper the route tag, one factor, the stored product,
and route-specific evidence. The verifier must reject all of them, including
the case where `verify_factorization(A, factors)` has been manually restored.

## Out of Scope

- No recursive column-peel support.
- No Quillen patch consumption.
- No public return-type change for `elementary_factorization(A)`.
- No new exported names or public API surface changes.

## Spec Self-Review

- The design keeps certificates internal and preserves the public API.
- The supported routes match issue 110 and exclude recursive and Quillen paths.
- Negative controls cover fake route metadata, factor tampering, product
  tampering, and evidence tampering.
- The verifier replays existing local and `SL_n` evidence instead of trusting
  decorative metadata.
