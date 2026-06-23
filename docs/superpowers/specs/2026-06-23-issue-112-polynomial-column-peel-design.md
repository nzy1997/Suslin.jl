# Issue 112 Polynomial Column-Peel Certificate Design

## Context

Issue 112 adds the ordinary-polynomial analogue of the existing Laurent
recursive column-peel replay. The repository now has the issue 109
Park-Woodburn ordinary-polynomial fixture catalog, the issue 110 internal
polynomial route certificate, and the issue 62 ECP column reducer exposed as
`reduce_unimodular_column(v, R)`. Public `elementary_factorization(A)` must not
route through this new path yet.

The new path is only for exact ordinary polynomial `SL_n` matrices whose last
column can be reduced by the ECP reducer and whose recursively peeled block is
eventually handled by an existing local route certificate.

## Approaches Considered

Recommended: add a focused internal `src/algorithm/polynomial_column_peel.jl`
file. This keeps the recursive peel evidence separate from the public driver,
mirrors the replay shape already used by `laurent_column_peel.jl`, and lets the
issue 110 route certificate consume the peel certificate as route-specific
evidence.

Alternative: put the implementation directly in `src/algorithm/factorization.jl`.
That would avoid a new include but would make an already broad file absorb the
recursive step verifier.

Alternative: only add an expert-test helper that composes existing APIs. That
would prove factor products for the fixtures but would not provide the reusable
certificate object that later Park-Woodburn routing issues need.

## Design

Add internal structs:

- `PolynomialColumnPeelStep`: one recursive peel step, storing dimension, input
  matrix, recorded last column, ECP left factors, left-normalized matrix,
  right-clearing factors, peeled block matrix, and the next block.
- `PolynomialColumnPeelCertificate`: full replay data, storing the original
  matrix, peel steps, final block, final route certificate, final factors,
  replayed factor sequence, product, and verification metadata.

The constructor `_polynomial_column_peel_certificate(A; final_route=nothing)`
will validate square exact ordinary-polynomial determinant-one input. It then
recurses from `A`. At each non-final step it records the last column, calls
`reduce_unimodular_column(last_column, R)`, verifies the left product sends the
recorded column to `e_d`, builds right-clearing elementary factors from the
bottom row, and extracts the upper-left `(d - 1) x (d - 1)` next block. It stops
only when `_polynomial_factorization_route_certificate(current;
allow_recursive_column_peel=false)` returns a supported final route.

The replayed factor sequence for a step is:

```text
inverse(left ECP factors), embedded recursive factors, inverse(right-clearing factors)
```

The verifier `_verify_polynomial_column_peel_certificate(cert)::Bool` fails
closed and recomputes all structural facts. It checks ordinary-polynomial
preconditions, exact step chain, every ECP left factor application, every
right-clearing factor, every next block, final route certificate verification,
factor replay equality, stored product equality, and
`verify_factorization(cert.original_matrix, cert.factors)`.

Extend the issue 110 route certificate with internal route
`:recursive_column_peel`. Its evidence is a `PolynomialColumnPeelCertificate`;
the route verifier must reject tampered peel evidence even if factors or product
are manually restored.

## Fixtures And Tests

Add `test/expert/park_woodburn_polynomial_column_peel.jl` and register it in
the `expert` group. The focused test will include the issue 109 catalog and
construct two determinant-one ordinary-polynomial inputs backed by catalog
entries:

- a direct catalog recursive column-peel entry over `GF(2)[x, y]`;
- a wrapped peel input whose final block is the catalog fast-local `SL_3`
  matrix over `QQ[X]`.

Each positive case must have at least one real peel step. The test checks that
the certificate verifier returns `true`, every step sends its recorded last
column to `e_n`, right-clearing isolates the next block, replayed factors
multiply exactly to the original matrix, and `verify_factorization(A,
certificate.factors)` returns `true`.

Negative controls create separate tampered certificates for one recorded last
column, one ECP left factor, one right-clearing factor, and one next block.
Each tampered certificate must fail verification even when its final product
field is manually restored.

## Files

- Create `src/algorithm/polynomial_column_peel.jl`: internal certificate
  structs, constructor, recursion, factor replay, and verifier.
- Modify `src/Suslin.jl`: include the focused implementation file after the
  route certificate and `SL_n` reduction helpers are available.
- Modify `src/algorithm/factorization.jl`: add the recursive column-peel route
  tag and route-specific evidence verification while keeping the public driver
  unchanged.
- Create `test/expert/park_woodburn_polynomial_column_peel.jl`: focused
  acceptance and tamper tests.
- Modify `test/runtests.jl`: register the expert test.

## Verification

Focused command:

```bash
julia --project=. -e 'include("test/expert/park_woodburn_polynomial_column_peel.jl")'
```

Required package command:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Out Of Scope

- No public `elementary_factorization(A)` routing change.
- No Quillen patching.
- No Laurent matrix handling in the polynomial certificate.
- No new exported names.

## Spec Self-Review

- The design keeps public API behavior unchanged.
- The ECP reducer is the source of the peel left factors.
- The final recursive block is handled through existing issue 110 route
  certificate paths with recursive routing disabled.
- Negative controls prove the verifier does not trust stored product metadata.
