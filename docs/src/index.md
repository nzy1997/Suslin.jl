```@meta
CurrentModule = Suslin
```

# Suslin

Constructive elementary-matrix factorizations for staged supported slices over polynomial and Laurent polynomial rings.

## Example

```julia
using Suslin, Oscar

R, (X,) = Oscar.polynomial_ring(QQ, ["X"])
A = matrix(R, [
    one(R)      one(R) + X      zero(R);
    X           one(R) + X + X^2 zero(R);
    zero(R)     zero(R)         one(R)
])

factors = elementary_factorization(A)
verify_factorization(A, factors)
```

## Scope

- `elementary_factorization(A)` is staged but now supports the evidence-backed
  ordinary-polynomial `SL_3` driver path for exact field-backed `3 x 3`
  determinant-one inputs whose #235 checked context, #236 local-form or
  variable-change witness, #237 ordinary Quillen local evidence, and #183/#220
  global patch evidence all replay before factors are returned. This includes
  the implemented Park-Woodburn `SL_3` special-form route and preserves the
  existing univariate local `SL_3`, selected `n > 3` ordinary-polynomial
  block-local/column-peel fixtures, and #183 Quillen patch gates. The #184
  closeout does not claim arbitrary determinant-one multivariate `SL_3`
  realization unless the corresponding local-form, variable-change, or
  normality/conjugation witness is implemented and replayed.
- `verify_factorization(A, factors)` checks exact multiplication against `A`.
- `laurent_gl_factorization_certificate(A)` defaults to the eager Laurent
  normalization/core certificate. With `determinant_strategy = :lazy`, it
  records the supported monomial-unit deferred determinant correction path for
  original Laurent `GL_n` inputs such as the issue #38 fixture.
- The supported ordinary-polynomial normality/conjugation certificates cover
  the completed #181 layer: Cohn-type realization certificates, rank-one
  normality certificates, and conjugated-elementary normality certificates. The
  staged ECP induction/normality adapter replays a nested
  conjugated-elementary certificate for its normality step.
- The Murthy local `SL_3` solver (#182) is supported for the proven
  ordinary/local-witness contract: ordinary factor vectors are exposed only
  when the certificate can materialize them over the base ring, while
  nontrivial local-witness cases are verified through localized
  denominator-cleared certificate replay.
- Staged ordinary-polynomial `SL_3` inputs include determinant-one matrices with
  no supported local-form, variable-change, normality/conjugation, Murthy, or
  Quillen evidence path. Outside the evidence-backed #184 slice above,
  Quillen automatic patching (#183), general `SL_3` (#184), the general ECP
  reducer (#185), recursive `SL_n` (#186), full public Park-Woodburn
  acceptance (#187), coefficient-ring support beyond exact field-backed
  ordinary polynomial rings, arbitrary Laurent `GL_n` determinant correction,
  Laurent/ToricBuilder mainline acceptance, and Steinberg factor-count
  optimization remain staged boundaries.

See [ToricBuilder Integration Contract](@ref) for the first recorded
consumer-boundary fixture contract.

## References

- Park and Woodburn, *An algorithmic proof of Suslin's stability theorem for polynomial rings*.
- Logar and Sturmfels, *Algorithms for the Quillen-Suslin theorem*.
- Fitchas and Galligo, *Nullstellensatz effectif et conjecture de Serre*.

```@index
```

```@autodocs
Modules = [Suslin]
```
