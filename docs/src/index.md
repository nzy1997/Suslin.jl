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
- The ordinary-polynomial ECP unimodular-column reducer (#185) is accepted for
  exact field-backed ordinary polynomial rings through
  `reduce_unimodular_column`, `ecp_column_reduction_certificate`, and replaying
  ECP certificate verifiers. The route covers the checked input context (#243),
  monicity normalization (#244), link witness extraction (#245), link-step
  replay (#246), induction/normality composition (#247), and public reducer
  dispatch (#248). Polynomial column peel records the verified ECP certificate
  used for each last-column peel step.
- The recursive ordinary-polynomial `SL_n` driver (#186) is supported for exact
  field-backed ordinary-polynomial `SL_n`, `n > 3`, inputs whose recursive peel
  steps verify #185 ECP evidence, whose final `SL_3` block verifies #184 route
  evidence, and whose public route certificate carries #186 mainline
  provenance; determinant-one `SL_n` inputs missing ECP, final `SL_3`,
  variable, local-form, Quillen/Murthy, or unsupported-ring evidence remain
  staged with stable reason codes such as `:missing_ecp_evidence` and
  `:missing_final_sl3_route`; legacy fast-local/disjoint-block examples may
  still verify factors but do not count as #186 mainline support by themselves.
- Staged ordinary-polynomial `SL_3` inputs include determinant-one matrices with
  no supported local-form, variable-change, normality/conjugation, Murthy, or
  Quillen evidence path. Outside the evidence-backed #184 and #186 slices
  above, Quillen automatic patching (#183), general `SL_3` (#184), recursive
  `SL_n` (#186), full public Park-Woodburn acceptance (#187),
  coefficient-ring support beyond exact field-backed ordinary polynomial rings,
  arbitrary Laurent `GL_n` determinant correction, Laurent/ToricBuilder
  mainline acceptance, and Steinberg factor-count optimization remain staged
  boundaries.

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
