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
- The final ordinary-polynomial Park-Woodburn public contract (#187) is
  supported for exact field-backed ordinary-polynomial determinant-one `SL_3`
  and exact field-backed ordinary-polynomial determinant-one `SL_n`, `n > 3`,
  inputs through the implemented evidence-backed route. Public acceptance
  requires the factors returned by `elementary_factorization(A)` to satisfy
  `verify_factorization(A, factors)`, and the route certificate must replay the
  required #184 `SL_3` evidence plus #185/#186 ECP-backed recursive provenance
  for larger sizes. The example above is the README-style ordinary-polynomial
  public call shape covered by the #187 acceptance tests. The #187 closeout
  coverage audit maps the accepted public cases and upstream gates in
  `docs/audits/2026-07-04-issue-187-park-woodburn-mainline-acceptance.md`.
- `verify_factorization(A, factors)` checks exact multiplication against `A`.
- The optional Steinberg factor-count optimizer (#188) is available only
  through `optimize_elementary_factor_sequence(factors; rules = :safe)`. It is
  not enabled by default, and #188 does not change the correctness contract of
  `elementary_factorization(A)`: public factorization still returns the
  evidence-backed sequence for `A`, and callers may separately optimize that
  sequence only by making the opt-in call. The safe rewrite set is the
  Park-Woodburn Section 6 subset covered by the #288 catalog:
  `:identity_removal`, `:same_position_merge`, `:inverse_cancellation`,
  `:commutator_forward`, `:commutator_reverse`, and
  `:disjoint_commutator_identity`. Every optimized sequence is accepted only
  through exact product verification by
  `verify_steinberg_optimization_certificate`; #188 does not claim global
  minimum factor counts and does not add Laurent `GL_n` or ToricBuilder
  support. The closeout audit is
  `docs/audits/2026-07-04-issue-188-steinberg-optimization.md`.
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
- Staged ordinary-polynomial inputs include determinant-one matrices missing the
  required local-form, variable-change, normality/conjugation, Murthy,
  Quillen, ECP, or final `SL_3` evidence path. Those inputs fail before public
  factors are returned, with stable reason codes such as
  `:missing_ecp_evidence` and `:missing_final_sl3_route` for recursive staged
  failures. Unsupported coefficient rings remain out of scope. Arbitrary
  Laurent `GL_n`, ToricBuilder mainline acceptance, and Steinberg factor-count
  optimization (#188) remain separate from #187.

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
