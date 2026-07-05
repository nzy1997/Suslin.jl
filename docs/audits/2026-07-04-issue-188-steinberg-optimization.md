# Issue 188 Steinberg Optimization Acceptance Audit

Date: 2026-07-04

Issue #294 is the parent-level closeout gate for #188. It records that the
optional Steinberg factor-count optimizer is accepted as an opt-in factor
sequence rewrite pass only. The optional Steinberg factor-count optimizer
(#188) is available only through
`optimize_elementary_factor_sequence(factors; rules = :safe)`. It is not
enabled by default, and #188 does not change the correctness contract of
`elementary_factorization(A)`.

The accepted rewrite set cites Park-Woodburn Section 6, "Eliminating
Redundancies", which lists the Steinberg relations for identity elementary
factors, same-position products, the two nontrivial commutators, and the
disjoint commutator identity.

## Upstream Layer Map

| Parent layer | Upstream closeout | Evidence boundary | #188 consumers |
| --- | --- | --- | --- |
| #288 fixture catalog | #288 | Ordinary-polynomial Section 6 Steinberg positive and negative catalog entries over exact field-backed rings | `steinberg-identity-removal-qq`, `steinberg-same-position-merge-qq`, `steinberg-inverse-cancellation-qq`, `steinberg-commutator-forward-qq`, `steinberg-commutator-reverse-qq`, and `steinberg-disjoint-commutator-identity-qq` |
| #289 canonical factors | #289 | Private canonical elementary-factor records for identity and single off-diagonal elementary matrices | #290 certificate validation and #291/#292 rewrite matching |
| #290 certificate replay | #290 | `SteinbergOptimizationCertificate` with exact original and optimized products, before/after counts, metrics, and verifier status | #291/#292 private rewrite certificates and #293 public optimizer certificates |
| #291 adjacent rewrites | #291 | Exact adjacent identity removal, same-position merge, and inverse cancellation | #293 `rules = :safe` adjacent pass |
| #292 commutator rewrites | #292 | Exact four-factor forward, reverse, and disjoint commutator windows with local product checks | #293 `rules = :safe` commutator pass |
| #293 public optimizer | #293 | Exported `optimize_elementary_factor_sequence` and `verify_steinberg_optimization_certificate` APIs | #294 closeout documentation and acceptance smoke |

## Accepted Safe Rewrite Set

The public safe rule set is exactly:

- `:identity_removal`
- `:same_position_merge`
- `:inverse_cancellation`
- `:commutator_forward`
- `:commutator_reverse`
- `:disjoint_commutator_identity`

Every optimized sequence is accepted only through exact product verification by
`verify_steinberg_optimization_certificate`. A certificate is accepted only
when the original and optimized products replay exactly, the summary agrees
with the replayed products and metrics, and `verification_status = true`.

## Documented Example

The compact example is the #288 `steinberg-commutator-forward-qq` sequence:

```julia
certificate = optimize_elementary_factor_sequence(factors; rules = :safe)
verify_steinberg_optimization_certificate(certificate)
```

| Metric | Original | Optimized | Delta |
| --- | --- | --- | --- |
| Factor count | 4 | 1 | -3 |
| Max monomial degree | 1 | 2 | 1 |
| Total off-diagonal monomial count | 6 | 2 | -4 |
| Applied rewrites | `:commutator_forward` | `:commutator_forward` | accepted safe rewrite |

For this documented example, `products_equal = true` and
`verification_status = true`.

## Closeout Note

#188 does not change the correctness contract of `elementary_factorization(A)`.
The public factorization call still returns the evidence-backed sequence for
`A` and remains accepted by `verify_factorization(A, factors)`. Users who want
factor-count cleanup must call `optimize_elementary_factor_sequence` explicitly
on the returned or supplied factor sequence and verify the returned certificate.

## Non-Claims

- #188 does not claim global minimum factor counts.
- #188 does not add Laurent `GL_n` or ToricBuilder support.
- Laurent/ToricBuilder mainline support remains separate from #188.
- This audit does not add algorithmic support beyond the landed #288-#293
  ordinary-polynomial Steinberg optimizer chain.

## Verification

The closeout gate is the focused optimizer/docs checks plus the full suite:

```bash
julia --project=. -e 'include("test/expert/steinberg_factor_count_optimization.jl")'
julia --project=. -e 'include("test/public/api_surface.jl")'
julia --project=. -e 'include("test/expert/documentation_smoke.jl")'
julia --project=. test/runtests.jl all
julia --project=. -e 'using Pkg; Pkg.test()'
```
