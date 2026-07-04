# Issue 187 Park-Woodburn Mainline Acceptance Audit

Date: 2026-07-04

Issue #274 is the parent-level closeout gate for #187. It records that the
final ordinary-polynomial Park-Woodburn public contract is accepted only for
exact field-backed ordinary-polynomial determinant-one `SL_3` and `SL_n`,
`n > 3`, inputs through the implemented evidence-backed route.

## Upstream Layer Map

| Parent layer | Upstream closeout | Evidence boundary | #187 consumers |
| --- | --- | --- | --- |
| #181 normality/conjugation | #195 | Cohn-type, rank-one, conjugated-elementary, and ECP nested-normality certificates | #270 catalog metadata for `pw-mainline-sl3-multivariate-issue184-qq`, `pw-mainline-readme-ordinary-polynomial-qq`, and `pw-mainline-sln-recursive-issue185-186-gf2` through the #184/#185/#186 route chain |
| #182 Murthy local `SL_3` | #212 | Local denominator and ordinary factor replay for supported `SL_3` slices | #271 public `SL_3` success coverage for `pw-mainline-sl3-multivariate-issue184-qq` and README-style public coverage for `pw-mainline-readme-ordinary-polynomial-qq` |
| #183 Quillen patching | #220 | Supplied-evidence global patch replay and base-term controls | #271 `SL_3` and README-style success cases; #272 `pw-mainline-negative-missing-sl3-quillen-evidence` |
| #184 `SL_3` mainline | #239 | Evidence-backed final `SL_3` route with #235/#236/#237/#183/#220 replay | #271 `pw-mainline-sl3-multivariate-issue184-qq`, `pw-mainline-readme-ordinary-polynomial-qq`, and recursive final-block evidence consumed by `pw-mainline-sln-recursive-issue185-186-gf2`; #272 `pw-mainline-negative-missing-sl3-local-form-evidence`, `pw-mainline-negative-missing-sl3-quillen-evidence`, `pw-mainline-negative-missing-evidence`, and `pw-mainline-negative-missing-final-sl3-evidence` |
| #185 ECP reducer | #249 | Verified ordinary-polynomial ECP peel certificates and public reducer dispatch | #271 recursive `pw-mainline-sln-recursive-issue185-186-gf2`; #272 `pw-mainline-negative-missing-ecp-evidence` |
| #186 recursive `SL_n` | #266 | Recursive column-peel route with #185 ECP steps, #184 final `SL_3` route, and #186 provenance | #271 recursive `pw-mainline-sln-recursive-issue185-186-gf2`; #272 recursive staged controls `pw-mainline-negative-missing-ecp-evidence`, `pw-mainline-negative-missing-evidence`, and `pw-mainline-negative-missing-final-sl3-evidence` |

## Final Acceptance Inventory

- #270 catalog entries accepted by #187:
  `pw-mainline-sl3-multivariate-issue184-qq`,
  `pw-mainline-sln-recursive-issue185-186-gf2`, and
  `pw-mainline-readme-ordinary-polynomial-qq`.
- #271 public success coverage proves these accepted entries through
  `test/public/park_woodburn_polynomial_factorization.jl` and
  `test/public/factorization_driver_shell.jl`.
- #272 negative controls remain part of the closeout gate:
  `pw-mainline-negative-det-not-one`,
  `pw-mainline-negative-unsupported-coefficient-ring`,
  `pw-mainline-negative-missing-sl3-local-form-evidence`,
  `pw-mainline-negative-missing-sl3-quillen-evidence`,
  `pw-mainline-negative-missing-ecp-evidence`,
  `pw-mainline-negative-missing-evidence`,
  `pw-mainline-negative-missing-final-sl3-evidence`, and
  `pw-mainline-negative-laurent-boundary`.

## Non-Claims

- Laurent/ToricBuilder mainline support remains separate from #187.
- Unsupported coefficient rings remain negative controls and out of scope.
- Steinberg factor-count optimization (#188) remains separate from #187.
- This audit does not add algorithmic support beyond the landed #181-#186
  ordinary-polynomial chain.

## Verification

The closeout gate is the full suite plus package test:

```bash
julia --project=. test/runtests.jl all
julia --project=. -e 'using Pkg; Pkg.test()'
```
