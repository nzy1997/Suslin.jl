# Issue 186 Recursive SLn Acceptance Coverage

Date: 2026-07-03

Issue #266 is the parent-level closeout gate for #186. It records that the
ordinary-polynomial recursive `SL_n` driver is accepted for the implemented
Park-Woodburn path when #185 ECP peel evidence and #184 final `SL_3` route
evidence both replay, without closing #187.

## Stage Map

| Issue | Stage | Evidence boundary |
| --- | --- | --- |
| #260 | Driver catalog | Representative ordinary-polynomial `SL_n`, `n > 3`, fixtures for Park-Woodburn Section 3 recursion |
| #261 | Driver context | Checked determinant-one, exact field-backed ordinary-polynomial input context and staged reason vocabulary |
| #262 | ECP peel step | Park-Woodburn Section 4 last-column peel records verified #185 ECP certificates |
| #263 | Final route | Park-Woodburn Section 5 final `SL_3` block routes through verified #184 evidence |
| #264 | Recursive certificate | Full recursive column-peel certificate records descent, factors, nested ECP evidence, final route evidence, and #186 metadata |
| #265 | Public route | `elementary_factorization` accepts only #186-mainline recursive route certificates and reports staged boundaries otherwise |
| #266 | Parent acceptance gate | Public/expert tests and docs prove the accepted #186 boundary |

## Acceptance Evidence

- `test/expert/park_woodburn_sln_recursive_driver.jl` checks one-step and
  two-step recursive certificates, final #184 route evidence, verified #185 ECP
  peel evidence, tampered nested evidence rejection, tampered factor rejection,
  and tampered #186 provenance rejection.
- `test/expert/park_woodburn_route_certificate.jl` checks public route
  certificate verification, legacy route rejection, staged reason-code mapping,
  and forged provenance rejection.
- `test/public/factorization_driver_shell.jl` and
  `test/public/park_woodburn_polynomial_factorization.jl` check that
  `elementary_factorization(A)` returns factors for representative `SL_4` and
  two-step `SL_5` #186 mainline fixtures, and that corrupted returned factors
  fail `verify_factorization`.

## Staged Boundaries

- `:missing_ecp_evidence` covers determinant-one recursive candidates whose
  peel step cannot replay verified #185 ECP evidence.
- `:missing_final_sl3_route` covers determinant-one recursive candidates whose
  final `SL_3` block lacks verified #184 route evidence.
- Determinant-not-one and unsupported coefficient-ring inputs fail before
  factors are returned.
- Legacy fast-local or disjoint-block examples can remain regression fixtures,
  but they do not count as #186 mainline support by themselves.

## Non-Claims

- This gate does not close #187 final public Park-Woodburn acceptance.
- It does not add arbitrary Laurent `GL_n`, ToricBuilder, or unsupported
  coefficient-ring support.
- It does not optimize Steinberg factor counts (#188).
