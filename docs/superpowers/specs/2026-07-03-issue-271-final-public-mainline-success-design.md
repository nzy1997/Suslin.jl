# Issue 271 Final Public Mainline Success Design

Issue #271 turns the #270 Park-Woodburn final acceptance catalog into public
behavioral evidence. The tests must prove that `elementary_factorization(A)`
and `verify_factorization(A, factors)` succeed on the final ordinary-polynomial
mainline cases, and that the route certificates reject corrupted nested
evidence even when the stored factor product still matches the input.

## Context

There is no repository `AGENTS.md`, `CLAUDE.md`, or `GEMINI.md`; the README,
existing Superpowers specs, and current test layout are the applicable local
instructions. GitHub issue #271 has no comments. The recent dependency PR
comments for #282 and #283 are Codecov reports only, so the supplied issue body
and the merged #270 catalog are the source of task context.

The implementation already exposes the public routes needed for this issue:

- #184 evidence-backed `SL_3` routes return `:quillen_patch` certificates with
  nested SL3/Murthy/Quillen evidence.
- #185 ECP evidence is recorded in polynomial column-peel steps.
- #186 recursive `SL_n` routes return `:polynomial_column_peel` certificates
  only when every peel step verifies ECP evidence and the final `SL_3` block
  verifies #184 provenance.
- #270 records the final #187 public acceptance entries in
  `test/fixtures/park_woodburn_mainline_acceptance_cases.jl`.

## Chosen Approach

Extend the existing public and expert tests, without changing production
mathematics.

1. In `test/public/park_woodburn_polynomial_factorization.jl`, include the
   #270 final acceptance catalog and exercise its three supported public cases:
   multivariate `SL_3`, recursive `SL_n` with `n > 3`, and the README-style
   ordinary-polynomial example. For each case call `elementary_factorization`
   directly on the matrix, check `verify_factorization`, and inspect the route
   certificate returned from the input matrix.
2. Keep legacy fast-local, older Quillen, and disjoint-block checks as
   regression coverage, but do not let them stand in for the new #187
   mainline success cases.
3. In `test/public/factorization_driver_shell.jl`, add a compact shell check
   against the #270 catalog so the public driver path proves the recursive
   catalog case rather than only the older SLn driver catalog.
4. In `test/expert/park_woodburn_route_certificate.jl`, add a focused #270
   route-certificate check that corrupts nested ECP or final `SL_3` evidence
   while preserving the stored factors and product. The verifier must reject
   the forged certificate.

This is preferable to adding new fixture entries because #270 already defined
the final acceptance surface. It is also preferable to production changes
because #271 is an evidence gate and explicitly rules out new mathematical
support.

## Acceptance Details

The public `SL_3` catalog case must show the #184 route by requiring a
`:quillen_patch` certificate whose evidence is one of the evidence-backed
`SL_3` route types and whose verifier passes.

The recursive `SL_n` catalog case must show the #186 route by requiring:

- `cert.route == :polynomial_column_peel`;
- `mainline_support_metadata.issue_id == "#186"`;
- `mainline_support_metadata.marker == :issue186_mainline`;
- verified #185 ECP evidence on every peel step;
- `final_route_provenance == :issue184_evidence_backed_sl3`;
- a final `SL_3` certificate accepted by the #184 route verifier.

The README-style case must be sourced from the #270 catalog entry and use the
same public API calls as a reader would use in the README example.

## Tests

Focused verification commands:

```bash
julia --project=. -e 'include("test/public/park_woodburn_polynomial_factorization.jl")'
julia --project=. -e 'include("test/expert/park_woodburn_route_certificate.jl")'
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Self-Review

The design is test-only and scoped to the #270 positive acceptance catalog plus
route-certificate corruption gates. It does not add public APIs, expand the
supported input class, claim Laurent/ToricBuilder support, change coefficient
ring boundaries, or optimize factor counts.
