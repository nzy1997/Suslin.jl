# Issue 265 Public SLn Recursive Route Design

## Context

Issue #264 added a recursive ordinary-polynomial `SL_n` column-peel certificate for
`n > 3`. That certificate already records exact ECP replay for each peel step,
the final `SL_3` route certificate, final-route provenance, and
`mainline_support_metadata` with the `#186` marker.

The public factorization route still needs to treat that metadata as the proof
boundary. A legacy column peel whose factors multiply, but whose final route is
`:fast_local_sl3`, `:disjoint_local_blocks`, or any non-`#184` route, must stay
available for internal regression tests without becoming a public `#186`
factorization certificate.

## Chosen Approach

Reuse the existing `:polynomial_column_peel` route tag. Public support is
recognized only when the nested `PolynomialColumnPeelCertificate` verifies and
its `mainline_support_metadata` has:

- `issue_id == "#186"`
- `marker == :issue186_mainline`
- `supported == true`
- `final_route_provenance == :issue184_evidence_backed_sl3`

This keeps API churn low while making the proof provenance explicit in the
route evidence, as requested by #265. The older `:recursive_column_peel` alias
continues to parse as a route tag, but public route-certificate verification
must reject any column-peel evidence that lacks the #186 mainline marker.

## Routing

For automatic ordinary-polynomial routing:

1. Keep `SL_3` fast-local and `SL_3` supplied/evidence-backed Quillen routes as
   the first choices for `n == 3`.
2. For `n > 3`, try the recursive column-peel route before legacy disjoint
   local blocks.
3. Return a recursive certificate only if it is #186-mainline supported.
4. If the recursive path is staged, preserve its staged reason instead of
   falling through to a legacy route that could hide missing ECP or missing
   final `SL_3` evidence.

Explicit `route = :disjoint_local_blocks` remains available for legacy expert
tests, but it is not accepted as #186 public proof.

## Diagnostics

Staged route evidence for ordinary-polynomial `n > 3` inputs carries stable
recursive reason information. The public staged evidence should expose:

- `reason_code = :missing_ecp_evidence` when a peel lacks verified ECP evidence.
- `reason_code = :missing_final_sl3_route` when the final `SL_3` block lacks
  accepted #184 route evidence.
- existing determinant and coefficient-ring messages for determinant and ring
  precondition failures.

The public throwing path continues to throw `ArgumentError`, but the message
includes the stable recursive reason code so tests can distinguish the staged
failure.

## Tests

Add public coverage using `test/fixtures/park_woodburn_sln_driver_cases.jl`:

- A multivariate `SL_4` or `SL_5` #260/#264 mainline fixture returns factors
  from `elementary_factorization`, verifies with `verify_factorization`, and
  selects `:polynomial_column_peel` with #186 mainline metadata.
- Legacy recursive/disjoint fixtures are not accepted as #186 public recursive
  proof even when their factor products multiply correctly.
- Tampering the nested mainline marker or final-route provenance makes
  `_verify_polynomial_factorization_route_certificate` return `false`.
- Missing final `SL_3` evidence and missing ECP evidence surface distinct
  staged reason codes/messages before public factors are returned.

## Out Of Scope

Do not add a new route tag unless tests reveal that reusing
`:polynomial_column_peel` cannot express the proof boundary. Do not broaden
Laurent, ToricBuilder, or arbitrary coefficient-ring support.
