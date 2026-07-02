# Issue 264 Recursive SLn Column-Peel Certificate Design

Issue #264 closes the internal #186 ordinary-polynomial recursive `SL_n`
certificate boundary. The existing column-peel implementation already records
per-step #262 ECP evidence and can finish selected size-3 blocks through the
#263/#184 route. The missing layer is a top-level certificate that binds every
recursive peel into one replayable proof with strict dimension descent,
final-route evidence, factor replay, and #186 mainline-support metadata.

## Context

- #261 added `SLnRecursiveDriverInputContext`, which classifies supported,
  staged, and legacy recursive-driver inputs.
- #262 added `PolynomialColumnPeelStep` evidence fields and verification that
  every peel step replays through `verify_ecp_column_reduction`.
- #263 added `final_route_provenance` and rejected raw adapter-only
  `:quillen_patch` finals for column peel.
- #260 fixture catalog contains the multivariate `SL_4` and `SL_5` driver
  cases required here. Their final `SL_3` blocks can be supported by supplied
  Quillen evidence, but the current column-peel final selector rejects raw
  `PolynomialQuillenPatchRouteAdapter` evidence because it is adapter-only.

## Chosen Approach

Upgrade `PolynomialColumnPeelCertificate` in place and add one narrow SL3 route
evidence wrapper for supplied Quillen patch finals.

The upgraded certificate will keep the existing fields and constructors, and
will add:

- `descent_metadata`: ordered dimensions, expected peel count, strict
  dimension-descent status, final dimension, final-size-3 status, and route
  provenance.
- `mainline_support_metadata`: #186 support marker, whether support is proved,
  blocking reason codes when it is not, per-step ECP evidence status, final
  #184/#263 route status, replay status, and reconstruction status.

The new route evidence record will wrap an already verified supplied-Quillen
`SL_3` final route. It will not treat a raw `PolynomialQuillenPatchRouteAdapter`
as enough. Column peel can therefore accept final `SL_3` blocks backed by this
wrapper while preserving the #263 adapter-only negative control.

## Alternatives Considered

1. Create a separate successor certificate type.
   This would avoid touching `PolynomialColumnPeelCertificate`, but
   `PolynomialFactorizationRouteCertificate` already uses that type as recursive
   column-peel evidence. A successor would add conversion glue without improving
   verification.

2. Broaden generic public `SL_3` route discovery.
   This would make more matrices route automatically, but the issue explicitly
   keeps public `elementary_factorization` dispatch and broad `SL_3` acceptance
   out of scope.

3. Accept `PolynomialQuillenPatchRouteAdapter` directly as a final route.
   This would be smaller, but it would undo #263's adapter-only rejection. The
   wrapper keeps adapter replay as a subproof instead of the final proof object.

## Verification Design

Certificate verification will recompute both new metadata records. It will
reject tampering when:

- peel steps are reordered;
- a peel step is duplicated or descent stalls;
- a final route certificate is replaced or its evidence is corrupted;
- the #186 mainline-support metadata is edited;
- the replayed factor sequence is edited;
- the stored product is changed.

The #186 mainline marker is true only when all of these are true:

- the original matrix is ordinary-polynomial determinant-one with size at least
  four;
- every peel step verifies through the #262 ECP certificate boundary;
- every step drops dimension by exactly one and the final block is exactly
  `3 x 3`;
- the final certificate is a supported `:quillen_patch` route with verified
  #184/#263 SL3 route evidence, including the new supplied-Quillen wrapper;
- replayed factors match the stored factor sequence, the stored product, and
  `verify_factorization(original_matrix, factors)`.

Legacy recursive/disjoint-block certificates continue to verify as regression
coverage, but their mainline-support marker remains false.

## Tests

Add `test/expert/park_woodburn_sln_recursive_driver.jl` with:

- the #260 `sln-driver-sl4-gf2-ecp-mainline` multivariate `SL_4` case;
- the #260 `sln-driver-sl5-gf2-two-step` case with two real peel steps;
- a legacy route check proving regression support does not set #186 mainline
  support;
- negative controls for reordered steps, duplicated steps, final-route evidence
  tampering, mainline metadata tampering, embedded factor tampering, and product
  tampering.

Update `test/runtests.jl` so the new expert file is part of the expert suite.

## Scope Boundaries

This design does not change public `elementary_factorization` dispatch, add
Laurent/ToricBuilder support, update public docs, or optimize factor counts.
