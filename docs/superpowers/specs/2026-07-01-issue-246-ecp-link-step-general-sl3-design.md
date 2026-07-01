# Issue 246 ECP Link Step General SL3 Design

## Context

Issue #246 follows the merged #245 link-witness extraction and #238/#239
ordinary-polynomial `SL_3` route work. The current ECP link-step replay can
verify supplied GF(2)/QQ fixture families, but each segment still treats the
transport as an identity `SL_2` endpoint shortcut. That leaves the #242
`ecp-mainline-sl3-route-qq` case staged even though the #184 route can verify
ordinary-polynomial `SL_3` factors for the elementary obligations needed inside
the segment transport.

There is no `AGENTS.md` in this checkout. The GitHub issue has no comments.
Relevant PR context: #256 added verified ECP link witnesses for #245, #258
connected evidence-backed `SL_3` contexts to route certificates for #238, #259
kept #184 supported only at the evidence-backed route boundary, and #252 added
the #242 ECP mainline fixture catalog with the staged `ecp-mainline-sl3-route-qq`
case.

## Approaches Considered

1. Route each endpoint-transport elementary obligation through
   `_polynomial_factorization_route_certificate`. This is the chosen approach.
   It keeps the ECP layer from implementing `SL_3` factorization, stores a
   route certificate and factor group for every embedded obligation, and lets
   the verifier check the #184 certificates and products directly.
2. Ask #184 to factor the full endpoint-transport matrix for each link
   segment. This would be a cleaner one-certificate-per-segment shape, but the
   existing route correctly stages the full transport matrices from the #242
   case. The issue asks to use the supported route where applicable, not to
   broaden #184.
3. Preserve only the existing fixture recognizers and add another family tag.
   This is rejected because the acceptance requirements explicitly disallow
   using fixture ids or GF(2)/QQ family recognition as the new support reason.

## Design

Extend `ecp_link_step_certificate` with an internal `route_mode` keyword:

- `:auto` preserves existing behavior by using the legacy fixture path for
  recognized GF(2)/QQ witnesses and the general route otherwise;
- `:legacy_fixture` forces the old identity endpoint transport and remains
  available for regression coverage;
- `:polynomial_sl3` forces the general #184-backed route and rejects unsupported
  route obligations.

Add `route_mode` to `ECPLinkStepCertificate` so replay can recompute the same
segment strategy. In the general route, a segment still derives the exact
endpoint transport from verified endpoint reductions, but it does not trust
those raw elementary factors as final evidence. Instead it builds a
`PolynomialFactorizationRouteCertificate` for each elementary transport matrix,
requires every route certificate to verify and be `:supported`, stores the
route matrix, route certificate, route factor group, and route metadata, and
uses the concatenated route-certificate factors as the segment's
`forward_factors`.

For reviewer-visible embedded-block evidence, the segment records the first
non-identity route matrix that is an embedded `SL_2` block in coordinates
`(1, 2)` as `sl2_block` and `sl2_embedding`. The verifier no longer requires
`sl2_block == identity_matrix(R, 2)` for general routed segments. It verifies
that the block has determinant one, that the embedding matches the stored
matrix, that the embedded block is one of the routed obligations, and that the
route-certificate products compose to the stored endpoint transport matrix.

The legacy verifier path remains strict: fixture segments must still have the
identity `SL_2` block, no stored `SL_3` route certificates, and the existing
endpoint certificates. This keeps the old GF(2)/QQ tests meaningful while the
new #242 route-mode test proves the same QQ mainline witness can be realized
without relying on the fixture family tag.

## Tests

Add `test/expert/ecp_link_step_general.jl` and register it in the expert group.
The focused command is:

```bash
julia --project=. -e 'include("test/expert/ecp_link_step_general.jl")'
```

The new test will:

- build the #242 `ecp-mainline-sl3-route-qq` link witness and call
  `ecp_link_step_certificate(...; route_mode = :polynomial_sl3)`;
- assert the certificate has multiple segments, no segment uses
  `:supplied_fixture_identity_sl2_endpoint_transport`, each segment stores
  supported `PolynomialFactorizationRouteCertificate`s, and at least one
  segment records a non-identity embedded `SL_2` block;
- verify each stored `SL_3` route certificate independently and verify the
  concatenated link factors send the lower endpoint to the recorded next
  endpoint and back;
- reject negative controls that corrupt a stored route factor, an embedded
  block, a path endpoint, and the link factor order.

Existing `test/expert/ecp_link_step.jl` remains as legacy fixture regression
coverage. Required verification commands:

```bash
julia --project=. -e 'include("test/expert/ecp_link_step_general.jl")'
julia --project=. -e 'include("test/expert/ecp_link_step.jl")'
julia --project=. -e 'include("test/expert/park_woodburn_route_certificate.jl")'
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Out Of Scope

This change does not assemble the full recursive ECP reducer, implement #186
matrix peeling, broaden Laurent/ToricBuilder support, or teach #184 to factor
the full endpoint-transport matrices that are currently staged.

## Automatic Approval

This Agent Desk run is non-interactive. Under the Standing Answer Policy, the
design is approved automatically because it chooses the conservative additive
route, uses the existing #184 certificate API as the authority, keeps fixture
regressions available, and rejects unsupported route obligations instead of
guessing.

## Self-Review

- No incomplete markers remain.
- The chosen route does not use fixture ids or family recognizers for the new
  #242 general test.
- The verifier checks route certificates, route products, embedded-block
  metadata, endpoint maps, and factor ordering.
- The out-of-scope recursive ECP, #186, and Laurent/ToricBuilder work remains
  excluded.
