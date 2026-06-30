# Issue 245 ECP Link Witness Extraction Design

## Context

Issue #245 follows #244. The merged #244 layer already builds an internal
`ECPMonicityNormalization` record for ordinary-polynomial ECP input contexts,
checks the selected variable, moves the selected monic entry to coordinate 1,
and verifies the exact reduction replay. The current `ecp_link_witness` stage
still requires supplied witness metadata and only verifies that metadata after
construction.

Issue #245 asks for a bounded extraction stage that turns a #244-normalized
context into checked Park-Woodburn Section 4 link evidence: residue probes, tail
reductions, resultants, Bezout identities, cover multipliers, and path points.
It must not realize elementary link factors, consume the #184 `SL_3` route, or
change the public reducer.

## Approaches Considered

1. Add a bounded exact extractor on top of `ECPMonicityNormalization`. This is
   the chosen approach. It keeps #244 normalization as the input contract, uses
   Oscar resultants and ideal membership to compute the algebraic witnesses, and
   returns a diagnostic when the bounded search cannot prove a cover.
2. Keep requiring supplied link metadata and only broaden verification. This
   preserves old behavior but fails the objective because extraction would still
   be hidden in fixtures.
3. Implement the full Park-Woodburn link realization now. This is out of scope:
   #245 stops at witness extraction and verified storage, before elementary link
   factors or the #184 route.

## Design

Add `ECPLinkWitnessExtractionFailure` beside `ECPLinkWitnessRecord`. The failure
record stores the normalized column, selected-variable data, search bounds,
attempt counts, valid resultants found, and a staged diagnostic message.

Extend `ecp_link_witness` in two compatible ways:

- keep the existing vector-and-ring entry point, including supplied metadata;
- add an `ECPMonicityNormalization` entry point that extracts from
  `normalization.normalized_column` with selected monic index 1.

When `supplied_link_witness` is present, construction still verifies and stores
the supplied equations. When it is absent, the extractor:

1. verifies the normalized column is ordinary-polynomial, unimodular, and has
   first entry monic in the selected variable;
2. generates bounded tail-combination candidates from the tail entries using
   monomial coefficients of total degree at most 1 and at most 2 nonzero terms
   by default;
3. computes each candidate `G = sum a_j v_{j+1}` and the resultant
   `Res_X(v1, G)` with Oscar;
4. obtains Bezout coefficients by asking Oscar for coordinates of the resultant
   in the ideal `(v1, G)`, then verifies `f*v1 + h*G == resultant`;
5. searches small subsets of valid candidates and asks Oscar whether their
   resultants generate the unit ideal; if so, it stores cover multipliers from
   exact ideal coordinates and verifies the coverage equation;
6. constructs path points from the cover terms so that consecutive differences
   equal `coverage_multiplier * resultant * selected_variable`;
7. builds an `ECPLinkWitnessRecord` with `metadata.source =
   :extracted_link_witness` and verifies it with exact replay before returning.

The existing verifier remains the authority for stored records. It will accept
both `:supplied_link_witness` and `:extracted_link_witness` metadata, recompute
tail reductions, resultants, Bezout identities, coverage, and path points, and
reject tampered equations. The extractor does not use fixture ids.

If no bounded subset proves a cover, `ecp_link_witness` returns
`ECPLinkWitnessExtractionFailure` rather than guessing. Downstream link-step
construction still requires a verified `ECPLinkWitnessRecord`, so no link
factors are produced from a failed extraction.

## Tests

Add `test/expert/ecp_link_witness_general.jl` and register it in the expert
group. The focused command is:

```bash
julia --project=. -e 'include("test/expert/ecp_link_witness_general.jl")'
```

The test covers:

- extraction from the #242 `ecp-mainline-qq-link-bezout` case after #244
  normalization;
- extraction from the #242 length-4 multivariate ordinary case after #244
  normalization;
- exact recomputation of each stored tail lift, resultant, Bezout identity,
  cover multiplier, and path point;
- supplied metadata still verifies through the existing constructor path;
- a staged diagnostic is returned when the bound is too small to prove a cover;
- negative controls that corrupt a resultant, a Bezout coefficient, a path
  point, and a tail-reduction lift and are rejected by the verifier or
  constructor before link factors are produced.

Required verification commands:

```bash
julia --project=. -e 'include("test/expert/ecp_link_witness_general.jl")'
julia --project=. -e 'include("test/internal/ecp_mainline_fixtures.jl")'
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Out Of Scope

This change does not materialize link elementary factors, consume the #184
`SL_3` route, broaden the public reducer, or attempt an unbounded/general
Park-Woodburn search. The first extractor is deliberately exact and bounded.

## Automatic Approval

This Agent Desk run is non-interactive. Under the Standing Answer Policy, the
design is approved automatically because it is additive, uses the #244 verified
normalization contract, proves equations exactly with Oscar, returns diagnostics
instead of guessing, and avoids all out-of-scope realization work.

## Self-Review

- No incomplete markers remain.
- The chosen design stores replayable equations rather than fixture ids.
- Failure is explicit through a diagnostic record.
- The tests map directly to the issue acceptance and negative controls.
