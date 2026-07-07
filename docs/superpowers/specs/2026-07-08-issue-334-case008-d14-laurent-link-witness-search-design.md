# Issue 334 Case008 d14 Laurent Link-Witness Search Design

## Context

Issue #334 is the bounded evidence step after #329's replay-derived
post-descent `case_008 d=14` column and #333's link-witness context. The
deliverable is expert-only test/report code. It must not create a production
Laurent ECP stage, widen the search beyond radius `1`, or claim endpoint
reduction unless a real replay-verified d14 candidate is recorded.

The repository pattern for this batch is to keep staged evidence in
`test/expert/*.jl`, register the file in `test/runtests.jl`, and expose stable
NamedTuple reports plus validators from the expert file itself. The new file
will follow the #329 and #333 shape.

## Approach

Use a narrow test-only report in
`test/expert/case008_d14_laurent_link_witness_search.jl`.

The report consumes:

- `case008_d14_laurent_link_witness_context()` and
  `validate_case008_d14_laurent_link_witness_context(...) == :ok`;
- `_case008_d14_replayed_post_descent_data(...)` from the #329 expert helper
  to recompute the post-descent column by replay, not by trusting stored
  after-data.

Alternatives rejected:

- Moving the search into `src/` would create public/internal algorithm surface
  before the issue asks for it.
- Reusing the older all-pair elementary move search would scan the wrong family
  and count reverse roles under the same name.
- Widening exponent or coefficient bounds would violate the explicit bounded
  report contract.

## Report Shape

`case008_d14_laurent_link_witness_search_report()` returns a NamedTuple with
the required stable fields:

- fixed context fields: `case_id`, `dimension`, `context_status`,
  `witness_families`, `witness_semantics`, `pivot_index`, `partner_indices`,
  `exponent_radius`, `exponent_vectors`, `coefficient_family`;
- scan accounting: `checked_candidate_count == 117`, `candidate_count`,
  `replay_verified_count`;
- boundary status: `status in (:candidate_found, :exhausted)`;
- `next_boundary == :laurent_link_witness_certificate` when candidates exist,
  otherwise `:laurent_link_witness_search_expansion`;
- `candidates`, where each candidate contains `witness`, `source_endpoint`,
  and `target_endpoint`.

The candidate witness fields are exactly the #333 required witness schema:
`(:family, :pivot_index, :partner_index, :coefficient, :exponent,
:ring_generators)`.

## Replay Semantics

The only family is `:two_entry_laurent_combination`. Its replay convention is:

```text
pivot_entry + coefficient * u^a * v^b * partner_entry
```

where `(a, b)` follows the declared `ring_generators` order. This is encoded as
an entry-addition replay with `target_index = pivot_index` and
`source_index = partner_index`. A reverse-direction experiment is invalid under
this family and should fail the verifier instead of being counted.

Endpoint metadata is recomputed from replay. It records the endpoint entry
index, ring generators, term count, support bounds, leading exponent, and
measure before/after. The metadata is intentionally a compact summary, not a
new trusted copy of the Laurent polynomial entry.

## Validation

Add `verify_laurent_link_witness_candidate(column, R, candidate)::Bool` for the
test-only verifier. It rejects:

- missing/malformed witness fields;
- `pivot_index == partner_index` or out-of-range indices;
- wrong `ring_generators`;
- any witness family other than `:two_entry_laurent_combination`;
- stale `source_endpoint` or `target_endpoint` metadata;
- accidental reverse pivot/partner replay under the same family;
- synthetic non-witnesses whose replay does not satisfy the endpoint identity.

Add `validate_case008_d14_laurent_link_witness_search_report(report, context,
source_report)::Symbol` to make the report self-checking. It validates the
source context, fixed bounds, checked count, candidate/replay counts, boundary,
and each candidate through exact replay from the recomputed post-descent
column.

## Tests

The new expert file includes:

- a d14 report testset asserting all fixed issue fields and the conditional
  `:candidate_found`/`:exhausted` invariants;
- synthetic calibration with a known-good pivot-plus-shifted-partner witness;
- synthetic controls for malformed fields, equal pivot/partner indices, wrong
  ring generators, stale endpoint metadata, reversed pivot/partner replay, and
  a non-witness endpoint identity.

`test/runtests.jl` registers the file immediately after
`expert/case008_d14_laurent_link_witness_context.jl`.

## Verification

Required focused command:

```bash
julia --project=. -e 'include("test/expert/case008_d14_laurent_link_witness_search.jl")'
```

Required full command from the Agent Desk run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

The expected d14 outcome may be either a candidate-backed certificate boundary
or an exhausted search-expansion boundary, but the declared 117 candidates must
be fully accounted for either way.
