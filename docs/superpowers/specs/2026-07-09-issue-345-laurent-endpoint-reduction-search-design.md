# Issue 345 Laurent Endpoint-Reduction Search Design

## Goal

Add an expert-only bounded search report for the certified `case_008 d=14`
Laurent endpoint-reduction context. The report must rebuild the issue 341
context from replay-certified link-witness data, record the exact small search
space, and return either replay-verified endpoint-reduction candidates or a
structured exhausted boundary without promoting production Laurent ECP support.

## Approaches Considered

1. Add a narrow expert-only search report in
   `test/expert/case008_d14_laurent_endpoint_reduction_search.jl`.
   This follows the recent d14 expert pattern, consumes the cached endpoint
   context helper from issue 341, keeps all endpoint-search logic test-only, and
   records a deterministic bounded result for issue 346 to consume later.

2. Promote endpoint-reduction replay helpers into `src/algorithm/column_reduction.jl`.
   This is premature because issue 345 is evidence-only and explicitly out of
   scope for production endpoint-reduction support.

3. Read endpoint metadata from the certified diagnostic detail added for issue
   340. This is rejected because the issue requires replay-derived source of
   truth from the issue 341 context, not copied diagnostic details.

The selected approach is option 1.

## Search Report

Add `case008_d14_laurent_endpoint_reduction_search_report(...)::NamedTuple` in
the new expert file. It will include
`case008_d14_laurent_endpoint_reduction_context.jl`, obtain the cached replay
source with `_case008_d14_endpoint_reduction_replay_source`, validate the
context, and build source and target endpoint columns from the certified
link-witness replay.

The bounded search will use one paired Laurent elementary endpoint-addition
family:

- endpoint indices: `(10,)`;
- source indices: `(1,)`, matching the certified link-witness partner;
- exponent bounds: `((-1, -1), (1, 1))`;
- exponent vectors: all nine pairs in that box;
- coefficient family: `(1,)`;
- operation family: `(:paired_laurent_endpoint_entry_addition,)`.

Each endpoint operation wraps a Laurent `:entry_addition` operation targeting
the endpoint index. The same operation is replayed against the source endpoint
column and the target endpoint column. A candidate is accepted only when replay
is valid and the operation strictly improves both endpoint column measures.

The report will expose the issue fields:

- `case_id`, `dimension`, `ring_generators`, `source_boundary`,
  `context_status`, `boundary`, `pivot_index`, `partner_index`,
  `witness_exponent`;
- top-level `source_endpoint` and `target_endpoint` from the replay-derived
  context;
- `required_endpoint_reduction_fields`;
- `endpoint_indices`, `source_indices`, `operation_families`,
  `operation_semantics`, `exponent_bounds`, `exponent_vectors`,
  `coefficient_family`, and `checked_candidate_count`;
- `status`, `candidate_count`, `replay_verified_count`, `next_boundary`, and
  `candidates`.

If candidates are found, `next_boundary` is
`:laurent_endpoint_reduction_certificate`. If the bounded space is exhausted,
`next_boundary` is `:laurent_endpoint_reduction_search_expansion`.

## Validation

Add `validate_case008_d14_laurent_endpoint_reduction_search_report(report,
fixture = ...)::Symbol`. The validator will recompute the context and replayed
endpoint columns, reject missing or stale report fields, check the exact bounded
search parameters and checked count, and replay every candidate from its
operation instead of trusting copied endpoint metadata.

Candidate validation will return stable rejection symbols for malformed
endpoint operations, wrong endpoint indices, wrong ring generators, stale
candidate source endpoint metadata, stale candidate target endpoint metadata,
and operations that replay but do not satisfy the endpoint-reduction criterion.

## Tests

The new expert test will assert the fixed issue 345 fields, the exact bounded
space, the conditional candidate/exhausted boundary, and registration in
`test/runtests.jl`.

Negative controls will verify rejection of:

- swapped ring generators;
- stale top-level source endpoint metadata;
- stale top-level target endpoint metadata;
- malformed endpoint operations;
- wrong endpoint indices;
- missing required endpoint-reduction fields;
- candidate endpoint metadata that was copied instead of recomputed from the
  replay source.

`test/runtests.jl` will include the new expert file immediately after
`expert/case008_d14_laurent_endpoint_reduction_context.jl`.

## Out of Scope

No production endpoint-reduction helper, endpoint-reduction certificate,
Laurent normality/conjugation replay, recursive peel integration, determinant
normalization, public API, or full `case_008` support claim is added.
