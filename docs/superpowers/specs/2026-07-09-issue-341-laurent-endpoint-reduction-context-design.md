# Issue 341 Laurent Endpoint-Reduction Context Design

## Goal

Define an expert-only `case_008 d=14` Laurent endpoint-reduction context from replay-certified link-witness data. The context must expose the fixed witness metadata, source and target endpoint metadata, ring, dimension, required endpoint-reduction fields, and unresolved `:laurent_endpoint_reduction` boundary for later search and certificate issues.

The design must not search for endpoint reductions, implement endpoint certificates, or consume the production diagnostic detail added for issue 340.

## Approaches Considered

1. Add a narrow expert-only context in `test/expert/case008_d14_laurent_endpoint_reduction_context.jl`.
   This follows the recent d14 context/search/certificate pattern, derives the post-descent source column through the existing replay path, validates the certified link-witness certificate through internal helpers, and keeps the artifact scoped to upcoming expert endpoint work.

2. Promote a generic endpoint-context helper into `src/algorithm/column_reduction.jl`.
   This would be premature because endpoint reductions are not implemented yet and the issue asks only for a stable expert context from the certified d14 certificate.

3. Read endpoint metadata from the issue 340 diagnostic detail.
   This is explicitly out of scope. Diagnostics prove production can see the boundary, but the endpoint context needs an independent replay source of truth.

The selected approach is option 1.

## Context Object

Add `case008_d14_laurent_endpoint_reduction_context(...)` in the new expert test file. It will:

- include the existing expert link-witness certificate chain when needed;
- rebuild the d14 post-descent source data with `_case008_d14_link_witness_source_data`;
- build the d14 link-witness certificate summary with `case008_d14_laurent_link_witness_certificate_summary`;
- require `d14_status == :link_witness_certificate`;
- validate the certificate with `validate_laurent_link_witness_certificate`;
- return a NamedTuple with the exact issue fields:
  `case_id`, `dimension`, `ring_generators`, `source_boundary`, `boundary`, `link_witness_status`, `witness_family`, `pivot_index`, `partner_index`, `witness_exponent`, `source_endpoint`, `target_endpoint`, `required_endpoint_reduction_fields`, and `status`.

The source boundary will be `:case008_d14_link_witness_certificate`, and the unresolved boundary will remain `:laurent_endpoint_reduction`.

## Validation

Add `validate_case008_d14_laurent_endpoint_reduction_context(context, fixture = ...)::Symbol`.

The validator will:

- check all required context fields are present;
- require the fixed d14 identity fields from the issue body;
- recompute the post-descent source data and certified link-witness certificate;
- assert the context source and target endpoint metadata exactly match the recomputed certificate endpoints;
- assert `required_endpoint_reduction_fields == (:family, :endpoint_index, :operation, :ring_generators)`;
- reject any stale or tampered certificate data before accepting the context.

The validator returns stable symbols for negative controls, following the neighboring expert tests.

## Tests

The new expert test will assert the positive contract and negative controls required by the issue:

- stale link-witness certificate;
- swapped ring generators;
- tampered source endpoint metadata;
- tampered target endpoint metadata;
- wrong witness exponent;
- omitted required endpoint-reduction field.

It will also assert that `test/runtests.jl` includes the new expert file.

`test/runtests.jl` will include the new test immediately after `expert/laurent_link_witness_certificate.jl`, keeping the d14 endpoint context downstream of the certificate shell.

## Out of Scope

No endpoint-reduction search, endpoint-reduction certificate, Laurent normality/conjugation replay, determinant normalization, recursive peel integration, public API, or full `case_008` success claim is added.
