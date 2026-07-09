# Issue 346 Laurent Endpoint-Reduction Certificate Design

## Context

Issue #346 follows the replay-derived endpoint context from #341 and the bounded
endpoint-reduction search report from #345. The reusable unit for this issue is
an expert-only certificate shell: it can certify one replayable endpoint
operation, but it must not promote helpers into production code or claim real
`case_008 d=14` progress from a synthetic example.

## Design Choices

Use a test-only shell in `test/expert/laurent_endpoint_reduction_certificate.jl`.
This matches the expert-only pattern used by the Laurent descent-step and
link-witness certificate shells, keeps the public/internal API unchanged, and
leaves helper promotion to issue #347.

Build the certificate from replay, not supplied endpoint metadata. The shell
will accept a context, an endpoint operation, the source column, the linked
target column, and the ring. It will replay the operation on both source and
target columns, recompute endpoint metadata from the replayed columns, and
return a NamedTuple with stable fields:
`case_id`, `dimension`, `ring_generators`, `context_status`, `operation`,
`source_endpoint`, `target_endpoint`, `replay_status`,
`identity_status`, `next_boundary`, and `status`.

Validate by recomputation. `validate_laurent_endpoint_reduction_certificate(cert,
column, R)::Symbol` will derive the paired target column from the source
context, replay `cert.operation`, recompute source and target endpoint metadata,
and reject stale endpoint fields, malformed fields, wrong ring generators, wrong
endpoint index, and operations that do not replay as strict endpoint reductions.

Consume the #345 report conditionally. A summary helper will call
`case008_d14_laurent_endpoint_reduction_search_report()`. When the report is
`:candidate_found`, it will certify the first replay-verified candidate and set
`d14_status = :endpoint_reduction_certificate` with
`next_boundary = :laurent_normality_replay`. When the report is `:exhausted`, it
will set `d14_status = :endpoint_reduction_search_expansion` and will not expose
a d14 certificate.

## Alternatives Considered

1. Promote endpoint-reduction certificate helpers into `src/algorithm`. Rejected
   because #346 explicitly keeps helper promotion and production diagnostics out
   of scope.
2. Add only a synthetic certificate shell. Rejected because the issue requires
   conditional consumption of the #345 `case_008 d=14` search report.
3. Split d14 behavior into a second expert file only when a candidate exists.
   Rejected because the current #345 contract already has a conditional report;
   one focused file can verify both candidate and exhausted outcomes without
   letting synthetic evidence count as real d14 progress.

## Test Strategy

The focused verification command is:

```bash
julia --project=. -e 'include("test/expert/laurent_endpoint_reduction_certificate.jl")'
```

The test will create a fast synthetic Laurent endpoint-reduction case over
`GF(2)[u^+-1, v^+-1]`, construct a certificate from replay, validate it, and
cover negative controls for tampered operations, swapped ring generators, stale
source endpoint metadata, stale target endpoint metadata, malformed certificate
fields, wrong endpoint index, and non-replaying endpoint operations.

The same file will consume #345. If the d14 search report has a candidate, the
first candidate must validate through the certificate shell and move the next
boundary to `:laurent_normality_replay`. If the report is exhausted, d14 remains
at `:endpoint_reduction_search_expansion`.

## Non-Interactive Approval

This Agent Desk run is non-interactive. Under the standing answer policy, the
design is approved because it follows the issue requirements, uses the
recommended test-only scope, and avoids irreversible or public API changes.
