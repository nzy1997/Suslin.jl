# Issue 347 Laurent Endpoint-Reduction Internals Design

## Context

Issue #347 follows the expert-only endpoint-reduction search report from #345
and the certificate shell from #346. Those files already have replay,
candidate, and certificate-validation mechanics, but the logic is test-only and
cannot be reused by later diagnostics. This issue promotes the generic pieces
into internal algorithm code without exposing a public API or adding a
case-specific diagnostic gate.

The existing model to follow is the issue #339 promotion for Laurent
link-witness helpers: internal helpers live in `src/algorithm/column_reduction.jl`,
use underscore-prefixed names, and are exercised by a focused internal test.

## Design Choices

Promote a narrow internal helper layer in `src/algorithm/column_reduction.jl`.
This keeps the helper code next to the Laurent descent and link-witness helpers
it reuses, avoids a new include boundary for one small unit, and follows the
current repository style. The helper layer will be generic over two-generator
Laurent columns and will not encode `case_008 d=14` bounds, source indices, or
fingerprints.

The internal surface will include:

- `_laurent_endpoint_reduction_status(endpoint_operation, n, R)::Symbol`
- `_replay_laurent_endpoint_reduction(column, R, endpoint_operation; case_id, require_strict = true)`
- `_laurent_endpoint_reduction_candidate_from_replay(source_column, target_column, R, endpoint_operation; case_id, require_strict = true)`
- `_verify_laurent_endpoint_reduction_candidate(source_column, target_column, R, candidate)::Bool`
- `_laurent_endpoint_reduction_certificate_from_replay(context, endpoint_operation, column, R; source_endpoint = nothing, require_strict = true)`
- `_validate_laurent_endpoint_reduction_certificate(cert, column, R)::Symbol`

Validation is by replay, not supplied metadata. Operation status checks required
fields, the two-generator Laurent ring, endpoint index bounds, nested
entry-addition shape, nested ring generators, and target/endpoint agreement.
Replay recomputes endpoint metadata before and after applying the operation.
Candidate verification replays one operation on both source and linked target
columns and requires both endpoint measures to strictly decrease. Certificate
validation rebuilds the expected certificate from the input column and rejects
stale source or target endpoint metadata.

Update the expert files to delegate to the promoted internals. The #345 search
test will keep its case-specific bounds and counters, but use the internal
status, replay, and candidate verifier for generic mechanics. The #346
certificate test will keep its test-facing constructor names for compatibility,
but make them thin wrappers around the internal certificate helper and validator.

## Alternatives Considered

1. Keep endpoint-reduction helpers test-only. Rejected because #347 exists to
   make these mechanics reusable by later diagnostics.
2. Add a new internal helper file included from `column_reduction.jl`. Rejected
   because the existing Laurent descent/link-witness helpers are already local
   to `column_reduction.jl`, and the promoted endpoint layer is small.
3. Promote `case_008 d=14` search bounds with the helper layer. Rejected because
   the issue explicitly requires generic production helpers and leaves fixture
   constants to later diagnostic integration.

## Test Strategy

Add `test/internal/laurent_endpoint_reduction_helpers.jl`. The internal test
will first use a synthetic two-column Laurent example over `GF(2)[u^+-1,
v^+-1]`, verify status/replay/candidate/certificate helpers, and cover negative
controls for missing fields, swapped generators, wrong endpoint index, malformed
nested operations, stale source metadata, stale target metadata, operation
tampering, and non-strict endpoint replay.

The same internal test will conditionally consume the #345 `case_008 d=14`
search report. If the report found a replay-verified candidate, the first
candidate must validate through the promoted internals. If the report is
exhausted, the internal helper test only asserts that no d14 endpoint
certificate is claimed.

Focused verification commands:

```bash
julia --project=. -e 'include("test/internal/laurent_endpoint_reduction_helpers.jl")'
julia --project=. -e 'include("test/expert/laurent_endpoint_reduction_certificate.jl")'
```

The final branch verification also runs:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Non-Interactive Approval

This Agent Desk run is non-interactive. Under the standing answer policy, the
recommended narrow internal-helper approach is approved because it follows the
issue requirements, mirrors #339, avoids public API changes, and keeps
case-specific constants in expert tests.
