# Issue 335 Laurent Link-Witness Certificate Design

## Context

Issue #335 follows the expert-only Laurent link-witness context from #333 and
the bounded post-descent d14 link-witness search report from #334. The new
deliverable is a replayable certificate shell, not a production Laurent ECP
stage. It must accept a synthetic known-good witness, consume the d14 search
report conditionally, and avoid claiming endpoint reduction unless the d14
search actually produced a replay-verified candidate.

This repository keeps this staged Laurent evidence in `test/expert/*.jl` as
NamedTuple reports plus validators, registered in `test/runtests.jl`. The
certificate shell should follow that pattern and should not export a public API.

## Chosen Approach

Recommended and approved non-interactively: add
`test/expert/laurent_link_witness_certificate.jl`.

The new file will include the #334 search helper and reuse its replay helpers
for witness schema checks and endpoint metadata. It will add a certificate
constructor named `laurent_link_witness_certificate(context, witness, column,
R)` and a validator named
`validate_laurent_link_witness_certificate(cert, column, R)::Symbol`.

Rejected alternatives:

- Moving this into `src/` would expose an internal/public API before the issue
  asks for one.
- Recomputing or widening the bounded d14 search would violate #334's fixed
  radius and checked-count contract.
- Treating a synthetic certificate as d14 progress would contradict #335's
  out-of-scope section.

## Certificate Shape

`laurent_link_witness_certificate(...)` returns a NamedTuple with stable fields:

- `case_id`;
- `dimension`;
- `ring_generators`;
- `context_status`;
- `witness`;
- `source_endpoint`;
- `target_endpoint`;
- `replay_status == :ok`;
- `identity_status == :verified`;
- `next_boundary == :laurent_endpoint_reduction`;
- `status == :link_witness_certificate`.

The constructor validates the witness schema, checks that the context has
`status == :link_witness_context`, checks that the witness ring generators match
the ring, replays the entry-addition identity over the supplied source column,
and stores endpoint metadata from the replay result.

## Validation

`validate_laurent_link_witness_certificate(cert, column, R)::Symbol` recomputes
all replay-derived endpoint data from `cert.witness` and the supplied `column`.
It rejects:

- missing certificate fields;
- malformed witness fields;
- wrong ring generators;
- wrong context/status/boundary fields;
- stale `source_endpoint`;
- stale `target_endpoint`;
- tampered coefficient or exponent data;
- pivot/partner equality or out-of-range indices;
- non-identity or non-decreasing replay.

The validator returns `:ok` only when the supplied certificate exactly matches
the replay-derived certificate and the accepted candidate verifier also accepts
the witness replay.

## D14 Consumption

Add `case008_d14_laurent_link_witness_certificate_summary()`.

It consumes `case008_d14_laurent_link_witness_search_report()`:

- when the search report has `status == :candidate_found`, it validates the
  first replay-verified candidate through the certificate shell and returns
  `d14_status == :link_witness_certificate`;
- when the search report has `status == :exhausted`, it returns
  `d14_status == :no_link_witness_candidate`;
- it never treats the synthetic certificate as real d14 endpoint-reduction
  progress.

## Tests

The new expert file will contain:

- a known-good synthetic Laurent certificate over a small two-entry column;
- validation of the synthetic certificate, including `identity_status`,
  `next_boundary`, and `status`;
- negative controls for tampered coefficients, wrong ring generators, stale
  endpoints, malformed witness fields, and non-identity/non-decreasing replay;
- conditional d14 summary assertions matching #335.

`test/runtests.jl` will register the file immediately after
`expert/case008_d14_laurent_link_witness_search.jl`.

## Verification

Required focused command:

```bash
julia --project=. -e 'include("test/expert/laurent_link_witness_certificate.jl")'
```

Required Agent Desk command:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

The focused command must pass whether the current bounded d14 search is
`candidate_found` or `exhausted`.
