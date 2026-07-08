# Issue 339 Laurent Link-Witness Internal Helpers Design

## Context

Issue #339 follows the expert-only Laurent link-witness context, bounded d14
search report, and certificate shell from #333, #334, and #335. Those tests
proved that the replay-derived `case_008 d=14` post-descent column has a real
candidate and that the first candidate validates through a certificate shell.
The reusable mechanics still live in expert tests, so production diagnostics
cannot depend on them yet.

This issue promotes generic mechanics only. It must not expose a public API,
must not claim full `case_008 d=14` support, and must not bake d14 fixture
fingerprints or fixture-specific helper names into internal code.

No repository-level instruction files were present in this worktree. The
related merged PRs establish the local pattern: add a Superpowers design and
plan, use TDD, keep helper names internal, register focused tests in
`test/runtests.jl`, and verify with focused Julia commands plus `Pkg.test()`.

## Approaches

Recommended and approved non-interactively: promote the generic helper layer
into `src/algorithm/column_reduction.jl`, and update expert tests to call the
new `Suslin._...` helpers. This keeps the future diagnostic integration issue
able to reuse the same primitives without introducing public API.

Rejected alternatives:

- Leave the helpers in expert tests and add an internal wrapper. That would
  preserve the current dependency direction problem and leave production code
  unable to reuse the mechanics.
- Move d14 search/report construction into `src/`. That would promote fixture
  constants and a bounded expert search into generic algorithm code, which is
  explicitly out of scope.
- Update `diagnose_unimodular_column_reduction` now. That belongs to the next
  integration issue after these primitives exist.

## Internal Interfaces

Add internal constants and helpers with leading underscores:

- `_LAURENT_LINK_WITNESS_FIELDS`;
- `_LAURENT_LINK_WITNESS_CANDIDATE_FIELDS`;
- `_LAURENT_LINK_WITNESS_CERTIFICATE_FIELDS`;
- `_laurent_link_witness_status(witness, n::Int, R)::Symbol`;
- `_laurent_link_witness_operation(witness)::NamedTuple`;
- `_laurent_link_endpoint_metadata(entry, R, entry_index::Int, column_measure; case_id)`;
- `_laurent_link_witness_candidate_from_replay(column, R, witness; case_id, require_strict = true)`;
- `_verify_laurent_link_witness_candidate(column, R, candidate)::Bool`;
- `_laurent_link_witness_certificate_from_replay(context, witness, column, R; require_strict = true)`;
- `_validate_laurent_link_witness_certificate(cert, column, R)::Symbol`.

These helpers reuse existing Laurent descent primitives:

- `_laurent_descent_has_fields`;
- `_laurent_descent_ring_generators`;
- `_laurent_descent_exponent_tuple`;
- `_laurent_descent_entry_support`;
- `_laurent_descent_support_bounds`;
- `_laurent_descent_checked_entry_index`;
- `_replay_laurent_elementary_entry_addition`;
- `_laurent_descent_measure_from_column`;
- `_strictly_decreases_laurent_measure`.

The only supported witness family in this issue is
`:two_entry_laurent_combination`. Its replay semantics are exactly the expert
contract from #334: replace the pivot entry with
`pivot_entry + coefficient * u^a * v^b * partner_entry`, where `(a, b)` is in
the declared `ring_generators` order.

## Validation Behavior

`_laurent_link_witness_status` rejects:

- missing witness fields;
- unsupported witness families;
- wrong ring generators;
- malformed pivot or partner indices;
- pivot/partner equality;
- malformed two-generator exponents;
- coefficients that cannot be coerced into the ring.

`_verify_laurent_link_witness_candidate` requires the candidate to match a
fresh replay-derived candidate exactly. It recomputes source and target
endpoint metadata and requires `measure_relation == :strict_decrease`.

`_validate_laurent_link_witness_certificate` recomputes the candidate and
certificate from the supplied column and `cert.witness`. It rejects stale source
endpoint metadata, stale target endpoint metadata, coefficient or exponent
tampering, wrong status or boundary fields, wrong context status, wrong ring
generators, missing fields, and any replay that does not satisfy the witness
identity.

## Test Plan

Add `test/internal/laurent_link_witness_helpers.jl` and register it in the
internal group immediately after `internal/laurent_descent_measure_helpers.jl`.

The internal test covers:

- a small synthetic positive witness over `GF(2)[u, v, u^-1, v^-1]`;
- all negative controls required by the issue;
- the first d14 candidate by consuming the existing expert d14 search report
  and asserting `case_id == "case_008"`, `pivot_index == 10`,
  `partner_index == 1`, `coefficient == 1`, `exponent == (1, -1)`, and
  `next_boundary == :laurent_endpoint_reduction`.

Update expert tests where practical:

- `test/expert/case008_d14_laurent_link_witness_search.jl` keeps d14 search
  constants and reports, but delegates generic witness status, replay,
  metadata, and candidate verification to `Suslin._...` helpers.
- `test/expert/laurent_link_witness_certificate.jl` keeps expert summary
  helpers, but delegates certificate construction and validation to
  `Suslin._...` helpers.

## Verification

Required issue commands:

```bash
julia --project=. -e 'include("test/internal/laurent_link_witness_helpers.jl")'
julia --project=. -e 'include("test/expert/laurent_link_witness_certificate.jl")'
```

Required Agent Desk command:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Run `git diff --check` before committing.
