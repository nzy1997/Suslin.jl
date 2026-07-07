# Issue 328 Certified Laurent Descent Diagnostic Design

## Context

Issue #328 advances the Laurent-native ECP diagnostic boundary for the checked-in
`case_008 d=14` column. The existing diagnostic reaches
`:laurent_native_ecp_boundary` and reports `requires_descent_measure = true`.
Merged prerequisite work added the d14 boundary fixture, the descent measure
contract, a bounded entry-addition search, a replayable descent-step
certificate, and internal helpers that recompute measures and validate replayed
certificates.

No repository `AGENTS.md`, `CLAUDE.md`, or `CONVENTIONS.md` file is present in
this worktree. Issue #328 has no comments. Relevant merged PR context:

- #319 exposed the terminal Laurent-native ECP boundary diagnostic.
- #324 defined the replayable d14 descent measure contract.
- #325 recorded the bounded d14 entry-addition evidence.
- #326 added the replayable descent-step certificate shell.
- #330 promoted the measure, replay, strict-decrease, and certificate
  validation helpers into internal production code.

## Approach Options

Recommended: add a narrow diagnostic helper that recognizes only the recorded
entry-addition operation shape on the recorded two-generator Laurent column,
checks compact before/after support fingerprints, constructs the candidate
certificate by recomputing before and after measures through the internal
helpers, validates that certificate with
`_validate_laurent_descent_step_certificate`, and emits
`:laurent_descent_step_certificate` only when validation returns `:ok`.

Alternative: key the diagnostic on the `case_008` fixture id or the full d14
baseline constants. This would be easy to target, but it would put fixture
identity into production code and make the diagnostic less algebraic.

Alternative: run the full bounded descent search inside the diagnostic. This is
out of scope and would imply a partial Laurent descent algorithm rather than one
certified staged evidence step.

Chosen approach: replay the one recorded operation conservatively and expose the
stage only when the supplied column matches the recorded support evidence and
the certificate validates. This keeps the public diagnostic truthful: one
descent step is certified, while Laurent link witness, endpoint reduction,
Laurent normality replay, and recursive peel integration remain required.

## Diagnostic Shape

When the recorded operation validates, `diagnose_unimodular_column_reduction`
adds a `:laurent_descent_step_certificate` attempted stage after
`:laurent_elementary_row_preconditioning` declines and before
`:laurent_native_ecp_boundary`.

The stage detail records only plain data:

- `outcome = :certified_descent_step`;
- `descent_scope = :single_certified_step`;
- `operation_family = :entry_addition`;
- `target_index = 1`;
- `source_index = 2`;
- `coefficient = 1`;
- `exponent = (-1, 1)`;
- `before_measure`;
- `after_measure`;
- `measure_relation = :strict_decrease`;
- `replay_status = :ok`;
- `next_boundary = :laurent_link_witness`.

The terminal `:laurent_native_ecp_boundary` detail accepts an optional
`certified_descent_step` flag. For the certified d14 path it sets
`requires_descent_measure = false`, `certified_descent_scope =
:single_certified_step`, and keeps `requires_link_witness`,
`requires_endpoint_reduction`, `requires_laurent_normality_replay`, and
`requires_recursive_peel_integration` true. Generic unsupported Laurent columns
still receive the old `requires_descent_measure = true` boundary.

The support fingerprint is a compact recorded-evidence guard, not a public
algorithm. It records per-entry support counts plus a stable support checksum
for the certified before and after columns. This prevents the diagnostic from
self-certifying an arbitrary column whose coarse measure happens to decrease
under the same operation.

## Controls

Tests extend `test/expert/laurent_native_ecp_boundary_diagnostics.jl` because it
already owns the d14 boundary diagnostic. Coverage asserts:

- the certified stage appears before the terminal boundary for the d14 fixture;
- the stage detail records the exact operation and replay status;
- the boundary no longer treats the descent-measure contract as missing for that
  certified path;
- d15 supported diagnostics do not include the certified stage;
- non-unimodular columns do not include the certified stage;
- a tampered copy of the d14 column does not include the certified stage and
  still reaches the generic boundary.

## Out Of Scope

Do not implement Laurent link witnesses, endpoint reductions,
normality/conjugation replay, determinant normalization, recursive peel
integration, full `case_008` success, or a general Laurent descent algorithm.

## Self Review

This design has no placeholders, preserves public API compatibility, avoids
fixture-id checks in production code, and maps every issue acceptance criterion
to a focused diagnostic assertion.
