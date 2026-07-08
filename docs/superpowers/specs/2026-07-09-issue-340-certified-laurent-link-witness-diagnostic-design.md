# Issue 340 Certified Laurent Link-Witness Diagnostic Design

## Context

Issue #340 advances the observable production diagnostic for the validated
`case_008 d=14` Laurent column path. The current diagnostic can expose a
certified single Laurent descent step, but the terminal Laurent-native ECP
boundary still reports `next_boundary = :laurent_link_witness` and
`requires_link_witness = true`. After issues #333-#335 and #339, the
post-descent d14 column has internal replay helpers that validate a
`two_entry_laurent_combination` link-witness certificate:

- `pivot_index = 10`
- `partner_index = 1`
- `coefficient = 1`
- `exponent = (1, -1)`
- `next_boundary = :laurent_endpoint_reduction`

This issue should surface that validated certificate as diagnostic evidence
only. It must not implement endpoint reduction, Laurent normality replay,
recursive peel integration, or a general Laurent link-witness solver.

## Recommended Approach

Use a narrow production diagnostic bridge that replays the already-certified
d14 path through internal helpers.

The diagnostic should first obtain the existing
`:laurent_descent_step_certificate` for the original d14 fixture column. It
should replay that descent operation to obtain the source column for the
link-witness certificate. A new helper should then require the exact
post-descent d14 support fingerprint, build the fixed d14 witness, create a
certificate with `_laurent_link_witness_certificate_from_replay`, and accept it
only when `_validate_laurent_link_witness_certificate` returns `:ok`.

The helper should return `nothing` for tampered columns, unsupported Laurent
columns, wrong rings, failed replay, or failed certificate validation. That
keeps the stage conservative and fixture-bound.

## Alternatives Considered

1. Add the narrow certified diagnostic bridge now. This is the recommended
   option because it reuses the internal certificate validator, advances the
   observable boundary, and keeps all new behavior tied to the exact d14
   fingerprint.
2. Recreate the expert search/context logic inside production diagnostics.
   This would duplicate test-only report code and make the diagnostic look
   more general than it is.
3. Wait for endpoint-reduction context support before changing diagnostics.
   This would leave diagnostics behind the already-validated link-witness
   certificate and would not satisfy #340.

## Diagnostic Shape

For the validated `case_008 d=14` fixture, `attempted_stages` should include
`:laurent_link_witness_certificate` after
`:laurent_descent_step_certificate` and before
`:laurent_native_ecp_boundary`.

The stage detail should be plain data derived from the internal certificate:

- `outcome = :certified_link_witness`
- `witness_family = :two_entry_laurent_combination`
- `pivot_index = 10`
- `partner_index = 1`
- `coefficient = 1`
- `exponent = (1, -1)`
- `replay_status = :ok`
- `identity_status = :verified`
- `source_endpoint`
- `target_endpoint`
- `next_boundary = :laurent_endpoint_reduction`

The later `:laurent_native_ecp_boundary` detail should set
`requires_link_witness = false` when this certified link witness is present.
It should still report:

- `requires_endpoint_reduction = true`
- `requires_laurent_normality_replay = true`
- `requires_recursive_peel_integration = true`

When no certified link witness is present, the existing Laurent boundary
behavior remains unchanged.

## Components

`src/algorithm/column_reduction.jl`

- Add fixed d14 link-witness constants beside the existing certified descent
  constants.
- Add `_laurent_link_witness_diagnostic_certificate(column, R)` that validates
  the exact post-descent d14 fingerprint and internal certificate replay.
- Add `_laurent_link_witness_certificate_stage_detail(R, certificate)` to
  expose only plain diagnostic fields.
- Extend `_laurent_native_ecp_boundary_stage_detail` with a
  `certified_link_witness` keyword so terminal details can clear
  `requires_link_witness` and advance `next_boundary` to
  `:laurent_endpoint_reduction` for the certified d14 path.
- Update the Laurent diagnostic flow to insert the link-witness stage only
  after the descent certificate stage validates.

`test/expert/laurent_native_ecp_boundary_diagnostics.jl`

- Update the d14 positive diagnostic assertions to require the new stage,
  ordering, witness fields, certificate statuses, endpoint metadata, and
  terminal boundary flags.
- Keep negative controls for d15, non-unimodular Laurent columns, and a
  tampered d14 column.

No new public exports are required.

## Testing

The primary focused verification is:

```bash
julia --project=. -e 'include("test/expert/laurent_native_ecp_boundary_diagnostics.jl")'
```

Expected result: the command exits 0, the d14 diagnostic contains the new
certified link-witness stage with the exact witness data, and negative
controls do not report the new stage.

The final verification gate also runs:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Scope Guard

This design intentionally does not implement Laurent endpoint reductions,
endpoint-reduction contexts, Laurent normality/conjugation replay, determinant
normalization, recursive peel integration, or full `case_008` success.
