# Issue 100 Quillen Local Realization Certificates Design

## Context

Issue 100 builds on the Issue 99 Quillen patch fixture catalog. The current
Quillen scaffolding can construct `QuillenPatch` objects from supplied
`QuillenLocalContribution` data, and Issue 99 records deterministic fixture
entries with denominators, coverage multipliers, local factors, expected local
corrections, and optional patched-substitution witness metadata.

The missing piece is a replayable local certificate that proves a supplied
local realization was consumed as algebraic data. Patch assembly must be able to
consume the exact local correction recorded by the certificate without trusting
a stage label or a manually supplied product.

## Design Choice

Add an expert/internal certificate type and constructor in
`src/algorithm/quillen_induction.jl`:

- `QuillenLocalRealizationCertificate`
- `quillen_local_realization_certificate(...)`
- `verify_quillen_local_certificate(cert)::Bool`

These names remain unexported. Tests use qualified `Suslin.<name>` access, so
the public API surface remains unchanged.

The certificate records:

- original input matrix or correction,
- target ring and matrix size,
- selected substitution variable,
- denominator and coverage multiplier,
- the preserved `LocalCertificate`,
- local correction row, column, and entry,
- supplied local factors,
- replayed local product,
- recorded local correction,
- optional patched-substitution witness,
- auxiliary witness metadata,
- replay verification summary.

The constructor accepts supplied factors or a supplied local correction. It
coerces denominator, coverage multiplier, and correction data into the target
ring, verifies the selected variable is a generator, rebuilds the expected
elementary local correction, multiplies supplied factors exactly, replays
patched-substitution witness fields when present, and stores the verified
result. The verifier recomputes all of that from the certificate fields and
returns `false` on tampering.

## Alternatives Considered

- Export the new type and verifier. This would force a public compatibility
  decision before Issue 63 finalizes the patching API.
- Extend `LocalCertificate` directly. The issue asks to preserve the existing
  constructor unless migration is necessary, and existing patch tests already
  depend on the small `LocalCertificate` shape.
- Add generic helpers in `src/core/groebner_tools.jl`. The new replay path is
  Quillen-specific and consumes #99 fixture metadata, so it belongs in
  `src/algorithm/quillen_induction.jl`.

The selected option is the smallest compatible path: keep `LocalCertificate`
unchanged, add a Quillen-specific replay certificate, and leave public exports
unchanged.

## Replay Rules

Construction and verification check:

- the original input is a supported exact ordinary polynomial-ring or Laurent
  polynomial-ring matrix/correction accepted by existing Quillen helpers;
- the selected substitution variable is a generator of the target ring;
- denominator, coverage multiplier, local certificate denominators, and
  correction entry are coercible into the target ring;
- the denominator paired with each correction index in `LocalCertificate`
  matches the recorded denominator;
- every supplied factor is a square matrix over the target ring with the
  certificate size;
- the product of supplied factors equals the recorded local correction;
- the recorded local correction equals the expected elementary matrix
  `elementary_matrix(n, row, col, coverage_multiplier * denominator * entry, R)`;
- patched-substitution witness metadata, when present, has the same schema as
  #99 and replays through `patched_substitution`;
- witness metadata stored in the verification summary is stable and must match
  on replay.

The constructor throws `ArgumentError` for invalid input. The verifier catches
ordinary replay failures and returns `false`, rethrowing only interrupts.

## Fixture Coverage

Add `test/expert/quillen_local_certificate.jl` and register it in the expert
group in `test/runtests.jl`.

The expert test includes the Issue 99 fixture catalog and builds replayable
certificates for at least two positive entries. It checks:

- `verify_quillen_local_certificate(cert) == true`;
- recorded local factors multiply exactly to the recorded local correction;
- denominators are coercible into the target ring;
- patched-substitution witness data is replayed when present;
- tampering with one local factor, one denominator, one selected variable, and
  one patched-substitution witness field separately is rejected by verification
  or construction.

The tests use qualified `Suslin.<name>` access for new expert/internal names and
do not update `src/Suslin.jl` exports or `test/public/api_surface.jl`.

## Files

- Modify `src/algorithm/quillen_induction.jl`: add certificate structs,
  constructor, replay summary, and verifier.
- Create `test/expert/quillen_local_certificate.jl`: fixture-backed positive
  and negative replay tests.
- Modify `test/runtests.jl`: register the expert test.
- Add `docs/superpowers/plans/2026-06-23-issue-100-quillen-local-realization-certificates.md`.

## Verification

Focused test:

```bash
julia --project=. -e 'include("test/expert/quillen_local_certificate.jl")'
```

Expert group:

```bash
julia --project=. test/runtests.jl expert
```

Full package verification:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Spec Self-Review

- No placeholders or incomplete sections remain.
- Scope excludes denominator cover choice, global factor assembly, and automatic
  local `SL_3` solving.
- #99 patched-substitution witness fields are reused unchanged.
- Public API compatibility is explicit: no new exports.
