# Issue 234 Park-Woodburn SL3 Driver Catalog Design

## Context

Issue #234 asks for a source-grounded fixture catalog for the parent #184
ordinary-polynomial `SL_3` driver. The existing
`test/fixtures/park_woodburn_polynomial_cases.jl` catalog mixes older route
shells, recursive `SL_n` examples, and fixture-backed Quillen cases. The new
catalog should be narrower: exact ordinary-polynomial `3 x 3` matrices, an
explicit selected induction variable, support-boundary metadata, source
references, and negative controls that prove unsupported determinant-one
matrices stay staged.

No `AGENTS.md`, `CLAUDE.md`, or `GEMINI.md` file is present in this checkout.
GitHub issue comments for #234, #184, #212, and #220 are empty through the
available connector. The local `gh issue view` path is blocked by the sandboxed
proxy, so the supplied Agent Desk issue body is the authoritative issue text.

## Approaches Considered

1. Add a new test-only driver catalog and validator. This is the chosen
   approach. It keeps #184 inputs reusable by later issues without changing
   public APIs or implementing the driver.
2. Extend `park_woodburn_polynomial_cases.jl`. This would keep all
   Park-Woodburn fixtures together, but it would also preserve the mixed old
   route-shell schema that #234 is trying to separate from the #184 driver
   boundary.
3. Add production driver data types now. This is out of scope because #234 is
   only the source-grounded catalog layer and explicitly excludes local solving,
   Quillen assembly, public dispatch, ECP, and recursive `SL_n`.

## Design

Create `test/fixtures/park_woodburn_sl3_driver_cases.jl` with module
`ParkWoodburnSL3DriverFixtureCatalog`. Each positive or staged-positive entry
records:

- `id`, `role`, and `expected_status` (`:supported` or `:staged`);
- `ring_constructor` and `ring` metadata for an exact field-backed ordinary
  polynomial ring;
- a `3 x 3` determinant-one matrix over that ring;
- a `selected_variable` record with name, generator, index, and status;
- explicit support-boundary fields:
  `local_form_status`, `selected_variable_status`,
  `supplied_witness_status`, `upstream_evidence_status`, and a replayable
  support or staged reason;
- optional `local_form_witness`, `supplied_witness`, and `upstream_evidence`;
- `source_refs` citing `refs/arXiv-alg-geom9405003v1 Section 3` or Section 5;
- `consumer_issue_ids`, which must include `#184`;
- optional `downstream_issue_ids` for later #219/#220 consumers.

Supported status is valid only when the validator can replay one of these
evidence routes:

- an explicit local-form witness with a valid selected generator and a monic
  entry in that variable;
- upstream Quillen mainline evidence whose denominator cover, local factors,
  and target product replay exactly.

Staged status is valid for determinant-one inputs that intentionally lack one
of those support proofs, including the legacy patched-substitution fixture and
an arbitrary multivariate determinant-one matrix with no selected monic or
local-form witness.

## Required Cases

The catalog will contain at least these positive or staged-positive entries:

- `sl3-driver-univariate-fast-local-qq`: supported univariate fast-local
  `SL_3` input with explicit local-form metadata.
- `sl3-driver-multivariate-monic-special-form-qq`: supported multivariate
  special-form `SL_3` input whose selected variable is monic in the top-left
  entry.
- `sl3-driver-quillen-mainline-evidence-gf2`: supported upstream
  evidence-backed case reusing the replayable Quillen mainline constructive
  fixture for later #219/#220 consumption.
- `sl3-driver-legacy-quillen-patched-substitution-qq`: staged legacy coverage
  for the old patched-substitution case. It is intentionally not considered
  supported by fixture-id equality alone.
- `sl3-driver-det-one-no-witness-staged-qq`: staged determinant-one
  multivariate `SL_3` input with no replayable local-form, variable-change, or
  upstream evidence.

Negative controls will cover determinant not one, unsupported coefficient ring,
selected variable not a generator, metadata claiming local evidence that is not
present, and a determinant-one multivariate `SL_3` input marked supported
without a replayable witness.

## Files

- Create `test/fixtures/park_woodburn_sl3_driver_cases.jl`.
- Create `test/internal/park_woodburn_sl3_driver_fixtures.jl`.
- Modify `test/runtests.jl` to include the new internal validator.
- Add `docs/superpowers/plans/2026-06-30-issue-234-park-woodburn-sl3-driver-catalog.md`.

## Verification

Focused command required by the issue:

```bash
julia --project=. -e 'include("test/internal/park_woodburn_sl3_driver_fixtures.jl")'
```

Package command required by Agent Desk:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Diff hygiene:

```bash
git diff --check
```

## Automatic Approval

This Agent Desk run is non-interactive. Under the Standing Answer Policy, the
design is approved automatically because it is additive, keeps the public API
unchanged, distinguishes supported from staged inputs with machine-checkable
evidence, and avoids all out-of-scope driver implementation work.

## Spec Self-Review

- No placeholder markers remain.
- The supported/staged boundary is explicit and machine-checkable.
- The legacy Quillen patched-substitution case is staged, not supported by an
  old fixture-id shortcut.
- Negative controls cover every rejection named in the issue.
