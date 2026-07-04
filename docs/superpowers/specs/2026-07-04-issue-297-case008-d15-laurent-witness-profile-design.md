# Issue 297 Case008 D15 Laurent Witness Profile Design

## Context

There is no `AGENTS.md`, `CLAUDE.md`, `GEMINI.md`, or `CODEX.md` file in
this checkout. Repository guidance comes from `README.md`: expert algorithm
checks live under `test/expert`, default package tests run through
`julia --project=. -e 'using Pkg; Pkg.test()'`, and expert files can be run
directly with `include(...)`.

Issue #295 is already present on `main` through the
`ToricBuilderCase008D15ColumnBoundary` fixture and default internal validator.
The fixture records the `case_008` Laurent column at peel dimension `d=15` and
already validates the current unsupported reducer diagnostic. Issue #297 asks
for a narrower expert-only profile test that documents exactly where the
current Laurent column reducer stops before any repair is added.

GitHub issue lookup was unavailable in this Agent Desk sandbox because the
configured local proxy rejected the API request. The issue body supplied by
Agent Desk is the source of issue context for this non-interactive run.

## Clarifying Decisions

No interactive questions were asked because this is a non-interactive Agent
Desk run. Under the Standing Answer Policy, the conservative choices are:

- mirror the existing `case008_d16_laurent_witness_profile.jl` style;
- keep profile helpers local to the new expert test;
- use `Suslin.diagnose_unimodular_column_reduction(column, R)` as the primary
  source of truth;
- avoid new public APIs and avoid reducer implementation changes;
- include a cheap supported Laurent direct-unit negative control so the helper
  cannot classify every Laurent column as the unsupported `case_008 d=15`
  family.

## Approaches Considered

1. Add a standalone expert diagnostic-profile test. This is the recommended
   approach because it is focused, reviewable, and follows the issue's
   requested pre-repair evidence shape.
2. Extend the broad `laurent_column_reduction_diagnostics.jl` expert test. This
   would reuse helper patterns but would make an already broad diagnostic file
   noisier and less focused on the `d=15` boundary.
3. Add fixture-level validation to the default internal test. This would put
   pre-repair unsupported evidence into the default suite, which the issue
   explicitly says to avoid unless later work needs it.

## Design

Create `test/expert/case008_d15_laurent_witness_profile.jl`. It will include
the #295 fixture, define local helper functions for extracting diagnostic stage
details and classifying the exact unsupported profile, and assert the ordered
stage boundary for the fixture.

The test will assert:

- `diagnostic.status == :unsupported`;
- `diagnostic.failure_code == :unsupported_laurent_column_family`;
- `diagnostic.column_length == 15`;
- attempted stages include `:unit_entry`, `:laurent_unit_creation`,
  `:laurent_witness_unit`, `:laurent_normalization`, and
  `:laurent_elementary_row_preconditioning`;
- no diagnostic detail is marked `:supported`;
- the witness stage reports `:witness_without_unit` and has no witness-unit
  index;
- the normalization stage records a normalized column length of `15` and a
  normalized precondition failure such as `:not_unimodular`;
- the row-preconditioning stage reports `:no_row_preconditioning_candidate`;
- the local helper recognizes the fixture profile as `true` and rejects a
  synthetic Laurent column that is supported by direct unit entry.

Any optional witness-level assertions will stay local to the expert test and
will avoid changing `src/`.

## Testing

Required focused command:

```bash
julia --project=. -e 'include("test/expert/case008_d15_laurent_witness_profile.jl")'
```

Required Agent Desk package command:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Spec Self-Review

- Placeholder scan: no unfinished placeholder markers remain.
- Consistency check: the design creates expert-only coverage and does not add
  default-suite registration.
- Scope check: the design is limited to one test file plus Superpowers docs.
- Ambiguity check: the accepted helper classification is explicitly tied to the
  unsupported status, failure code, length, and required stage detail values.
