# Issue 316 Laurent Native ECP Boundary Design

## Goal

Expose a structured `:laurent_native_ecp_boundary` diagnostic for validated
Laurent unimodular columns that pass or exhaust the existing Laurent reduction
stages but still lack a Laurent-native ECP route. This is a diagnostic boundary
only; it must not implement Laurent link witnesses, endpoint reductions,
normality replay, or recursive peel integration.

## Context

No `AGENTS.md`, `CLAUDE.md`, or `GEMINI.md` file is present in this worktree.
GitHub issue #316 is open, issue #314 remains the Laurent-native ECP roadmap,
and issue #315 was completed by merged PR #318. Current `main` already contains
the `case_008 d=14` boundary fixture from #315.

The current reducer diagnostic in `src/algorithm/column_reduction.jl` validates
the column first, then runs Laurent stages in this order:

1. `:unit_entry`
2. `:laurent_unit_creation`
3. `:laurent_witness_unit`
4. `:laurent_normalization`
5. delegated ordinary-polynomial diagnostic stages with
   `allow_general_ecp_pipeline = false`
6. `:laurent_elementary_row_preconditioning`

When every stage declines, the top-level diagnostic currently reports
`failure_code == :unsupported_laurent_column_family`. That generic family hides
the architectural boundary described by #314: the missing Laurent-native ECP
pieces are a descent measure, Laurent link witness, endpoint reduction, Laurent
normality/conjugation replay, and recursive peel integration.

## Approach Options

Recommended: add a terminal diagnostic-only Laurent stage after
`:laurent_elementary_row_preconditioning` declines. The stage will be named
`:laurent_native_ecp_boundary` and its `stage_details` entry will carry plain
data fields for the missing route pieces. The overall diagnostic remains
`status == :unsupported` and keeps the established
`:unsupported_laurent_column_family` failure code for compatibility.

Alternative: replace the top-level Laurent failure code with
`:laurent_native_ecp_boundary`. This would make the boundary visible, but it
would break existing callers and tests that treat
`:unsupported_laurent_column_family` as the stable family code.

Alternative: add a stub Laurent ECP reducer that returns `nothing`. That would
look more like an algorithmic route, but it risks implying support and invites
the ordinary ECP pipeline to be used for Laurent inputs.

Chosen approach: terminal diagnostic stage only. It preserves existing
supported stages, leaves certificate reduction unchanged, and makes the
architectural boundary visible without claiming algorithmic progress.

## Diagnostic Shape

Add a small helper in `src/algorithm/column_reduction.jl` that returns the
terminal detail:

```julia
_column_reduction_stage_detail(
    :laurent_native_ecp_boundary,
    R,
    :staged_boundary;
    boundary = :laurent_native_ecp,
    requires_descent_measure = true,
    requires_link_witness = true,
    requires_endpoint_reduction = true,
    requires_laurent_normality_replay = true,
    requires_recursive_peel_integration = true,
    fallback_policy = :diagnostic_only,
)
```

Append this detail only from the Laurent diagnostic path after
`:laurent_elementary_row_preconditioning` records
`:no_row_preconditioning_candidate`. Do not append it when:

- preconditions fail, including `:not_unimodular`;
- any earlier Laurent stage succeeds;
- Laurent normalization cannot produce a normalized ordinary column at all;
- row preconditioning succeeds;
- the ordinary-polynomial diagnostic path is running.

Do not route Laurent columns into
`_reduce_via_general_ecp_pipeline_certificate`. The existing Laurent
normalization diagnostic may continue to inspect the normalized ordinary
polynomial column with `allow_general_ecp_pipeline = false`.

The #315 d14 fixture is already validated by its fixture helper, but running
the generic ideal-membership check and later large-support witness checks is
too expensive for this boundary diagnostic. Add an internal diagnostic keyword
`assume_unimodular = true` so tests can reuse the fixture validator and skip
only the precondition unimodularity proof. Keep the default at `false` so
non-unimodular inputs still fail with `:not_unimodular` before any stage
attempts.

For the already-validated d14 fixture, opt into
`laurent_large_support_diagnostic_decline = true` so entries exceeding the
current staged diagnostic support size decline the witness-unit and normalized
ordinary-delegation probes with plain data instead of invoking slow
solver/ideal-membership calls. The d14 boundary records
`max_entry_term_count = 3734` against a conservative `support_term_limit =
1000`, then continues to row preconditioning and the terminal Laurent-native
ECP boundary. This large-support decline is off by default so ordinary
diagnostics continue to run existing supported stages, including large columns
that still have a Laurent witness-unit route.

## Tests

Create `test/expert/laurent_native_ecp_boundary_diagnostics.jl` for the new
boundary contract. It should include:

- the #315 `case_008 d=14` fixture, asserting the diagnostic reaches
  `:laurent_native_ecp_boundary` when called with
  `assume_unimodular = true` and
  `laurent_large_support_diagnostic_decline = true` after fixture validation;
- d14 witness and normalization details asserting the large-support diagnostic
  declines instead of entering expensive solver/ideal-membership paths;
- field-level checks for `outcome == :staged_boundary`,
  `boundary == :laurent_native_ecp`, every `requires_* == true`, and
  `fallback_policy == :diagnostic_only`;
- a stage-order check that the boundary is attempted after
  `:laurent_elementary_row_preconditioning`;
- the supported `case_008 d=15` negative control, asserting it remains
  supported by `:laurent_elementary_row_preconditioning` and does not report
  `:laurent_native_ecp_boundary`;
- the d14 non-unimodular negative control, asserting it fails with
  `:not_unimodular` before any stage attempts.
- a large-support Laurent column with an explicit unit witness, asserting the
  default diagnostic still supports `:laurent_witness_unit` and does not report
  `:laurent_native_ecp_boundary`.
- a runner-registration assertion so the new expert file remains in
  `TEST_GROUP_FILES["expert"]` and is covered when CI runs the expert group.

Update `test/expert/laurent_column_reduction_diagnostics.jl` so its existing
unsupported Laurent example expects the new terminal stage and detail. Existing
supported fixture checks should explicitly assert that the boundary stage is
absent.

## Verification

Run:

```bash
julia --project=. -e 'include("test/expert/laurent_native_ecp_boundary_diagnostics.jl")'
julia --project=. -e 'include("test/expert/laurent_column_reduction_diagnostics.jl")'
julia --project=. test/runtests.jl expert
julia --project=. -e 'using Pkg; Pkg.test()'
git diff --check
```

Expected result: all commands exit 0.

## Out Of Scope

Do not implement Laurent ECP. Do not add Laurent link witnesses, endpoint
reductions, Laurent normality/conjugation replay, or recursive Laurent peel
integration. Do not claim arbitrary Laurent `GL_n` support. Do not make
diagonal monomial balancing or polynomialization the primary route.

## Automatic Decisions

- Visual companion skipped because this is a structured diagnostic and test
  change, not a visual design problem.
- Clarifying questions skipped because Agent Desk is non-interactive and the
  issue body gives exact diagnostic fields and verification commands.
- Recommended approach selected: append a terminal diagnostic-only stage after
  row preconditioning declines, because it preserves existing stage behavior
  and exposes the #314 boundary without implementing Laurent ECP.
- Existing `:unsupported_laurent_column_family` failure code retained for
  compatibility; the more specific boundary is exposed through
  `attempted_stages` and `stage_details`.
- Design approval auto-approved under the Standing Answer Policy.
- Added `assume_unimodular = true` only for the already-validated d14 fixture,
  because the default diagnostic path must still prove non-unimodular negative
  controls before stage attempts.
- Added an explicit opt-in large-support diagnostic decline at 1000 terms for
  the validated d14 boundary, because the default path must preserve existing
  supported witness-unit routes even for large entries.
- Added a regression test for a large-support Laurent column with a unit
  witness after code review showed that a default term-count cutoff could hide
  supported `:laurent_witness_unit` diagnostics.
- Registered the new expert file in `test/runtests.jl` after final review
  found that direct issue verification did not automatically protect the new
  regression in the static CI runner.
