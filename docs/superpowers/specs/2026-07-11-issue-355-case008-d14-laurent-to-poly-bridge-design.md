# Issue 355 Case008 d14 LaurentToPoly Bridge Design

## Context

Issue #355 is the first full-column acceptance check for the LaurentToPoly
delivery track. The existing #353 conversion certificate and #354 bridge are
certificate-first, but the current conversion implementation only handles
already-polynomial Noether outputs or outputs with a Laurent unit entry. The
validated `ToricBuilderCase008D14ColumnBoundary` source has no unit entry after
the recorded Noether variable change, so the bounded report must exercise the
general Park LaurentToPoly preprocessing path before the ordinary ECP bridge.

No `AGENTS.md` or override instruction file is present in this worktree. The
relevant issue context is #355, its dependency #354 / PR #364, #353, #351, and
the #314 LaurentToPoly delivery-track note.

## Approach

Extend the internal generic Laurent-to-polynomial conversion certificate, not
production reducer dispatch. For a verified Laurent Noether certificate over a
two-generator exact field-backed Laurent ring, keep the existing fast paths for
already-polynomial and unit-entry transformed columns. Add the Park Algorithm
6.1 preprocessing path for the remaining supported case:

- scale the selected Noether entry by an elementary diagonal pair so its lowest
  selected-variable coefficient becomes the constant `1`;
- add one computed multiple of the selected row to each other row so every
  other entry has no selected-variable terms of degree `<= 0`;
- choose a nonnegative substitution power `l` for
  `x_selected = y_selected * y_other^l` so all entries become ordinary
  polynomials in the two recorded polynomial generators.

The conversion certificate verifier recomputes the preprocessing factors, the
substitution power, the forward map, the inverse lift map, the polynomialized
column, and the inverse-lift replay. The existing #354 bridge then lifts the
ordinary ECP child factors through the recorded inverse lift and Noether inverse
substitution.

## Report Shape

Create `test/expert/case008_d14_laurent_to_polynomial_bridge.jl` as the
machine-checkable report and replay test. It returns a NamedTuple with:

- `case_id == "case_008"`;
- `source_dimension == 14`;
- `ring_generators == ("u", "v")`;
- `source_boundary == :ToricBuilderCase008D14ColumnBoundary`;
- `conversion_certificate_status`;
- `ordinary_ecp_child_status`;
- `lifted_factor_count`;
- term-count and degree-growth metrics before and after conversion;
- `elapsed_seconds`;
- `final_target_column`;
- `verification_status == :laurent_to_polynomial_bridge_pass` on success.

Create `scripts/report_case008_d14_laurent_to_polynomial.jl` as a bounded CLI
wrapper around the expert report. Timeouts return a nonpassing explicit status
and never masquerade as a successful certificate.

## Validation And Negative Controls

The expert test verifies that the complete ordered factor sequence sends the
original d14 column exactly to `e_14`, that the recomputed factor product
matches the stored product, and that every lifted factor is elementary over the
original Laurent ring.

Negative controls reject:

- the existing non-unimodular d14 fixture;
- a modified source-column entry;
- a modified Noether inverse substitution map;
- a reordered lifted factor sequence;
- a modified lifted factor;
- a forged successful ordinary ECP child certificate.

The test also scans generic source paths to ensure no `case_008` fingerprint
check is added outside fixtures, tests, scripts, docs, or plans.

## Out Of Scope

Do not add production fallback dispatch, recursive column-peel integration,
full Laurent matrix decomposition, a general Laurent `GL_n` correction, public
API surface, or factor-count optimization beyond recording bounded metrics.

## Automatic Decisions

- Clarifying questions: answered from the issue body because this is a
  non-interactive Agent Desk run.
- Visual companion: not offered because the design has no visual UI decision.
- Approach: choose the generic Park LaurentToPoly preprocessing extension plus
  focused report/test, because the current fast-path conversion rejects the d14
  Noether output and fixture-specific shortcuts are forbidden.
- Design approval: accepted under the Standing Answer Policy so the run can
  proceed non-interactively.
