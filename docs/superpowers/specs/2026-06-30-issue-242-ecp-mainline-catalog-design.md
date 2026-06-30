# Issue 242 ECP Mainline Catalog Design

## Context

Issue #242 asks for a source-grounded fixture catalog for parent #185, the
general ordinary-polynomial Elementary Column Property reducer. The existing
`test/fixtures/ecp_column_cases.jl` catalog covers older reducer families for
#62, #85, #87, and #88, but it does not record the full Park-Woodburn Section 4
mainline obligations as reusable #185 cases. The new catalog should be narrower
than production ECP support: it records exact inputs and staged evidence
requirements without implementing monicity search, link extraction, lower-ring
induction, normality, public dispatch, or recursive matrix factorization.

No `AGENTS.md`, `CLAUDE.md`, or `GEMINI.md` file is present in this checkout.
GitHub issue comments for #242, #185, #184, and #181 are empty through the
available connector. PR #241 is the nearest merged catalog PR and establishes
the current pattern: a test-only fixture module, an internal validator, explicit
source refs and consumers, staged support boundaries, negative controls, and no
public API changes.

## Approaches Considered

1. Add a new test-only ECP mainline catalog and validator. This is the chosen
   approach. It keeps the #185 schema separate from older ECP fixtures while
   allowing later child issues to cite stable case ids.
2. Extend `ecp_column_cases.jl`. This would reuse old cases directly, but it
   would blur #185 mainline metadata with the older reducer contract and would
   make it harder to distinguish staged support from fixture-only support.
3. Add production ECP data types now. This is out of scope because #242 is a
   validated catalog issue and explicitly excludes the algorithms that would
   consume the catalog.

## Design

Create `test/fixtures/ecp_mainline_cases.jl` with module
`ECPMainlineFixtureCatalog`. Each positive or staged-positive entry records:

- `id`, `role`, and `expected_status` (`:supported` or `:staged`);
- `ring_constructor` and `ring` metadata for an exact field-backed ordinary
  polynomial ring;
- named column entries and `column_order`, with length at least three;
- selected-variable metadata naming an actual ring generator;
- monicity metadata, including selected entry and transformed entry;
- a concrete unimodularity witness whose coefficients reconstruct one;
- explicit Section 4 stage metadata for monicity, link witness, link step,
  lower-variable reduction, normality, and `SL_3` realization expectations;
- `source_refs` containing `refs/arXiv-alg-geom9405003v1 Section 4`;
- `consumer_issue_ids = ("#185",)` or a tuple containing `#185`;
- optional `reused_fixture_id` for cases derived from the older ECP catalog.

Supported status is valid only when required mainline evidence is present and
replayable. Staged status must record which required evidence is missing. The
validator checks ordinary exact field-backed polynomial rings, selected
generator metadata, length, unimodularity for positive cases, exact witness
equations, Park-Woodburn Section 4 source refs, #185 consumers, unique ids, and
negative-control rejection.

## Required Cases

The catalog will contain at least these positive or staged-positive entries:

- `ecp-mainline-gf2-hard-slice`: reuses the existing GF(2) variable-change hard
  slice with supplied link and lower-variable expectations.
- `ecp-mainline-qq-link-bezout`: reuses the existing QQ monic-first/link-Bezout
  family with exact supplied link witness metadata.
- `ecp-mainline-length4-coupled-qq`: a length-four column whose witness uses all
  entries, so it is not just a supported three-entry subcolumn.
- `ecp-mainline-monicity-change-gf2`: records a monicity-changing substitution
  as a first-class mainline case.
- `ecp-mainline-sl3-route-qq`: records a case whose link realization is shaped
  to consume the #184 `SL_3` route rather than fixture-only endpoint transport.

Negative controls cover the four issue-required rejections: non-unimodular
column, corrupted Bezout/resultant/link witness, selected variable that is not
a ring generator, and a supported-status claim without required link or
lower-variable evidence.

## Files

- Create `test/fixtures/ecp_mainline_cases.jl`.
- Create `test/internal/ecp_mainline_fixtures.jl`.
- Modify `test/runtests.jl` to include the new internal validator.
- Add `docs/superpowers/plans/2026-06-30-issue-242-ecp-mainline-catalog.md`.

## Verification

Focused command required by the issue:

```bash
julia --project=. -e 'include("test/internal/ecp_mainline_fixtures.jl")'
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
design is approved automatically because it is additive, follows the supplied
issue schema, keeps public APIs unchanged, and avoids all out-of-scope
algorithm implementation.

## Spec Self-Review

- No placeholder markers remain.
- The supported/staged boundary is explicit and machine-checkable.
- The new schema is separate from `ecp_column_cases.jl` but can reuse old case
  ids as provenance.
- Negative controls cover every rejection named in the issue.
