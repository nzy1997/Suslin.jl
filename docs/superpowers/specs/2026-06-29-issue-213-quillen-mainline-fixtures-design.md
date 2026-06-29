# Issue 213 Quillen Mainline Fixtures Design

## Context

Issue #213 prepares fixture data for the broader #183 Quillen local-to-global
patching route. The repository already has the #99 Quillen fixture catalog in
`test/fixtures/quillen_patch_cases.jl`, and later issues #100-#105, #114, and
#115 already consume those ids for local certificates, denominator covers,
normalization, deterministic global assembly, and the narrow Park-Woodburn
Quillen route. PR context confirms this lineage: PR #106 introduced the #99
catalog, PRs #107/#108/#120/#124 built the replay layers on top, PR #118 added
the Park-Woodburn polynomial catalog, and PR #127 wired the current narrow
Quillen fixture route.

Park-Woodburn Section 3 describes local realizability over localizations
`R_M[X]`, raw local denominators `d_j`, their product denominator `r_i`, a
cover identity `sum r_i^l g_i = 1`, the telescoping substitutions
`A(X - sum X r_i^l g_i)`, and an `A(0)` base term. The existing #99 catalog has
exact ordinary-polynomial examples but does not expose those Section 3
mainline concepts directly enough for later #183 children to consume cold.

## Approaches Considered

Recommended: add a thin `test/fixtures/quillen_mainline_cases.jl` wrapper that
loads and reuses entries from `quillen_patch_cases.jl` by id, augmenting them
with Park-Woodburn mainline metadata. Add
`test/internal/quillen_mainline_fixtures.jl` as the validator. This preserves
#99 ids, avoids duplicating algebraic matrices, and makes the new fields
explicit.

Alternative: extend every #99 fixture in place. That would avoid one extra
fixture file, but it would churn an older schema that existing issues already
consume.

Alternative: create a standalone Park-Woodburn Quillen catalog with copied
matrices. That risks divergence from #99 and is explicitly discouraged by the
issue.

The chosen design is the recommended wrapper because it is test-only,
backward-compatible, and source-grounded while avoiding new solver behavior.

## Catalog Shape

Create module `QuillenMainlineFixtureCatalog` with:

- `catalog()` returning `cases` and `negative_controls`.
- `cases_by_id()` returning positive entries by stable id.
- Each positive entry contains the reused #99 fixture as `patch_case` and keeps
  the same `id`.

Each positive entry records:

- ordinary exact ring metadata and determinant-one target through the reused
  `patch_case`;
- `raw_denominator_provenance`, separate from any Park-Woodburn exponent;
- `denominator_cover` with denominators, multipliers, coverage terms, and
  exact coverage sum;
- `local_evidence` records that mirror local factors today and can later grow
  into #214 factor-sequence records;
- `patched_substitution_chain` using the Park-Woodburn sign convention
  `X -> X - X * r_i^l * g_i`;
- `base_term_evidence`, describing whether `A(0) = I`, whether base-term
  factors are supplied, or whether the case is staged until base-term evidence
  exists;
- `expected_global_product`, `source_refs`, and `consumer_issue_ids`.

Positive entries include:

- `quillen-two-open-cover-qq`: the existing two-open `QQ[X,r,g]` shape with
  `A(0) = I`.
- `quillen-nontrivial-multipliers-qq`: a nonconstant multiplier cover, also
  with `A(0) = I`.
- `quillen-patched-substitution-witness-qq`: a patched-substitution chain case
  whose base term is intentionally staged.
- `quillen-constructive-acceptance-gf2`: a denominator-one local-evidence case
  shaped for the future #211 Murthy adapter handoff, without invoking Murthy
  solving.

Negative controls reuse positive entries and corrupt:

- the denominator coverage sum;
- local evidence replay/product.

## Validation

The validator includes the #99 catalog and the new mainline wrapper. It checks:

- unique positive and negative ids;
- ordinary exact polynomial-ring metadata through the reused patch fixture;
- determinant-one expected global product and exact equality with the reused
  target matrix;
- denominator cover terms and multipliers sum exactly to one;
- raw denominator provenance is present and separate from `exponent_l`;
- local evidence factors equal the source local factors and their product
  equals the expected global product;
- patched-substitution chain metadata has a Park-Woodburn sign convention and
  either replays through `Suslin.patched_substitution` where a source witness is
  available or records explicit staged metadata;
- base-term status is one of `:assumes_identity`, `:supplied_factors`, or
  `:staged`;
- every `source_refs` includes `refs/arXiv-alg-geom9405003v1 Section 3`;
- every positive entry includes `#183` in `consumer_issue_ids`.

The validator must reject the two negative controls with `ArgumentError`.

## Files

- Create `test/fixtures/quillen_mainline_cases.jl`.
- Create `test/internal/quillen_mainline_fixtures.jl`.
- Modify `test/runtests.jl` to include the new internal validator.
- Add `docs/superpowers/plans/2026-06-29-issue-213-quillen-mainline-fixtures.md`.

No public API, solver, route-selection, Murthy, ECP, recursive `SL_n`, Laurent,
or ToricBuilder behavior changes are included.

## Verification

Focused validator:

```bash
julia --project=. -e 'include("test/internal/quillen_mainline_fixtures.jl")'
```

Full package verification:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Automatic Approval

This Agent Desk run is non-interactive. Under the Standing Answer Policy, the
design is approved automatically because it chooses the conservative wrapper
over duplicating or reshaping the #99 schema, preserves existing fixture ids,
keeps all changes in test/support files, and leaves broad #183 automation out
of scope.

## Spec Self-Review

- No placeholders remain.
- Scope is limited to fixture catalog data and validation.
- The design preserves #99 ids and reuses existing algebraic entries.
- Park-Woodburn Section 3 source grounding and #183 consumer metadata are
  explicit.
- Negative controls cover denominator coverage and local evidence replay.
