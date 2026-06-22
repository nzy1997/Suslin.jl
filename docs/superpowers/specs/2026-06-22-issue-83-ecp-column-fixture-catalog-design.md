# Issue 83 ECP Column Fixture Catalog Design

## Context

Issue 62 will expand the Elementary Column Property reducer for unimodular
columns over exact polynomial rings. The current reducer already handles unit
entries, witness-unit shortcuts, bounded monicity normalization, embedded
three-entry blocks, and Laurent normalization, but later issues need stable
ordinary polynomial examples with exact metadata instead of redefining examples
inside each solver test.

Issue 83 asks for a shared catalog over `GF(2)[x,y]` and/or `QQ[x,y]` with
named column entries, variable order, optional monicity and ideal-witness
records, and current status tags. The issue comment narrows the catalog toward
later link-theorem work: include monic-first-entry cases, link-witness cases
with Bezout/resultant/coverage/path metadata, and staged-failure cases where
that metadata is missing.

## Design Choice

Add a test-only fixture module at `test/fixtures/ecp_column_cases.jl` and a
focused internal validator at `test/internal/ecp_column_fixtures.jl`.

The fixture module returns structured named tuples. Each valid case records:

- `id`
- `kind`
- `stage_coverage`
- `ring_constructor`
- `ring`
- `variable_order`
- `entries`
- `column_order`
- `monicity`
- `witnesses`
- `expected`
- `source_refs`
- `consumer_issue_ids`

The catalog also carries `negative_controls`: deliberately corrupted fixture
entries with their own ids. The validator proves these controls are rejected,
so corrupt witness or monicity metadata cannot be decorative.

This remains internal test support. No public Suslin API is added, and
`reduce_unimodular_column` is not changed.

Alternatives considered:

- Export a public fixture API from `Suslin`. That would create compatibility
  burden before the ECP reducer API stabilizes.
- Keep examples inside expert reducer tests. That repeats examples across
  later issues and loses the issue's single-source-of-truth goal.
- Build a broad generic example bank. The issue comment rejects examples whose
  only value is metadata coverage, so the seed catalog stays small and tied to
  current and planned stages.

## Fixture Scope

The initial catalog contains at least eight valid entries:

- `ecp-unit-entry-gf2`: a current-pass unit-entry column over `GF(2)[x,y]`.
- `ecp-witness-unit-gf2`: a current-pass nonunit-entry column whose ideal
  membership witness has a unit coefficient.
- `ecp-variable-change-monic-gf2`: the existing variable-change monicity
  example with exact substitution metadata.
- `ecp-link-bezout-nonunit-witness-qq`: a `QQ[x,y]` link-witness case whose
  supplied ideal witness has no unit coefficients and records
  Bezout/resultant/coverage/path metadata for later link-theorem issues.
- `ecp-longer-embedded-block-gf2`: a longer ordinary column supported through
  the embedded three-entry block path.
- `ecp-unsupported-unimodular-gf2`: a unimodular staged-failure case where the
  current reducer has no supported stage and the link witness is explicitly
  missing.
- `ecp-non-unimodular-gf2`: a negative non-unimodular ordinary column.
- `ecp-monic-first-entry-qq`: a `QQ[x,y]` monic-first-entry case consumable by
  future link-theorem work.

The negative controls include a corrupt witness fixture and a corrupt monicity
fixture. They are not valid catalog cases; they are expected to be rejected by
the same validator used for valid cases.

## Validation Rules

The validator checks required metadata, unique ids, ring constructor metadata,
variable order, and column entry reconstruction from named entries. It calls
`Suslin.is_unimodular_column(column, R)` and checks the exact expected
unimodularity result.

For current reducer status:

- `:passes` cases must run `Suslin.reduce_unimodular_column`, multiply the
  returned factors into the column, and equal `e_n` exactly.
- `:staged_fail` cases must be unimodular and throw an `ArgumentError` whose
  text contains the expected staged-failure substring.
- `:rejects_non_unimodular` cases must be non-unimodular and throw the current
  non-unimodular validation error.

For monicity records, direct claims check that the selected named entry is
monic in the selected ring variable. Variable-change claims reconstruct the
substitution exactly, compare the transformed entry to the recorded value, and
then check monicity in the selected variable.

For witness records, ideal-membership witnesses must satisfy
`sum(w[i] * v[i]) == 1` exactly. Nonunit witness records also assert that none
of the supplied witness coefficients are units. Link-Bezout witness records
must carry a unit resultant plus nonempty coverage and path metadata; missing
link witnesses are only allowed on staged-failure entries.

The validator runs negative controls by expecting `ArgumentError`, and the
test also mutates at least one valid witness and one valid monicity claim at
runtime to prove the validator rejects corrupted metadata.

## Files

- Create `test/fixtures/ecp_column_cases.jl`: exact catalog construction.
- Create `test/internal/ecp_column_fixtures.jl`: validator, focused tests,
  current reducer checks, witness and monicity negative controls.
- Modify `test/runtests.jl`: include the validator in the internal group.

## Verification

Focused validator:

```bash
julia --project=. -e 'include("test/internal/ecp_column_fixtures.jl")'
```

Package entry point:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Spec Self-Review

- No incomplete markers remain.
- Scope is limited to test fixtures and internal validation.
- The catalog includes all fixture categories named by Issue 83 and the issue
  comment's link-theorem guardrail.
- Current reducer behavior is observed, not extended.
