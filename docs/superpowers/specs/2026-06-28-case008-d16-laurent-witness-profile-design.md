# case008 d16 Laurent Witness Profile Design

Issue: #167

## Objective

Add an expert-only regression/profile test for the `case_008` d=16 Laurent
column witness returned by `Suslin._laurent_unimodular_witness(column, R)`.
The test records stable algebraic facts that explain why the existing
witness-unit reducer cannot apply and preserves a small positive control for
the reducer route.

## Context

The repository already has a generated d=16 boundary fixture in
`test/fixtures/toricbuilder_case008_d16_column_boundary.jl` and an internal
diagnostic test in `test/internal/toricbuilder_case008_d16_column_boundary.jl`.
That diagnostic confirms the d=16 Laurent column is unimodular but unsupported
because the Laurent witness stage returns a witness without a unit entry and
normalization does not leave a unimodular polynomial column.

The new expert test should inspect the raw witness family, not just the
diagnostic status. It must remain exact and fixture-derived, with no
production reducer changes and no preconditioning search.

## Approaches Considered

1. Local expert-test profile helpers. This creates
   `test/expert/case008_d16_laurent_witness_profile.jl` with helper functions
   that compute a `NamedTuple` profile from a column, ring, and witness. It is
   narrow, exact, and avoids committing an exploratory public API.

2. Production diagnostic API. This would add a reusable exported or internal
   profile function under `src/algorithm`. It might be useful in future
   reducer work, but would prematurely stabilize internal diagnostic semantics for an
   exploratory issue.

3. Extend the existing d=16 internal fixture test. This would keep related
   checks together, but the issue explicitly asks for an expert profile and
   the witness solve is more appropriate for the expert group.

Chosen approach: local expert-test profile helpers. It has the smallest API
surface, matches the issue guidance, and produces a stable profile that
follow-up reducer work can cite.

## Profile Shape

The profile is a `NamedTuple` with stable field names:

- `witness_length`
- `witness_unit_entry_count`
- `witness_unit_indices`
- `witness_nonunit_obstructions`
- `witness_entry_term_counts`
- `witness_entry_support_bounds`
- `witness_entry_is_unit`
- `witness_gcd_is_unit`
- `column_witness_entry_gcd_is_unit`
- `column_dot_witness`
- `column_dot_witness_is_one`
- `existing_witness_unit_stage_applicable`
- `elementary_pair_unit_attempts`
- `unit_obstruction_reason`

Support bounds are per-entry tuples of exponent minima and maxima over the
Laurent generators. GCD checks are best-effort exact checks: when Oscar can
compute them, the profile records whether the gcd of all witness entries is a
unit and whether each same-index column/witness entry gcd is a unit. If an
exact gcd call is unavailable, the corresponding field records `nothing`.

The elementary pair attempts are deliberately simple. For every pair of
witness entries, the profile checks the two combinations `a + b` and
`a - b` and records whether any pair creates a unit. This is not a
preconditioning search; it is a bounded algebraic probe that records whether
the most basic elementary combinations explain the obstruction.

## Tests

The d=16 test uses the fixture and asserts:

- `witness_length == 16`
- `witness_unit_entry_count == 0`
- `column_dot_witness == 1`
- `column_dot_witness_is_one`
- every witness entry fails `is_unit`
- per-entry obstruction records explain why the existing witness-unit stage
  cannot apply
- the profile carries `unit_obstruction_reason == :no_witness_unit_entry`

The negative control in the issue is implemented as a positive synthetic
control: a small Laurent column with no unit entries and a known witness that
has a unit entry. The profile must report at least one unit witness entry, and
`Suslin._reduce_via_laurent_witness_unit_certificate(column, R)` must produce
a valid witness-unit stage whose factors reduce the column to the target.

## Runtime

The d=16 witness solve is fixture-heavy and should stay in the expert group.
The test file is intended to be run directly with:

```bash
julia --project=. -e 'include("test/expert/case008_d16_laurent_witness_profile.jl")'
```

The default package test entry point does not include expert tests.

## Out of Scope

No production reducer stage, preconditioning search, public API, or fixture
rewriting is included.
