# Issue 165 Stage-Level Laurent Column Reduction Diagnostics Design

## Goal

Extend `diagnose_unimodular_column_reduction(column, R)` with additive
stage-level metadata so Laurent column-reduction boundaries can distinguish why
an attempted stage did not support a validated unimodular column.

The immediate target is the `case_008` `d=16` fixture from issue #164. Its
diagnostic must show that the Laurent witness solve succeeds but the witness has
no unit entry, and that Laurent normalization delegates to a length-16 ordinary
polynomial column whose reducer profile is `:not_unimodular`.

## Context

The current diagnostic returns a named tuple with these stable fields:

- `status`
- `failure_code`
- `ring_profile`
- `column_length`
- `attempted_stages`
- `message`

Those fields must remain valid for existing callers. The missing information is
inside `attempted_stages`: `:laurent_witness_unit` currently covers both "no
witness exists" and "a witness exists but contains no unit entry." The d=16
fixture needs the second case pinned exactly. It also needs the Laurent
normalization fallback to describe the normalized ordinary polynomial column as
not unimodular, not as a generic ordinary reducer unsupported case.

## Approach Options

Recommended: add a `stage_details` tuple of named tuples alongside the existing
fields. Each diagnostic stage appends a compact record with `stage`,
`ring_kind`, `outcome`, and relevant metadata such as `pivot_index`,
`witness_unit_index`, `normalized_column_length`, `normalized_ring_kind`,
`normalized_status`, and `normalized_failure_code`. This follows the existing
named-tuple diagnostic style and is fully backward-compatible.

Alternative: replace `attempted_stages` with richer records. That would be more
compact internally, but it would break callers and tests that already consume
the stage symbols.

Alternative: introduce a typed diagnostic struct hierarchy. That may be useful
if diagnostics become public API, but it is too much surface area for this
internal/test-support function and would add compatibility questions outside the
issue.

## Chosen Design

Keep `_column_reduction_diagnostic` as the single diagnostic constructor and add
a final `stage_details` argument. The returned named tuple will include the
existing fields unchanged plus:

```julia
stage_details = tuple(stage_details...)
```

Precondition failures return `stage_details == ()`, just as they return
`attempted_stages == ()`.

Reducer diagnosis helpers thread both collections:

- `attempted::Vector{Symbol}` preserves the existing ordered stage list.
- `details::Vector` records exactly one named tuple per attempted stage.

For successful stages, the stage detail uses `outcome = :supported` and records
the selected index when available:

- `:unit_entry`: `pivot_index`
- `:laurent_unit_creation`: `pivot_index` when the certificate exposes one,
  otherwise no pivot field
- `:laurent_witness_unit`: `witness_unit_index`
- `:witness_unit`: `witness_unit_index`
- `:monicity_normalization`: `normalized_column_length`
- `:three_entry_block`: `normalized_column_length`

For failed Laurent witness-unit attempts:

- if `_laurent_unimodular_witness(column, R)` returns `nothing`, record
  `outcome = :witness_unavailable` and `witness_unit_index = nothing`;
- if a witness exists but `findfirst(is_unit, witness) === nothing`, record
  `outcome = :witness_without_unit` and `witness_unit_index = nothing`;
- if a unit witness entry exists, record `outcome = :supported` and the
  one-based `witness_unit_index`.

For Laurent normalization attempts, normalize the Laurent column and compute the
ordinary polynomial profile before delegating to ordinary reducer diagnosis:

- if `is_unimodular_column(poly_column, P) == false`, record
  `outcome = :normalized_not_unimodular`, `normalized_status =
  :precondition_failed`, `normalized_failure_code = :not_unimodular`, and
  `normalized_column_length = length(poly_column)`;
- if the unimodularity check throws, record `outcome =
  :normalized_unimodularity_check_failed`, `normalized_status =
  :precondition_failed`, and `normalized_failure_code =
  :unimodularity_check_failed`;
- if the normalized ordinary column is unimodular, delegate to the ordinary
  diagnostic path and then record `outcome = :delegated_to_polynomial`,
  `normalized_status`, and `normalized_failure_code` from that ordinary result.

The d=16 fixture should hit the first normalization case: a length-16 ordinary
polynomial column with `normalized_status = :precondition_failed` and
`normalized_failure_code = :not_unimodular`.

## Testing

Update `test/expert/laurent_column_reduction_diagnostics.jl` to assert:

- all diagnostic named tuples have `stage_details`;
- a non-unimodular Laurent negative control still reports
  `status == :precondition_failed`, `failure_code == :not_unimodular`,
  `attempted_stages == ()`, and `stage_details == ()`;
- the supported d=21 fixture still reports `status == :supported` and has a
  `:laurent_witness_unit` detail with `outcome == :supported` and a non-nothing
  `witness_unit_index`;
- the unsupported Laurent fixture still records failed-stage details without
  inventing a successful stage.

Update `test/internal/toricbuilder_case008_d16_column_boundary.jl` to assert:

- the fixture diagnostic contains a `:laurent_witness_unit` detail with
  `outcome == :witness_without_unit` and `witness_unit_index === nothing`;
- the diagnostic contains a `:laurent_normalization` detail with
  `normalized_column_length == 16`, `normalized_status ==
  :precondition_failed`, and `normalized_failure_code == :not_unimodular`;
- the fixture validator still returns `:ok`;
- corrupted non-unimodular fixture data remains rejected as `:not_unimodular`.

## Verification

Run:

```bash
julia --project=. -e 'include("test/expert/laurent_column_reduction_diagnostics.jl")'
julia --project=. -e 'include("test/internal/toricbuilder_case008_d16_column_boundary.jl")'
julia --project=. test/runtests.jl
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected result: all commands exit 0.

## Out Of Scope

Do not add a new reducer stage. Do not change `reduce_unimodular_column`
behavior. Do not export a new public diagnostic type. Do not require `case_008`
to pass beyond the current d=16 boundary.

## Automatic Decisions

- Visual companion: skipped because this is an API/test diagnostic change, not a
  visual design task.
- Clarifying questions: skipped because Agent Desk requires a non-interactive
  run and the issue body defines the interface, acceptance criteria, and out of
  scope boundary.
- Approach: selected additive `stage_details` named tuples because they preserve
  every existing diagnostic field.
- Laurent normalization outcome: record normalized ordinary-column
  precondition failures before ordinary reducer stage attempts, because the
  issue specifically requires `:not_unimodular` rather than an ordinary
  unsupported reducer case.
