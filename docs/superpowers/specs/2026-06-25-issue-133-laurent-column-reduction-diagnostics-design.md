# Issue 133 Laurent Column Reduction Diagnostics Design

## Goal

Add an internal diagnostic entry point for exact unimodular-column reduction
failures so tests can inspect unsupported Laurent column boundaries without
parsing staged error text.

## Chosen Approach

Use an unexported function, `diagnose_unimodular_column_reduction(column, R)`,
implemented near `reduce_unimodular_column` in
`src/algorithm/column_reduction.jl`. The function returns a named tuple with
stable fields:

- `status`: `:supported`, `:unsupported`, or `:precondition_failed`;
- `failure_code`: `nothing` for supported columns, or a stable symbol for the
  failure family;
- `ring_profile`: a named tuple describing whether the ring is Laurent or
  ordinary polynomial, with generator names and coefficient-ring text;
- `column_length`: the input column length;
- `attempted_stages`: a tuple of stage symbols tried in reducer order;
- `message`: a human-readable diagnostic.

This keeps the API internal and test-support-only. It is narrower than adding a
public exported diagnostic type and less invasive than changing the reducer to
throw structured exceptions. `reduce_unimodular_column` may keep throwing the
current `ArgumentError`.

## Stage Model

The diagnostic will share the existing reducer predicates and helpers rather
than reimplement algebraic decisions. For ordinary polynomial columns it will
record attempts in the reducer order:

1. `:unit_entry`
2. `:witness_unit`
3. `:monicity_normalization`
4. `:three_entry_block`

For Laurent columns it will first try the existing Laurent unit-entry fast path,
then record `:laurent_normalization` and diagnose the normalized polynomial
column with the same ordinary stage sequence. A validated Laurent column that
cannot be reduced by the current exact reducer reports
`failure_code == :unsupported_laurent_column_family`.

Precondition failures stop before stage attempts. A non-unimodular column
therefore reports `status == :precondition_failed` and
`failure_code == :not_unimodular`, never
`:unsupported_laurent_column_family`.

## Testing

Add `test/expert/laurent_column_reduction_diagnostics.jl` with three checks:

- the #132 `case_010` boundary fixture reports `status == :unsupported`,
  `failure_code == :unsupported_laurent_column_family`,
  `column_length == 5`, and attempted stages including `:unit_entry`,
  `:witness_unit`, `:monicity_normalization`, and `:three_entry_block`;
- the first successful Laurent column-peel column from the same fixture reports
  `status == :supported`;
- the fixture's non-unimodular negative control reports a precondition failure,
  not `:unsupported_laurent_column_family`.

The issue verification command is:

```bash
julia --project=. -e 'include("test/expert/laurent_column_reduction_diagnostics.jl")'
```

The full repository verification remains:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Out Of Scope

Do not add a new reduction algorithm. Do not export the diagnostic function.
Do not change Laurent column-peel behavior or ToricBuilder report statuses.
