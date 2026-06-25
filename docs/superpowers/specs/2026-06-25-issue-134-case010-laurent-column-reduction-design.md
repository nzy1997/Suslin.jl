# Issue 134 Case 010 Laurent Column Reduction Design

## Goal

Reduce the extracted `case_010` length-5 Laurent unimodular column with a
small certified stage in the existing column-reduction certificate pipeline.

## Context

Issue #132 added `test/fixtures/toricbuilder_case010_column_boundary.jl`, which
extracts the blocking column over `GF(2)[u^+/-1, v^+/-1]`:

```julia
[
    u^-1*v^-2 + u^-1*v^-3,
    v^2 + v + 1 + v^-1 + v^-2,
    v^2 + v,
    1 + v^-3 + v^-4,
    v + 1,
]
```

Issue #133 added internal diagnostics showing that the current reducer tries
unit-entry, Laurent normalization, witness-unit, monicity normalization, and
three-entry block stages before reporting the column as unsupported.

The useful exact relation in the fixture is:

```julia
(v^2 + v + 1 + v^-1 + v^-2) + (v + v^-2) * (v + 1) == 1
```

So one elementary row operation can create a unit entry before the existing
unit-entry stage reduces the whole column.

## Approach Options

Recommended: add a certified Laurent one-step unit-creation stage. The stage
searches for a row pair `(pivot, source)` where `pivot_entry + coeff *
source_entry == 1`, using exact Laurent division to compute `coeff`. It then
records the elementary creation factor, the created column, the nested
unit-entry stage, and the final factors. Replay recomputes the coefficient and
factors from the stored metadata, so tampering is rejected.

Alternative: hard-code the `case_010` column and its factor sequence. This is
smaller but would be an escape hatch rather than an algorithm stage and would
not fit the existing certificate model.

Alternative: implement a general Euclidean two-entry Laurent reduction. This is
more mathematically general but broader than the issue asks for and would need
more coverage before claiming support.

## Chosen Design

Add `_reduce_via_laurent_unit_creation_certificate(column, R)` in a focused
column-reduction extension included immediately after
`src/algorithm/column_reduction.jl`. It specializes the existing Laurent
certificate and replay hooks after the direct unit-entry fast path and before
Laurent normalization.

The stage is deliberately narrow:

- It only runs for Laurent polynomial rings.
- It only runs for length-5 columns, matching the extracted boundary target.
- It only targets creation of the literal unit `one(R)`.
- It only accepts exact division from `divexact(one(R) - column[pivot],
  column[source])`.
- It only succeeds after `_checked_reduction_factors` proves the resulting
  elementary factors send the input column to `[0, ..., 0, 1]`.

This supports the `case_010` boundary without claiming arbitrary Laurent
unimodular-column support.

## Certificate And Replay

The new stage kind is `:laurent_unit_creation`. Its metadata records:

- the original input column;
- `pivot_index` and `source_index`;
- `target_unit == one(R)`;
- the exact `creation_coefficient`;
- the single `creation_factors` row operation;
- the `created_column`;
- the nested `unit_stage`;
- the complete `factors` and `output_column`.

`_ecp_replay_stage` recomputes the creation factor and nested unit-stage
reduction from the current input column and stored indices. Replay checks the
stored coefficient, created column, nested stage, factor sequence, and output
column exactly.

## Diagnostics

`diagnose_unimodular_column_reduction(column, R)` should report the `case_010`
boundary as supported and include `:laurent_unit_creation` in
`attempted_stages`. The existing unsupported Laurent diagnostic coverage should
move to a different still-unsupported Laurent column so #133's failure-code
contract remains covered.

## Tests

Add `test/expert/case010_laurent_column_reduction.jl`. The test loads the #132
fixture, calls `Suslin.reduce_unimodular_column(column, R)`, multiplies the
returned factors by the column, and asserts the exact target column
`[0, 0, 0, 0, 1]`.

The test also:

- obtains `Suslin.ecp_column_reduction_certificate(column, R)`;
- asserts `Suslin.verify_ecp_column_reduction(certificate)`;
- checks that the certificate contains a `:laurent_unit_creation` stage;
- replaces the first returned factor with the identity and asserts the mutated
  factor product does not reduce the column;
- replaces the first certificate factor with the identity and asserts
  `verify_ecp_column_reduction` rejects the tampered certificate.

Update the #132 boundary fixture validator so it still validates the extracted
column metadata and unimodularity after the algorithm becomes supported. Update
the #133 diagnostics test so the old `case_010` unsupported assertion becomes a
supported-stage assertion, while a separate Laurent column continues to cover
`:unsupported_laurent_column_family`.

The issue verification command is:

```bash
julia --project=. -e 'include("test/expert/case010_laurent_column_reduction.jl")'
```

The full verification command remains:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Out Of Scope

Do not update the `case_010` report status. Do not claim support for arbitrary
Laurent unimodular columns. Do not export new public APIs.

## Automatic Decisions

- Clarifying question auto-answer: treat the #134 objective as authorization to
  change old #132/#133 tests that asserted the former unsupported state, because
  those tests otherwise contradict the requested new reducer behavior.
- Approach auto-answer: choose the certified one-step Laurent unit-creation
  stage, because it is the smallest reusable extension and matches the issue's
  preference for certificate/replay machinery.
- Visual companion auto-answer: skip it, because this is algebraic code and no
  visual decision would clarify the design.
