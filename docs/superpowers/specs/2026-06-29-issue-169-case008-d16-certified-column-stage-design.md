# Issue 169 Case008 D16 Certified Column Stage Design

## Context

Issue #169 asks for `reduce_unimodular_column(fixture.failing_column, fixture.ring)`
and `ecp_column_reduction_certificate(fixture.failing_column, fixture.ring)` to
reduce the `case_008` `d=16` Laurent column to `e_16` with replayable
elementary factors and certificate verification. The stage must be algebraic,
not a `case_008` fixture-id special case.

The dependency from issue #168 is merged. Its bounded search finds that one
right-column addition on the full `d=16` matrix, target column 16 from source
column 8 with coefficient `1`, reaches an existing supported Laurent
witness-unit family. The public reducer does not receive that matrix source
column, so the implementation cannot depend on it directly.

## Approaches

1. Certified row-preconditioned Laurent column stage. Search bounded elementary
   row additions on the input column itself, using coefficient `1`, and accept a
   candidate only when the transformed column reduces through the existing
   certified Laurent reducer. This keeps the stage algebraic and replayable.
   This is the chosen approach.
2. Full matrix preconditioning replay. Reconstruct or pass the source matrix
   column found by issue #168. This does not fit the requested public interface,
   which only passes the column and ring.
3. Fixture recognition. Detect the stored `case_008` column and replay a stored
   factor list. This is out of scope because it special-cases the fixture id.

## Design

Add a Laurent-only stage named `:laurent_elementary_row_preconditioning`.
It runs after the current direct Laurent stages fail and before normalization
returns unsupported. The stage tries elementary factors
`elementary_matrix(n, target, source, one(R), R)` for `target != source`.
For each transformed column it invokes the existing Laurent base reducer that
does not recursively include this new stage. A candidate is accepted only when
the transformed column has a verified base certificate.

The final factor sequence is:

```julia
vcat(transformed_certificate.factors, [precondition_factor])
```

This order replays as `S * P * column`, where `P` is the row-preconditioning
factor and `S` is the certified reducer for the transformed column.

## Certificate Replay

The stage stores the input column, target/source indices, coefficient,
precondition factor, transformed column, transformed certificate, composed
factors, and output column. `_ecp_replay_stage` recomputes the elementary
precondition factor, replays the transformed certificate, recomposes the factors,
and checks exact equality against the stored fields. Tampering with either a
factor or a stored stage field must make `verify_ecp_column_reduction` return
`false`.

## Diagnostics

`diagnose_unimodular_column_reduction` records the new attempted stage as
`:laurent_elementary_row_preconditioning`. On success it reports
`status == :supported`, `failure_code === nothing`, and a stage detail with
`outcome == :supported`, `target_index`, `source_index`, `coefficient`, and the
transformed certificate stage kind. Non-unimodular columns still stop in
precondition validation with `status == :precondition_failed` and do not attempt
the new stage.

## Tests

Add `test/expert/case008_d16_laurent_column_reduction.jl` covering:

- factor replay reduces the fixture column exactly to `e_16`;
- the certificate verifies and contains
  `:laurent_elementary_row_preconditioning`;
- diagnostics report `status == :supported`;
- tampering with a factor and with the stage coefficient makes certificate
  verification fail;
- the fixture's non-unimodular negative control still fails preconditions.

Update diagnostics and fixture metadata for the new supported boundary. Keep
the stored stage name algebraic and do not claim full `case_008` route success
unless the bounded report verifies it.

## Verification

Run:

```bash
julia --project=. -e 'include("test/expert/case008_d16_laurent_column_reduction.jl")'
julia --project=. scripts/report_toricbuilder_cache_q_blocks.jl --exercise=case_008 --timeout-seconds=180 --output=/tmp/qblock-case008-after-d16.md
julia --project=. -e 'using Pkg; Pkg.test()'
```
