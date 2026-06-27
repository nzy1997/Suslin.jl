# Issue 149 Case008 D21 Laurent Witness Unit Design

## Goal

Reduce the extracted `case_008` `d=21` Laurent unimodular column through a
small certified stage in the existing column-reduction certificate pipeline.

The repair target is the fixture from issue #148:

```julia
test/fixtures/toricbuilder_case008_d21_column_boundary.jl
```

`Suslin.reduce_unimodular_column(column, R)` should return elementary factors
that send the fixture column exactly to `e_21`. The ECP column-reduction
certificate should verify, and diagnostics should report the column as
supported instead of `:unsupported_laurent_column_family`.

## Context

The `case_008` Q-block is larger than the already repaired `case_010` boundary.
Issue #148 isolated the first remaining Laurent column-reduction boundary as a
deterministic offline fixture over `GF(2)[u^+/-1, v^+/-1]`. The fixture column
has length 21 and is Laurent-unimodular, but the current reducer cannot reduce
it:

- it has no direct unit entry;
- it does not satisfy the existing one-step Laurent unit-creation predicate;
- normalizing the column as an ordinary polynomial vector does not preserve a
  supported unimodular ordinary-polynomial column for the current reducer.

The useful algebraic predicate is a Laurent Bezout witness. Solving

```julia
[column[1] column[2] ... column[21]] * witness == [1]
```

with the existing `solve_laurent_linear` path produces a sparse witness with a
unit coefficient at index 16. That gives a certified route into the existing
`:witness_unit` reduction stage.

## Approach Options

Recommended: add a Laurent witness-unit route. The route solves for a Laurent
Bezout witness, requires at least one unit witness coefficient, and delegates
factor construction to the existing `_witness_unit_reduction_certificate_stage`.
This is narrow enough for issue #149, but it is still described by an algebraic
predicate rather than by a fixture id.

Alternative: hard-code the `case_008` `d=21` factor sequence. This would be the
smallest local edit, but it would behave like a bypass and would not fit the
certificate/replay style used by prior repairs.

Alternative: implement general Laurent unimodular-column reduction. This is
broader than issue #149, would need substantially more theory and test coverage,
and would blur the staged support boundary.

## Chosen Design

Add a private helper near the existing column-reduction helpers:

```julia
_laurent_unimodular_witness(column, R)
```

It builds the `1 x n` row matrix from the column and calls
`solve_laurent_linear` against a `1 x 1` matrix containing `one(R)`. It returns
the witness vector when the exact solve succeeds and `nothing` for expected
no-solution or known solver capability failures.

Add a reducer route:

```julia
_reduce_via_laurent_witness_unit_certificate(column, R)
```

The route only runs for Laurent polynomial rings. It asks for a witness, finds a
unit witness coefficient, and calls:

```julia
_witness_unit_reduction_certificate_stage(column, witness, pivot_idx, R)
```

The stored certificate stage remains `kind = :witness_unit`; no new replay
stage kind is needed. The new helper is just a Laurent-specific route into
already replayable machinery.

Wire the route into `_reduce_laurent_unimodular_column_certificate` after the
cheap direct unit and Laurent unit-creation stages, before the current Laurent
normalization fallback. Wire diagnostics into
`_diagnose_laurent_unimodular_column_reduction` with an attempted stage marker
named `:laurent_witness_unit`.

## Certificate And Replay

The returned certificate should continue to use the existing
`ECPColumnReductionCertificate` shape:

- validation stage;
- `:witness_unit` stage;
- factor sequence;
- final column;
- replay summary.

Replay already recomputes the witness-unit creation factors, checks that the
witness combines with the input column to `one(R)`, verifies the nested unit
stage, compares stored factors exactly, and confirms the final reduced column.
Tampering with returned factors or stored certificate factors should make exact
verification fail.

## Error Handling

The Laurent witness route should be quiet for unsupported inputs:

- if `solve_laurent_linear` reports no exact solution, return `nothing`;
- if known Laurent solver capability errors occur, return `nothing`;
- if the witness has no unit coefficient, return `nothing`;
- if a stage is selected, `_checked_reduction_factors` remains the final guard
  before factors escape.

Unexpected exceptions should rethrow so real internal defects remain visible.

## Tests

Add `test/expert/case008_d21_laurent_column_reduction.jl`. The test should load
the issue #148 fixture and assert:

- `Suslin.reduce_unimodular_column(column, R)` sends the column exactly to
  `e_21`;
- replacing the first returned factor with the identity no longer reduces the
  column to `e_21`;
- `Suslin.ecp_column_reduction_certificate(column, R)` verifies;
- the certificate contains a `:witness_unit` stage;
- tampering with the certificate factors makes
  `Suslin.verify_ecp_column_reduction` return `false`;
- `Suslin.diagnose_unimodular_column_reduction(column, R)` reports
  `status == :supported`, `failure_code === nothing`, and includes
  `:laurent_witness_unit` in `attempted_stages`.

Update the issue #148 internal boundary validator and fixture expectation so
the same extracted boundary remains valid after it becomes supported. Preserve
the non-unimodular negative control as a precondition failure.

Update `test/expert/laurent_column_reduction_diagnostics.jl` if needed so the
old unsupported Laurent diagnostic contract remains covered by a still
unsupported Laurent column.

## Verification

Run the issue expert test:

```bash
julia --project=. -e 'include("test/expert/case008_d21_laurent_column_reduction.jl")'
```

Run the bounded report command from issue #149:

```bash
julia --project=. scripts/report_toricbuilder_cache_q_blocks.jl --exercise=case_008 --timeout-seconds=180 --output=/tmp/qblock-case008-after-d21.md
```

Expected result: the expert test passes, and the bounded `case_008` report no
longer fails at the old `current d=21` `unsupported_laurent_column_family`
boundary. The report may pass the Laurent `GL_n` certificate route or advance to
a later, smaller structured boundary.

Run the focused diagnostics and fixture checks touched by the change:

```bash
julia --project=. -e 'include("test/internal/toricbuilder_case008_d21_column_boundary.jl")'
julia --project=. -e 'include("test/expert/laurent_column_reduction_diagnostics.jl")'
```

## Out Of Scope

Do not claim arbitrary Laurent unimodular-column support. Do not require full
`case_008` certificate success if the repair exposes a later boundary. Do not
change the public `elementary_factorization` contract for original Laurent
`GL_n` inputs. Do not export new public APIs.

## Automatic Decisions

- Visual companion: skipped because the design choice is algebraic and no visual
  layout or diagram decision would clarify it.
- Scope: keep the repair tied to the issue #149 algebraic predicate, not to the
  fixture id and not to full Laurent column reduction.
- Stage naming: use `:laurent_witness_unit` only in diagnostics; keep the stored
  replay stage as `:witness_unit`.
