# Issue 43 Side-Aware Elementary Preconditioning Design

## Context

Issue #43 adds the elementary preconditioning layer needed before later Laurent
reduction-certificate work can replay local obligations. Existing helpers can
construct elementary matrices and compose ordinary left-to-right factor
sequences, but they do not record whether a factor was applied on the left or on
the right. That distinction is essential: left multiplication changes rows,
right multiplication changes columns, and replay must fail when side metadata is
wrong.

There is no `AGENTS.md` in this worktree. The repository guidance comes from
`README.md`, the existing helper style in `src/core/elementary_matrices.jl`, and
nearby expert tests such as `test/expert/block_embeddings.jl`.

## Approach Options

1. Add focused side-aware preconditioning helpers beside
   `elementary_matrix`. This keeps row/column elementary primitives near the
   existing matrix helper, reuses exact Oscar matrix operations, and avoids
   touching algorithm reducers. This is the selected approach.
2. Add a new algorithm module for preconditioning. That may become useful once a
   real preconditioning search exists, but issue #43 explicitly excludes
   searching for a useful sequence.
3. Introduce a certificate or optimizer abstraction now. That would be
   premature because the issue only asks for primitives and exact replay checks.

## Selected Design

Extend `src/core/elementary_matrices.jl` with three public helpers:

- `elementary_preconditioning_step(A, side::Symbol, target::Integer, source::Integer, coefficient)`
- `replay_elementary_preconditioning(A, steps)`
- `verify_elementary_preconditioning(A, steps, expected)`

`elementary_preconditioning_step` accepts `side == :left` or `side == :right`.
For `:left`, it constructs the elementary factor `E[target, source] =
coefficient` and returns `factor * A`; this records the row operation
`row[target] += coefficient * row[source]`. For `:right`, it constructs
`E[source, target] = coefficient` and returns `A * factor`; this records the
column operation `column[target] += coefficient * column[source]`.

The return value is a named tuple with these fields:

- `side`
- `target`
- `source`
- `coefficient`
- `factor`
- `transformed_matrix`

The helper coerces `coefficient` into `base_ring(A)` with the existing
`_coerce_into_ring` behavior before constructing the elementary matrix. It
validates side, integer coordinates, in-range coordinates, and distinct target
and source indices. The API works for square matrices and rectangular matrices:
left factors have size `nrows(A)`, and right factors have size `ncols(A)`.

`replay_elementary_preconditioning` starts from `A`, applies each step in order,
and uses each step's recorded `side` and `factor` to decide whether to multiply
on the left or right. It validates each factor's dimensions and base ring
against the current matrix before applying it.

`verify_elementary_preconditioning` catches validation and multiplication
failures and returns `false`; otherwise it returns whether replay reconstructs
`expected` exactly. This follows the existing repository pattern for boolean
verification helpers such as `verify_factorization` and
`verify_laurent_gl_normalization`.

Export all three helpers from `src/Suslin.jl` and add public API surface checks.
Register the new expert test under `test/runtests.jl`'s expert group.

## Validation Contract

The helpers throw `ArgumentError` for unsupported sides, non-integer indices,
repeated target/source coordinates, out-of-range coordinates, missing step
metadata, wrong factor dimensions for the recorded side, and mismatched factor
base rings.

No optimizer, sequence search, determinant normalization, or certificate type is
introduced. The helper layer only expresses elementary row additions, elementary
column additions, ordered replay, and exact final reconstruction checks.

## Tests

Add `test/expert/elementary_preconditioning.jl`. The focused tests must verify:

- a left operation returns the expected factor and transformed matrix,
- a right operation returns the expected factor and transformed matrix with the
  documented column-addition ordering,
- replaying a known sequence over a Laurent matrix reconstructs the final matrix
  exactly,
- `verify_elementary_preconditioning` returns `true` for the recorded sequence,
- swapping each step's side metadata makes verification return `false`, proving
  side metadata matters,
- invalid side, repeated coordinates, out-of-range coordinates, mixed base
  rings, and wrong factor dimensions are rejected.

## Verification

Run the issue-specific command:

```bash
julia --project=. -e 'include("test/expert/elementary_preconditioning.jl")'
```

Run the package entry point required by Agent Desk:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Run the documented expert/full coverage command:

```bash
julia --project=. test/runtests.jl all
```

## Automatic Approval

This Agent Desk run is non-interactive. Under the standing answer policy, the
design and written spec are approved automatically because the selected approach
is the narrowest change that satisfies issue #43, follows existing repository
patterns, and preserves the issue's out-of-scope boundary around searching for
preconditioning sequences.
