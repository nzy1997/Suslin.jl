# Issue 11 Block Embedding and Factor Sequence Helpers Design

## Context

Later SL_n-to-local-SL_3 reduction needs a small, exact helper layer for placing
constructive SL_2 and SL_3 blocks into larger matrices. The current package has
only `elementary_matrix(n, i, j, a, R)` in `src/core/elementary_matrices.jl` and
algorithm code builds larger factors ad hoc.

Issue #11 depends on the Laurent parent validators from #7 and the test-command
contract from #21. Both dependencies are present on `main`: Laurent validators
live in `src/core/rings.jl`, and the documented full-suite command is
`julia --project=. test/runtests.jl all`.

## Approach Options

1. Add focused block and factor-sequence helpers beside
   `elementary_matrix`. This keeps the helper layer near the existing matrix
   primitive, shares exact ring coercion and parent checks, and avoids changing
   algorithm modules. This is the selected approach.
2. Add the helpers inside the SL_3 factorization module. That would keep the
   first consumer nearby, but the helpers are matrix primitives needed by
   multiple later reductions, not part of the current narrow SL_3 factorizer.
3. Introduce a new embedding abstraction or wrapper type. That is premature:
   the issue asks for exact embeddings and sequence composition checks, not a
   new factor-sequence domain model.

## Selected Design

Extend `src/core/elementary_matrices.jl` with three public helpers:

- `block_embedding(block, n::Int, indices)`
- `embed_factor_sequence(factors, n::Int, indices)`
- `compose_factor_sequences(sequences...)`

`block_embedding` validates that `block` is square, `n` can contain the block,
and `indices` are distinct in-range coordinates. It returns `identity_matrix(R,
n)` over `R = base_ring(block)` with the block copied into `indices[i],
indices[j]`. Untouched entries remain exact identity entries in the same parent.

`embed_factor_sequence` validates a nonempty factor collection, verifies all
small factors have the same square size and base ring, checks that the index
count matches that size, and returns each factor embedded with `block_embedding`.
It must preserve order so multiplying the returned factors equals embedding the
small product.

`compose_factor_sequences` accepts one or more factor collections, validates that
all factors are square matrices with the same size and base ring, and returns a
single vector in the same order. Empty input collections are allowed only when at
least one nonempty collection provides the matrix type; all-empty composition is
ambiguous and throws `ArgumentError`.

Export the three helpers from `src/Suslin.jl` and cover them in the public API
surface test. Keep implementation minimal: no optimization, no normalization,
and no algorithm-specific assumptions about determinant or elementary shape.

## Validation Contract

The helpers throw:

- `DimensionMismatch` for nonsquare blocks, nonsquare factors, blocks larger
  than `n`, and index-count mismatches.
- `ArgumentError` for repeated indices, out-of-range indices, empty factor
  sequences where a matrix parent cannot be inferred, mixed dimensions, and
  mixed base rings.

Base-ring comparison should accept both equality and object identity, matching
the existing `verify_factorization` behavior. Coercion is intentionally not
performed between different matrix parents.

The Laurent example should use `suslin_laurent_polynomial_ring` and values such
as `x^-1 * y`, so the helpers exercise the supported Laurent ring path from #7.
The ordinary polynomial example should use `suslin_polynomial_ring` or Oscar's
ordinary polynomial constructor.

## Tests

Add `test/expert/block_embeddings.jl` and register it in the expert group.
The focused tests must verify:

- embedding a 2x2 ordinary polynomial block into a larger identity preserves all
  untouched entries and places the block at the requested coordinates,
- embedding a 3x3 Laurent block works over a Laurent parent and preserves
  Laurent entries such as `x^-1 * y`,
- embedded factor sequences multiply to the same result as embedding the small
  product directly,
- composed embedded factor sequences preserve order and product,
- repeated indices, out-of-range indices, nonsquare blocks, index-count
  mismatches, nonsquare factors, and mixed parent/dimension sequences throw the
  required error classes.

## Verification

Run the issue-specific command:

```bash
julia --project=. -e 'include("test/expert/block_embeddings.jl")'
```

Run the package entry point required by Agent Desk:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Run the full suite documented by #21:

```bash
julia --project=. test/runtests.jl all
```

## Automatic Approval

This Agent Desk run is non-interactive. Under the standing answer policy, the
design and written spec are approved automatically because the selected approach
is the narrowest option that satisfies issue #11, follows existing repository
patterns, and keeps public API compatibility limited to the requested helpers.
