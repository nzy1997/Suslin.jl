# Issue 11 Block Embedding and Factor Sequence Helpers Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add exact block-embedding and factor-sequence helpers for placing 2x2 and 3x3 constructive blocks into larger Suslin matrices.

**Architecture:** Keep the helper layer in `src/core/elementary_matrices.jl` beside `elementary_matrix`, because it is a matrix primitive rather than an algorithm-specific routine. Add expert tests that exercise ordinary polynomial and Laurent parents and update the public API surface for the exported helpers.

**Tech Stack:** Julia, Oscar exact matrices, Suslin ring constructors, `Test`.

## Global Constraints

- Work in the existing linked worktree on branch `agent/issue-11-add-block-embedding-and-factor-sequence-helpers-run-1`; do not create another worktree.
- Follow TDD: write the failing helper tests first, run them and confirm the expected failure, then write production code.
- Keep the helpers exact and minimal: no optimization, no normalization, and no algorithm-specific determinant or elementary-shape assumptions.
- Preserve matrix base rings and dimensions; do not coerce factors between different matrix parents.
- Throw `DimensionMismatch` for nonsquare blocks, nonsquare factors, blocks larger than the target size, and index-count mismatches.
- Throw `ArgumentError` for repeated indices, out-of-range indices, empty factor sequences where no matrix parent can be inferred, mixed factor dimensions, and mixed base rings.
- The issue-specific verification command is `julia --project=. -e 'include("test/expert/block_embeddings.jl")'`.
- The Agent Desk package verification command is `julia --project=. -e 'using Pkg; Pkg.test()'`.
- The full suite documented by #21 is `julia --project=. test/runtests.jl all`.

---

## File Structure

- Modify `src/Suslin.jl` to export `block_embedding`, `embed_factor_sequence`, and `compose_factor_sequences`.
- Modify `src/core/elementary_matrices.jl` to add validation helpers plus the three public helper routines.
- Create `test/expert/block_embeddings.jl` for focused issue #11 coverage.
- Modify `test/expert/elementary_matrices.jl` only if the new tests reveal a local helper belongs with the existing elementary test. The planned implementation should not need this.
- Modify `test/runtests.jl` to register `expert/block_embeddings.jl`.
- Modify `test/public/api_surface.jl` to assert the new exports exist and are the same bindings as the public imports.

### Task 1: Block Embedding and Factor Sequence Helpers

**Files:**
- Modify: `src/Suslin.jl`
- Modify: `src/core/elementary_matrices.jl`
- Create: `test/expert/block_embeddings.jl`
- Modify: `test/runtests.jl`
- Modify: `test/public/api_surface.jl`

**Interfaces:**
- Consumes: existing `elementary_matrix(n::Int, i::Int, j::Int, a, R)`, `identity_matrix`, `base_ring`, `nrows`, `ncols`, and Suslin ring constructors.
- Produces: `block_embedding(block, n::Int, indices)`, `embed_factor_sequence(factors, n::Int, indices)`, and `compose_factor_sequences(sequences...)`.

- [ ] **Step 1: Write the failing expert tests and public API assertions**

Create `test/expert/block_embeddings.jl` with this content:

```julia
using Test
using Suslin
using Oscar

function block_embedding_product(factors, R, n::Int)
    product = identity_matrix(R, n)
    for factor in factors
        product *= factor
    end
    return product
end

function assert_identity_outside_block(embedded, indices, R)
    index_set = Set(indices)
    for i in 1:nrows(embedded), j in 1:ncols(embedded)
        if !(i in index_set && j in index_set)
            @test embedded[i, j] == (i == j ? one(R) : zero(R))
        end
    end
end

@testset "block embeddings" begin
    R, (x,) = suslin_polynomial_ring(QQ, ["x"])
    block2 = matrix(R, [
        one(R) + x  x;
        x^2         one(R) - x
    ])

    embedded2 = block_embedding(block2, 4, [1, 3])
    expected2 = identity_matrix(R, 4)
    expected2[1, 1] = block2[1, 1]
    expected2[1, 3] = block2[1, 2]
    expected2[3, 1] = block2[2, 1]
    expected2[3, 3] = block2[2, 2]

    @test size(embedded2) == (4, 4)
    @test base_ring(embedded2) == R || base_ring(embedded2) === R
    @test embedded2 == expected2
    assert_identity_outside_block(embedded2, [1, 3], R)

    L, (u, v) = suslin_laurent_polynomial_ring(GF(2), ["u", "v"])
    block3 = matrix(L, [
        one(L)      u^-1          zero(L);
        v           one(L) + u^-1 * v  u;
        zero(L)     v^-1          one(L)
    ])

    embedded3 = block_embedding(block3, 5, [2, 4, 5])
    expected3 = identity_matrix(L, 5)
    for local_i in 1:3, local_j in 1:3
        expected3[[2, 4, 5][local_i], [2, 4, 5][local_j]] = block3[local_i, local_j]
    end

    @test size(embedded3) == (5, 5)
    @test base_ring(embedded3) == L || base_ring(embedded3) === L
    @test embedded3 == expected3
    @test embedded3[4, 4] == one(L) + u^-1 * v
    assert_identity_outside_block(embedded3, [2, 4, 5], L)
end

@testset "embedded factor sequences" begin
    R, (x,) = suslin_polynomial_ring(QQ, ["x"])
    small_factors = [
        elementary_matrix(2, 1, 2, x, R),
        elementary_matrix(2, 2, 1, x + 1, R),
    ]

    small_product = block_embedding_product(small_factors, R, 2)
    embedded_factors = embed_factor_sequence(small_factors, 5, [2, 5])

    @test length(embedded_factors) == length(small_factors)
    @test all(size(factor) == (5, 5) for factor in embedded_factors)
    @test block_embedding_product(embedded_factors, R, 5) == block_embedding(small_product, 5, [2, 5])

    final_factor = embed_factor_sequence([elementary_matrix(2, 1, 2, one(R), R)], 5, [2, 5])
    composed = compose_factor_sequences(embedded_factors[1:1], embedded_factors[2:end], final_factor)
    expected_sequence = vcat(embedded_factors[1:1], embedded_factors[2:end], final_factor)

    @test composed == expected_sequence
    @test block_embedding_product(composed, R, 5) == block_embedding_product(expected_sequence, R, 5)
    @test compose_factor_sequences(typeof(embedded_factors)(), embedded_factors) == embedded_factors
end

@testset "block embedding validation" begin
    R, (x,) = suslin_polynomial_ring(QQ, ["x"])
    S, (y,) = suslin_polynomial_ring(QQ, ["y"])

    block2 = matrix(R, [
        one(R)  x;
        zero(R) one(R)
    ])
    nonsquare = matrix(R, 2, 3, [
        one(R), zero(R), zero(R),
        zero(R), one(R), zero(R),
    ])

    @test_throws ArgumentError block_embedding(block2, 4, [1, 1])
    @test_throws ArgumentError block_embedding(block2, 4, [0, 2])
    @test_throws ArgumentError block_embedding(block2, 4, [1, 5])
    @test_throws DimensionMismatch block_embedding(nonsquare, 4, [1, 2])
    @test_throws DimensionMismatch block_embedding(block2, 4, [1])
    @test_throws DimensionMismatch block_embedding(block2, 1, [1, 2])

    factors = [elementary_matrix(2, 1, 2, x, R)]
    mixed_dimension = [elementary_matrix(2, 1, 2, x, R), elementary_matrix(3, 1, 2, x, R)]
    mixed_parent = [elementary_matrix(2, 1, 2, x, R), elementary_matrix(2, 1, 2, y, S)]

    @test_throws ArgumentError embed_factor_sequence([], 4, [1, 2])
    @test_throws DimensionMismatch embed_factor_sequence([nonsquare], 4, [1, 2])
    @test_throws DimensionMismatch embed_factor_sequence(factors, 4, [1])
    @test_throws ArgumentError embed_factor_sequence(mixed_dimension, 4, [1, 2])
    @test_throws ArgumentError embed_factor_sequence(mixed_parent, 4, [1, 2])

    @test_throws ArgumentError compose_factor_sequences()
    @test_throws ArgumentError compose_factor_sequences(typeof(factors)())
    @test_throws ArgumentError compose_factor_sequences(factors, mixed_dimension)
    @test_throws ArgumentError compose_factor_sequences(factors, mixed_parent)
end
```

Modify `test/runtests.jl` by adding the new expert file after
`expert/elementary_matrices.jl`:

```julia
"expert/block_embeddings.jl",
```

Modify `test/public/api_surface.jl` in the public API testset by adding:

```julia
@test isdefined(Suslin, :block_embedding)
@test isdefined(Suslin, :embed_factor_sequence)
@test isdefined(Suslin, :compose_factor_sequences)
@test Suslin.block_embedding === block_embedding
@test Suslin.embed_factor_sequence === embed_factor_sequence
@test Suslin.compose_factor_sequences === compose_factor_sequences
```

- [ ] **Step 2: Run the new expert test and confirm the expected red failure**

Run:

```bash
julia --project=. -e 'include("test/expert/block_embeddings.jl")'
```

Expected: FAIL/UndefVarError because `block_embedding` is not defined yet. This
is the required TDD red step.

- [ ] **Step 3: Implement the helper layer and exports**

Modify `src/Suslin.jl` exports:

```julia
export block_embedding
export embed_factor_sequence
export compose_factor_sequences
```

Replace `src/core/elementary_matrices.jl` with:

```julia
function elementary_matrix(n::Int, i::Int, j::Int, a, R)
    i == j && throw(ArgumentError("i and j must differ"))

    E = identity_matrix(R, n)
    coerced_a = _coerce_into_ring(R, a, "a")
    E[i, j] = coerced_a
    return E
end

function _require_square_matrix(M, label::AbstractString)
    nrows(M) == ncols(M) || throw(DimensionMismatch("$label must be square"))
    return nrows(M)
end

function _same_base_ring(left, right)::Bool
    return left == right || left === right
end

function _embedding_indices(n::Int, block_size::Int, indices)
    block_size <= n || throw(DimensionMismatch("target size must be at least the block size"))
    length(indices) == block_size || throw(DimensionMismatch("number of indices must match the block size"))

    result = Int[]
    seen = Set{Int}()
    for index in indices
        index isa Integer || throw(ArgumentError("indices must be integers"))
        idx = Int(index)
        1 <= idx <= n || throw(ArgumentError("indices must be between 1 and the target size"))
        idx in seen && throw(ArgumentError("indices must be distinct"))
        push!(seen, idx)
        push!(result, idx)
    end
    return result
end

function block_embedding(block, n::Int, indices)
    block_size = _require_square_matrix(block, "block")
    target_indices = _embedding_indices(n, block_size, indices)

    R = base_ring(block)
    embedded = identity_matrix(R, n)
    for local_i in 1:block_size, local_j in 1:block_size
        embedded[target_indices[local_i], target_indices[local_j]] = block[local_i, local_j]
    end
    return embedded
end

function _factor_shape_and_ring(factor, label::AbstractString)
    factor_size = _require_square_matrix(factor, label)
    return factor_size, base_ring(factor)
end

function _require_matching_factor(factor, expected_size::Int, expected_ring, label::AbstractString)
    factor_size, factor_ring = _factor_shape_and_ring(factor, label)
    factor_size == expected_size || throw(ArgumentError("$label must have the same size as the first factor"))
    _same_base_ring(factor_ring, expected_ring) || throw(ArgumentError("$label must have the same base ring as the first factor"))
    return factor
end

function embed_factor_sequence(factors, n::Int, indices)
    collected = collect(factors)
    isempty(collected) && throw(ArgumentError("factor sequence must be nonempty"))

    factor_size, factor_ring = _factor_shape_and_ring(first(collected), "factor")
    _embedding_indices(n, factor_size, indices)
    for factor in Iterators.drop(collected, 1)
        _require_matching_factor(factor, factor_size, factor_ring, "factor")
    end

    return [block_embedding(factor, n, indices) for factor in collected]
end

function compose_factor_sequences(sequences...)
    isempty(sequences) && throw(ArgumentError("at least one factor sequence is required"))

    composed = nothing
    expected_size = nothing
    expected_ring = nothing

    for sequence in sequences
        for factor in sequence
            if composed === nothing
                expected_size, expected_ring = _factor_shape_and_ring(factor, "factor")
                composed = typeof(factor)[factor]
            else
                _require_matching_factor(factor, expected_size, expected_ring, "factor")
                push!(composed, factor)
            end
        end
    end

    composed === nothing && throw(ArgumentError("at least one factor is required"))
    return composed
end
```

- [ ] **Step 4: Run the focused test and confirm green**

Run:

```bash
julia --project=. -e 'include("test/expert/block_embeddings.jl")'
```

Expected: PASS with all block-embedding tests passing.

- [ ] **Step 5: Run public API coverage**

Run:

```bash
julia --project=. -e 'include("test/public/api_surface.jl")'
```

Expected: PASS, including the new exported helper checks.

- [ ] **Step 6: Run package and full-suite verification**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
julia --project=. test/runtests.jl all
```

Expected: `Pkg.test()` passes the default public/internal groups, and the full
suite passes public, internal, and expert groups including
`expert/block_embeddings.jl`.

- [ ] **Step 7: Commit**

Run:

```bash
git add src/Suslin.jl src/core/elementary_matrices.jl test/expert/block_embeddings.jl test/runtests.jl test/public/api_surface.jl
git commit -m "feat: add block embedding helpers"
```

Expected: one feature commit containing helper implementation and tests.

## Plan Self-Review

- Spec coverage: the plan adds exact 2x2/3x3 block embeddings, embedded factor
  sequence composition, polynomial and Laurent examples, and all negative
  controls requested by issue #11.
- Placeholder scan: no placeholder text remains.
- Type consistency: the produced public helpers match the design names and the
  public API assertions.

## Automatic Execution Choice

This Agent Desk run is non-interactive. When the writing-plans workflow asks for
an execution approach, choose **Subagent-Driven (recommended)** because it is the
option marked recommended. With one tightly coupled implementation task, dispatch
one worker subagent for Task 1, then run the required task review and final
branch review before finishing.
