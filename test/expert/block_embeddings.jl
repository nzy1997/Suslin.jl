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
