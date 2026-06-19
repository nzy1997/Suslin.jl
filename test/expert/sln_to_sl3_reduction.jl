using Test
using Suslin
using Oscar

function _issue15_product(factors, R, n::Int)
    product = identity_matrix(R, n)
    for factor in factors
        product *= factor
    end
    return product
end

function _issue15_local_block(p, q, r, s, R)
    return matrix(R, [
        p q zero(R);
        r s zero(R);
        zero(R) zero(R) one(R)
    ])
end

function _issue15_supported_matrix(R, blocks, n::Int)
    product = identity_matrix(R, n)
    for (block, indices) in blocks
        product *= block_embedding(block, n, indices)
    end
    return product
end

function _issue15_captured_error(f)
    try
        f()
        return nothing
    catch err
        return err
    end
end

function _issue15_assert_reduction(A, expected_obligations::Int)
    reduction = reduce_sln_to_sl3(A)
    @test reduction isa SLNToSL3Reduction
    @test length(reduction.obligations) == expected_obligations
    @test verify_sln_to_sl3_reduction(reduction)
    @test verify_factorization(A, reduction.factors)
    @test _issue15_product(reduction.factors, base_ring(A), nrows(A)) == A
    @test all(obligation -> obligation isa SL3LocalObligation, reduction.obligations)
    @test all(obligation -> obligation.reassembly_data.embedded_product_ok, reduction.obligations)
    @test all(obligation -> obligation.reassembly_data.local_product_ok, reduction.obligations)
    return reduction
end

@testset "SL_n to local SL3 reduction supported examples" begin
    R, (X,) = Oscar.polynomial_ring(QQ, ["X"])

    block_a = _issue15_local_block(one(R) + X, one(R), X, one(R), R)
    block_b = _issue15_local_block(one(R), one(R) + X, X, one(R) + X + X^2, R)
    block_c = _issue15_local_block(one(R) + X, X, one(R), one(R), R)

    matrix6 = _issue15_supported_matrix(R, [(block_a, [1, 2, 3]), (block_b, [4, 5, 6])], 6)
    reduction6 = _issue15_assert_reduction(matrix6, 2)
    @test reduction6.obligations[1].block_location == [1, 2, 3]
    @test reduction6.obligations[2].block_location == [4, 5, 6]
    @test elementary_factorization(matrix6) == reduction6.factors

    dropped6 = reduction6.obligations[1].embedded_factors
    @test !verify_factorization(matrix6, dropped6)

    matrix8 = _issue15_supported_matrix(R, [(block_b, [1, 2, 3]), (block_c, [4, 5, 6])], 8)
    reduction8 = _issue15_assert_reduction(matrix8, 2)
    @test reduction8.obligations[1].block_location == [1, 2, 3]
    @test reduction8.obligations[2].block_location == [4, 5, 6]
    @test matrix8[7, 7] == one(R)
    @test matrix8[8, 8] == one(R)

    custom = _issue15_supported_matrix(R, [(block_a, [2, 4, 6])], 6)
    custom_reduction = reduce_sln_to_sl3(custom; block_locations = [[2, 4, 6]])
    @test verify_sln_to_sl3_reduction(custom_reduction)
    @test verify_factorization(custom, custom_reduction.factors)
end

@testset "SL_n to local SL3 reduction staged failures" begin
    R, (X,) = Oscar.polynomial_ring(QQ, ["X"])

    unsupported = identity_matrix(R, 6)
    unsupported[1, 4] = X
    unsupported[4, 1] = zero(R)
    unsupported_err = _issue15_captured_error(() -> reduce_sln_to_sl3(unsupported))
    @test unsupported_err isa ArgumentError
    @test occursin("staged SL_n to local SL_3 reduction failure", sprint(showerror, unsupported_err))
    @test !occursin("local SL_3 special-form recognition failed", sprint(showerror, unsupported_err))

    bad_locations_err = _issue15_captured_error(() -> reduce_sln_to_sl3(identity_matrix(R, 6); block_locations = [[1, 2, 2]]))
    @test bad_locations_err isa ArgumentError
    @test occursin("block locations", sprint(showerror, bad_locations_err))

    S, (X, Y) = Oscar.polynomial_ring(QQ, ["X", "Y"])
    multivariate = identity_matrix(S, 6)
    multivariate[1:3, 1:3] = _issue15_local_block(one(S) + X, one(S), X, one(S), S)
    multivariate_err = _issue15_captured_error(() -> reduce_sln_to_sl3(multivariate))
    @test multivariate_err isa ArgumentError
    @test occursin("ordinary polynomial reduction currently requires a univariate base ring", sprint(showerror, multivariate_err))
end
