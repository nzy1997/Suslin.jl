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

function _issue15_wrap_column_peel_matrix(final_block, tail_entries)
    R = base_ring(final_block)
    n = nrows(final_block) + 1
    length(tail_entries) == n - 1 || throw(ArgumentError("tail_entries must match final block size"))
    wrapped = block_embedding(final_block, n, collect(1:(n - 1)))
    for row in 1:(n - 1)
        wrapped[row, n] = tail_entries[row]
    end
    return wrapped
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

function _issue15_replace_reduction(
        reduction;
        ring = reduction.ring,
        size = reduction.size,
        original_matrix = reduction.original_matrix,
        normalized_matrix = reduction.normalized_matrix,
        normalization = reduction.normalization,
        obligations = reduction.obligations,
        factors = reduction.factors,
        product = reduction.product,
        verification = reduction.verification)
    return SLNToSL3Reduction(
        ring,
        size,
        original_matrix,
        normalized_matrix,
        normalization,
        obligations,
        factors,
        product,
        verification,
    )
end

function _issue15_replace_obligation(
        obligation;
        block_location = obligation.block_location,
        ring = obligation.ring,
        target_local_matrix = obligation.target_local_matrix,
        required_assumptions = obligation.required_assumptions,
        embedded_target = obligation.embedded_target,
        local_factors = obligation.local_factors,
        embedded_factors = obligation.embedded_factors,
        reassembly_data = obligation.reassembly_data)
    return SL3LocalObligation(
        block_location,
        ring,
        target_local_matrix,
        required_assumptions,
        embedded_target,
        local_factors,
        embedded_factors,
        reassembly_data,
    )
end

struct Issue15ExplodingEq end

Base.:(==)(::Issue15ExplodingEq, _) = throw(ArgumentError("issue 15 equality sentinel"))

struct Issue15InterruptMatrix end

Oscar.nrows(::Issue15InterruptMatrix) = throw(InterruptException())
Oscar.base_ring(::Issue15InterruptMatrix) = throw(InterruptException())

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

    tampered_obligation_factors = reduce_sln_to_sl3(matrix6)
    deleteat!(tampered_obligation_factors.obligations[1].embedded_factors, 1)
    @test !verify_sln_to_sl3_reduction(tampered_obligation_factors)

    tampered_obligations = reduce_sln_to_sl3(matrix6)
    deleteat!(tampered_obligations.obligations, 1)
    @test !verify_sln_to_sl3_reduction(tampered_obligations)

    tampered_metadata = reduce_sln_to_sl3(matrix6)
    tampered_metadata.obligations[1].block_location[1] = 2
    @test !verify_sln_to_sl3_reduction(tampered_metadata)

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

    peel_matrix = _issue15_wrap_column_peel_matrix(
        block_a,
        [X, X + one(R), X^2 + X],
    )
    @test nrows(peel_matrix) > 3
    peel_cert = Suslin._polynomial_column_peel_certificate(peel_matrix)
    @test Suslin._verify_polynomial_column_peel_certificate(peel_cert)
    @test length(peel_cert.peel_steps) == 1
    first_step = only(peel_cert.peel_steps)
    @test first_step.dimension == nrows(peel_matrix)
    @test first_step.left_certificate isa Suslin.ECPColumnReductionCertificate
    @test Suslin.verify_ecp_column_reduction(first_step.left_certificate)
    @test first_step.left_certificate.original_column == first_step.last_column
    @test first_step.left_certificate.factors == first_step.left_factors
    @test first_step.left_certificate.final_column ==
          Suslin._column_peel_target_column(R, first_step.dimension)
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

    det_not_one = identity_matrix(R, 4)
    det_not_one[1, 1] = one(R) + X
    det_err = _issue15_captured_error(() -> Suslin._polynomial_column_peel_certificate(det_not_one))
    @test det_err isa ArgumentError
    @test occursin("determinant-one input", sprint(showerror, det_err))

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

@testset "SL_n to local SL3 reduction defensive verification branches" begin
    R, (X,) = Oscar.polynomial_ring(QQ, ["X"])

    block = _issue15_local_block(one(R) + X, one(R), X, one(R), R)
    matrix6 = _issue15_supported_matrix(R, [(block, [1, 2, 3])], 6)
    reduction = reduce_sln_to_sl3(matrix6)

    L, _ = suslin_laurent_polynomial_ring(QQ, ["t"])
    laurent_identity = identity_matrix(L, 3)
    laurent_reduction = reduce_sln_to_sl3(laurent_identity)
    @test laurent_reduction isa SLNToSL3Reduction
    @test length(laurent_reduction.obligations) == 0
    @test verify_sln_to_sl3_reduction(laurent_reduction)
    @test laurent_reduction.verification.overall_ok
    @test laurent_reduction.product == laurent_identity

    t = gen(L, 1)
    laurent_corrected = matrix(L, [
        t zero(L) zero(L);
        zero(L) one(L) zero(L);
        zero(L) zero(L) one(L)
    ])
    corrected_err = _issue15_captured_error(() -> reduce_sln_to_sl3(laurent_corrected))
    @test corrected_err isa ArgumentError
    @test occursin("Laurent determinant correction", sprint(showerror, corrected_err))

    malformed_local = identity_matrix(R, 6)
    malformed_local[1, 3] = X
    malformed_local_err = _issue15_captured_error(() -> reduce_sln_to_sl3(malformed_local))
    @test malformed_local_err isa ArgumentError
    @test occursin("failed to solve local SL_3 obligation on block [1, 2, 3]", sprint(showerror, malformed_local_err))

    wrong_size_factor = identity_matrix(R, 2)
    wrong_factor_reduction = _issue15_replace_reduction(reduction; factors = [wrong_size_factor])
    wrong_factor_verification = Suslin._sln_to_sl3_reduction_verification(wrong_factor_reduction)
    @test !wrong_factor_verification.factors_ok
    @test !wrong_factor_verification.obligation_factors_ok
    @test !verify_sln_to_sl3_reduction(wrong_factor_reduction)

    exploding_product_reduction = _issue15_replace_reduction(reduction; product = Issue15ExplodingEq())
    exploding_product_verification = Suslin._sln_to_sl3_reduction_verification(exploding_product_reduction)
    @test !exploding_product_verification.product_ok
    @test !exploding_product_verification.original_reconstruction_ok

    bad_location_obligation = _issue15_replace_obligation(
        reduction.obligations[1];
        block_location = [1, 2, 7],
    )
    bad_location_reduction = _issue15_replace_reduction(reduction; obligations = [bad_location_obligation])
    bad_location_verification = Suslin._sln_to_sl3_reduction_verification(bad_location_reduction)
    @test !bad_location_verification.obligation_locations_ok
    @test !bad_location_verification.obligations_ok

    empty_embedded_obligation = _issue15_replace_obligation(
        reduction.obligations[1];
        embedded_factors = typeof(identity_matrix(R, 6))[],
    )
    empty_embedded_reduction = _issue15_replace_reduction(reduction; obligations = [empty_embedded_obligation])
    empty_embedded_verification = Suslin._sln_to_sl3_reduction_verification(empty_embedded_reduction)
    @test !empty_embedded_verification.obligation_factors_ok

    interrupt_obligation_reduction = _issue15_replace_reduction(
        reduction;
        normalized_matrix = Issue15InterruptMatrix(),
    )
    @test_throws InterruptException Suslin._sln_to_sl3_reduction_verification(interrupt_obligation_reduction)

    normalization_interrupt_reduction = _issue15_replace_reduction(
        reduction;
        original_matrix = Issue15InterruptMatrix(),
        normalization = (;),
        obligations = SL3LocalObligation[],
        factors = typeof(identity_matrix(R, 6))[],
        product = identity_matrix(R, 6),
    )
    @test_throws InterruptException Suslin._sln_to_sl3_reduction_verification(normalization_interrupt_reduction)
end
