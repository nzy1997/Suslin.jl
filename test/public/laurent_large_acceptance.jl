using Suslin
using Test
using Oscar

include("../fixtures/laurent_large_acceptance_cases.jl")
using .LaurentLargeAcceptanceCases

function _acceptance_product(factors, R, n::Int)
    product = identity_matrix(R, n)
    for factor in factors
        product *= factor
    end
    return product
end

function _assert_toricbuilder_contract_case(case)
    R = base_ring(case.matrix)
    normalization = normalize_laurent_gl_matrix(case.matrix)
    @test normalization.determinant_classification == :one
    @test normalization.normalized_matrix == case.matrix
    @test verify_laurent_gl_normalization(case.matrix, normalization)
    @test case.source_matrix * case.matrix == identity_matrix(R, case.size[1])

    corrupted = copy(case.matrix)
    corrupted[case.negative_control.row, case.negative_control.col] += one(R)
    @test case.source_matrix * corrupted != identity_matrix(R, case.size[1])
end

function _assert_large_factorization_case(case)
    R = base_ring(case.matrix)
    n = case.size[1]
    normalization = normalize_laurent_gl_matrix(case.matrix)
    @test normalization.determinant_classification == :one
    @test normalization.normalized_matrix == case.matrix
    @test verify_laurent_gl_normalization(case.matrix, normalization)

    reduction = reduce_sln_to_sl3(case.matrix; block_locations = case.block_locations)
    @test verify_sln_to_sl3_reduction(reduction)
    @test length(reduction.obligations) == length(case.block_locations)
    @test reduction.product == case.matrix

    factors = elementary_factorization(case.matrix)
    @test !isempty(factors)
    @test _acceptance_product(factors, R, n) == case.matrix
    @test verify_factorization(case.matrix, factors)

    corrupted_factors = copy(factors)
    corrupted_factors[1] = identity_matrix(R, n)
    @test !verify_factorization(case.matrix, corrupted_factors)
end

@testset "large Laurent acceptance catalog" begin
    catalog = LaurentLargeAcceptanceCases.acceptance_catalog()
    cases = catalog.cases

    @test any(case -> case.kind == :toricbuilder_normalized_contract, cases)
    @test any(case -> case.size == (40, 40), cases)
    @test any(case -> case.size == (48, 48), cases)

    for case in cases
        if case.expected_path == :normalized_contract
            _assert_toricbuilder_contract_case(case)
        elseif case.expected_path == :elementary_factorization
            _assert_large_factorization_case(case)
        else
            error("unknown acceptance path $(case.expected_path) for $(case.id)")
        end
    end
end
