using Suslin
using Test
using Oscar

include("../fixtures/toricbuilder_factor_toric_block_3.jl")

function _issue58_product(factors, R, n::Int)
    product = identity_matrix(R, n)
    for factor in factors
        product *= factor
    end
    return product
end

function _issue58_entry_by_role(fixture, role::AbstractString)
    return only(filter(entry -> entry.toricbuilder_role == role, fixture.cases))
end

function _issue58_assert_step_replay(step)
    R = base_ring(step.input_matrix)
    left_product = _issue58_product(step.left_factors, R, step.dimension)
    right_product = _issue58_product(step.right_factors, R, step.dimension)
    @test left_product * step.input_matrix * right_product == step.peeled_matrix
    @test step.peeled_matrix[step.dimension, step.dimension] == one(R)
    @test all(step.peeled_matrix[row, step.dimension] == zero(R) for row in 1:(step.dimension - 1))
    @test all(step.peeled_matrix[step.dimension, col] == zero(R) for col in 1:(step.dimension - 1))
    @test step.next_block == matrix(R, [
        step.peeled_matrix[row, col]
        for row in 1:(step.dimension - 1), col in 1:(step.dimension - 1)
    ])
end

function _issue58_corrupt_source_relation(entry)
    corrupted = copy(entry.source_matrix)
    corrupted[1, 1] += one(base_ring(corrupted))
    return corrupted
end

function _issue58_replace_first_factor_with_identity(factors, R, n::Int)
    corrupted = copy(factors)
    corrupted[1] = identity_matrix(R, n)
    return corrupted
end

function _issue58_assert_column_peel_entry(entry, expected_dimensions, expected_final_block)
    R = base_ring(entry.matrix)
    n = entry.size[1]

    @test entry.expected_suslin_status == :supported_column_peel
    @test entry.expected_suslin_path == :laurent_column_peel
    @test entry.source_matrix * entry.matrix == identity_matrix(R, n)

    certificate = Suslin._factor_laurent_sl_column_peel(entry.matrix)
    @test Suslin._verify_laurent_column_peel_replay(certificate)
    @test [step.dimension for step in certificate.peel_steps] == expected_dimensions
    @test certificate.final_block == expected_final_block
    for step in certificate.peel_steps
        _issue58_assert_step_replay(step)
    end

    factors = elementary_factorization(entry.matrix)
    @test !isempty(factors)
    @test verify_factorization(entry.matrix, factors)

    corrupted_factors = _issue58_replace_first_factor_with_identity(factors, R, n)
    @test !verify_factorization(entry.matrix, corrupted_factors)

    corrupted_source = _issue58_corrupt_source_relation(entry)
    @test corrupted_source * entry.matrix != identity_matrix(R, n)
end

@testset "ToricBuilder factor_toric_block column peel acceptance" begin
    fixture = ToricBuilderFactorToricBlock3Fixture.fixture()
    pinv = _issue58_entry_by_role(fixture, "Pinv")
    qinv = _issue58_entry_by_role(fixture, "Qinv")
    R = base_ring(pinv.matrix)
    x, y = gens(R)

    pinv_expected_final = identity_matrix(R, 2)
    qinv_expected_final = matrix(R, [
        y                  zero(R);
        x * y^-1 + y^-1   y^-1
    ])

    _issue58_assert_column_peel_entry(pinv, collect(8:-1:3), pinv_expected_final)
    _issue58_assert_column_peel_entry(qinv, collect(16:-1:3), qinv_expected_final)
end
