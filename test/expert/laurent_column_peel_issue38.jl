using Test
using Suslin
using Oscar

include("../fixtures/toricbuilder_issue38_cases.jl")

struct _Issue57BadFactorSequence
    factor
end

Base.iterate(sequence::_Issue57BadFactorSequence, state = 1) =
    state == 1 ? (sequence.factor, 2) : nothing
Base.length(::_Issue57BadFactorSequence) = 1
Base.:(==)(::_Issue57BadFactorSequence, ::Vector) = true

struct _Issue57FallbackProbe end

struct _Issue57FallbackCertificate
    factors
end

function Suslin.reduce_sln_to_sl3(::_Issue57FallbackProbe)
    throw(ArgumentError("probe reduction unsupported"))
end

function Suslin._factor_laurent_sl_column_peel(::_Issue57FallbackProbe)
    return _Issue57FallbackCertificate([])
end

function Suslin.verify_factorization(::_Issue57FallbackProbe, factors)::Bool
    return false
end

function _issue57_product(factors, R, n::Int)
    product = identity_matrix(R, n)
    for factor in factors
        product *= factor
    end
    return product
end

function _issue57_assert_step(step)
    R = base_ring(step.input_matrix)
    left_product = _issue57_product(step.left_factors, R, step.dimension)
    right_product = _issue57_product(step.right_factors, R, step.dimension)
    @test left_product * step.input_matrix * right_product == step.peeled_matrix
    @test step.peeled_matrix[step.dimension, step.dimension] == one(R)
    @test all(step.peeled_matrix[row, step.dimension] == zero(R) for row in 1:(step.dimension - 1))
    @test all(step.peeled_matrix[step.dimension, col] == zero(R) for col in 1:(step.dimension - 1))
    @test step.next_block == matrix(R, [
        step.peeled_matrix[row, col]
        for row in 1:(step.dimension - 1), col in 1:(step.dimension - 1)
    ])
end

function _issue57_corrupt_left_factor(certificate)
    corrupted = collect(certificate.peel_steps)
    first_step = first(corrupted)
    bad_left = copy(first_step.left_factors)
    R = base_ring(first_step.input_matrix)
    bad_left[1] = elementary_matrix(first_step.dimension, 1, 2, one(R), R) * bad_left[1]
    corrupted[1] = Suslin.LaurentColumnPeelStep(
        first_step.dimension,
        first_step.input_matrix,
        first_step.last_column,
        bad_left,
        first_step.after_left_matrix,
        first_step.right_factors,
        first_step.peeled_matrix,
        first_step.next_block,
    )
    return Suslin.LaurentColumnPeelFactorization(
        certificate.original_matrix,
        certificate.final_block,
        certificate.final_local_target,
        certificate.final_local_factors,
        certificate.final_factors,
        certificate.factors,
        certificate.product,
        corrupted,
        certificate.verification,
    )
end

function _issue57_delete_clearing_factor(certificate)
    corrupted = collect(certificate.peel_steps)
    first_step = first(corrupted)
    shortened = first_step.right_factors[1:(end - 1)]
    corrupted[1] = Suslin.LaurentColumnPeelStep(
        first_step.dimension,
        first_step.input_matrix,
        first_step.last_column,
        first_step.left_factors,
        first_step.after_left_matrix,
        shortened,
        first_step.peeled_matrix,
        first_step.next_block,
    )
    return Suslin.LaurentColumnPeelFactorization(
        certificate.original_matrix,
        certificate.final_block,
        certificate.final_local_target,
        certificate.final_local_factors,
        certificate.final_factors,
        certificate.factors[1:(end - 1)],
        certificate.product,
        corrupted,
        certificate.verification,
    )
end

function _issue57_corrupt_last_column(certificate)
    corrupted = collect(certificate.peel_steps)
    first_step = first(corrupted)
    bad_column = copy(first_step.last_column)
    R = base_ring(first_step.input_matrix)
    bad_column[1] += one(R)
    bad_left = copy(first_step.left_factors)
    bad_left[1] = elementary_matrix(first_step.dimension, 1, 2, one(R), R) * bad_left[1]
    corrupted[1] = Suslin.LaurentColumnPeelStep(
        first_step.dimension,
        first_step.input_matrix,
        bad_column,
        bad_left,
        first_step.after_left_matrix,
        first_step.right_factors,
        first_step.peeled_matrix,
        first_step.next_block,
    )
    return Suslin.LaurentColumnPeelFactorization(
        certificate.original_matrix,
        certificate.final_block,
        certificate.final_local_target,
        certificate.final_local_factors,
        certificate.final_factors,
        certificate.factors,
        certificate.product,
        corrupted,
        certificate.verification,
    )
end

function _issue57_corrupt_final_local_metadata(certificate)
    R = base_ring(certificate.original_matrix)
    return Suslin.LaurentColumnPeelFactorization(
        certificate.original_matrix,
        certificate.final_block,
        identity_matrix(R, 3),
        typeof(certificate.final_local_factors)(),
        typeof(certificate.final_factors)(),
        typeof(certificate.factors)(),
        certificate.product,
        certificate.peel_steps,
        certificate.verification,
    )
end

function _issue57_append_cancelling_factor_pair(certificate)
    mutated = Suslin.LaurentColumnPeelFactorization(
        certificate.original_matrix,
        certificate.final_block,
        certificate.final_local_target,
        certificate.final_local_factors,
        certificate.final_factors,
        certificate.factors,
        certificate.product,
        certificate.peel_steps,
        certificate.verification,
    )
    R = base_ring(mutated.original_matrix)
    n = nrows(mutated.original_matrix)
    push!(mutated.factors, elementary_matrix(n, 1, 2, one(R), R))
    push!(mutated.factors, elementary_matrix(n, 1, 2, -one(R), R))
    return mutated
end

function _issue57_certificate_tuple(certificate; overrides...)
    base = (
        original_matrix = certificate.original_matrix,
        final_block = certificate.final_block,
        final_local_target = certificate.final_local_target,
        final_local_factors = certificate.final_local_factors,
        final_factors = certificate.final_factors,
        factors = certificate.factors,
        product = certificate.product,
        peel_steps = certificate.peel_steps,
    )
    return merge(base, (; overrides...))
end

function _issue57_assert_malformed_replay_guards(certificate)
    first_step = first(certificate.peel_steps)
    R = base_ring(first_step.input_matrix)
    bad_factor = identity_matrix(R, 1)

    @test !Suslin._is_valid_laurent_column_peel_step_data(
        first_step.dimension,
        first_step.input_matrix,
        first_step.last_column,
        [bad_factor],
        first_step.after_left_matrix,
        first_step.right_factors,
        first_step.peeled_matrix,
        first_step.next_block,
    )

    @test !Suslin._is_valid_laurent_column_peel_step_data(
        first_step.dimension,
        first_step.input_matrix,
        first_step.last_column,
        first_step.left_factors,
        first_step.after_left_matrix,
        _Issue57BadFactorSequence(bad_factor),
        first_step.peeled_matrix,
        first_step.next_block,
    )

    invalid_step = Suslin.LaurentColumnPeelStep(
        first_step.dimension,
        nothing,
        first_step.last_column,
        first_step.left_factors,
        first_step.after_left_matrix,
        first_step.right_factors,
        first_step.peeled_matrix,
        first_step.next_block,
    )
    invalid_certificate = _issue57_certificate_tuple(
        certificate;
        peel_steps = Suslin.LaurentColumnPeelStep[invalid_step],
        final_block = nothing,
        final_local_target = nothing,
        final_factors = nothing,
        factors = [bad_factor],
    )
    invalid_verification = Suslin._laurent_column_peel_verification(invalid_certificate)
    @test !invalid_verification.overall_ok
    @test !invalid_verification.steps_ok
    @test !invalid_verification.final_metadata_ok
    @test !invalid_verification.final_local_ok
    @test !invalid_verification.final_factors_ok
    @test !invalid_verification.factor_sequence_ok
    @test !invalid_verification.product_ok
    @test !invalid_verification.factors_ok
end

function _issue57_assert_core(core, expected_final_block)
    certificate = Suslin._factor_laurent_sl_column_peel(core)

    @test certificate.final_block == expected_final_block
    @test length(certificate.factors) > 0
    @test verify_factorization(core, certificate.factors)
    @test verify_factorization(core, elementary_factorization(core))
    @test Suslin._verify_laurent_column_peel_replay(certificate)
    @test [step.dimension for step in certificate.peel_steps] == [6, 5, 4, 3]
    for step in certificate.peel_steps
        _issue57_assert_step(step)
    end

    corrupted_left = _issue57_corrupt_left_factor(certificate)
    @test !Suslin._verify_laurent_column_peel_replay(corrupted_left)
    @test !verify_factorization(core, corrupted_left.factors)

    deleted_clear = _issue57_delete_clearing_factor(certificate)
    @test !Suslin._verify_laurent_column_peel_replay(deleted_clear)
    @test !verify_factorization(core, deleted_clear.factors)

    corrupted_column = _issue57_corrupt_last_column(certificate)
    @test !Suslin._verify_laurent_column_peel_replay(corrupted_column)
    @test !verify_factorization(core, corrupted_column.factors)

    corrupted_final_local = _issue57_corrupt_final_local_metadata(certificate)
    @test !Suslin._verify_laurent_column_peel_replay(corrupted_final_local)
    @test !verify_factorization(core, corrupted_final_local.factors)

    appended_cancel = _issue57_append_cancelling_factor_pair(certificate)
    @test verify_factorization(core, appended_cancel.factors)
    @test !Suslin._verify_laurent_column_peel_replay(appended_cancel)

    _issue57_assert_malformed_replay_guards(certificate)
end

@testset "Issue 38 Laurent column peel" begin
    entry = only(ToricBuilderIssue38Cases.catalog().cases)
    R = entry.ring.object
    u, v = entry.ring.generators

    row_expected = matrix(R, [
        u      u * v;
        zero(R) u^-1
    ])
    column_expected = matrix(R, [
        v^-1   u * v;
        zero(R) v
    ])

    _issue57_assert_core(entry.normalizations.row.core, row_expected)
    _issue57_assert_core(entry.normalizations.column.core, column_expected)

    original_err = try
        elementary_factorization(entry.inputs.matrix)
        nothing
    catch err
        err
    end
    @test original_err isa ArgumentError
    @test occursin("Laurent GL_n normalization boundary", sprint(showerror, original_err))

    @test_throws ErrorException Suslin._laurent_sl_fallback_factorization(_Issue57FallbackProbe())
end
