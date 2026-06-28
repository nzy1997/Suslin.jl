using Test
using Suslin
using Oscar

include("../fixtures/laurent_lazy_determinant_cases.jl")

function _issue156_fixture(id::AbstractString)
    catalog = LaurentLazyDeterminantCases.catalog()
    return only(filter(entry -> entry.id == id, catalog.cases))
end

function _issue156_blockdiag_deferred(deferred, original_dimension::Int)
    return block_embedding(deferred, original_dimension, collect(1:nrows(deferred)))
end

function _issue156_replace_first_factor_with_identity(factors, R, n::Int)
    tampered = copy(factors)
    tampered[1] = identity_matrix(R, n)
    return tampered
end

function _issue156_tamper_first_step_factor(certificate)
    steps = collect(certificate.peel_steps)
    first_step = first(steps)
    R = base_ring(first_step.input_matrix)

    if !isempty(first_step.left_factors)
        bad_left = _issue156_replace_first_factor_with_identity(
            first_step.left_factors,
            R,
            first_step.dimension,
        )
        steps[1] = Suslin.LaurentColumnPeelStep(
            first_step.dimension,
            first_step.input_matrix,
            first_step.last_column,
            bad_left,
            first_step.after_left_matrix,
            first_step.right_factors,
            first_step.peeled_matrix,
            first_step.next_block,
        )
    elseif !isempty(first_step.right_factors)
        bad_right = _issue156_replace_first_factor_with_identity(
            first_step.right_factors,
            R,
            first_step.dimension,
        )
        steps[1] = Suslin.LaurentColumnPeelStep(
            first_step.dimension,
            first_step.input_matrix,
            first_step.last_column,
            first_step.left_factors,
            first_step.after_left_matrix,
            bad_right,
            first_step.peeled_matrix,
            first_step.next_block,
        )
    else
        error("fixture did not produce a tamperable peel factor")
    end

    return Suslin.LaurentDeterminantDeferredPeelCertificate(
        certificate.original_matrix,
        steps,
        certificate.deferred_submatrix,
        certificate.determinant_source,
        certificate.verification,
    )
end

@testset "determinant-deferred lazy Laurent peel certificate" begin
    entry = _issue156_fixture("determinant-one-triangular")
    A = entry.inputs.matrix
    n = nrows(A)

    deferred = Suslin._laurent_determinant_deferred_peel_certificate(A)
    sl_certificate = Suslin._factor_laurent_sl_column_peel(A)

    @test deferred.original_matrix == A
    @test deferred.determinant_source == :deferred_submatrix
    @test !isempty(deferred.peel_steps)
    @test nrows(deferred.deferred_submatrix) < n
    @test ncols(deferred.deferred_submatrix) < ncols(A)
    @test deferred.deferred_submatrix == first(sl_certificate.peel_steps).next_block
    @test first(deferred.peel_steps).left_factors == first(sl_certificate.peel_steps).left_factors
    @test first(deferred.peel_steps).right_factors == first(sl_certificate.peel_steps).right_factors
    @test Suslin._verify_laurent_determinant_deferred_peel_replay(deferred)
    @test deferred.left_product * A * deferred.right_product == deferred.target_matrix
    @test deferred.target_matrix == _issue156_blockdiag_deferred(deferred.deferred_submatrix, n)
    @test deferred.verification.overall_ok
    @test deferred.verification.relation_ok
    @test deferred.verification.replay_metadata_ok

    first_step = first(deferred.peel_steps)
    @test !isempty(first_step.left_factors) || !isempty(first_step.right_factors)
    tampered = _issue156_tamper_first_step_factor(deferred)
    replay_ok = try
        Suslin._verify_laurent_determinant_deferred_peel_replay(tampered)
    catch err
        err isa InterruptException && rethrow()
        false
    end
    @test !replay_ok
end

@testset "determinant-deferred lazy Laurent peel certificate min_steps" begin
    entry = _issue156_fixture("issue-38-q-block-lazy-determinant")
    A = entry.normalizations.row.core
    n = nrows(A)

    deferred = Suslin._laurent_determinant_deferred_peel_certificate(A; min_steps = 2)

    @test deferred.original_matrix == A
    @test deferred.determinant_source == :deferred_submatrix
    @test length(deferred.peel_steps) == 2
    @test nrows(deferred.deferred_submatrix) < n
    @test ncols(deferred.deferred_submatrix) < ncols(A)
    @test deferred.deferred_submatrix == last(deferred.peel_steps).next_block
    @test deferred.left_product * A * deferred.right_product == deferred.target_matrix
    @test deferred.target_matrix ==
        block_embedding(deferred.deferred_submatrix, n, collect(1:nrows(deferred.deferred_submatrix)))
    @test Suslin._verify_laurent_determinant_deferred_peel_replay(deferred)
    @test deferred.verification.overall_ok
    @test deferred.verification.relation_ok
    @test deferred.verification.replay_metadata_ok
end
