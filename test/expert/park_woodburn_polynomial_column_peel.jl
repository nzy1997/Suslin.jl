using Test
using Suslin
using Oscar

const PARK_WOODBURN_POLY_PEEL_CATALOG_PATH =
    joinpath(@__DIR__, "..", "fixtures", "park_woodburn_polynomial_cases.jl")

function _pw_poly_peel_product(factors, R, n::Int)
    product = identity_matrix(R, n)
    for factor in factors
        product *= factor
    end
    return product
end

function _pw_poly_peel_target_column(R, n::Int)
    target = zero_matrix(R, n, 1)
    target[n, 1] = one(R)
    return target
end

function _pw_poly_replace_step(
        step;
        dimension = step.dimension,
        input_matrix = step.input_matrix,
        last_column = step.last_column,
        left_factors = step.left_factors,
        after_left_matrix = step.after_left_matrix,
        right_factors = step.right_factors,
        peeled_matrix = step.peeled_matrix,
        next_block = step.next_block)
    return Suslin.PolynomialColumnPeelStep(
        dimension,
        input_matrix,
        last_column,
        left_factors,
        after_left_matrix,
        right_factors,
        peeled_matrix,
        next_block,
    )
end

function _pw_poly_replace_certificate(
        cert;
        original_matrix = cert.original_matrix,
        peel_steps = cert.peel_steps,
        final_block = cert.final_block,
        final_certificate = cert.final_certificate,
        final_factors = cert.final_factors,
        factors = cert.factors,
        product = cert.product,
        verification = cert.verification)
    return Suslin.PolynomialColumnPeelCertificate(
        original_matrix,
        peel_steps,
        final_block,
        final_certificate,
        final_factors,
        factors,
        product,
        verification,
    )
end

function _pw_poly_assert_step(step)
    R = base_ring(step.input_matrix)
    left_product = _pw_poly_peel_product(step.left_factors, R, step.dimension)
    right_product = _pw_poly_peel_product(step.right_factors, R, step.dimension)
    recorded_column = matrix(R, step.dimension, 1, step.last_column)
    @test left_product * recorded_column == _pw_poly_peel_target_column(R, step.dimension)
    @test left_product * step.input_matrix == step.after_left_matrix
    @test step.after_left_matrix * right_product == step.peeled_matrix
    @test step.peeled_matrix == block_embedding(step.next_block, step.dimension, collect(1:(step.dimension - 1)))
end

function _pw_poly_assert_real_peel_certificate(cert, A)
    R = base_ring(A)
    @test !isempty(cert.peel_steps)
    @test Suslin._verify_polynomial_column_peel_certificate(cert)
    @test cert.product == A
    @test _pw_poly_peel_product(cert.factors, R, nrows(A)) == A
    @test verify_factorization(A, cert.factors)
    for step in cert.peel_steps
        _pw_poly_assert_step(step)
    end
end

function _pw_poly_corrupt_last_column(cert)
    corrupted = collect(cert.peel_steps)
    step = first(corrupted)
    bad_column = copy(step.last_column)
    bad_column[1] += one(base_ring(step.input_matrix))
    corrupted[1] = _pw_poly_replace_step(step; last_column = bad_column)
    return _pw_poly_replace_certificate(cert; peel_steps = corrupted, product = cert.product)
end

function _pw_poly_corrupt_left_factor(cert)
    corrupted = collect(cert.peel_steps)
    step = first(corrupted)
    R = base_ring(step.input_matrix)
    bad_left = vcat(
        typeof(identity_matrix(R, step.dimension))[elementary_matrix(step.dimension, 1, step.dimension, one(R), R)],
        copy(step.left_factors),
    )
    corrupted[1] = _pw_poly_replace_step(step; left_factors = bad_left)
    return _pw_poly_replace_certificate(cert; peel_steps = corrupted, product = cert.product)
end

function _pw_poly_corrupt_right_factor(cert)
    corrupted = collect(cert.peel_steps)
    step = first(corrupted)
    R = base_ring(step.input_matrix)
    bad_right = vcat(
        copy(step.right_factors),
        typeof(identity_matrix(R, step.dimension))[elementary_matrix(step.dimension, step.dimension, 1, one(R), R)],
    )
    corrupted[1] = _pw_poly_replace_step(step; right_factors = bad_right)
    return _pw_poly_replace_certificate(cert; peel_steps = corrupted, product = cert.product)
end

function _pw_poly_corrupt_next_block(cert)
    corrupted = collect(cert.peel_steps)
    step = first(corrupted)
    bad_block = copy(step.next_block)
    bad_block[1, 1] += one(base_ring(step.input_matrix))
    corrupted[1] = _pw_poly_replace_step(step; next_block = bad_block)
    return _pw_poly_replace_certificate(cert; peel_steps = corrupted, product = cert.product)
end

@testset "Park-Woodburn ordinary polynomial column-peel certificates" begin
    if !isdefined(Main, :ParkWoodburnPolynomialFixtureCatalog)
        include(PARK_WOODBURN_POLY_PEEL_CATALOG_PATH)
    end
    entries = ParkWoodburnPolynomialFixtureCatalog.cases_by_id()

    recursive_entry = entries["pw-poly-recursive-column-peel-sl3-qq"]
    recursive_cert = Suslin._polynomial_column_peel_certificate(recursive_entry.matrix)
    @test recursive_cert.final_block == entries[recursive_entry.provenance.final_case_id].matrix
    @test recursive_cert.final_certificate.route == :fast_local_sl3
    _pw_poly_assert_real_peel_certificate(recursive_cert, recursive_entry.matrix)

    block_recursive_entry = entries["pw-poly-recursive-column-peel-sln-block-qq"]
    block_recursive_cert = Suslin._polynomial_column_peel_certificate(block_recursive_entry.matrix)
    @test block_recursive_cert.final_block == entries[block_recursive_entry.provenance.final_case_id].matrix
    @test block_recursive_cert.final_certificate.route == :disjoint_local_blocks
    _pw_poly_assert_real_peel_certificate(block_recursive_cert, block_recursive_entry.matrix)

    route_cert = Suslin._polynomial_factorization_route_certificate(
        recursive_entry.matrix;
        route = :recursive_column_peel,
    )
    @test route_cert.route == :recursive_column_peel
    @test route_cert.evidence isa Suslin.PolynomialColumnPeelCertificate
    @test Suslin._verify_polynomial_factorization_route_certificate(route_cert)

    staged_entry = entries["pw-poly-recursive-column-peel-gf2"]
    @test_throws ArgumentError Suslin._polynomial_column_peel_certificate(staged_entry.matrix)

    bad_last_column = _pw_poly_corrupt_last_column(recursive_cert)
    @test !Suslin._verify_polynomial_column_peel_certificate(bad_last_column)

    bad_left_factor = _pw_poly_corrupt_left_factor(recursive_cert)
    @test !Suslin._verify_polynomial_column_peel_certificate(bad_left_factor)

    bad_right_factor = _pw_poly_corrupt_right_factor(recursive_cert)
    @test !Suslin._verify_polynomial_column_peel_certificate(bad_right_factor)

    bad_next_block = _pw_poly_corrupt_next_block(recursive_cert)
    @test !Suslin._verify_polynomial_column_peel_certificate(bad_next_block)

    R = base_ring(recursive_entry.matrix)
    @test_throws ArgumentError Suslin._polynomial_column_peel_certificate(identity_matrix(R, 2))
    @test_throws ArgumentError Suslin._polynomial_column_peel_certificate(last(recursive_cert.peel_steps).peeled_matrix)
end
