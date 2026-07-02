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
        left_certificate = step.left_certificate,
        ecp_evidence = step.ecp_evidence,
        ecp_route_provenance = step.ecp_route_provenance,
        after_left_matrix = step.after_left_matrix,
        right_factors = step.right_factors,
        right_clearing_coefficients = step.right_clearing_coefficients,
        peeled_matrix = step.peeled_matrix,
        next_block = step.next_block,
        block_embedding_indices = step.block_embedding_indices,
        determinant_metadata = step.determinant_metadata,
        descent_metadata = step.descent_metadata,
        verification = step.verification)
    return Suslin.PolynomialColumnPeelStep(
        dimension,
        input_matrix,
        last_column,
        left_factors,
        left_certificate,
        ecp_evidence,
        ecp_route_provenance,
        after_left_matrix,
        right_factors,
        right_clearing_coefficients,
        peeled_matrix,
        next_block,
        block_embedding_indices,
        determinant_metadata,
        descent_metadata,
        verification,
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

function _pw_poly_replace_route_certificate(
        cert;
        matrix = cert.matrix,
        route = cert.route,
        factors = cert.factors,
        product = cert.product,
        evidence = cert.evidence,
        status = cert.status,
        verification = cert.verification)
    return Suslin.PolynomialFactorizationRouteCertificate(
        matrix,
        route,
        factors,
        product,
        evidence,
        status,
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
    @test step.left_certificate isa Suslin.ECPColumnReductionCertificate
    @test step.ecp_evidence == step.left_certificate
    @test Suslin.verify_ecp_column_reduction(step.left_certificate)
    @test step.ecp_route_provenance.verifier == :verify_ecp_column_reduction
    @test step.ecp_route_provenance.status == :verified
    @test step.left_certificate.original_column == step.last_column
    @test step.left_certificate.factors == step.left_factors
    @test step.left_certificate.final_column == _pw_poly_peel_target_column(R, step.dimension)
    @test step.right_clearing_coefficients ==
          tuple((step.after_left_matrix[step.dimension, col] for col in 1:(step.dimension - 1))...)
    @test det(step.next_block) == one(R)
    @test step.determinant_metadata.next_block_determinant == one(R)
    @test step.descent_metadata.next_dimension == step.dimension - 1
    @test Suslin._polynomial_column_peel_step_verification(step).overall_ok
    left_stage = step.left_certificate.stages[end]
    if hasproperty(left_stage, :route_metadata)
        @test hasproperty(left_stage.route_metadata, :route)
    end
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

function _pw_poly_tamper_ecp_certificate(ecp_certificate)
    factors = copy(ecp_certificate.factors)
    factors[1] = identity_matrix(ecp_certificate.ring, length(ecp_certificate.original_column))
    return Suslin.ECPColumnReductionCertificate(
        ecp_certificate.original_column,
        ecp_certificate.ring,
        ecp_certificate.stages,
        factors,
        ecp_certificate.final_column,
        ecp_certificate.verification,
    )
end

function _pw_poly_corrupt_ecp_certificate(cert)
    corrupted = collect(cert.peel_steps)
    step = first(corrupted)
    bad_ecp = _pw_poly_tamper_ecp_certificate(step.ecp_evidence)
    corrupted[1] = _pw_poly_replace_step(step; ecp_evidence = bad_ecp, left_certificate = bad_ecp)
    return _pw_poly_replace_certificate(cert; peel_steps = corrupted, product = cert.product)
end

function _pw_poly_corrupt_ecp_route_provenance(cert)
    corrupted = collect(cert.peel_steps)
    step = first(corrupted)
    bad_provenance = merge(step.ecp_route_provenance, (; route = :tampered_route))
    corrupted[1] = _pw_poly_replace_step(step; ecp_route_provenance = bad_provenance)
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

function _pw_poly_corrupt_right_clearing_coefficient(cert)
    corrupted = collect(cert.peel_steps)
    step = first(corrupted)
    coefficients = collect(step.right_clearing_coefficients)
    coefficients[1] += one(base_ring(step.input_matrix))
    corrupted[1] = _pw_poly_replace_step(step; right_clearing_coefficients = tuple(coefficients...))
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

struct _PWPolyBadFactorList
end

Base.iterate(::_PWPolyBadFactorList, state = 1) = state == 1 ? ((;), 2) : nothing
Base.:(==)(::_PWPolyBadFactorList, other) = true
Base.:(==)(other, ::_PWPolyBadFactorList) = true

@testset "Park-Woodburn ordinary polynomial column-peel certificates" begin
    if !isdefined(Main, :ParkWoodburnPolynomialFixtureCatalog)
        include(PARK_WOODBURN_POLY_PEEL_CATALOG_PATH)
    end
    entries = ParkWoodburnPolynomialFixtureCatalog.cases_by_id()

    recursive_entry = entries["pw-poly-recursive-column-peel-sl3-qq"]
    recursive_cert = Suslin._polynomial_column_peel_certificate(recursive_entry.matrix)
    explicit_recursive_cert = Suslin._polynomial_column_peel_certificate(
        recursive_entry.matrix;
        final_route = :fast_local_sl3,
    )
    @test explicit_recursive_cert.final_certificate.route == :fast_local_sl3
    @test recursive_cert.final_block == entries[recursive_entry.provenance.final_case_id].matrix
    @test recursive_cert.final_certificate.route == :fast_local_sl3
    @test recursive_cert.verification.left_certificates_ok
    _pw_poly_assert_real_peel_certificate(recursive_cert, recursive_entry.matrix)
    first_recursive_step = first(recursive_cert.peel_steps)
    legacy_step = Suslin.PolynomialColumnPeelStep(
        first_recursive_step.dimension,
        first_recursive_step.input_matrix,
        first_recursive_step.last_column,
        first_recursive_step.left_factors,
        first_recursive_step.after_left_matrix,
        first_recursive_step.right_factors,
        first_recursive_step.peeled_matrix,
        first_recursive_step.next_block,
    )
    @test legacy_step.left_certificate === nothing
    stripped_steps = copy(recursive_cert.peel_steps)
    stripped_steps[1] = legacy_step
    stripped_cert = _pw_poly_replace_certificate(recursive_cert; peel_steps = stripped_steps)
    @test !Suslin._polynomial_column_peel_core_verification(stripped_cert).left_certificates_ok
    @test !Suslin._verify_polynomial_column_peel_certificate(stripped_cert)

    certified_legacy_step = Suslin.PolynomialColumnPeelStep(
        first_recursive_step.dimension,
        first_recursive_step.input_matrix,
        first_recursive_step.last_column,
        first_recursive_step.left_factors,
        first_recursive_step.left_certificate,
        first_recursive_step.after_left_matrix,
        first_recursive_step.right_factors,
        first_recursive_step.peeled_matrix,
        first_recursive_step.next_block,
    )
    @test certified_legacy_step.left_certificate == first_recursive_step.left_certificate
    @test certified_legacy_step.ecp_evidence == first_recursive_step.left_certificate
    @test certified_legacy_step.ecp_route_provenance == first_recursive_step.ecp_route_provenance
    @test Suslin._polynomial_column_peel_step_verification(certified_legacy_step).overall_ok

    route_R, (route_x, route_y) = Oscar.polynomial_ring(GF(2), ["route_x", "route_y"])
    route_column = [
        route_x + route_y^2,
        route_x * route_y + route_x + one(route_R),
        route_x^2 + route_x * route_y + route_y + one(route_R),
    ]
    route_metadata_cert = Suslin.ecp_column_reduction_certificate(route_column, route_R)
    @test Suslin.verify_ecp_column_reduction(route_metadata_cert)
    @test hasproperty(route_metadata_cert.stages[end], :route_metadata)
    @test route_metadata_cert.stages[end].route_metadata.route == :general_ecp_pipeline
    @test Suslin._polynomial_column_peel_ecp_route(route_metadata_cert) == :general_ecp_pipeline

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
    bad_route_evidence = _pw_poly_corrupt_last_column(route_cert.evidence)
    bad_route_cert = _pw_poly_replace_route_certificate(
        route_cert;
        evidence = bad_route_evidence,
        product = route_cert.product,
    )
    @test verify_factorization(bad_route_cert.matrix, bad_route_cert.factors)
    @test !Suslin._verify_polynomial_factorization_route_certificate(bad_route_cert)

    staged_entry = entries["pw-poly-recursive-column-peel-gf2"]
    @test_throws ArgumentError Suslin._polynomial_column_peel_certificate(staged_entry.matrix)

    bad_last_column = _pw_poly_corrupt_last_column(recursive_cert)
    @test !Suslin._verify_polynomial_column_peel_certificate(bad_last_column)

    bad_left_factor = _pw_poly_corrupt_left_factor(recursive_cert)
    @test !Suslin._verify_polynomial_column_peel_certificate(bad_left_factor)

    bad_ecp_certificate = _pw_poly_corrupt_ecp_certificate(recursive_cert)
    @test !Suslin._verify_polynomial_column_peel_certificate(bad_ecp_certificate)

    bad_ecp_route_provenance = _pw_poly_corrupt_ecp_route_provenance(recursive_cert)
    @test !Suslin._verify_polynomial_column_peel_certificate(bad_ecp_route_provenance)

    bad_right_factor = _pw_poly_corrupt_right_factor(recursive_cert)
    @test !Suslin._verify_polynomial_column_peel_certificate(bad_right_factor)

    bad_right_clearing_coefficient = _pw_poly_corrupt_right_clearing_coefficient(recursive_cert)
    @test !Suslin._verify_polynomial_column_peel_certificate(bad_right_clearing_coefficient)

    bad_next_block = _pw_poly_corrupt_next_block(recursive_cert)
    @test !Suslin._verify_polynomial_column_peel_certificate(bad_next_block)

    R = base_ring(recursive_entry.matrix)
    @test_throws ArgumentError Suslin._polynomial_column_peel_certificate(identity_matrix(R, 2))
    @test_throws ArgumentError Suslin._polynomial_column_peel_certificate(last(recursive_cert.peel_steps).peeled_matrix)
    @test_throws ArgumentError Suslin._polynomial_column_peel_certificate(
        recursive_entry.matrix;
        final_route = "fast_local_sl3",
    )
    @test_throws ArgumentError Suslin._polynomial_column_peel_certificate(
        recursive_entry.matrix;
        final_route = :staged_failure,
    )
    @test Suslin._polynomial_column_peel_try_final_route(identity_matrix(R, 2)) === nothing
    @test !Suslin._verify_polynomial_column_peel_certificate((;))

    first_step = first(recursive_cert.peel_steps)
    @test !Suslin._is_valid_polynomial_column_peel_step_data(
        first_step.dimension,
        first_step.input_matrix,
        first_step.last_column,
        Any[(;)],
        first_step.after_left_matrix,
        first_step.right_factors,
        first_step.peeled_matrix,
        first_step.next_block,
    )
    @test !Suslin._is_valid_polynomial_column_peel_step_data(
        first_step.dimension,
        first_step.input_matrix,
        first_step.last_column,
        first_step.left_factors,
        first_step.after_left_matrix,
        _PWPolyBadFactorList(),
        first_step.peeled_matrix,
        first_step.next_block,
    )
    malformed_step = (;
        input_matrix = recursive_entry.matrix,
        dimension = nrows(recursive_entry.matrix),
        next_block = recursive_cert.final_block,
    )
    malformed_cert = (;
        original_matrix = recursive_entry.matrix,
        peel_steps = Any[malformed_step],
        final_block = recursive_cert.final_block,
        final_certificate = recursive_cert.final_certificate,
        final_factors = recursive_cert.final_factors,
        factors = Any[(;)],
        product = recursive_cert.product,
    )
    malformed_verification = Suslin._polynomial_column_peel_core_verification(malformed_cert)
    @test !malformed_verification.steps_ok
    @test !malformed_verification.factor_sequence_ok
    @test !malformed_verification.product_ok
    @test !malformed_verification.factors_ok
    @test !Suslin._polynomial_column_peel_preconditions_ok((; original_matrix = 1, peel_steps = []))
    @test !Suslin._polynomial_column_peel_final_certificate_ok((;))
end
