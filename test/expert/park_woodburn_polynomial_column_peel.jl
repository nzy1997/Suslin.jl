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
        after_left_matrix = step.after_left_matrix,
        right_factors = step.right_factors,
        peeled_matrix = step.peeled_matrix,
        next_block = step.next_block)
    return Suslin.PolynomialColumnPeelStep(
        dimension,
        input_matrix,
        last_column,
        left_factors,
        left_certificate,
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

function _pw_poly_rebuild(record; kwargs...)
    overrides = Dict{Symbol,Any}()
    for pair in kwargs
        overrides[pair.first] = pair.second
    end
    values = map(fieldnames(typeof(record))) do name
        get(overrides, name, getproperty(record, name))
    end
    return typeof(record)(values...)
end

function _pw_poly_issue184_sl3_route_case()
    R, (X, r, g) = Oscar.polynomial_ring(QQ, ["X", "r", "g"])
    p = X + r * g + one(R)
    q = one(R)
    s = one(R)
    lower = X + r * g
    A = matrix(R, [
        p q zero(R);
        lower s zero(R);
        zero(R) zero(R) one(R)
    ])
    @assert det(A) == one(R)
    return (; R, X, r, g, p, q, s, lower, A)
end

function _pw_poly_wrap_sl4_final_block(final_block, tail_entries)
    R = base_ring(final_block)
    length(tail_entries) == 3 || throw(ArgumentError("SL4 wrapper needs three tail entries"))
    wrapped = block_embedding(final_block, 4, [1, 2, 3])
    for row in 1:3
        wrapped[row, 4] = tail_entries[row]
    end
    return wrapped
end

function _pw_poly_certificate_with_provenance(
        cert,
        final_route_provenance)
    return Suslin.PolynomialColumnPeelCertificate(
        cert.original_matrix,
        cert.peel_steps,
        cert.final_block,
        cert.final_certificate,
        cert.final_factors,
        cert.factors,
        cert.product,
        cert.verification,
        final_route_provenance,
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
    @test Suslin.verify_ecp_column_reduction(step.left_certificate)
    @test step.left_certificate.original_column == step.last_column
    @test step.left_certificate.factors == step.left_factors
    @test step.left_certificate.final_column == _pw_poly_peel_target_column(R, step.dimension)
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

    bad_right_factor = _pw_poly_corrupt_right_factor(recursive_cert)
    @test !Suslin._verify_polynomial_column_peel_certificate(bad_right_factor)

    bad_next_block = _pw_poly_corrupt_next_block(recursive_cert)
    @test !Suslin._verify_polynomial_column_peel_certificate(bad_next_block)

    issue184 = _pw_poly_issue184_sl3_route_case()
    issue184_wrapped = _pw_poly_wrap_sl4_final_block(
        issue184.A,
        [issue184.X + issue184.r, issue184.g + one(issue184.R), issue184.r * issue184.g],
    )
    issue184_cert = Suslin._polynomial_column_peel_certificate(issue184_wrapped)
    explicit_issue184_cert = Suslin._polynomial_column_peel_certificate(
        issue184_wrapped;
        final_route = :quillen_patch,
    )
    for cert in (issue184_cert, explicit_issue184_cert)
        @test cert.final_block == issue184.A
        @test cert.final_certificate.route == :quillen_patch
        @test cert.final_certificate.evidence isa Suslin.PolynomialSL3QuillenMurthyRouteEvidence
        @test cert.final_route_provenance == :issue184_evidence_backed_sl3
        @test cert.verification.final_route_provenance_ok
        _pw_poly_assert_real_peel_certificate(cert, issue184_wrapped)
    end

    adapter_only_final = _pw_poly_replace_route_certificate(
        issue184_cert.final_certificate;
        evidence = issue184_cert.final_certificate.evidence.quillen_route_adapter,
    )
    @test Suslin._verify_polynomial_factorization_route_certificate(adapter_only_final)
    adapter_only_peel = _pw_poly_replace_certificate(
        issue184_cert;
        final_certificate = adapter_only_final,
    )
    @test verify_factorization(adapter_only_peel.original_matrix, adapter_only_peel.factors)
    @test !Suslin._verify_polynomial_column_peel_certificate(adapter_only_peel)

    tampered_provenance =
        _pw_poly_certificate_with_provenance(issue184_cert, :tampered_quillen_patch)
    @test verify_factorization(tampered_provenance.original_matrix, tampered_provenance.factors)
    @test !Suslin._verify_polynomial_column_peel_certificate(tampered_provenance)

    tampered_witness = merge(
        issue184_cert.final_certificate.evidence.context.local_form_witness,
        (; monic_entry_position = (1, 2)),
    )
    tampered_context = _pw_poly_rebuild(
        issue184_cert.final_certificate.evidence.context;
        local_form_witness = tampered_witness,
    )
    tampered_evidence = _pw_poly_rebuild(
        issue184_cert.final_certificate.evidence;
        context = tampered_context,
    )
    tampered_final_certificate = _pw_poly_replace_route_certificate(
        issue184_cert.final_certificate;
        evidence = tampered_evidence,
    )
    tampered_final_evidence_peel = _pw_poly_replace_certificate(
        issue184_cert;
        final_certificate = tampered_final_certificate,
    )
    @test verify_factorization(
        tampered_final_evidence_peel.original_matrix,
        tampered_final_evidence_peel.factors,
    )
    @test !Suslin._verify_polynomial_column_peel_certificate(tampered_final_evidence_peel)

    tampered_final_factors = copy(issue184_cert.final_factors)
    tampered_final_factors[1] =
        tampered_final_factors[1] *
        elementary_matrix(3, 1, 2, one(issue184.R), issue184.R)
    tampered_final_factor_peel = _pw_poly_replace_certificate(
        issue184_cert;
        final_factors = tampered_final_factors,
    )
    @test verify_factorization(tampered_final_factor_peel.original_matrix, tampered_final_factor_peel.factors)
    @test !Suslin._verify_polynomial_column_peel_certificate(tampered_final_factor_peel)

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
