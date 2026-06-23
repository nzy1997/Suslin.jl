using Test
using Suslin
using Oscar

const PARK_WOODBURN_ROUTE_CATALOG_PATH =
    joinpath(@__DIR__, "..", "fixtures", "park_woodburn_polynomial_cases.jl")

struct PWRouteExplodingEq end

Base.:(==)(::PWRouteExplodingEq, _) = throw(ArgumentError("route evidence equality sentinel"))

function _pw_route_product(factors, R, n::Int)
    product = identity_matrix(R, n)
    for factor in factors
        product *= factor
    end
    return product
end

function _pw_captured_error(f)
    try
        f()
        return nothing
    catch err
        return err
    end
end

function _pw_replace_certificate(
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

function _pw_replace_reduction(
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
    return Suslin.SLNToSL3Reduction(
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

function _pw_route_assert_success(cert, A)
    R = base_ring(A)
    @test Suslin._verify_polynomial_factorization_route_certificate(cert)
    @test _pw_route_product(cert.factors, R, nrows(A)) == A
    @test cert.product == A
    @test verify_factorization(A, cert.factors)
    return nothing
end

@testset "Park-Woodburn polynomial route certificates" begin
    if !isdefined(Main, :ParkWoodburnPolynomialFixtureCatalog)
        include(PARK_WOODBURN_ROUTE_CATALOG_PATH)
    end
    entries = ParkWoodburnPolynomialFixtureCatalog.cases_by_id()

    fast_entry = entries["pw-poly-univariate-sl3-fast-local-qq"]
    fast_cert = Suslin._polynomial_factorization_route_certificate(
        fast_entry.matrix;
        route = fast_entry.route,
    )
    @test fast_cert.route == :fast_local_sl3
    @test fast_cert.status == :supported
    @test Suslin.verify_sl3_local_realization(fast_cert.evidence)
    @test fast_cert.evidence.target == fast_entry.matrix
    @test fast_cert.factors == fast_cert.evidence.factors
    _pw_route_assert_success(fast_cert, fast_entry.matrix)

    auto_fast_cert = Suslin._polynomial_factorization_route_certificate(fast_entry.matrix)
    @test auto_fast_cert.route == :fast_local_sl3
    @test Suslin._verify_polynomial_factorization_route_certificate(auto_fast_cert)

    block_entry = entries["pw-poly-univariate-sln-disjoint-blocks-qq"]
    block_cert = Suslin._polynomial_factorization_route_certificate(
        block_entry.matrix;
        route = block_entry.route,
    )
    @test block_cert.route == :disjoint_local_blocks
    @test block_cert.status == :supported
    @test Suslin.verify_sln_to_sl3_reduction(block_cert.evidence)
    @test block_cert.evidence.original_matrix == block_entry.matrix
    @test block_cert.factors == block_cert.evidence.factors
    _pw_route_assert_success(block_cert, block_entry.matrix)

    auto_block_cert = Suslin._polynomial_factorization_route_certificate(block_entry.matrix)
    @test auto_block_cert.route == :disjoint_local_blocks
    @test Suslin._verify_polynomial_factorization_route_certificate(auto_block_cert)
    @test Suslin._polynomial_staged_failure_evidence(block_entry.matrix).error_type == :none

    recursive_entry = entries["pw-poly-recursive-column-peel-gf2"]
    staged_cert = Suslin._polynomial_factorization_route_certificate(
        recursive_entry.matrix;
        route = :staged_failure,
    )
    @test staged_cert.route == :staged_failure
    @test staged_cert.status == :staged
    @test isempty(staged_cert.factors)
    @test hasproperty(staged_cert.evidence, :message)
    @test !isempty(staged_cert.evidence.message)
    @test Suslin._verify_polynomial_factorization_route_certificate(staged_cert)

    auto_staged_cert = Suslin._polynomial_factorization_route_certificate(recursive_entry.matrix)
    @test auto_staged_cert.route == :staged_failure
    @test Suslin._verify_polynomial_factorization_route_certificate(auto_staged_cert)
    @test Suslin._polynomial_staged_failure_evidence(fast_entry.matrix).error_type == :none

    R = base_ring(fast_cert.matrix)
    n = nrows(fast_cert.matrix)

    staged_n_gt_3_err = try
        Suslin._throw_staged_factorization_failure(identity_matrix(R, 4), :polynomial, nothing)
        nothing
    catch err
        err
    end
    @test staged_n_gt_3_err isa ArgumentError
    @test occursin("SL_n reduction layer for n > 3", sprint(showerror, staged_n_gt_3_err))

    @test_throws ArgumentError Suslin._polynomial_factorization_route_certificate(
        fast_entry.matrix;
        route = "fast_local_sl3",
    )
    @test_throws ArgumentError Suslin._polynomial_factorization_route_certificate(
        fast_entry.matrix;
        route = :quillen_patched_substitution,
    )

    bad_route = _pw_replace_certificate(fast_cert; route = :quillen_patched_substitution)
    @test !Suslin._verify_polynomial_factorization_route_certificate(bad_route)

    bad_factors = copy(fast_cert.factors)
    bad_factors[1] = identity_matrix(R, n)
    bad_factor_cert = _pw_replace_certificate(fast_cert; factors = bad_factors)
    @test !Suslin._verify_polynomial_factorization_route_certificate(bad_factor_cert)

    wrong_size_factors = copy(fast_cert.factors)
    wrong_size_factors[1] = identity_matrix(R, 2)
    wrong_size_factor_cert = _pw_replace_certificate(fast_cert; factors = wrong_size_factors)
    @test !Suslin._verify_polynomial_factorization_route_certificate(wrong_size_factor_cert)

    bad_product = identity_matrix(R, n)
    bad_product_cert = _pw_replace_certificate(fast_cert; product = bad_product)
    @test !Suslin._verify_polynomial_factorization_route_certificate(bad_product_cert)

    bad_evidence = Suslin.SL3LocalRealizationCertificate(
        fast_cert.evidence.target,
        fast_cert.evidence.branch,
        fast_cert.evidence.factors,
        fast_cert.evidence.selected_variable,
        merge(fast_cert.evidence.witness, (; q = fast_cert.evidence.witness.q + one(R))),
    )
    bad_evidence_cert = _pw_replace_certificate(fast_cert; evidence = bad_evidence)
    @test Suslin.verify_factorization(bad_evidence_cert.matrix, bad_evidence_cert.factors)
    @test !Suslin._verify_polynomial_factorization_route_certificate(bad_evidence_cert)

    exploding_evidence = Suslin.SL3LocalRealizationCertificate(
        PWRouteExplodingEq(),
        fast_cert.evidence.branch,
        fast_cert.evidence.factors,
        fast_cert.evidence.selected_variable,
        fast_cert.evidence.witness,
    )
    exploding_evidence_cert = _pw_replace_certificate(fast_cert; evidence = exploding_evidence)
    @test !Suslin._verify_polynomial_factorization_route_certificate(exploding_evidence_cert)

    fake_staged_evidence = (;
        error_type = :ArgumentError,
        message = "fake staged failure for unrelated matrix",
    )
    bad_staged_evidence_cert = _pw_replace_certificate(
        staged_cert;
        evidence = fake_staged_evidence,
    )
    @test !Suslin._verify_polynomial_factorization_route_certificate(bad_staged_evidence_cert)

    staged_empty_message_cert = _pw_replace_certificate(
        staged_cert;
        evidence = (; error_type = :ArgumentError, message = ""),
    )
    staged_empty_message_err =
        _pw_captured_error(() -> Suslin._polynomial_verified_certificate_factors(staged_empty_message_cert))
    @test staged_empty_message_err isa ArgumentError
    @test occursin(
        "SL_n reduction layer for n > 3",
        sprint(showerror, staged_empty_message_err),
    )

    bad_block_evidence = _pw_replace_reduction(
        block_cert.evidence;
        obligations = Suslin.SL3LocalObligation[],
    )
    bad_block_evidence_cert = _pw_replace_certificate(block_cert; evidence = bad_block_evidence)
    @test Suslin.verify_factorization(bad_block_evidence_cert.matrix, bad_block_evidence_cert.factors)
    @test !Suslin._verify_polynomial_factorization_route_certificate(bad_block_evidence_cert)

    unsupported_status_cert = _pw_replace_certificate(fast_cert; status = :unknown)
    unsupported_status_err =
        _pw_captured_error(() -> Suslin._polynomial_verified_certificate_factors(unsupported_status_cert))
    @test unsupported_status_err isa ArgumentError
    @test occursin(
        "unsupported polynomial factorization route certificate status unknown",
        sprint(showerror, unsupported_status_err),
    )

    public_bad_factors = copy(fast_cert.factors)
    public_bad_factors[1] = identity_matrix(R, n)
    public_bad_cert = _pw_replace_certificate(fast_cert; factors = public_bad_factors)
    matrix_type = typeof(fast_cert.matrix)
    @eval Suslin function _polynomial_factorization_route_certificate(
            A::$matrix_type;
            route = nothing)
        return $public_bad_cert
    end
    injected_method = which(
        Suslin._polynomial_factorization_route_certificate,
        (matrix_type,),
    )
    try
        public_err = _pw_captured_error(() -> elementary_factorization(fast_cert.matrix))
        @test public_err isa ErrorException
        @test occursin(
            "internal polynomial factorization route certificate verification failed",
            sprint(showerror, public_err),
        )
    finally
        Base.delete_method(injected_method)
    end
end
