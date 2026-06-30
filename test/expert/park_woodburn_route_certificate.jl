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

function _pw_rebuild(record; kwargs...)
    overrides = Dict{Symbol,Any}()
    for pair in kwargs
        overrides[pair.first] = pair.second
    end
    values = map(fieldnames(typeof(record))) do name
        get(overrides, name, getproperty(record, name))
    end
    return typeof(record)(values...)
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

function _pw_replace_quillen_adapter(
        adapter;
        target = adapter.target,
        route = adapter.route,
        quillen_patch = adapter.quillen_patch,
        global_elementary_factors = adapter.global_elementary_factors,
        product = adapter.product,
        target_matrix = adapter.target_matrix,
        replay_metadata = adapter.replay_metadata,
        verification = adapter.verification)
    return Suslin.PolynomialQuillenPatchRouteAdapter(
        target,
        route,
        quillen_patch,
        global_elementary_factors,
        product,
        target_matrix,
        replay_metadata,
        verification,
    )
end

function _pw_corrupt_route_peel_evidence(cert)
    evidence = cert.evidence
    first_step = first(evidence.peel_steps)
    bad_last_column = copy(first_step.last_column)
    bad_last_column[1] += one(base_ring(evidence.original_matrix))
    bad_step = Suslin.PolynomialColumnPeelStep(
        first_step.dimension,
        first_step.input_matrix,
        bad_last_column,
        first_step.left_factors,
        first_step.after_left_matrix,
        first_step.right_factors,
        first_step.peeled_matrix,
        first_step.next_block,
    )
    bad_evidence = Suslin.PolynomialColumnPeelCertificate(
        evidence.original_matrix,
        Suslin.PolynomialColumnPeelStep[bad_step; evidence.peel_steps[2:end]],
        evidence.final_block,
        evidence.final_certificate,
        evidence.final_factors,
        evidence.factors,
        evidence.product,
        evidence.verification,
    )
    return _pw_replace_certificate(cert; evidence = bad_evidence)
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

    recursive_supported_entry = entries["pw-poly-recursive-column-peel-sl3-qq"]
    auto_peel_cert = Suslin._polynomial_factorization_route_certificate(
        recursive_supported_entry.matrix,
    )
    @test auto_peel_cert.route == :polynomial_column_peel
    @test auto_peel_cert.status == :supported
    @test auto_peel_cert.evidence isa Suslin.PolynomialColumnPeelCertificate
    @test Suslin._verify_polynomial_factorization_route_certificate(auto_peel_cert)
    @test verify_factorization(auto_peel_cert.matrix, auto_peel_cert.factors)

    alias_peel_cert = Suslin._polynomial_factorization_route_certificate(
        recursive_supported_entry.matrix;
        route = :recursive_column_peel,
    )
    @test alias_peel_cert.route == :recursive_column_peel
    @test Suslin._verify_polynomial_factorization_route_certificate(alias_peel_cert)
    @test Suslin._polynomial_staged_failure_evidence(recursive_supported_entry.matrix).error_type == :none
    @test_throws ErrorException Suslin._polynomial_factorization_route_certificate(
        recursive_supported_entry.matrix;
        route = :staged_failure,
    )

    bad_peel_route_cert = _pw_corrupt_route_peel_evidence(auto_peel_cert)
    @test verify_factorization(bad_peel_route_cert.matrix, bad_peel_route_cert.factors)
    @test !Suslin._verify_polynomial_factorization_route_certificate(bad_peel_route_cert)

    quillen_entry = entries["quillen-patched-substitution-witness-qq"]
    quillen_cert = Suslin._polynomial_factorization_route_certificate(quillen_entry.matrix)
    @test quillen_cert.route == :quillen_patch
    @test quillen_cert.status == :supported
    @test quillen_cert.evidence isa Suslin.PolynomialQuillenPatchRouteAdapter
    @test quillen_cert.evidence.quillen_patch isa Suslin.QuillenSuppliedEvidencePatchAssembly
    @test Suslin._verify_polynomial_quillen_patch_route_adapter(quillen_cert.evidence)
    @test verify_factorization(quillen_entry.matrix, quillen_cert.factors)
    @test Suslin._verify_polynomial_factorization_route_certificate(quillen_cert)
    quillen_staged_evidence = Suslin._polynomial_staged_failure_evidence(quillen_entry.matrix)
    @test quillen_staged_evidence.error_type == :none
    @test isempty(quillen_staged_evidence.message)

    if quillen_cert.evidence isa Suslin.PolynomialQuillenPatchRouteAdapter
        quillen_patch = quillen_cert.evidence.quillen_patch
        @test quillen_patch.base_term_policy == :supplied

        bad_quillen_factors = copy(quillen_cert.evidence.global_elementary_factors)
        bad_quillen_factors[1] =
            bad_quillen_factors[1] *
            elementary_matrix(
                nrows(quillen_entry.matrix),
                1,
                3,
                one(base_ring(quillen_entry.matrix)),
                base_ring(quillen_entry.matrix),
            )
        bad_quillen_evidence = _pw_replace_quillen_adapter(
            quillen_cert.evidence;
            global_elementary_factors = bad_quillen_factors,
        )
        bad_quillen_cert = _pw_replace_certificate(quillen_cert; evidence = bad_quillen_evidence)
        @test verify_factorization(bad_quillen_cert.matrix, bad_quillen_cert.factors)
        @test !Suslin._verify_polynomial_factorization_route_certificate(bad_quillen_cert)

        tampered_local_certificates = copy(quillen_patch.local_certificates)
        tampered_first_sequence = tampered_local_certificates[1]
        tampered_first_factors = copy(tampered_first_sequence.factors)
        tampered_first_factors[1] = _pw_rebuild(
            tampered_first_factors[1];
            numerator = tampered_first_factors[1].numerator + one(quillen_patch.ring),
        )
        tampered_local_certificates[1] = _pw_rebuild(
            tampered_first_sequence;
            factors = tampered_first_factors,
        )
        tampered_local_patch = _pw_rebuild(
            quillen_patch;
            local_certificates = tampered_local_certificates,
        )
        @test !Suslin.verify_quillen_patch(tampered_local_patch)
        @test_throws ArgumentError Suslin._polynomial_quillen_patch_route_adapter(
            quillen_entry.matrix,
            tampered_local_patch,
        )

        tampered_chain = _pw_rebuild(
            quillen_patch.substitution_chain;
            sign_convention = :park_woodburn_plus,
        )
        tampered_chain_patch = _pw_rebuild(
            quillen_patch;
            substitution_chain = tampered_chain,
        )
        @test !Suslin.verify_quillen_patch(tampered_chain_patch)
        @test_throws ArgumentError Suslin._polynomial_quillen_patch_route_adapter(
            quillen_entry.matrix,
            tampered_chain_patch,
        )
    end

    S = base_ring(quillen_entry.matrix)
    X, r, g = collect(gens(S))
    nonfixture_quillen = elementary_matrix(
        3,
        1,
        3,
        X * r + g + one(S),
        S,
    )
    nonfixture_quillen_cert =
        Suslin._polynomial_factorization_route_certificate(nonfixture_quillen)
    @test nonfixture_quillen_cert.route == :quillen_patch
    @test nonfixture_quillen_cert.evidence.quillen_patch isa
          Suslin.QuillenSuppliedEvidencePatchAssembly
    @test all(
        Suslin.verify_quillen_local_factor_sequence_certificate,
        nonfixture_quillen_cert.evidence.quillen_patch.local_certificates,
    )
    @test nonfixture_quillen_cert.evidence.quillen_patch.base_term_policy == :supplied
    @test nonfixture_quillen_cert.evidence.quillen_patch.base_term ==
          elementary_matrix(3, 1, 3, g + one(S), S)
    @test nonfixture_quillen_cert.evidence.quillen_patch.substitution_chain.verification.telescope_ok
    @test verify_factorization(nonfixture_quillen, nonfixture_quillen_cert.factors)
    nonfixture_quillen_data = Suslin._polynomial_quillen_supplied_evidence_data(nonfixture_quillen)
    @test nonfixture_quillen_data !== nothing
    @test all(
        Suslin.verify_quillen_local_factor_sequence_certificate,
        nonfixture_quillen_data.local_certificates,
    )
    missing_base_term_err = _pw_captured_error(() ->
        Suslin.assemble_quillen_patch_from_local_evidence(
            nonfixture_quillen,
            nonfixture_quillen_data.selected_variable,
            nonfixture_quillen_data.local_certificates;
            exponent = 1,
            coverage_multipliers = [one(S), one(S)],
        )
    )
    @test missing_base_term_err isa ArgumentError
    @test occursin("A(0)", sprint(showerror, missing_base_term_err))

    wrong_base_factor = elementary_matrix(3, 1, 3, g, S)
    wrong_base_term_err = _pw_captured_error(() ->
        Suslin.assemble_quillen_patch_from_local_evidence(
            nonfixture_quillen,
            nonfixture_quillen_data.selected_variable,
            nonfixture_quillen_data.local_certificates;
            exponent = 1,
            coverage_multipliers = [one(S), one(S)],
            base_term_policy = :supplied,
            base_term_factors = [wrong_base_factor],
        )
    )
    @test wrong_base_term_err isa ArgumentError
    @test occursin("base-term evidence", sprint(showerror, wrong_base_term_err))

    tampered_patch_certificate = _pw_rebuild(
        nonfixture_quillen_cert.evidence.quillen_patch;
        replay_metadata = (; source = :tampered_quillen_patch_certificate),
    )
    @test !Suslin.verify_quillen_patch(tampered_patch_certificate)
    @test_throws ArgumentError Suslin._polynomial_quillen_patch_route_adapter(
        nonfixture_quillen,
        tampered_patch_certificate,
    )

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
        "missing Quillen/local realizability witness",
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
    injected_route_override_active = Ref(true)
    @eval Suslin function _polynomial_factorization_route_certificate(
            A::$matrix_type;
            route = nothing,
            quillen_patch = nothing,
            allow_recursive_column_peel::Bool = false)
        if $injected_route_override_active[] &&
                A == $fast_cert.matrix &&
                route === nothing &&
                quillen_patch === nothing &&
                !allow_recursive_column_peel
            return $public_bad_cert
        end

        return invoke(
            Suslin._polynomial_factorization_route_certificate,
            Tuple{Any},
            A;
            route = route,
            quillen_patch = quillen_patch,
            allow_recursive_column_peel = allow_recursive_column_peel,
        )
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
        injected_route_override_active[] = false
        Base.delete_method(injected_method)
    end
end
