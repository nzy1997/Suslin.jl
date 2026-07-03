using Test
using Oscar
using Suslin

const PARK_WOODBURN_SLN_RECURSIVE_DRIVER_CATALOG_PATH =
    joinpath(@__DIR__, "..", "fixtures", "park_woodburn_sln_driver_cases.jl")
const PARK_WOODBURN_POLYNOMIAL_CATALOG_PATH =
    joinpath(@__DIR__, "..", "fixtures", "park_woodburn_polynomial_cases.jl")

function _sln_recursive_product(factors, R, n::Int)
    product = identity_matrix(R, n)
    for factor in factors
        product *= factor
    end
    return product
end

function _sln_recursive_replace_certificate(
        cert;
        original_matrix = cert.original_matrix,
        peel_steps = cert.peel_steps,
        final_block = cert.final_block,
        final_certificate = cert.final_certificate,
        final_factors = cert.final_factors,
        factors = cert.factors,
        product = cert.product,
        verification = cert.verification,
        descent_metadata = cert.descent_metadata,
        mainline_support_metadata = cert.mainline_support_metadata,
        final_route_provenance = cert.final_route_provenance)
    return Suslin.PolynomialColumnPeelCertificate(
        original_matrix,
        peel_steps,
        final_block,
        final_certificate,
        final_factors,
        factors,
        product,
        verification,
        descent_metadata,
        mainline_support_metadata,
        final_route_provenance,
    )
end

function _sln_recursive_replace_route_certificate(
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

function _sln_recursive_rebuild(record; kwargs...)
    overrides = Dict{Symbol,Any}(pair.first => pair.second for pair in kwargs)
    values = map(fieldnames(typeof(record))) do name
        get(overrides, name, getproperty(record, name))
    end
    if record isa NamedTuple
        return typeof(record)(Tuple(values))
    end
    return typeof(record)(values...)
end

function _sln_recursive_assert_mainline_certificate(cert, entry)
    @test Suslin._verify_polynomial_column_peel_certificate(cert)
    @test cert.original_matrix == entry.matrix
    @test length(cert.peel_steps) == entry.expected_peel_count
    @test tuple((step.dimension for step in cert.peel_steps)..., nrows(cert.final_block)) ==
          entry.descent_dimensions
    @test cert.final_block == entry.final_route.matrix
    @test nrows(cert.final_block) == 3
    @test ncols(cert.final_block) == 3
    @test cert.final_certificate.route == :quillen_patch
    @test cert.final_certificate.evidence isa Suslin.PolynomialSL3SuppliedQuillenRouteEvidence
    @test Suslin._verify_polynomial_sl3_supplied_quillen_route_evidence(
        cert.final_certificate.evidence,
    )
    supplied_evidence = cert.final_certificate.evidence.supplied_evidence
    tampered_local_certificates = copy(supplied_evidence.local_certificates)
    tampered_local_certificates[1] = _sln_recursive_rebuild(
        tampered_local_certificates[1];
        witness_metadata = merge(
            tampered_local_certificates[1].witness_metadata,
            (; tampered = true),
        ),
    )
    tampered_supplied_evidence = _sln_recursive_rebuild(
        supplied_evidence;
        local_certificates = tampered_local_certificates,
    )
    tampered_evidence = _sln_recursive_rebuild(
        cert.final_certificate.evidence;
        supplied_evidence = tampered_supplied_evidence,
    )
    @test !Suslin._verify_polynomial_sl3_supplied_quillen_route_evidence(tampered_evidence)
    quillen_route_adapter = cert.final_certificate.evidence.quillen_route_adapter
    tampered_quillen_patch = _sln_recursive_rebuild(
        quillen_route_adapter.quillen_patch;
        replay_metadata = merge(
            quillen_route_adapter.quillen_patch.replay_metadata,
            (; tampered = true),
        ),
    )
    tampered_adapter = _sln_recursive_rebuild(
        quillen_route_adapter;
        quillen_patch = tampered_quillen_patch,
    )
    tampered_patch_evidence = _sln_recursive_rebuild(
        cert.final_certificate.evidence;
        quillen_route_adapter = tampered_adapter,
    )
    @test !Suslin._verify_polynomial_sl3_supplied_quillen_route_evidence(
        tampered_patch_evidence,
    )
    @test cert.final_route_provenance == :issue184_evidence_backed_sl3
    @test cert.descent_metadata.strict_dimension_descent
    @test cert.descent_metadata.final_block_is_sl3
    @test cert.mainline_support_metadata.supported
    @test cert.mainline_support_metadata.marker == :issue186_mainline
    @test cert.mainline_support_metadata.final_route_issue184_ok
    @test cert.mainline_support_metadata.peel_steps_ecp_verified
    @test cert.mainline_support_metadata.factor_replay_ok
    @test cert.mainline_support_metadata.reconstruction_ok
    @test cert.product == entry.matrix
    @test _sln_recursive_product(cert.factors, base_ring(entry.matrix), nrows(entry.matrix)) == entry.matrix
    @test verify_factorization(entry.matrix, cert.factors)
    for step in cert.peel_steps
        @test Suslin._polynomial_column_peel_step_verification(step).overall_ok
        @test Suslin.verify_ecp_column_reduction(step.ecp_evidence)
    end
end

@testset "Park-Woodburn recursive SLn column-peel certificate" begin
    if !isdefined(Main, :ParkWoodburnSLnDriverFixtureCatalog)
        include(PARK_WOODBURN_SLN_RECURSIVE_DRIVER_CATALOG_PATH)
    end
    if !isdefined(Main, :ParkWoodburnPolynomialFixtureCatalog)
        include(PARK_WOODBURN_POLYNOMIAL_CATALOG_PATH)
    end
    entries = ParkWoodburnSLnDriverFixtureCatalog.cases_by_id()
    polynomial_entries = ParkWoodburnPolynomialFixtureCatalog.cases_by_id()

    sl4 = entries["sln-driver-sl4-gf2-ecp-mainline"]
    sl4_cert = Suslin._polynomial_column_peel_certificate(sl4.matrix)
    _sln_recursive_assert_mainline_certificate(sl4_cert, sl4)

    sl5 = entries["sln-driver-sl5-gf2-two-step"]
    sl5_cert = Suslin._polynomial_column_peel_certificate(sl5.matrix)
    _sln_recursive_assert_mainline_certificate(sl5_cert, sl5)
    @test length(sl5_cert.peel_steps) >= 2

    legacy = entries["sln-driver-legacy-recursive-column-peel-qq"]
    legacy_cert = Suslin._polynomial_column_peel_certificate(legacy.matrix)
    @test Suslin._verify_polynomial_column_peel_certificate(legacy_cert)
    @test !legacy_cert.mainline_support_metadata.supported
    @test legacy_cert.mainline_support_metadata.marker == :not_issue186_mainline
    @test :missing_issue184_final_sl3_route in legacy_cert.mainline_support_metadata.reason_codes

    block_recursive = polynomial_entries["pw-poly-recursive-column-peel-sln-block-qq"]
    block_recursive_cert = Suslin._polynomial_column_peel_certificate(block_recursive.matrix)
    @test Suslin._verify_polynomial_column_peel_certificate(block_recursive_cert)
    @test block_recursive_cert.final_certificate.route == :disjoint_local_blocks
    @test block_recursive_cert.final_route_provenance == :disjoint_local_blocks
    @test !block_recursive_cert.mainline_support_metadata.supported
    @test block_recursive_cert.mainline_support_metadata.marker == :not_issue186_mainline
    @test :missing_issue184_final_sl3_route in block_recursive_cert.mainline_support_metadata.reason_codes

    reordered_steps = reverse(copy(sl5_cert.peel_steps))
    reordered = _sln_recursive_replace_certificate(sl5_cert; peel_steps = reordered_steps)
    @test !Suslin._verify_polynomial_column_peel_certificate(reordered)
    @test !Suslin._polynomial_column_peel_core_verification(reordered).descent_metadata_ok

    duplicated_steps = vcat(copy(sl5_cert.peel_steps), [last(sl5_cert.peel_steps)])
    duplicated = _sln_recursive_replace_certificate(sl5_cert; peel_steps = duplicated_steps)
    @test !Suslin._verify_polynomial_column_peel_certificate(duplicated)

    bad_final_evidence = _sln_recursive_rebuild(
        sl5_cert.final_certificate.evidence;
        replay_metadata = merge(
            sl5_cert.final_certificate.evidence.replay_metadata,
            (; source = :tampered_recursive_final_route),
        ),
    )
    bad_final_certificate = _sln_recursive_replace_route_certificate(
        sl5_cert.final_certificate;
        evidence = bad_final_evidence,
    )
    bad_final = _sln_recursive_replace_certificate(sl5_cert; final_certificate = bad_final_certificate)
    @test !Suslin._verify_polynomial_column_peel_certificate(bad_final)

    bad_mainline_metadata = merge(
        sl5_cert.mainline_support_metadata,
        (; marker = :not_issue186_mainline),
    )
    bad_metadata = _sln_recursive_replace_certificate(
        sl5_cert;
        mainline_support_metadata = bad_mainline_metadata,
    )
    @test verify_factorization(bad_metadata.original_matrix, bad_metadata.factors)
    @test !Suslin._verify_polynomial_column_peel_certificate(bad_metadata)

    bad_factors = copy(sl5_cert.factors)
    bad_factors[1] =
        bad_factors[1] *
        elementary_matrix(nrows(sl5.matrix), 1, 2, one(base_ring(sl5.matrix)), base_ring(sl5.matrix))
    bad_factor_cert = _sln_recursive_replace_certificate(sl5_cert; factors = bad_factors)
    @test !Suslin._verify_polynomial_column_peel_certificate(bad_factor_cert)

    bad_product = _sln_recursive_replace_certificate(
        sl5_cert;
        product = identity_matrix(base_ring(sl5.matrix), nrows(sl5.matrix)),
    )
    @test !Suslin._verify_polynomial_column_peel_certificate(bad_product)
end
