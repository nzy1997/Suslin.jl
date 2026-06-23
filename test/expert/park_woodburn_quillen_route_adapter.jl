using Test
using Suslin
using Oscar

const PW_QUILLEN_ROUTE_CATALOG_PATH =
    joinpath(@__DIR__, "..", "fixtures", "park_woodburn_polynomial_cases.jl")

if !isdefined(Main, :ParkWoodburnPolynomialFixtureCatalog)
    include(PW_QUILLEN_ROUTE_CATALOG_PATH)
end
if !isdefined(Main, :constructive_patch)
    include(joinpath(@__DIR__, "quillen_induction_constructive.jl"))
end

function _pwq_product(factors, R, n::Int)
    product = identity_matrix(R, n)
    for factor in factors
        product *= factor
    end
    return product
end

function _pwq_rebuild(record; kwargs...)
    overrides = Dict{Symbol,Any}()
    for pair in kwargs
        overrides[pair.first] = pair.second
    end
    values = map(fieldnames(typeof(record))) do name
        get(overrides, name, getproperty(record, name))
    end
    return typeof(record)(values...)
end

function _pwq_replace_adapter(
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

function _pwq_replace_route_certificate(
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

function _pwq_adapter_accepts(target, patch)::Bool
    try
        adapter = Suslin._polynomial_quillen_patch_route_adapter(target, patch)
        return Suslin._verify_polynomial_quillen_patch_route_adapter(adapter)
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end

function _pwq_route_accepts(target, patch)::Bool
    try
        cert = Suslin._polynomial_factorization_route_certificate(
            target;
            route = :quillen_patch,
            quillen_patch = patch,
        )
        return Suslin._verify_polynomial_factorization_route_certificate(cert)
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end

@testset "Park-Woodburn Quillen patch route adapter" begin
    route_entries = ParkWoodburnPolynomialFixtureCatalog.cases_by_id()
    route_entry = route_entries["quillen-patched-substitution-witness-qq"]
    quillen_entries = Main.QuillenPatchFixtureCatalog.cases_by_id()
    quillen_entry = quillen_entries[route_entry.provenance.quillen_fixture_id]

    _, _, _, patch = constructive_patch(quillen_entry)
    @test Suslin.verify_quillen_patch(patch)
    @test patch.target == route_entry.matrix

    adapter = Suslin._polynomial_quillen_patch_route_adapter(route_entry.matrix, patch)
    @test adapter isa Suslin.PolynomialQuillenPatchRouteAdapter
    @test adapter.route == :quillen_patch
    @test adapter.target_matrix == route_entry.matrix
    @test adapter.product == route_entry.matrix
    @test adapter.global_elementary_factors == patch.global_elementary_factors
    @test _pwq_product(adapter.global_elementary_factors, base_ring(route_entry.matrix), nrows(route_entry.matrix)) == route_entry.matrix
    @test Suslin._verify_polynomial_quillen_patch_route_adapter(adapter)

    auto_cert = Suslin._polynomial_factorization_route_certificate(route_entry.matrix)
    @test auto_cert.route == :quillen_patch
    @test auto_cert.status == :supported
    @test auto_cert.evidence isa Suslin.PolynomialQuillenPatchRouteAdapter
    @test verify_factorization(route_entry.matrix, auto_cert.factors)
    @test Suslin._verify_polynomial_factorization_route_certificate(auto_cert)

    cert = Suslin._polynomial_factorization_route_certificate(
        route_entry.matrix;
        route = :quillen_patch,
        quillen_patch = patch,
    )
    @test cert.route == :quillen_patch
    @test cert.status == :supported
    @test cert.evidence isa Suslin.PolynomialQuillenPatchRouteAdapter
    @test cert.factors == adapter.global_elementary_factors
    @test cert.product == route_entry.matrix
    @test verify_factorization(route_entry.matrix, cert.factors)
    @test Suslin._verify_polynomial_factorization_route_certificate(cert)

    bad_cover = Suslin.QuillenDenominatorCoverCertificate(
        patch.cover_certificate.ring,
        patch.cover_certificate.denominators,
        [
            patch.cover_certificate.coverage_multipliers[1] + one(patch.ring),
            patch.cover_certificate.coverage_multipliers[2],
        ],
        patch.cover_certificate.coverage_sum + patch.cover_certificate.denominators[1],
        patch.cover_certificate.verification,
    )
    tampered_cover_patch = _pwq_rebuild(patch; cover_certificate = bad_cover)
    @test !Suslin.verify_quillen_patch(tampered_cover_patch)
    @test !_pwq_adapter_accepts(route_entry.matrix, tampered_cover_patch)
    @test !_pwq_route_accepts(route_entry.matrix, tampered_cover_patch)

    tampered_local_certificates = copy(patch.local_certificates)
    tampered_factors = copy(tampered_local_certificates[1].factors)
    tampered_factors[1] =
        tampered_factors[1] *
        elementary_matrix(patch.size, 1, 3, one(patch.ring), patch.ring)
    tampered_local_certificates[1] = _pwq_rebuild(
        tampered_local_certificates[1];
        factors = tampered_factors,
    )
    tampered_local_patch = _pwq_rebuild(
        patch;
        local_certificates = tampered_local_certificates,
    )
    @test !Suslin.verify_quillen_patch(tampered_local_patch)
    @test !_pwq_adapter_accepts(route_entry.matrix, tampered_local_patch)
    @test !_pwq_route_accepts(route_entry.matrix, tampered_local_patch)

    overwritten_patch = _pwq_rebuild(
        tampered_local_patch;
        patched_product = route_entry.matrix,
        target = route_entry.matrix,
    )
    @test overwritten_patch.patched_product == route_entry.matrix
    @test overwritten_patch.target == route_entry.matrix
    @test !Suslin.verify_quillen_patch(overwritten_patch)
    @test !_pwq_adapter_accepts(route_entry.matrix, overwritten_patch)

    bad_adapter_factors = copy(adapter.global_elementary_factors)
    bad_adapter_factors[1] =
        bad_adapter_factors[1] *
        elementary_matrix(patch.size, 1, 3, one(patch.ring), patch.ring)
    bad_adapter = _pwq_replace_adapter(
        adapter;
        global_elementary_factors = bad_adapter_factors,
    )
    @test !Suslin._verify_polynomial_quillen_patch_route_adapter(bad_adapter)

    malformed_adapter = _pwq_replace_adapter(adapter; quillen_patch = nothing)
    @test !Suslin._verify_polynomial_quillen_patch_route_adapter(malformed_adapter)

    bad_route_cert = _pwq_replace_route_certificate(
        cert;
        factors = bad_adapter_factors,
        evidence = bad_adapter,
    )
    @test !Suslin._verify_polynomial_factorization_route_certificate(bad_route_cert)

    @test_throws ArgumentError Suslin._polynomial_factorization_route_certificate(
        route_entry.matrix;
        route = :quillen_patch,
    )
end
