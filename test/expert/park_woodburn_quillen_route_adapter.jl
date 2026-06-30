using Test
using Suslin
using Oscar

const PW_QUILLEN_ROUTE_CATALOG_PATH =
    joinpath(@__DIR__, "..", "fixtures", "park_woodburn_polynomial_cases.jl")
const SL3_MURTHY_QUILLEN_ADAPTER_FIXTURE_PATH =
    joinpath(@__DIR__, "..", "fixtures", "sl3_murthy_gupta_cases.jl")

if !isdefined(Main, :ParkWoodburnPolynomialFixtureCatalog)
    include(PW_QUILLEN_ROUTE_CATALOG_PATH)
end
if !isdefined(Main, :SL3MurthyGuptaFixtureCatalog)
    include(SL3_MURTHY_QUILLEN_ADAPTER_FIXTURE_PATH)
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

function _pwq_rebuild_sl3_certificate(
        cert;
        target = cert.target,
        branch = cert.branch,
        factors = cert.factors,
        selected_variable = cert.selected_variable,
        witness = cert.witness)
    return Suslin.SL3LocalRealizationCertificate(
        target,
        branch,
        factors,
        selected_variable,
        witness,
    )
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

@testset "Murthy local SL3 Quillen handoff adapter" begin
    catalog = SL3MurthyGuptaFixtureCatalog.catalog()
    by_id = Dict(entry.id => entry for entry in catalog.cases)

    ordinary_fixture = by_id["mg-q0-unit-recursion"]
    ordinary_cert = Suslin.realize_sl3_local_certificate(
        ordinary_fixture.entries.p,
        ordinary_fixture.entries.q,
        ordinary_fixture.entries.r,
        ordinary_fixture.entries.s,
        ordinary_fixture.variable,
    )
    @test ordinary_cert isa Suslin.SL3LocalRealizationCertificate
    @test ordinary_cert.branch == :murthy_q0_unit
    @test Suslin.verify_sl3_local_realization(ordinary_cert)

    ordinary_adapter = Suslin._murthy_quillen_local_adapter(
        ordinary_cert,
        ordinary_fixture.target,
        ordinary_fixture.variable;
        witness_metadata = (;
            fixture_id = ordinary_fixture.id,
            consumer_issue = 211,
        ),
    )
    @test ordinary_adapter isa Suslin.MurthyQuillenLocalAdapter
    @test ordinary_adapter.mode == :ordinary_quillen_factor_sequence
    @test ordinary_adapter.selected_variable == ordinary_fixture.variable
    @test ordinary_adapter.local_product == ordinary_fixture.target
    @test ordinary_adapter.local_correction == ordinary_fixture.target
    @test ordinary_adapter.local_factor_replay.mode == :ordinary
    @test ordinary_adapter.local_factor_replay.materialized_factors == ordinary_cert.factors
    @test ordinary_adapter.witness_metadata.fixture_id == ordinary_fixture.id
    @test ordinary_adapter.replay_metadata.murthy_branch == ordinary_cert.branch
    @test ordinary_adapter.replay_metadata.denominator_product == one(base_ring(ordinary_fixture.target))
    @test Suslin._verify_murthy_quillen_local_adapter(ordinary_adapter)

    ordinary_sequence = Suslin._murthy_quillen_local_factor_sequence_certificate(ordinary_adapter)
    @test ordinary_sequence isa Suslin.QuillenLocalFactorSequenceCertificate
    @test Suslin.verify_quillen_local_factor_sequence_certificate(ordinary_sequence)
    @test ordinary_sequence.selected_variable == ordinary_fixture.variable
    @test ordinary_sequence.local_product == ordinary_fixture.target
    @test ordinary_sequence.local_correction == ordinary_fixture.target
    @test ordinary_sequence.witness_metadata.fixture_id == ordinary_fixture.id
    @test ordinary_sequence.verification.local_product == ordinary_fixture.target

    local_fixture = by_id["mg-local-q0-nonunit-bezout-at-u"]
    local_context = Suslin.sl3_local_murthy_input_context(
        local_fixture.target,
        local_fixture.variable;
        witness = first(local_fixture.witnesses),
    )
    localized_cert = Suslin.realize_sl3_local_certificate(local_context)
    @test localized_cert.branch == :murthy_q0_nonunit_bezout_resultant
    @test Suslin.verify_sl3_local_realization(localized_cert)

    localized_adapter = Suslin._murthy_quillen_local_adapter(
        localized_cert,
        local_fixture.target,
        local_fixture.variable;
        witness_metadata = (;
            fixture_id = local_fixture.id,
            consumer_issue = 211,
        ),
    )
    @test localized_adapter.mode == :localized_replay_handoff
    @test localized_adapter.local_factor_replay.mode == :denominator_cleared
    @test localized_adapter.quillen_factor_sequence === nothing
    @test localized_adapter.quillen_local_certificate === nothing
    @test localized_adapter.local_product === nothing
    @test localized_adapter.local_correction == local_fixture.target
    @test localized_adapter.replay_metadata.denominator_product != one(base_ring(local_fixture.target))
    @test localized_adapter.replay_metadata.cleared_product ==
          localized_adapter.local_factor_replay.cleared_product
    @test Suslin._verify_murthy_quillen_local_adapter(localized_adapter)
    @test_throws ArgumentError Suslin._murthy_quillen_local_factor_sequence_certificate(localized_adapter)
    @test_throws ArgumentError Suslin._murthy_quillen_local_realization_certificate(localized_adapter)

    R = base_ring(ordinary_fixture.target)
    tampered_factors = copy(ordinary_cert.factors)
    tampered_factors[1] =
        tampered_factors[1] * elementary_matrix(3, 1, 3, one(R), R)
    tampered_cert = _pwq_rebuild_sl3_certificate(
        ordinary_cert;
        factors = tampered_factors,
    )
    @test !Suslin.verify_sl3_local_realization(tampered_cert)
    @test_throws ArgumentError Suslin._murthy_quillen_local_adapter(
        tampered_cert,
        ordinary_fixture.target,
        ordinary_fixture.variable,
    )

    local_generators = collect(gens(base_ring(local_fixture.target)))
    mismatched_variable = first(filter(gen -> gen != local_fixture.variable, local_generators))
    @test_throws ArgumentError Suslin._murthy_quillen_local_adapter(
        localized_cert,
        local_fixture.target,
        mismatched_variable,
    )

    @test_throws MethodError Suslin._murthy_quillen_local_adapter(
        ordinary_cert.factors,
        ordinary_fixture.target,
        ordinary_fixture.variable,
    )
end
