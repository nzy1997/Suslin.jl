using Test
using Suslin
using Oscar

if !(@isdefined case008_d14_laurent_link_witness_search_report)
    include(joinpath(@__DIR__, "case008_d14_laurent_link_witness_search.jl"))
end

const LAURENT_LINK_WITNESS_CERTIFICATE_STATUS = :link_witness_certificate
const LAURENT_LINK_WITNESS_CERTIFICATE_NEXT_BOUNDARY =
    :laurent_endpoint_reduction

const LAURENT_LINK_WITNESS_CERTIFICATE_FIELDS = (
    :case_id,
    :dimension,
    :ring_generators,
    :context_status,
    :witness,
    :source_endpoint,
    :target_endpoint,
    :replay_status,
    :identity_status,
    :next_boundary,
    :status,
)

function _laurent_link_witness_certificate_has_fields(value, fields)::Bool
    return all(field -> hasproperty(value, field), fields)
end

function _laurent_link_witness_certificate_from_replay(
    context,
    witness,
    column,
    R;
    require_strict::Bool = true,
)
    _laurent_link_witness_certificate_has_fields(
        context,
        (:case_id, :dimension, :ring_generators, :status),
    ) || throw(ArgumentError("context must include case_id, dimension, ring_generators, and status"))
    context.status == :link_witness_context ||
        throw(ArgumentError("context must have status :link_witness_context"))
    context.dimension == length(column) ||
        throw(ArgumentError("context dimension does not match source column"))
    context.ring_generators == _ring_generator_names(R) ||
        throw(ArgumentError("context ring_generators do not match the ring"))

    candidate = _laurent_link_witness_candidate_from_replay(
        column,
        R,
        witness;
        case_id = context.case_id,
        require_strict,
    )
    verify_laurent_link_witness_candidate(column, R, candidate) ||
        throw(ArgumentError("witness replay did not verify as a link witness"))

    return (;
        case_id = context.case_id,
        dimension = context.dimension,
        ring_generators = context.ring_generators,
        context_status = context.status,
        witness = candidate.witness,
        source_endpoint = candidate.source_endpoint,
        target_endpoint = candidate.target_endpoint,
        replay_status = :ok,
        identity_status = :verified,
        next_boundary = LAURENT_LINK_WITNESS_CERTIFICATE_NEXT_BOUNDARY,
        status = LAURENT_LINK_WITNESS_CERTIFICATE_STATUS,
    )
end

function laurent_link_witness_certificate(context, witness, column, R)
    return _laurent_link_witness_certificate_from_replay(
        context,
        witness,
        column,
        R;
        require_strict = true,
    )
end

function validate_laurent_link_witness_certificate(cert, column, R)::Symbol
    try
        _laurent_link_witness_certificate_has_fields(
            cert,
            LAURENT_LINK_WITNESS_CERTIFICATE_FIELDS,
        ) || return :missing_certificate_fields
        cert.status == LAURENT_LINK_WITNESS_CERTIFICATE_STATUS ||
            return :wrong_status
        cert.replay_status == :ok || return :wrong_replay_status
        cert.identity_status == :verified || return :wrong_identity_status
        cert.next_boundary == LAURENT_LINK_WITNESS_CERTIFICATE_NEXT_BOUNDARY ||
            return :wrong_next_boundary
        cert.dimension == length(column) || return :wrong_dimension
        cert.ring_generators == _ring_generator_names(R) ||
            return :wrong_ring_generators
        cert.context_status == :link_witness_context ||
            return :wrong_context_status

        context = (;
            case_id = cert.case_id,
            dimension = cert.dimension,
            ring_generators = cert.ring_generators,
            status = cert.context_status,
        )
        candidate = _laurent_link_witness_candidate_from_replay(
            column,
            R,
            cert.witness;
            case_id = cert.case_id,
            require_strict = true,
        )
        candidate.source_endpoint == cert.source_endpoint ||
            return :stale_source_endpoint
        candidate.target_endpoint == cert.target_endpoint ||
            return :stale_target_endpoint
        verify_laurent_link_witness_candidate(column, R, candidate) ||
            return :identity_replay_failed

        expected = _laurent_link_witness_certificate_from_replay(
            context,
            cert.witness,
            column,
            R;
            require_strict = true,
        )
        cert == expected || return :stale_certificate
        return :ok
    catch err
        err isa InterruptException && rethrow()
        return :identity_replay_failed
    end
end

function case008_d14_laurent_link_witness_certificate_summary(
    fixture = ToricBuilderCase008D14ColumnBoundary.boundary_fixture(),
)
    source = _case008_d14_link_witness_source_data(fixture)
    search_report = case008_d14_laurent_link_witness_search_report(fixture)
    search_validation =
        validate_case008_d14_laurent_link_witness_search_report(
            search_report,
            fixture,
        )
    search_validation == :ok ||
        throw(ArgumentError("d14 link-witness search report must validate; got $(search_validation)"))

    common = (;
        case_id = search_report.case_id,
        dimension = search_report.dimension,
        search_status = search_report.status,
        search_next_boundary = search_report.next_boundary,
        candidate_count = search_report.candidate_count,
        replay_verified_count = search_report.replay_verified_count,
        source_column = source.replay.after_column,
        ring = fixture.ring,
    )

    if search_report.status == :candidate_found
        candidate = first(search_report.candidates)
        cert = laurent_link_witness_certificate(
            source.context,
            candidate.witness,
            source.replay.after_column,
            fixture.ring,
        )
        validate_laurent_link_witness_certificate(
            cert,
            source.replay.after_column,
            fixture.ring,
        ) == :ok ||
            throw(ArgumentError("first d14 candidate did not validate through the certificate shell"))
        return merge(
            common,
            (;
                d14_status = :link_witness_certificate,
                certificate = cert,
            ),
        )
    end

    return merge(common, (; d14_status = :no_link_witness_candidate))
end

@testset "Laurent link-witness certificate shell" begin
    runtests = read(joinpath(@__DIR__, "..", "runtests.jl"), String)
    @test occursin(
        "\"expert/laurent_link_witness_certificate.jl\"",
        runtests,
    )

    R, (u, v) = Suslin.suslin_laurent_polynomial_ring(GF(2), ["u", "v"])
    column = [one(R), u + v]
    context = (;
        case_id = "synthetic",
        dimension = 2,
        ring_generators = ("u", "v"),
        status = :link_witness_context,
    )
    witness = (;
        family = :two_entry_laurent_combination,
        pivot_index = 2,
        partner_index = 1,
        coefficient = 1,
        exponent = (1, 0),
        ring_generators = ("u", "v"),
    )

    cert = laurent_link_witness_certificate(context, witness, column, R)
    @test cert.case_id == "synthetic"
    @test cert.dimension == 2
    @test cert.ring_generators == ("u", "v")
    @test cert.context_status == :link_witness_context
    @test cert.witness == witness
    @test cert.source_endpoint.leading_exponent == (1, 0)
    @test cert.target_endpoint.leading_exponent == (0, 1)
    @test cert.replay_status == :ok
    @test cert.identity_status == :verified
    @test cert.next_boundary == :laurent_endpoint_reduction
    @test cert.status == :link_witness_certificate
    @test validate_laurent_link_witness_certificate(cert, column, R) == :ok

    tampered_coefficient = merge(
        cert,
        (; witness = merge(witness, (; coefficient = 0))),
    )
    @test validate_laurent_link_witness_certificate(
        tampered_coefficient,
        column,
        R,
    ) == :identity_replay_failed

    wrong_generators = merge(cert, (; ring_generators = ("v", "u")))
    @test validate_laurent_link_witness_certificate(
        wrong_generators,
        column,
        R,
    ) == :wrong_ring_generators

    stale_source = merge(
        cert,
        (;
            source_endpoint = merge(
                cert.source_endpoint,
                (; term_count = cert.source_endpoint.term_count + 1),
            ),
        ),
    )
    @test validate_laurent_link_witness_certificate(
        stale_source,
        column,
        R,
    ) == :stale_source_endpoint

    stale_target = merge(
        cert,
        (;
            target_endpoint = merge(
                cert.target_endpoint,
                (; term_count = cert.target_endpoint.term_count + 1),
            ),
        ),
    )
    @test validate_laurent_link_witness_certificate(
        stale_target,
        column,
        R,
    ) == :stale_target_endpoint

    malformed_witness = merge(
        cert,
        (; witness = _laurent_link_without_field(witness, :family)),
    )
    @test validate_laurent_link_witness_certificate(
        malformed_witness,
        column,
        R,
    ) == :identity_replay_failed

    nonidentity_cert = merge(
        cert,
        (; witness = merge(witness, (; exponent = (0, 0)))),
    )
    @test validate_laurent_link_witness_certificate(
        nonidentity_cert,
        column,
        R,
    ) == :identity_replay_failed
end

@testset "case_008 d=14 Laurent link-witness certificate summary" begin
    summary = case008_d14_laurent_link_witness_certificate_summary()
    @test summary.case_id == "case_008"
    @test summary.dimension == 14
    @test summary.search_status in (:candidate_found, :exhausted)
    @test summary.d14_status in (:link_witness_certificate, :no_link_witness_candidate)

    if summary.search_status == :candidate_found
        @test summary.d14_status == :link_witness_certificate
        @test summary.search_next_boundary == :laurent_link_witness_certificate
        @test summary.candidate_count > 0
        @test summary.replay_verified_count == summary.candidate_count
        @test summary.certificate.case_id == "case_008"
        @test summary.certificate.next_boundary == :laurent_endpoint_reduction
        @test summary.certificate.status == :link_witness_certificate
        @test validate_laurent_link_witness_certificate(
            summary.certificate,
            summary.source_column,
            summary.ring,
        ) == :ok
    else
        @test summary.d14_status == :no_link_witness_candidate
        @test summary.search_next_boundary == :laurent_link_witness_search_expansion
        @test summary.candidate_count == 0
        @test summary.replay_verified_count == 0
        @test !hasproperty(summary, :certificate)
    end
end
