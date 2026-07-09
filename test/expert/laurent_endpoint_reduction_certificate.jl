using Test
using Suslin
using Oscar

if !(@isdefined case008_d14_laurent_endpoint_reduction_search_report)
    include(joinpath(@__DIR__, "case008_d14_laurent_endpoint_reduction_search.jl"))
end

const LAURENT_ENDPOINT_REDUCTION_CERTIFICATE_STATUS =
    :endpoint_reduction_certificate
const LAURENT_ENDPOINT_REDUCTION_CERTIFICATE_NEXT_BOUNDARY =
    :laurent_normality_replay
const LAURENT_ENDPOINT_REDUCTION_OPERATION_FIELDS = (
    :family,
    :endpoint_index,
    :operation,
    :ring_generators,
)
const LAURENT_ENDPOINT_REDUCTION_CERTIFICATE_FIELDS = (
    :case_id,
    :dimension,
    :ring_generators,
    :context_status,
    :operation,
    :source_endpoint,
    :target_endpoint,
    :replay_status,
    :identity_status,
    :next_boundary,
    :status,
)

function _laurent_endpoint_reduction_has_fields(value, fields)::Bool
    return all(field -> hasproperty(value, field), fields)
end

function _laurent_endpoint_reduction_without_field(
    value::NamedTuple,
    field::Symbol,
)
    kept = tuple((name for name in keys(value) if name != field)...)
    return NamedTuple{kept}(tuple((getproperty(value, name) for name in kept)...))
end

function _laurent_endpoint_metadata_from_column(
    column,
    R,
    endpoint_index::Int;
    case_id,
)
    measure = Suslin._laurent_descent_measure_from_column(column, R; case_id)
    return Suslin._laurent_link_endpoint_metadata(
        column[endpoint_index],
        R,
        endpoint_index,
        measure;
        case_id,
    )
end

function _laurent_endpoint_reduction_operation_status(
    endpoint_operation,
    n::Int,
    R,
)::Symbol
    _laurent_endpoint_reduction_has_fields(
        endpoint_operation,
        LAURENT_ENDPOINT_REDUCTION_OPERATION_FIELDS,
    ) || return :malformed_endpoint_operation
    endpoint_operation.family in (
        :laurent_endpoint_entry_addition,
        :paired_laurent_endpoint_entry_addition,
    ) || return :malformed_endpoint_operation
    endpoint_operation.ring_generators ==
        Suslin._laurent_descent_ring_generators(R) ||
        return :wrong_ring_generators
    _laurent_endpoint_reduction_has_fields(
        endpoint_operation.operation,
        (
            :family,
            :target_index,
            :source_index,
            :coefficient,
            :exponent,
            :ring_generators,
        ),
    ) || return :malformed_endpoint_operation
    endpoint_operation.operation.target_index == endpoint_operation.endpoint_index ||
        return :wrong_endpoint_index
    endpoint_operation.operation.ring_generators ==
        endpoint_operation.ring_generators || return :wrong_ring_generators
    Suslin._laurent_descent_operation_status(
        endpoint_operation.operation,
        n,
        R,
    ) == :ok || return :malformed_endpoint_operation
    1 <= Int(endpoint_operation.endpoint_index) <= n || return :wrong_endpoint_index
    return :ok
end

function _laurent_endpoint_reduction_replay(
    context,
    endpoint_operation,
    column,
    R;
    require_strict::Bool,
)
    status = _laurent_endpoint_reduction_operation_status(
        endpoint_operation,
        length(column),
        R,
    )
    status == :ok ||
        throw(ArgumentError("invalid Laurent endpoint operation: $(status)"))
    endpoint_index = Int(endpoint_operation.endpoint_index)
    source_endpoint = _laurent_endpoint_metadata_from_column(
        column,
        R,
        endpoint_index;
        case_id = context.case_id,
    )
    replayed_column = Suslin._replay_laurent_elementary_entry_addition(
        column,
        R,
        endpoint_operation.operation,
    )
    target_endpoint = _laurent_endpoint_metadata_from_column(
        replayed_column,
        R,
        endpoint_index;
        case_id = context.case_id,
    )
    relation = Suslin._strictly_decreases_laurent_measure(
        source_endpoint.column_measure,
        target_endpoint.column_measure,
    ) ? :strict_decrease : :not_strict_decrease
    require_strict && relation != :strict_decrease &&
        throw(
            ArgumentError(
                "endpoint operation does not strictly decrease the endpoint measure",
            ),
        )
    return (; endpoint_operation, source_endpoint, target_endpoint, relation)
end

function laurent_endpoint_reduction_certificate(
    context,
    endpoint_operation,
    column,
    R;
    source_endpoint = nothing,
    require_strict::Bool = true,
)
    _laurent_endpoint_reduction_has_fields(
        context,
        (:case_id, :dimension, :ring_generators, :status),
    ) || throw(
        ArgumentError(
            "context must include case_id, dimension, ring_generators, and status",
        ),
    )
    context.status == :endpoint_reduction_context ||
        throw(ArgumentError("context must have status :endpoint_reduction_context"))
    context.dimension == length(column) ||
        throw(ArgumentError("context dimension does not match source column"))
    context.ring_generators == Suslin._laurent_descent_ring_generators(R) ||
        throw(ArgumentError("context ring_generators do not match the ring"))
    replay = _laurent_endpoint_reduction_replay(
        context,
        endpoint_operation,
        column,
        R;
        require_strict,
    )
    source_endpoint === nothing || source_endpoint == replay.source_endpoint ||
        throw(ArgumentError("source endpoint is stale for the input column"))
    return (;
        case_id = context.case_id,
        dimension = context.dimension,
        ring_generators = context.ring_generators,
        context_status = context.status,
        operation = replay.endpoint_operation,
        source_endpoint = replay.source_endpoint,
        target_endpoint = replay.target_endpoint,
        endpoint_measure_relation = replay.relation,
        replay_status = :ok,
        identity_status = :verified,
        next_boundary = LAURENT_ENDPOINT_REDUCTION_CERTIFICATE_NEXT_BOUNDARY,
        status = LAURENT_ENDPOINT_REDUCTION_CERTIFICATE_STATUS,
    )
end

function validate_laurent_endpoint_reduction_certificate(cert, column, R)::Symbol
    try
        _laurent_endpoint_reduction_has_fields(
            cert,
            LAURENT_ENDPOINT_REDUCTION_CERTIFICATE_FIELDS,
        ) || return :missing_certificate_fields
        cert.status == LAURENT_ENDPOINT_REDUCTION_CERTIFICATE_STATUS ||
            return :wrong_status
        cert.replay_status == :ok || return :wrong_replay_status
        cert.identity_status == :verified || return :wrong_identity_status
        cert.next_boundary == LAURENT_ENDPOINT_REDUCTION_CERTIFICATE_NEXT_BOUNDARY ||
            return :wrong_next_boundary
        cert.dimension == length(column) || return :wrong_dimension
        cert.ring_generators == Suslin._laurent_descent_ring_generators(R) ||
            return :wrong_ring_generators
        cert.context_status == :endpoint_reduction_context ||
            return :wrong_context_status

        context = (;
            case_id = cert.case_id,
            dimension = cert.dimension,
            ring_generators = cert.ring_generators,
            status = cert.context_status,
        )
        operation_status = _laurent_endpoint_reduction_operation_status(
            cert.operation,
            length(column),
            R,
        )
        operation_status == :ok || return operation_status
        replay = _laurent_endpoint_reduction_replay(
            context,
            cert.operation,
            column,
            R;
            require_strict = false,
        )
        replay.source_endpoint == cert.source_endpoint ||
            return :stale_source_endpoint
        replay.target_endpoint == cert.target_endpoint ||
            return :stale_target_endpoint
        hasproperty(cert, :endpoint_measure_relation) ||
            return :missing_certificate_fields
        cert.endpoint_measure_relation == replay.relation ||
            return :wrong_endpoint_measure_relation
        replay.relation == :strict_decrease || return :not_endpoint_reduction
        expected = laurent_endpoint_reduction_certificate(
            context,
            cert.operation,
            column,
            R;
            source_endpoint = cert.source_endpoint,
        )
        cert == expected || return :stale_certificate
        return :ok
    catch err
        err isa InterruptException && rethrow()
        return :operation_replay_failed
    end
end

const CASE008_D14_ENDPOINT_REDUCTION_CERTIFICATE_SUMMARY_CACHE = Ref{Any}(nothing)

function case008_d14_laurent_endpoint_reduction_certificate_summary(
    fixture = ToricBuilderCase008D14ColumnBoundary.boundary_fixture(),
)
    cached = CASE008_D14_ENDPOINT_REDUCTION_CERTIFICATE_SUMMARY_CACHE[]
    cached !== nothing && cached.fixture == fixture && return cached.summary

    report = case008_d14_laurent_endpoint_reduction_search_report(fixture)
    report_validation =
        validate_case008_d14_laurent_endpoint_reduction_search_report(
            report,
            fixture,
        )
    report_validation == :ok ||
        throw(ArgumentError("d14 endpoint-reduction search report must validate; got $(report_validation)"))

    replay_source = _case008_d14_endpoint_reduction_replay_source(fixture)
    columns = _case008_d14_endpoint_reduction_columns(replay_source)
    common = (;
        case_id = report.case_id,
        dimension = report.dimension,
        search_status = report.status,
        search_next_boundary = report.next_boundary,
        candidate_count = report.candidate_count,
        replay_verified_count = report.replay_verified_count,
        source_column = columns.source_column,
        ring = columns.ring,
    )

    if report.status == :candidate_found
        context = case008_d14_laurent_endpoint_reduction_context(
            fixture;
            replay_source,
        )
        candidate = first(report.candidates)
        cert = laurent_endpoint_reduction_certificate(
            context,
            candidate.endpoint_operation,
            columns.source_column,
            columns.ring;
            source_endpoint = report.source_endpoint,
        )
        validate_laurent_endpoint_reduction_certificate(
            cert,
            columns.source_column,
            columns.ring,
        ) == :ok ||
            throw(ArgumentError("first d14 endpoint-reduction candidate did not validate through the certificate shell"))
        cert.target_endpoint == candidate.source_endpoint ||
            throw(ArgumentError("d14 certificate target endpoint must match the replayed source-side candidate endpoint"))
        summary = merge(
            common,
            (;
                d14_status = :endpoint_reduction_certificate,
                certificate = cert,
            ),
        )
        CASE008_D14_ENDPOINT_REDUCTION_CERTIFICATE_SUMMARY_CACHE[] = (;
            fixture,
            summary,
        )
        return summary
    end

    summary = merge(
        common,
        (; d14_status = :endpoint_reduction_search_expansion),
    )
    CASE008_D14_ENDPOINT_REDUCTION_CERTIFICATE_SUMMARY_CACHE[] = (;
        fixture,
        summary,
    )
    return summary
end

@testset "Laurent endpoint-reduction certificate shell" begin
    R, (u, v) = Suslin.suslin_laurent_polynomial_ring(GF(2), ["u", "v"])
    column = [one(R), u + one(R)]
    context = (;
        case_id = "synthetic",
        dimension = 2,
        ring_generators = ("u", "v"),
        status = :endpoint_reduction_context,
    )
    operation = (;
        family = :paired_laurent_endpoint_entry_addition,
        endpoint_index = 2,
        operation = (;
            family = :entry_addition,
            target_index = 2,
            source_index = 1,
            coefficient = 1,
            exponent = (1, 0),
            ring_generators = ("u", "v"),
        ),
        ring_generators = ("u", "v"),
    )

    source_endpoint = _laurent_endpoint_metadata_from_column(
        column,
        R,
        2;
        case_id = "synthetic",
    )
    cert = laurent_endpoint_reduction_certificate(
        context,
        operation,
        column,
        R;
        source_endpoint,
    )

    @test cert.case_id == "synthetic"
    @test cert.dimension == 2
    @test cert.ring_generators == ("u", "v")
    @test cert.context_status == :endpoint_reduction_context
    @test cert.operation == operation
    @test cert.source_endpoint.leading_exponent == (1, 0)
    @test cert.target_endpoint.leading_exponent == (0, 0)
    @test cert.replay_status == :ok
    @test cert.identity_status == :verified
    @test cert.next_boundary == :laurent_normality_replay
    @test cert.status == :endpoint_reduction_certificate
    @test validate_laurent_endpoint_reduction_certificate(cert, column, R) == :ok

    tampered_operation = merge(
        cert,
        (; operation = merge(operation, (; operation = merge(operation.operation, (; coefficient = 0))))),
    )
    @test validate_laurent_endpoint_reduction_certificate(
        tampered_operation,
        column,
        R,
    ) == :stale_target_endpoint

    wrong_generators = merge(cert, (; ring_generators = ("v", "u")))
    @test validate_laurent_endpoint_reduction_certificate(
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
    @test validate_laurent_endpoint_reduction_certificate(
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
    @test validate_laurent_endpoint_reduction_certificate(
        stale_target,
        column,
        R,
    ) == :stale_target_endpoint

    malformed_cert = _laurent_endpoint_reduction_without_field(cert, :status)
    @test validate_laurent_endpoint_reduction_certificate(
        malformed_cert,
        column,
        R,
    ) == :missing_certificate_fields

    wrong_endpoint_index = merge(
        cert,
        (;
            operation = merge(
                operation,
                (; endpoint_index = 1),
            ),
        ),
    )
    @test validate_laurent_endpoint_reduction_certificate(
        wrong_endpoint_index,
        column,
        R,
    ) == :wrong_endpoint_index

    nonstrict_operation = merge(
        operation,
        (; operation = merge(operation.operation, (; exponent = (0, -1)))),
    )
    nonstrict_cert = laurent_endpoint_reduction_certificate(
        context,
        nonstrict_operation,
        column,
        R;
        source_endpoint,
        require_strict = false,
    )
    @test validate_laurent_endpoint_reduction_certificate(
        nonstrict_cert,
        column,
        R,
    ) == :not_endpoint_reduction
end

@testset "case_008 d=14 Laurent endpoint-reduction certificate summary" begin
    summary = case008_d14_laurent_endpoint_reduction_certificate_summary()

    @test summary.case_id == "case_008"
    @test summary.dimension == 14
    @test summary.search_status in (:candidate_found, :exhausted)
    @test summary.d14_status in (
        :endpoint_reduction_certificate,
        :endpoint_reduction_search_expansion,
    )

    if summary.search_status == :candidate_found
        @test summary.d14_status == :endpoint_reduction_certificate
        @test summary.search_next_boundary ==
              :laurent_endpoint_reduction_certificate
        @test summary.candidate_count > 0
        @test summary.replay_verified_count == summary.candidate_count
        @test summary.certificate.case_id == "case_008"
        @test summary.certificate.status == :endpoint_reduction_certificate
        @test summary.certificate.next_boundary == :laurent_normality_replay
        @test validate_laurent_endpoint_reduction_certificate(
            summary.certificate,
            summary.source_column,
            summary.ring,
        ) == :ok
    else
        @test summary.d14_status == :endpoint_reduction_search_expansion
        @test summary.search_next_boundary ==
              :laurent_endpoint_reduction_search_expansion
        @test summary.candidate_count == 0
        @test summary.replay_verified_count == 0
        @test !hasproperty(summary, :certificate)
    end
end
