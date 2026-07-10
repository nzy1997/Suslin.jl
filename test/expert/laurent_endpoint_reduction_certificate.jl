using Test
using Suslin
using Oscar

if !(@isdefined case008_d14_laurent_endpoint_reduction_search_report)
    include(joinpath(@__DIR__, "case008_d14_laurent_endpoint_reduction_search.jl"))
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
    return Suslin._laurent_endpoint_metadata_from_column(
        column,
        R,
        endpoint_index,
        ;
        case_id,
    )
end

function _laurent_endpoint_reduction_operation_status(
    endpoint_operation,
    n::Int,
    R,
)::Symbol
    return Suslin._laurent_endpoint_reduction_status(endpoint_operation, n, R)
end

function _laurent_endpoint_reduction_replay(
    context,
    endpoint_operation,
    column,
    R;
    require_strict::Bool,
)
    return Suslin._replay_laurent_endpoint_reduction(
        column,
        R,
        endpoint_operation;
        case_id = context.case_id,
        require_strict,
    )
end

function laurent_endpoint_reduction_certificate(
    context,
    endpoint_operation,
    column,
    R;
    source_endpoint = nothing,
    require_strict::Bool = true,
)
    return Suslin._laurent_endpoint_reduction_certificate_from_replay(
        context,
        endpoint_operation,
        column,
        R;
        source_endpoint,
        require_strict,
    )
end

function validate_laurent_endpoint_reduction_certificate(cert, column, R)::Symbol
    return Suslin._validate_laurent_endpoint_reduction_certificate(cert, column, R)
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
        exact_source_endpoint = Suslin._laurent_endpoint_metadata_from_column(
            columns.source_column,
            columns.ring,
            candidate.endpoint_operation.endpoint_index;
            case_id = report.case_id,
        )
        cert = laurent_endpoint_reduction_certificate(
            context,
            candidate.endpoint_operation,
            columns.source_column,
            columns.ring;
            source_endpoint = exact_source_endpoint,
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
    @test cert.source_endpoint.entry == column[2]
    @test cert.target_endpoint.entry == one(R)
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

    exact_R, (exact_u, exact_v) =
        Suslin.suslin_laurent_polynomial_ring(QQ, ["u", "v"])
    exact_column = [one(exact_R), exact_u + one(exact_R)]
    exact_context = merge(context, (; case_id = "synthetic-exact"))
    exact_operation = merge(
        operation,
        (;
            operation = merge(
                operation.operation,
                (; coefficient = -exact_u, exponent = (0, 0)),
            ),
        ),
    )
    exact_cert = laurent_endpoint_reduction_certificate(
        exact_context,
        exact_operation,
        exact_column,
        exact_R,
    )
    support_preserving_tamper = merge(
        exact_cert,
        (;
            operation = merge(
                exact_operation,
                (;
                    operation = merge(
                        exact_operation.operation,
                        (; coefficient = -exact_u + one(exact_R)),
                    ),
                ),
            ),
        ),
    )
    @test validate_laurent_endpoint_reduction_certificate(
        support_preserving_tamper,
        exact_column,
        exact_R,
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
