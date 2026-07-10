using Test
using Suslin
using Oscar

if !(@isdefined case008_d14_laurent_endpoint_reduction_search_report)
    include(joinpath(@__DIR__, "..", "expert", "case008_d14_laurent_endpoint_reduction_search.jl"))
end

function _internal_endpoint_without_field(value::NamedTuple, field::Symbol)
    kept = tuple((name for name in keys(value) if name != field)...)
    return NamedTuple{kept}(tuple((getproperty(value, name) for name in kept)...))
end

@testset "internal Laurent endpoint-reduction helpers" begin
    runtests = read(joinpath(@__DIR__, "..", "runtests.jl"), String)
    @test occursin("\"internal/laurent_endpoint_reduction_helpers.jl\"", runtests)

    R, (u, v) = Suslin.suslin_laurent_polynomial_ring(GF(2), ["u", "v"])
    source_column = [one(R), u + one(R)]
    target_column = [one(R), u + v]
    context = (;
        case_id = "synthetic",
        dimension = 2,
        ring_generators = ("u", "v"),
        status = :endpoint_reduction_context,
    )
    operation = (;
        family = :paired_laurent_entry_addition,
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

    @test Suslin._laurent_endpoint_reduction_status(operation, length(source_column), R) == :ok
    replay = Suslin._replay_laurent_endpoint_reduction(
        source_column,
        R,
        operation;
        case_id = context.case_id,
    )
    @test replay.source_endpoint.leading_exponent == (1, 0)
    @test replay.target_endpoint.leading_exponent == (0, 0)
    @test replay.relation == :strict_decrease

    candidate = Suslin._laurent_endpoint_reduction_candidate_from_replay(
        source_column,
        target_column,
        R,
        operation;
        case_id = context.case_id,
    )
    @test candidate.source_endpoint.leading_exponent == (0, 0)
    @test candidate.target_endpoint.leading_exponent == (0, 1)
    @test candidate.source_measure_relation == :strict_decrease
    @test candidate.target_measure_relation == :strict_decrease
    @test candidate.identity_status == :verified
    @test candidate.status == :strict_endpoint_decrease
    @test Suslin._verify_laurent_endpoint_reduction_candidate(
        source_column,
        target_column,
        R,
        candidate,
    )

    cert = Suslin._laurent_endpoint_reduction_certificate_from_replay(
        context,
        operation,
        source_column,
        R;
        source_endpoint = replay.source_endpoint,
    )
    @test cert.case_id == "synthetic"
    @test cert.target_endpoint.leading_exponent == (0, 0)
    @test cert.endpoint_measure_relation == :strict_decrease
    @test cert.next_boundary == :laurent_normality_replay
    @test cert.status == :endpoint_reduction_certificate
    @test Suslin._validate_laurent_endpoint_reduction_certificate(
        cert,
        source_column,
        R,
    ) == :ok

    missing_field = _internal_endpoint_without_field(operation, :family)
    @test Suslin._laurent_endpoint_reduction_status(
        missing_field,
        length(source_column),
        R,
    ) == :malformed_endpoint_operation

    wrong_generators = merge(operation, (; ring_generators = ("v", "u")))
    @test Suslin._laurent_endpoint_reduction_status(
        wrong_generators,
        length(source_column),
        R,
    ) == :wrong_ring_generators

    wrong_endpoint_index = merge(operation, (; endpoint_index = 1))
    @test Suslin._laurent_endpoint_reduction_status(
        wrong_endpoint_index,
        length(source_column),
        R,
    ) == :wrong_endpoint_index

    malformed_nested = merge(
        operation,
        (; operation = _internal_endpoint_without_field(operation.operation, :coefficient)),
    )
    @test Suslin._laurent_endpoint_reduction_status(
        malformed_nested,
        length(source_column),
        R,
    ) == :malformed_endpoint_operation

    stale_source = merge(
        cert,
        (;
            source_endpoint = merge(
                cert.source_endpoint,
                (; term_count = cert.source_endpoint.term_count + 1),
            ),
        ),
    )
    @test Suslin._validate_laurent_endpoint_reduction_certificate(
        stale_source,
        source_column,
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
    @test Suslin._validate_laurent_endpoint_reduction_certificate(
        stale_target,
        source_column,
        R,
    ) == :stale_target_endpoint

    tampered_operation = merge(
        cert,
        (;
            operation = merge(
                operation,
                (; operation = merge(operation.operation, (; coefficient = 0))),
            ),
        ),
    )
    @test Suslin._validate_laurent_endpoint_reduction_certificate(
        tampered_operation,
        source_column,
        R,
    ) == :stale_target_endpoint

    nonstrict_operation = merge(
        operation,
        (; operation = merge(operation.operation, (; exponent = (0, -1)))),
    )
    nonstrict_cert = Suslin._laurent_endpoint_reduction_certificate_from_replay(
        context,
        nonstrict_operation,
        source_column,
        R;
        source_endpoint = replay.source_endpoint,
        require_strict = false,
    )
    @test Suslin._validate_laurent_endpoint_reduction_certificate(
        nonstrict_cert,
        source_column,
        R,
    ) == :not_endpoint_reduction
    nonstrict_candidate = Suslin._laurent_endpoint_reduction_candidate_from_replay(
        source_column,
        target_column,
        R,
        nonstrict_operation;
        case_id = context.case_id,
        require_strict = false,
    )
    @test !Suslin._verify_laurent_endpoint_reduction_candidate(
        source_column,
        target_column,
        R,
        nonstrict_candidate,
    )
end

@testset "internal d14 Laurent endpoint-reduction candidate" begin
    report = case008_d14_laurent_endpoint_reduction_search_report()
    replay_source = _case008_d14_endpoint_reduction_replay_source()
    columns = _case008_d14_endpoint_reduction_columns(replay_source)

    if report.status == :candidate_found
        candidate = first(report.candidates)
        @test Suslin._verify_laurent_endpoint_reduction_candidate(
            columns.source_column,
            columns.target_column,
            columns.ring,
            candidate,
        )
        expected = Suslin._laurent_endpoint_reduction_candidate_from_replay(
            columns.source_column,
            columns.target_column,
            columns.ring,
            candidate.endpoint_operation;
            case_id = report.case_id,
        )
        @test candidate == expected

        context = case008_d14_laurent_endpoint_reduction_context(
            replay_source.fixture;
            replay_source,
        )
        cert = Suslin._laurent_endpoint_reduction_certificate_from_replay(
            context,
            candidate.endpoint_operation,
            columns.source_column,
            columns.ring;
            source_endpoint = report.source_endpoint,
        )
        @test cert.target_endpoint == candidate.source_endpoint
        @test Suslin._validate_laurent_endpoint_reduction_certificate(
            cert,
            columns.source_column,
            columns.ring,
        ) == :ok
    else
        @test report.status == :exhausted
        @test report.candidate_count == 0
        @test report.replay_verified_count == 0
        @test report.next_boundary == :laurent_endpoint_reduction_search_expansion
    end
end
