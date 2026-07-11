using Test
using Suslin
using Oscar

if !(@isdefined case008_d14_laurent_link_witness_certificate_summary)
    include(joinpath(@__DIR__, "laurent_link_witness_certificate.jl"))
end

const CASE008_D14_ENDPOINT_REDUCTION_SOURCE_BOUNDARY =
    :case008_d14_link_witness_certificate
const CASE008_D14_ENDPOINT_REDUCTION_BOUNDARY =
    :laurent_endpoint_reduction
const CASE008_D14_ENDPOINT_REDUCTION_CONTEXT_STATUS =
    :endpoint_reduction_context
const CASE008_D14_ENDPOINT_REDUCTION_REQUIRED_FIELDS = (
    :family,
    :endpoint_index,
    :operation,
    :ring_generators,
)
const CASE008_D14_ENDPOINT_REDUCTION_CONTEXT_FIELDS = (
    :case_id,
    :dimension,
    :ring_generators,
    :source_boundary,
    :boundary,
    :link_witness_status,
    :witness_family,
    :pivot_index,
    :partner_index,
    :witness_exponent,
    :source_endpoint,
    :target_endpoint,
    :required_endpoint_reduction_fields,
    :status,
)
const CASE008_D14_ENDPOINT_REDUCTION_REPLAY_SOURCE_CACHE = IdDict{Any, Any}()

function _case008_d14_endpoint_context_has_required_fields(context)::Bool
    return all(
        field -> hasproperty(context, field),
        CASE008_D14_ENDPOINT_REDUCTION_CONTEXT_FIELDS,
    )
end

function _case008_d14_endpoint_without_field(value::NamedTuple, field::Symbol)
    kept = tuple((name for name in keys(value) if name != field)...)
    return NamedTuple{kept}(tuple((getproperty(value, name) for name in kept)...))
end

function _case008_d14_endpoint_reduction_replay_source()
    cached = _case008_d14_cached_laurent_link_witness_certificate_summary()
    cached !== nothing && return _case008_d14_endpoint_reduction_replay_source(
        cached.fixture;
        summary = cached.summary,
    )
    return _case008_d14_endpoint_reduction_replay_source(
        ToricBuilderCase008D14ColumnBoundary.boundary_fixture(),
    )
end

function _case008_d14_endpoint_reduction_replay_source(
    fixture;
    summary = nothing,
)
    haskey(CASE008_D14_ENDPOINT_REDUCTION_REPLAY_SOURCE_CACHE, fixture) &&
        return CASE008_D14_ENDPOINT_REDUCTION_REPLAY_SOURCE_CACHE[fixture]
    if summary === nothing
        summary = case008_d14_laurent_link_witness_certificate_summary(fixture)
    end
    summary.d14_status == :link_witness_certificate ||
        throw(ArgumentError("case_008 d14 requires a certified Laurent link-witness certificate before endpoint context construction"))
    certificate_validation = validate_laurent_link_witness_certificate(
        summary.certificate,
        summary.source_column,
        summary.ring,
    )
    certificate_validation == :ok ||
        throw(ArgumentError("case_008 d14 link-witness certificate must validate; got $(certificate_validation)"))
    replay_source = (; fixture, summary, certificate_validation)
    CASE008_D14_ENDPOINT_REDUCTION_REPLAY_SOURCE_CACHE[fixture] = replay_source
    return replay_source
end

function case008_d14_laurent_endpoint_reduction_context()
    replay_source = _case008_d14_endpoint_reduction_replay_source()
    return case008_d14_laurent_endpoint_reduction_context(
        replay_source.fixture;
        replay_source,
    )
end

function case008_d14_laurent_endpoint_reduction_context(
    fixture;
    replay_source = _case008_d14_endpoint_reduction_replay_source(fixture),
)
    summary = replay_source.summary
    replay_source.certificate_validation == :ok ||
        throw(ArgumentError("case_008 d14 link-witness certificate must validate; got $(replay_source.certificate_validation)"))
    witness = summary.certificate.witness
    return (;
        case_id = summary.certificate.case_id,
        dimension = summary.certificate.dimension,
        ring_generators = summary.certificate.ring_generators,
        source_boundary = CASE008_D14_ENDPOINT_REDUCTION_SOURCE_BOUNDARY,
        boundary = CASE008_D14_ENDPOINT_REDUCTION_BOUNDARY,
        link_witness_status = summary.certificate.status,
        witness_family = witness.family,
        pivot_index = witness.pivot_index,
        partner_index = witness.partner_index,
        witness_exponent = witness.exponent,
        source_endpoint = summary.certificate.source_endpoint,
        target_endpoint = summary.certificate.target_endpoint,
        required_endpoint_reduction_fields =
            CASE008_D14_ENDPOINT_REDUCTION_REQUIRED_FIELDS,
        status = CASE008_D14_ENDPOINT_REDUCTION_CONTEXT_STATUS,
    )
end

function validate_case008_d14_laurent_endpoint_reduction_context(context)::Symbol
    replay_source = _case008_d14_endpoint_reduction_replay_source()
    return validate_case008_d14_laurent_endpoint_reduction_context(
        context,
        replay_source.fixture,
    )
end

function validate_case008_d14_laurent_endpoint_reduction_context(
    context,
    fixture,
)::Symbol
    try
        _case008_d14_endpoint_context_has_required_fields(context) ||
            return :missing_context_fields
        fixture_validation =
            ToricBuilderCase008D14ColumnBoundary.validate_boundary_fixture(fixture)
        fixture_validation == :ok || return :invalid_fixture
        replay_source = _case008_d14_endpoint_reduction_replay_source(fixture)
        summary = replay_source.summary
        summary.d14_status == :link_witness_certificate ||
            return :missing_link_witness_certificate
        replay_source.certificate_validation == :ok ||
            return :invalid_link_witness_certificate

        context.case_id == "case_008" || return :wrong_case
        context.dimension == 14 || return :wrong_dimension
        context.ring_generators == ("u", "v") || return :wrong_ring_generators
        context.source_boundary ==
            CASE008_D14_ENDPOINT_REDUCTION_SOURCE_BOUNDARY ||
            return :wrong_source_boundary
        context.boundary == CASE008_D14_ENDPOINT_REDUCTION_BOUNDARY ||
            return :wrong_boundary
        context.link_witness_status == :link_witness_certificate ||
            return :wrong_link_witness_status
        context.witness_family == :two_entry_laurent_combination ||
            return :wrong_witness_family
        context.pivot_index == 10 || return :wrong_pivot_index
        context.partner_index == 1 || return :wrong_partner_index
        context.witness_exponent == (1, -1) ||
            return :wrong_witness_exponent
        context.source_endpoint == summary.certificate.source_endpoint ||
            return :stale_source_endpoint
        context.target_endpoint == summary.certificate.target_endpoint ||
            return :stale_target_endpoint
        context.required_endpoint_reduction_fields isa Tuple ||
            return :wrong_required_endpoint_reduction_fields
        all(
            field -> field in context.required_endpoint_reduction_fields,
            CASE008_D14_ENDPOINT_REDUCTION_REQUIRED_FIELDS,
        ) || return :missing_required_endpoint_reduction_field
        context.required_endpoint_reduction_fields ==
            CASE008_D14_ENDPOINT_REDUCTION_REQUIRED_FIELDS ||
            return :wrong_required_endpoint_reduction_fields
        context.status == CASE008_D14_ENDPOINT_REDUCTION_CONTEXT_STATUS ||
            return :wrong_status

        expected = case008_d14_laurent_endpoint_reduction_context(
            fixture;
            replay_source,
        )
        context == expected || return :stale_context
        return :ok
    catch err
        err isa InterruptException && rethrow()
        return :invalid_context
    end
end

@testset "case_008 d=14 Laurent endpoint reduction context" begin
    replay_source = _case008_d14_endpoint_reduction_replay_source()
    fixture = replay_source.fixture
    summary = replay_source.summary
    @test summary.d14_status == :link_witness_certificate
    context = case008_d14_laurent_endpoint_reduction_context(
        fixture;
        replay_source,
    )

    @test context.case_id == "case_008"
    @test context.dimension == 14
    @test context.ring_generators == ("u", "v")
    @test context.source_boundary ==
        CASE008_D14_ENDPOINT_REDUCTION_SOURCE_BOUNDARY
    @test context.boundary == CASE008_D14_ENDPOINT_REDUCTION_BOUNDARY
    @test context.link_witness_status == :link_witness_certificate
    @test context.witness_family == :two_entry_laurent_combination
    @test context.pivot_index == 10
    @test context.partner_index == 1
    @test context.witness_exponent == (1, -1)
    @test context.source_endpoint == summary.certificate.source_endpoint
    @test context.target_endpoint == summary.certificate.target_endpoint
    @test context.required_endpoint_reduction_fields ==
        CASE008_D14_ENDPOINT_REDUCTION_REQUIRED_FIELDS
    @test context.status == CASE008_D14_ENDPOINT_REDUCTION_CONTEXT_STATUS
    @test context == case008_d14_laurent_endpoint_reduction_context(
        fixture;
        replay_source,
    )
    @test validate_case008_d14_laurent_endpoint_reduction_context(
        context,
        fixture,
    ) == :ok
end

@testset "case_008 d=14 Laurent endpoint reduction context validator" begin
    replay_source = _case008_d14_endpoint_reduction_replay_source()
    fixture = replay_source.fixture
    context = case008_d14_laurent_endpoint_reduction_context(
        fixture;
        replay_source,
    )

    @test validate_case008_d14_laurent_endpoint_reduction_context(
        merge(context, (; link_witness_status = :stale_link_witness_certificate)),
        fixture,
    ) == :wrong_link_witness_status

    @test validate_case008_d14_laurent_endpoint_reduction_context(
        merge(context, (; ring_generators = ("v", "u"))),
        fixture,
    ) == :wrong_ring_generators

    @test validate_case008_d14_laurent_endpoint_reduction_context(
        merge(context, (; pivot_index = 9)),
        fixture,
    ) == :wrong_pivot_index

    @test validate_case008_d14_laurent_endpoint_reduction_context(
        merge(context, (; partner_index = 2)),
        fixture,
    ) == :wrong_partner_index

    @test validate_case008_d14_laurent_endpoint_reduction_context(
        merge(
            context,
            (;
                source_endpoint = merge(
                    context.source_endpoint,
                    (; term_count = context.source_endpoint.term_count + 1),
                ),
            ),
        ),
        fixture,
    ) == :stale_source_endpoint

    @test validate_case008_d14_laurent_endpoint_reduction_context(
        merge(
            context,
            (;
                target_endpoint = merge(
                    context.target_endpoint,
                    (; term_count = context.target_endpoint.term_count + 1),
                ),
            ),
        ),
        fixture,
    ) == :stale_target_endpoint

    @test validate_case008_d14_laurent_endpoint_reduction_context(
        merge(context, (; witness_exponent = (0, 0))),
        fixture,
    ) == :wrong_witness_exponent

    @test validate_case008_d14_laurent_endpoint_reduction_context(
        merge(
            context,
            (;
                required_endpoint_reduction_fields = (
                    :family,
                    :endpoint_index,
                    :operation,
                ),
            ),
        ),
        fixture,
    ) == :missing_required_endpoint_reduction_field

    @test validate_case008_d14_laurent_endpoint_reduction_context(
        _case008_d14_endpoint_without_field(
            context,
            :required_endpoint_reduction_fields,
        ),
        fixture,
    ) == :missing_context_fields
end
