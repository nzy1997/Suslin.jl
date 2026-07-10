using Test
using Suslin
using Oscar

if !(@isdefined case008_d14_laurent_endpoint_reduction_context)
    include(joinpath(@__DIR__, "case008_d14_laurent_endpoint_reduction_context.jl"))
end

const CASE008_D14_ENDPOINT_REDUCTION_SEARCH_OPERATION_FAMILIES =
    (:paired_laurent_endpoint_entry_addition,)
const CASE008_D14_ENDPOINT_REDUCTION_OPERATION_SEMANTICS =
    :paired_source_target_endpoint_measure_decrease
const CASE008_D14_ENDPOINT_REDUCTION_ENDPOINT_INDICES = (10,)
const CASE008_D14_ENDPOINT_REDUCTION_SOURCE_INDICES = (1,)
const CASE008_D14_ENDPOINT_REDUCTION_EXPONENT_BOUNDS =
    ((-1, -1), (1, 1))
const CASE008_D14_ENDPOINT_REDUCTION_EXPONENT_VECTORS = Tuple(
    (a, b) for a in -1:1 for b in -1:1
)
const CASE008_D14_ENDPOINT_REDUCTION_COEFFICIENT_FAMILIES = (1,)
const CASE008_D14_ENDPOINT_REDUCTION_SEARCH_REPORT_FIELDS = (
    :case_id,
    :dimension,
    :ring_generators,
    :source_boundary,
    :context_status,
    :boundary,
    :pivot_index,
    :partner_index,
    :witness_exponent,
    :source_endpoint,
    :target_endpoint,
    :required_endpoint_reduction_fields,
    :endpoint_indices,
    :source_indices,
    :operation_families,
    :operation_semantics,
    :exponent_bounds,
    :exponent_vectors,
    :coefficient_families,
    :checked_candidate_count,
    :status,
    :candidate_count,
    :replay_verified_count,
    :next_boundary,
    :candidates,
)
const CASE008_D14_ENDPOINT_REDUCTION_CANDIDATE_FIELDS = (
    :endpoint_operation,
    :source_endpoint,
    :target_endpoint,
    :source_measure_relation,
    :target_measure_relation,
    :replay_status,
    :identity_status,
    :status,
)
const CASE008_D14_ENDPOINT_REDUCTION_NONSTRICT_ERROR =
    "endpoint operation does not strictly decrease both endpoint measures"

function _case008_d14_endpoint_reduction_has_fields(value, fields)::Bool
    return all(field -> hasproperty(value, field), fields)
end

function _case008_d14_endpoint_reduction_without_field(
    value::NamedTuple,
    field::Symbol,
)
    kept = tuple((name for name in keys(value) if name != field)...)
    return NamedTuple{kept}(tuple((getproperty(value, name) for name in kept)...))
end

function _case008_d14_endpoint_reduction_link_operation(witness)
    return (;
        family = :entry_addition,
        target_index = witness.pivot_index,
        source_index = witness.partner_index,
        coefficient = witness.coefficient,
        exponent = witness.exponent,
        ring_generators = witness.ring_generators,
    )
end

function _case008_d14_endpoint_reduction_in_exponent_bounds(exponent)::Bool
    lower, upper = CASE008_D14_ENDPOINT_REDUCTION_EXPONENT_BOUNDS
    return lower[1] <= exponent[1] <= upper[1] &&
           lower[2] <= exponent[2] <= upper[2]
end

function _case008_d14_endpoint_reduction_is_expected_noncandidate_error(err)::Bool
    return err isa ArgumentError &&
           occursin(
               CASE008_D14_ENDPOINT_REDUCTION_NONSTRICT_ERROR,
               sprint(showerror, err),
           )
end

function _case008_d14_endpoint_reduction_columns(replay_source)
    summary = replay_source.summary
    source_column = summary.source_column
    ring = summary.ring
    target_column = Suslin._replay_laurent_elementary_entry_addition(
        source_column,
        ring,
        _case008_d14_endpoint_reduction_link_operation(summary.certificate.witness),
    )
    return (; source_column, target_column, ring)
end

function _case008_d14_endpoint_reduction_operation(
    endpoint_index::Int,
    source_index::Int,
    exponent::Tuple{Int, Int},
    coefficient,
    ring_generators,
)
    return (;
        family = :paired_laurent_endpoint_entry_addition,
        endpoint_index,
        operation = (;
            family = :entry_addition,
            target_index = endpoint_index,
            source_index,
            coefficient,
            exponent,
            ring_generators = Tuple(ring_generators),
        ),
        ring_generators = Tuple(ring_generators),
    )
end

function _case008_d14_endpoint_reduction_operation_status(
    endpoint_operation,
    n::Int,
    R,
)::Symbol
    generic_status = Suslin._laurent_endpoint_reduction_status(
        endpoint_operation,
        n,
        R,
    )
    generic_status == :ok || return generic_status
    endpoint_operation.family == :paired_laurent_endpoint_entry_addition ||
        return :malformed_endpoint_operation
    endpoint_operation.endpoint_index in
        CASE008_D14_ENDPOINT_REDUCTION_ENDPOINT_INDICES ||
        return :wrong_endpoint_index
    endpoint_operation.operation.source_index in
        CASE008_D14_ENDPOINT_REDUCTION_SOURCE_INDICES ||
        return :wrong_source_index
    exponent = try
        Suslin._laurent_descent_exponent_tuple(endpoint_operation.operation.exponent)
    catch err
        err isa InterruptException && rethrow()
        return :malformed_endpoint_operation
    end
    exponent in
        CASE008_D14_ENDPOINT_REDUCTION_EXPONENT_VECTORS ||
        return :wrong_exponent
    _case008_d14_endpoint_reduction_in_exponent_bounds(
        exponent,
    ) || return :wrong_exponent
    endpoint_operation.operation.coefficient in
        CASE008_D14_ENDPOINT_REDUCTION_COEFFICIENT_FAMILIES ||
        return :wrong_coefficient
    return :ok
end

function _case008_d14_endpoint_metadata_from_column(
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

function _case008_d14_endpoint_relation(before_endpoint, after_endpoint)::Symbol
    return Suslin._strictly_decreases_laurent_measure(
        before_endpoint.column_measure,
        after_endpoint.column_measure,
    ) ? :strict_decrease : :not_strict_decrease
end

function _case008_d14_endpoint_reduction_candidate_from_replay(
    context,
    source_column,
    target_column,
    R,
    endpoint_operation;
    require_strict::Bool = true,
)
    status = _case008_d14_endpoint_reduction_operation_status(
        endpoint_operation,
        length(source_column),
        R,
    )
    status == :ok ||
        throw(ArgumentError("invalid Laurent endpoint operation: $(status)"))
    return Suslin._laurent_endpoint_reduction_candidate_from_replay(
        source_column,
        target_column,
        R,
        endpoint_operation;
        case_id = context.case_id,
        require_strict,
    )
end

function validate_case008_d14_laurent_endpoint_reduction_candidate(
    context,
    candidate,
    replay_source;
    require_strict::Bool = true,
)::Symbol
    try
        _case008_d14_endpoint_reduction_has_fields(
            candidate,
            CASE008_D14_ENDPOINT_REDUCTION_CANDIDATE_FIELDS,
        ) || return :missing_candidate_fields
        columns = _case008_d14_endpoint_reduction_columns(replay_source)
        operation_status = _case008_d14_endpoint_reduction_operation_status(
            candidate.endpoint_operation,
            length(columns.source_column),
            columns.ring,
        )
        operation_status == :ok && candidate.replay_status == :ok ||
            return operation_status == :ok ?
                   :wrong_candidate_replay_status : operation_status
        expected = _case008_d14_endpoint_reduction_candidate_from_replay(
            context,
            columns.source_column,
            columns.target_column,
            columns.ring,
            candidate.endpoint_operation;
            require_strict = false,
        )
        candidate.source_endpoint == expected.source_endpoint ||
            return :stale_candidate_source_endpoint
        candidate.target_endpoint == expected.target_endpoint ||
            return :stale_candidate_target_endpoint
        candidate.source_measure_relation == expected.source_measure_relation ||
            return :wrong_candidate_source_relation
        candidate.target_measure_relation == expected.target_measure_relation ||
            return :wrong_candidate_target_relation
        candidate.identity_status == expected.identity_status ||
            return :identity_replay_failed
        candidate.status == expected.status || return :wrong_candidate_status
        require_strict && expected.status != :strict_endpoint_decrease &&
            return :not_endpoint_reduction
        require_strict && !Suslin._verify_laurent_endpoint_reduction_candidate(
            columns.source_column,
            columns.target_column,
            columns.ring,
            candidate,
        ) && return :identity_replay_failed
        return candidate == expected ? :ok : :stale_candidate
    catch err
        err isa InterruptException && rethrow()
        return :invalid_candidate
    end
end

function case008_d14_laurent_endpoint_reduction_search_report(
    fixture = ToricBuilderCase008D14ColumnBoundary.boundary_fixture(),
)::NamedTuple
    replay_source = _case008_d14_endpoint_reduction_replay_source(fixture)
    context = case008_d14_laurent_endpoint_reduction_context(
        fixture;
        replay_source,
    )
    context_validation =
        validate_case008_d14_laurent_endpoint_reduction_context(
            context,
            fixture,
        )
    context_validation == :ok ||
        throw(
            ArgumentError(
                "endpoint-reduction context must validate; got $(context_validation)",
            ),
        )
    columns = _case008_d14_endpoint_reduction_columns(replay_source)

    candidates = NamedTuple[]
    checked_candidate_count = 0
    for endpoint_index in CASE008_D14_ENDPOINT_REDUCTION_ENDPOINT_INDICES
        for source_index in CASE008_D14_ENDPOINT_REDUCTION_SOURCE_INDICES
            for exponent in CASE008_D14_ENDPOINT_REDUCTION_EXPONENT_VECTORS
                for coefficient in CASE008_D14_ENDPOINT_REDUCTION_COEFFICIENT_FAMILIES
                    checked_candidate_count += 1
                    endpoint_operation = _case008_d14_endpoint_reduction_operation(
                        endpoint_index,
                        source_index,
                        exponent,
                        coefficient,
                        context.ring_generators,
                    )
                    candidate = try
                        _case008_d14_endpoint_reduction_candidate_from_replay(
                            context,
                            columns.source_column,
                            columns.target_column,
                            columns.ring,
                            endpoint_operation;
                            require_strict = true,
                        )
                    catch err
                        err isa InterruptException && rethrow()
                        _case008_d14_endpoint_reduction_is_expected_noncandidate_error(
                            err,
                        ) || rethrow()
                        nothing
                    end
                    candidate === nothing && continue
                    validation =
                        validate_case008_d14_laurent_endpoint_reduction_candidate(
                            context,
                            candidate,
                            replay_source,
                        )
                    validation == :ok ||
                        throw(
                            ArgumentError(
                                "endpoint-reduction candidate replay failed validation; got $(validation)",
                            ),
                        )
                    push!(candidates, candidate)
                end
            end
        end
    end

    status = isempty(candidates) ? :exhausted : :candidate_found
    return (;
        case_id = context.case_id,
        dimension = context.dimension,
        ring_generators = context.ring_generators,
        source_boundary = context.source_boundary,
        context_status = context.status,
        boundary = context.boundary,
        pivot_index = context.pivot_index,
        partner_index = context.partner_index,
        witness_exponent = context.witness_exponent,
        source_endpoint = context.source_endpoint,
        target_endpoint = context.target_endpoint,
        required_endpoint_reduction_fields =
            context.required_endpoint_reduction_fields,
        endpoint_indices = CASE008_D14_ENDPOINT_REDUCTION_ENDPOINT_INDICES,
        source_indices = CASE008_D14_ENDPOINT_REDUCTION_SOURCE_INDICES,
        operation_families =
            CASE008_D14_ENDPOINT_REDUCTION_SEARCH_OPERATION_FAMILIES,
        operation_semantics =
            CASE008_D14_ENDPOINT_REDUCTION_OPERATION_SEMANTICS,
        exponent_bounds = CASE008_D14_ENDPOINT_REDUCTION_EXPONENT_BOUNDS,
        exponent_vectors = CASE008_D14_ENDPOINT_REDUCTION_EXPONENT_VECTORS,
        coefficient_families =
            CASE008_D14_ENDPOINT_REDUCTION_COEFFICIENT_FAMILIES,
        checked_candidate_count,
        status,
        candidate_count = length(candidates),
        replay_verified_count = length(candidates),
        next_boundary = status == :candidate_found ?
            :laurent_endpoint_reduction_certificate :
            :laurent_endpoint_reduction_search_expansion,
        candidates = Tuple(candidates),
    )
end

function validate_case008_d14_laurent_endpoint_reduction_search_report(
    report,
    fixture = ToricBuilderCase008D14ColumnBoundary.boundary_fixture(),
)::Symbol
    try
        _case008_d14_endpoint_reduction_has_fields(
            report,
            CASE008_D14_ENDPOINT_REDUCTION_SEARCH_REPORT_FIELDS,
        ) || return :missing_report_fields
        replay_source = _case008_d14_endpoint_reduction_replay_source(fixture)
        context = case008_d14_laurent_endpoint_reduction_context(
            fixture;
            replay_source,
        )
        validate_case008_d14_laurent_endpoint_reduction_context(
            context,
            fixture,
        ) == :ok || return :invalid_context

        report.case_id == "case_008" || return :wrong_case
        report.dimension == 14 || return :wrong_dimension
        report.ring_generators == ("u", "v") || return :wrong_ring_generators
        report.source_boundary == :case008_d14_link_witness_certificate ||
            return :wrong_source_boundary
        report.context_status == :endpoint_reduction_context ||
            return :wrong_context_status
        report.boundary == :laurent_endpoint_reduction ||
            return :wrong_boundary
        report.pivot_index == 10 || return :wrong_pivot_index
        report.partner_index == 1 || return :wrong_partner_index
        report.witness_exponent == (1, -1) || return :wrong_witness_exponent
        report.source_endpoint == context.source_endpoint ||
            return :stale_source_endpoint
        report.target_endpoint == context.target_endpoint ||
            return :stale_target_endpoint
        report.required_endpoint_reduction_fields isa Tuple ||
            return :wrong_required_endpoint_reduction_fields
        all(
            field -> field in report.required_endpoint_reduction_fields,
            CASE008_D14_ENDPOINT_REDUCTION_REQUIRED_FIELDS,
        ) || return :missing_required_endpoint_reduction_field
        report.required_endpoint_reduction_fields ==
            CASE008_D14_ENDPOINT_REDUCTION_REQUIRED_FIELDS ||
            return :wrong_required_endpoint_reduction_fields
        report.endpoint_indices ==
            CASE008_D14_ENDPOINT_REDUCTION_ENDPOINT_INDICES ||
            return :wrong_endpoint_indices
        report.source_indices ==
            CASE008_D14_ENDPOINT_REDUCTION_SOURCE_INDICES ||
            return :wrong_source_indices
        report.operation_families ==
            CASE008_D14_ENDPOINT_REDUCTION_SEARCH_OPERATION_FAMILIES ||
            return :wrong_operation_families
        report.operation_semantics ==
            CASE008_D14_ENDPOINT_REDUCTION_OPERATION_SEMANTICS ||
            return :wrong_operation_semantics
        report.exponent_bounds ==
            CASE008_D14_ENDPOINT_REDUCTION_EXPONENT_BOUNDS ||
            return :wrong_exponent_bounds
        report.exponent_vectors ==
            CASE008_D14_ENDPOINT_REDUCTION_EXPONENT_VECTORS ||
            return :wrong_exponent_vectors
        report.coefficient_families ==
            CASE008_D14_ENDPOINT_REDUCTION_COEFFICIENT_FAMILIES ||
            return :wrong_coefficient_families
        report.checked_candidate_count ==
            length(report.endpoint_indices) *
            length(report.source_indices) *
            length(report.exponent_vectors) *
            length(report.coefficient_families) ||
            return :wrong_checked_candidate_count
        report.checked_candidate_count == 9 ||
            return :wrong_checked_candidate_count
        report.status in (:candidate_found, :exhausted) ||
            return :wrong_status
        report.candidate_count == length(report.candidates) ||
            return :wrong_candidate_count
        report.replay_verified_count == report.candidate_count ||
            return :wrong_replay_verified_count

        for candidate in report.candidates
            candidate_validation =
                validate_case008_d14_laurent_endpoint_reduction_candidate(
                    context,
                    candidate,
                    replay_source,
                )
            candidate_validation == :ok || return candidate_validation
        end
        if report.status == :candidate_found
            report.candidate_count > 0 || return :wrong_candidate_count
            report.next_boundary == :laurent_endpoint_reduction_certificate ||
                return :wrong_next_boundary
        else
            report.candidate_count == 0 || return :wrong_candidate_count
            report.replay_verified_count == 0 ||
                return :wrong_replay_verified_count
            report.next_boundary ==
                :laurent_endpoint_reduction_search_expansion ||
                return :wrong_next_boundary
        end

        expected = case008_d14_laurent_endpoint_reduction_search_report(fixture)
        report == expected || return :stale_report
        return :ok
    catch err
        err isa InterruptException && rethrow()
        return :invalid_report
    end
end

function _case008_d14_endpoint_reduction_report_with_candidate(report, candidate)
    return merge(
        report,
        (;
            status = :candidate_found,
            candidate_count = 1,
            replay_verified_count = 1,
            next_boundary = :laurent_endpoint_reduction_certificate,
            candidates = (candidate,),
        ),
    )
end

@testset "case_008 d=14 Laurent endpoint reduction search report" begin
    runtests = read(joinpath(@__DIR__, "..", "runtests.jl"), String)
    @test occursin(
        "\"expert/case008_d14_laurent_endpoint_reduction_search.jl\"",
        runtests,
    )

    replay_source = _case008_d14_endpoint_reduction_replay_source()
    fixture = replay_source.fixture
    @test hasmethod(
        case008_d14_laurent_endpoint_reduction_search_report,
        Tuple{},
    )

    report = case008_d14_laurent_endpoint_reduction_search_report(fixture)
    @test report.case_id == "case_008"
    @test report.dimension == 14
    @test report.ring_generators == ("u", "v")
    @test report.source_boundary == :case008_d14_link_witness_certificate
    @test report.context_status == :endpoint_reduction_context
    @test report.boundary == :laurent_endpoint_reduction
    @test report.pivot_index == 10
    @test report.partner_index == 1
    @test report.witness_exponent == (1, -1)
    @test report.source_endpoint.entry_index == 10
    @test report.target_endpoint.entry_index == 10
    @test report.required_endpoint_reduction_fields ==
          (:family, :endpoint_index, :operation, :ring_generators)
    @test report.endpoint_indices == (10,)
    @test report.source_indices == (1,)
    @test report.operation_families ==
          (:paired_laurent_endpoint_entry_addition,)
    @test report.operation_semantics ==
          :paired_source_target_endpoint_measure_decrease
    @test report.exponent_bounds == ((-1, -1), (1, 1))
    @test report.exponent_vectors ==
          ((-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 0), (0, 1), (1, -1), (1, 0), (1, 1))
    @test report.coefficient_families == (1,)
    @test report.checked_candidate_count == 9
    @test report.status in (:candidate_found, :exhausted)
    @test report.candidate_count == length(report.candidates)
    @test report.replay_verified_count == report.candidate_count
    @test validate_case008_d14_laurent_endpoint_reduction_search_report(
        report,
        fixture,
    ) == :ok

    if report.status == :candidate_found
        @test report.candidate_count > 0
        @test report.next_boundary == :laurent_endpoint_reduction_certificate
        @test all(candidate -> candidate.replay_status == :ok, report.candidates)
        @test all(
            candidate -> candidate.status == :strict_endpoint_decrease,
            report.candidates,
        )
    else
        @test report.candidate_count == 0
        @test report.replay_verified_count == 0
        @test isempty(report.candidates)
        @test report.next_boundary ==
              :laurent_endpoint_reduction_search_expansion
    end
end

@testset "case_008 d=14 Laurent endpoint reduction search validator controls" begin
    replay_source = _case008_d14_endpoint_reduction_replay_source()
    fixture = replay_source.fixture
    context = case008_d14_laurent_endpoint_reduction_context(
        fixture;
        replay_source,
    )
    report = case008_d14_laurent_endpoint_reduction_search_report(
        fixture,
    )
    columns = _case008_d14_endpoint_reduction_columns(replay_source)
    endpoint_operation = _case008_d14_endpoint_reduction_operation(
        10,
        1,
        (0, 0),
        1,
        ("u", "v"),
    )
    candidate = _case008_d14_endpoint_reduction_candidate_from_replay(
        context,
        columns.source_column,
        columns.target_column,
        columns.ring,
        endpoint_operation;
        require_strict = false,
    )

    @test validate_case008_d14_laurent_endpoint_reduction_search_report(
        merge(report, (; ring_generators = ("v", "u"))),
        fixture,
    ) == :wrong_ring_generators
    @test validate_case008_d14_laurent_endpoint_reduction_search_report(
        merge(
            report,
            (;
                source_endpoint = merge(
                    report.source_endpoint,
                    (; term_count = report.source_endpoint.term_count + 1),
                ),
            ),
        ),
        fixture,
    ) == :stale_source_endpoint
    @test validate_case008_d14_laurent_endpoint_reduction_search_report(
        merge(
            report,
            (;
                target_endpoint = merge(
                    report.target_endpoint,
                    (; term_count = report.target_endpoint.term_count + 1),
                ),
            ),
        ),
        fixture,
    ) == :stale_target_endpoint
    @test validate_case008_d14_laurent_endpoint_reduction_search_report(
        merge(
            report,
            (;
                required_endpoint_reduction_fields =
                    (:family, :endpoint_index, :ring_generators),
            ),
        ),
        fixture,
    ) == :missing_required_endpoint_reduction_field
    @test validate_case008_d14_laurent_endpoint_reduction_search_report(
        _case008_d14_endpoint_reduction_without_field(
            report,
            :required_endpoint_reduction_fields,
        ),
        fixture,
    ) == :missing_report_fields

    malformed_candidate = merge(
        candidate,
        (;
            endpoint_operation = _case008_d14_endpoint_reduction_without_field(
                candidate.endpoint_operation,
                :family,
            ),
        ),
    )
    @test validate_case008_d14_laurent_endpoint_reduction_search_report(
        _case008_d14_endpoint_reduction_report_with_candidate(
            report,
            malformed_candidate,
        ),
        fixture,
    ) == :malformed_endpoint_operation

    wrong_index_candidate = merge(
        candidate,
        (;
            endpoint_operation = merge(
                candidate.endpoint_operation,
                (; endpoint_index = 1),
            ),
        ),
    )
    @test validate_case008_d14_laurent_endpoint_reduction_search_report(
        _case008_d14_endpoint_reduction_report_with_candidate(
            report,
            wrong_index_candidate,
        ),
        fixture,
    ) == :wrong_endpoint_index

    wrong_source_candidate = merge(
        candidate,
        (;
            endpoint_operation = merge(
                candidate.endpoint_operation,
                (;
                    operation = merge(
                        candidate.endpoint_operation.operation,
                        (; source_index = 2),
                    ),
                ),
            ),
        ),
    )
    @test validate_case008_d14_laurent_endpoint_reduction_search_report(
        _case008_d14_endpoint_reduction_report_with_candidate(
            report,
            wrong_source_candidate,
        ),
        fixture,
    ) == :wrong_source_index

    wrong_exponent_candidate = merge(
        candidate,
        (;
            endpoint_operation = merge(
                candidate.endpoint_operation,
                (;
                    operation = merge(
                        candidate.endpoint_operation.operation,
                        (; exponent = (2, 0)),
                    ),
                ),
            ),
        ),
    )
    @test validate_case008_d14_laurent_endpoint_reduction_search_report(
        _case008_d14_endpoint_reduction_report_with_candidate(
            report,
            wrong_exponent_candidate,
        ),
        fixture,
    ) == :wrong_exponent

    wrong_coefficient_candidate = merge(
        candidate,
        (;
            endpoint_operation = merge(
                candidate.endpoint_operation,
                (;
                    operation = merge(
                        candidate.endpoint_operation.operation,
                        (; coefficient = 0),
                    ),
                ),
            ),
        ),
    )
    @test validate_case008_d14_laurent_endpoint_reduction_search_report(
        _case008_d14_endpoint_reduction_report_with_candidate(
            report,
            wrong_coefficient_candidate,
        ),
        fixture,
    ) == :wrong_coefficient

    copied_metadata_operation = _case008_d14_endpoint_reduction_operation(
        10,
        1,
        (-1, -1),
        1,
        ("u", "v"),
    )
    copied_metadata_candidate = merge(
        candidate,
        (; endpoint_operation = copied_metadata_operation),
    )
    @test validate_case008_d14_laurent_endpoint_reduction_search_report(
        _case008_d14_endpoint_reduction_report_with_candidate(
            report,
            copied_metadata_candidate,
        ),
        fixture,
    ) == :stale_candidate_source_endpoint

    stale_candidate = merge(
        candidate,
        (;
            target_endpoint = merge(
                candidate.target_endpoint,
                (; term_count = candidate.target_endpoint.term_count + 1),
            ),
        ),
    )
    @test validate_case008_d14_laurent_endpoint_reduction_search_report(
        _case008_d14_endpoint_reduction_report_with_candidate(
            report,
            stale_candidate,
        ),
        fixture,
    ) == :stale_candidate_target_endpoint
end
