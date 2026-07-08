using Test
using Suslin
using Oscar

if !(@isdefined case008_d14_laurent_link_witness_context)
    include(joinpath(@__DIR__, "case008_d14_laurent_link_witness_context.jl"))
end

const CASE008_D14_LINK_WITNESS_SEARCH_FAMILIES =
    (:two_entry_laurent_combination,)
const CASE008_D14_LINK_WITNESS_SEARCH_SEMANTICS =
    :pivot_plus_shifted_partner_endpoint
const CASE008_D14_LINK_WITNESS_SEARCH_EXPONENT_RADIUS = 1
const CASE008_D14_LINK_WITNESS_SEARCH_EXPONENT_VECTORS = Tuple(
    (a, b)
    for a in -CASE008_D14_LINK_WITNESS_SEARCH_EXPONENT_RADIUS:
        CASE008_D14_LINK_WITNESS_SEARCH_EXPONENT_RADIUS
    for b in -CASE008_D14_LINK_WITNESS_SEARCH_EXPONENT_RADIUS:
        CASE008_D14_LINK_WITNESS_SEARCH_EXPONENT_RADIUS
)
const CASE008_D14_LINK_WITNESS_SEARCH_COEFFICIENT_FAMILY = (1,)

const LAURENT_LINK_WITNESS_FIELDS = (
    :family,
    :pivot_index,
    :partner_index,
    :coefficient,
    :exponent,
    :ring_generators,
)

const LAURENT_LINK_WITNESS_CANDIDATE_FIELDS = (
    :witness,
    :source_endpoint,
    :target_endpoint,
    :replay_status,
    :identity_status,
    :measure_relation,
)

const CASE008_D14_LINK_WITNESS_SEARCH_REPORT_FIELDS = (
    :case_id,
    :dimension,
    :ring_generators,
    :source_boundary,
    :context_status,
    :source_measure,
    :witness_families,
    :witness_semantics,
    :pivot_index,
    :partner_indices,
    :exponent_radius,
    :exponent_vectors,
    :coefficient_family,
    :checked_candidate_count,
    :status,
    :candidate_count,
    :replay_verified_count,
    :next_boundary,
    :candidates,
)

function _laurent_link_witness_has_fields(value, fields)::Bool
    return all(field -> hasproperty(value, field), fields)
end

function _laurent_link_witness_checked_entry_index(
    index,
    n::Int,
    name::AbstractString,
)::Int
    index isa Integer ||
        throw(ArgumentError("$(name) must be an integer entry index"))
    checked = Int(index)
    1 <= checked <= n || throw(ArgumentError("$(name) must be between 1 and $(n)"))
    return checked
end

function _laurent_link_witness_checked_exponent(exponent)::Tuple{Int, Int}
    return _checked_exponent_pair(exponent)
end

function _laurent_link_witness_status(witness, n::Int, R)::Symbol
    _laurent_link_witness_has_fields(witness, LAURENT_LINK_WITNESS_FIELDS) ||
        return :malformed_witness
    witness.family == :two_entry_laurent_combination ||
        return :wrong_witness_family
    witness.ring_generators == _ring_generator_names(R) ||
        return :wrong_ring_generators
    try
        pivot = _laurent_link_witness_checked_entry_index(
            witness.pivot_index,
            n,
            "pivot_index",
        )
        partner = _laurent_link_witness_checked_entry_index(
            witness.partner_index,
            n,
            "partner_index",
        )
        pivot != partner || return :pivot_partner_index_equality
        _laurent_link_witness_checked_exponent(witness.exponent)
        R(witness.coefficient)
    catch err
        err isa InterruptException && rethrow()
        return :malformed_witness
    end
    return :ok
end

function _laurent_link_witness_operation(witness)
    return (;
        family = :entry_addition,
        target_index = witness.pivot_index,
        source_index = witness.partner_index,
        coefficient = witness.coefficient,
        exponent = witness.exponent,
        ring_generators = witness.ring_generators,
    )
end

function _laurent_link_endpoint_metadata(
    entry,
    R,
    entry_index::Int,
    column_measure;
    case_id,
)
    support = _entry_support(R(entry))
    return (;
        case_id,
        status = :link_witness_endpoint_metadata,
        entry_index,
        ring_generators = _ring_generator_names(R),
        is_zero = isempty(support),
        term_count = length(support),
        support_bounds = _support_bounds(support),
        leading_exponent = isempty(support) ? nothing : maximum(support),
        column_measure,
    )
end

function _laurent_link_witness_candidate_from_replay(
    column,
    R,
    witness;
    case_id,
    require_strict::Bool = true,
)
    status = _laurent_link_witness_status(witness, length(column), R)
    status == :ok ||
        throw(ArgumentError("invalid Laurent link witness: $(status)"))

    pivot = Int(witness.pivot_index)
    before_measure = Suslin._laurent_descent_measure_from_column(
        column,
        R;
        case_id,
    )
    after_column = Suslin._replay_laurent_elementary_entry_addition(
        column,
        R,
        _laurent_link_witness_operation(witness),
    )
    after_measure = Suslin._laurent_descent_measure_from_column(
        after_column,
        R;
        case_id,
    )
    relation = strictly_decreases_laurent_measure(before_measure, after_measure)
    require_strict && !relation &&
        throw(ArgumentError("witness does not strictly decrease the replay measure"))
    return (;
        witness,
        source_endpoint = _laurent_link_endpoint_metadata(
            column[pivot],
            R,
            pivot,
            before_measure;
            case_id,
        ),
        target_endpoint = _laurent_link_endpoint_metadata(
            after_column[pivot],
            R,
            pivot,
            after_measure;
            case_id,
        ),
        replay_status = :ok,
        identity_status = :verified,
        measure_relation = relation ? :strict_decrease : :not_strict_decrease,
    )
end

function verify_laurent_link_witness_candidate(column, R, candidate)::Bool
    try
        _laurent_link_witness_has_fields(
            candidate,
            LAURENT_LINK_WITNESS_CANDIDATE_FIELDS,
        ) || return false
        candidate.replay_status == :ok || return false
        candidate.identity_status == :verified || return false
        candidate.measure_relation == :strict_decrease || return false

        witness = candidate.witness
        _laurent_link_witness_status(witness, length(column), R) == :ok ||
            return false
        case_id = candidate.source_endpoint.case_id
        expected = _laurent_link_witness_candidate_from_replay(
            column,
            R,
            witness;
            case_id,
        )
        return candidate == expected
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end

function _case008_d14_link_witness_source_data(fixture)
    source_report = case008_d14_laurent_post_descent_profile_report(fixture)
    source_validation =
        validate_case008_d14_laurent_post_descent_profile_report(
            source_report,
            fixture,
        )
    source_validation == :ok ||
        throw(ArgumentError("post-descent report must validate; got $(source_validation)"))
    context = case008_d14_laurent_link_witness_context(
        source_report;
        source_validation,
    )
    context_validation = validate_case008_d14_laurent_link_witness_context(
        context,
        source_report;
        source_validation,
    )
    context_validation == :ok ||
        throw(ArgumentError("link-witness context must validate; got $(context_validation)"))
    replay = _case008_d14_replayed_post_descent_data(fixture)
    return (; source_report, source_validation, context, replay)
end

function _case008_d14_link_witness(
    pivot::Int,
    partner::Int,
    exponent::Tuple{Int, Int},
    coefficient,
    ring_generators,
)
    return (;
        family = :two_entry_laurent_combination,
        pivot_index = pivot,
        partner_index = partner,
        coefficient,
        exponent,
        ring_generators = Tuple(ring_generators),
    )
end

function case008_d14_laurent_link_witness_search_report(
    fixture = ToricBuilderCase008D14ColumnBoundary.boundary_fixture(),
)::NamedTuple
    source = _case008_d14_link_witness_source_data(fixture)
    context = source.context
    column = source.replay.after_column
    R = fixture.ring
    generator_names = _ring_generator_names(R)
    supports = _case008_d14_column_support_sets(column)
    baseline_measure = Suslin._laurent_descent_measure_from_column(
        column,
        R;
        case_id = context.case_id,
    )
    support_measure = _case008_d14_measure_from_supports(
        supports,
        generator_names;
        case_id = context.case_id,
    )
    support_measure == baseline_measure ||
        throw(ArgumentError("support-set baseline measure does not match replayed post-descent measure"))

    candidates = NamedTuple[]
    checked_candidate_count = 0
    for partner in context.candidate_partner_indices
        for exponent in CASE008_D14_LINK_WITNESS_SEARCH_EXPONENT_VECTORS
            for coefficient in CASE008_D14_LINK_WITNESS_SEARCH_COEFFICIENT_FAMILY
                checked_candidate_count += 1
                witness = _case008_d14_link_witness(
                    context.pivot_entry_index,
                    partner,
                    exponent,
                    coefficient,
                    generator_names,
                )
                after_supports = _case008_d14_update_entry_addition_supports(
                    supports,
                    R,
                    _laurent_link_witness_operation(witness),
                )
                after_measure = _case008_d14_measure_from_supports(
                    after_supports,
                    generator_names;
                    case_id = context.case_id,
                )
                strictly_decreases_laurent_measure(
                    baseline_measure,
                    after_measure,
                ) || continue
                candidate = _laurent_link_witness_candidate_from_replay(
                    column,
                    R,
                    witness;
                    case_id = context.case_id,
                )
                verify_laurent_link_witness_candidate(column, R, candidate) ||
                    continue
                push!(candidates, candidate)
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
        source_measure = context.source_measure,
        witness_families = CASE008_D14_LINK_WITNESS_SEARCH_FAMILIES,
        witness_semantics = CASE008_D14_LINK_WITNESS_SEARCH_SEMANTICS,
        pivot_index = context.pivot_entry_index,
        partner_indices = context.candidate_partner_indices,
        exponent_radius = CASE008_D14_LINK_WITNESS_SEARCH_EXPONENT_RADIUS,
        exponent_vectors = CASE008_D14_LINK_WITNESS_SEARCH_EXPONENT_VECTORS,
        coefficient_family = CASE008_D14_LINK_WITNESS_SEARCH_COEFFICIENT_FAMILY,
        checked_candidate_count,
        status,
        candidate_count = length(candidates),
        replay_verified_count = length(candidates),
        next_boundary = status == :candidate_found ?
            :laurent_link_witness_certificate :
            :laurent_link_witness_search_expansion,
        candidates = Tuple(candidates),
    )
end

function validate_case008_d14_laurent_link_witness_search_report(
    report,
    fixture = ToricBuilderCase008D14ColumnBoundary.boundary_fixture(),
)::Symbol
    try
        _laurent_link_witness_has_fields(
            report,
            CASE008_D14_LINK_WITNESS_SEARCH_REPORT_FIELDS,
        ) || return :missing_report_fields
        source = _case008_d14_link_witness_source_data(fixture)
        context = source.context
        report.case_id == "case_008" || return :wrong_case
        report.dimension == 14 || return :wrong_dimension
        report.ring_generators == ("u", "v") || return :wrong_ring_generators
        report.source_boundary == :case008_d14_post_descent ||
            return :wrong_source_boundary
        report.context_status == :link_witness_context ||
            return :wrong_context_status
        report.source_measure == context.source_measure ||
            return :stale_source_measure
        report.witness_families == CASE008_D14_LINK_WITNESS_SEARCH_FAMILIES ||
            return :wrong_witness_families
        report.witness_semantics == CASE008_D14_LINK_WITNESS_SEARCH_SEMANTICS ||
            return :wrong_witness_semantics
        report.pivot_index == 10 || return :wrong_pivot_index
        report.partner_indices == context.candidate_partner_indices ||
            return :wrong_partner_indices
        report.exponent_radius == CASE008_D14_LINK_WITNESS_SEARCH_EXPONENT_RADIUS ||
            return :wrong_exponent_radius
        report.exponent_vectors == CASE008_D14_LINK_WITNESS_SEARCH_EXPONENT_VECTORS ||
            return :wrong_exponent_vectors
        report.coefficient_family ==
            CASE008_D14_LINK_WITNESS_SEARCH_COEFFICIENT_FAMILY ||
            return :wrong_coefficient_family
        report.checked_candidate_count ==
            length(report.partner_indices) *
            length(report.exponent_vectors) *
            length(report.coefficient_family) ||
            return :wrong_checked_candidate_count
        report.checked_candidate_count == 117 ||
            return :wrong_checked_candidate_count
        report.status in (:candidate_found, :exhausted) || return :wrong_status
        report.candidate_count == length(report.candidates) ||
            return :wrong_candidate_count
        report.replay_verified_count == report.candidate_count ||
            return :wrong_replay_verified_count

        column = source.replay.after_column
        R = fixture.ring
        for candidate in report.candidates
            verify_laurent_link_witness_candidate(column, R, candidate) ||
                return :invalid_candidate
        end
        if report.status == :candidate_found
            report.candidate_count > 0 || return :wrong_candidate_count
            report.next_boundary == :laurent_link_witness_certificate ||
                return :wrong_next_boundary
        else
            report.candidate_count == 0 || return :wrong_candidate_count
            report.replay_verified_count == 0 || return :wrong_replay_verified_count
            report.next_boundary == :laurent_link_witness_search_expansion ||
                return :wrong_next_boundary
        end

        expected = case008_d14_laurent_link_witness_search_report(fixture)
        report == expected || return :stale_report
        return :ok
    catch err
        err isa InterruptException && rethrow()
        return :invalid_report
    end
end

function _laurent_link_replace_tuple_entry(values::Tuple, idx::Int, value)
    return ntuple(j -> j == idx ? value : values[j], length(values))
end

function _laurent_link_without_field(value::NamedTuple, field::Symbol)
    kept = tuple((name for name in keys(value) if name != field)...)
    return NamedTuple{kept}(tuple((getproperty(value, name) for name in kept)...))
end

@testset "case_008 d=14 Laurent link-witness search report" begin
    runtests = read(joinpath(@__DIR__, "..", "runtests.jl"), String)
    @test occursin(
        "\"expert/case008_d14_laurent_link_witness_search.jl\"",
        runtests,
    )

    report = case008_d14_laurent_link_witness_search_report()
    @test report.case_id == "case_008"
    @test report.dimension == 14
    @test report.context_status == :link_witness_context
    @test report.witness_families == (:two_entry_laurent_combination,)
    @test report.witness_semantics == :pivot_plus_shifted_partner_endpoint
    @test report.pivot_index == 10
    @test report.partner_indices ==
          (1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 14)
    @test report.exponent_radius == 1
    @test report.exponent_vectors ==
          ((-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 0), (0, 1), (1, -1), (1, 0), (1, 1))
    @test report.coefficient_family == (1,)
    @test report.checked_candidate_count == 117
    @test report.status in (:candidate_found, :exhausted)
    @test report.candidate_count == length(report.candidates)
    @test report.replay_verified_count == report.candidate_count
    @test validate_case008_d14_laurent_link_witness_search_report(report) == :ok

    if report.status == :candidate_found
        @test report.candidate_count > 0
        @test report.next_boundary == :laurent_link_witness_certificate
        @test all(candidate -> candidate.replay_status == :ok, report.candidates)
        @test all(candidate -> candidate.identity_status == :verified, report.candidates)
    else
        @test report.candidate_count == 0
        @test report.replay_verified_count == 0
        @test isempty(report.candidates)
        @test report.next_boundary == :laurent_link_witness_search_expansion
    end

    @test validate_case008_d14_laurent_link_witness_search_report(
        merge(report, (; checked_candidate_count = report.checked_candidate_count + 1)),
    ) == :wrong_checked_candidate_count
    @test validate_case008_d14_laurent_link_witness_search_report(
        merge(report, (; witness_semantics = :reversed_partner_plus_pivot_endpoint)),
    ) == :wrong_witness_semantics
    @test validate_case008_d14_laurent_link_witness_search_report(
        _laurent_link_without_field(report, :witness_semantics),
    ) == :missing_report_fields

    if !isempty(report.candidates)
        tampered_candidate = merge(
            first(report.candidates),
            (;
                target_endpoint = merge(
                    first(report.candidates).target_endpoint,
                    (; term_count = first(report.candidates).target_endpoint.term_count + 1),
                ),
            ),
        )
        tampered_report = merge(
            report,
            (; candidates = _laurent_link_replace_tuple_entry(report.candidates, 1, tampered_candidate)),
        )
        @test validate_case008_d14_laurent_link_witness_search_report(
            tampered_report,
        ) == :invalid_candidate
    end
end

@testset "Laurent link-witness candidate verifier controls" begin
    R, (u, v) = Suslin.suslin_laurent_polynomial_ring(GF(2), ["u", "v"])
    column = [one(R), u + v]
    good_witness = (;
        family = :two_entry_laurent_combination,
        pivot_index = 2,
        partner_index = 1,
        coefficient = 1,
        exponent = (1, 0),
        ring_generators = ("u", "v"),
    )
    good_candidate = _laurent_link_witness_candidate_from_replay(
        column,
        R,
        good_witness;
        case_id = "synthetic",
    )
    @test verify_laurent_link_witness_candidate(column, R, good_candidate)
    @test good_candidate.source_endpoint.leading_exponent == (1, 0)
    @test good_candidate.target_endpoint.leading_exponent == (0, 1)
    @test good_candidate.measure_relation == :strict_decrease

    missing_family = merge(
        good_candidate,
        (; witness = _laurent_link_without_field(good_witness, :family)),
    )
    @test !verify_laurent_link_witness_candidate(column, R, missing_family)

    equal_indices = merge(
        good_candidate,
        (; witness = merge(good_witness, (; partner_index = 2))),
    )
    @test !verify_laurent_link_witness_candidate(column, R, equal_indices)

    wrong_generators = merge(
        good_candidate,
        (; witness = merge(good_witness, (; ring_generators = ("v", "u")))),
    )
    @test !verify_laurent_link_witness_candidate(column, R, wrong_generators)

    stale_target = merge(
        good_candidate,
        (;
            target_endpoint = merge(
                good_candidate.target_endpoint,
                (; term_count = good_candidate.target_endpoint.term_count + 1),
            ),
        ),
    )
    @test !verify_laurent_link_witness_candidate(column, R, stale_target)

    reversed_roles = merge(
        good_candidate,
        (;
            witness = merge(
                good_witness,
                (; pivot_index = good_witness.partner_index, partner_index = good_witness.pivot_index),
            ),
        ),
    )
    @test !verify_laurent_link_witness_candidate(column, R, reversed_roles)

    nonwitness = _laurent_link_witness_candidate_from_replay(
        column,
        R,
        merge(good_witness, (; exponent = (0, 0)));
        case_id = "synthetic",
        require_strict = false,
    )
    @test nonwitness.measure_relation == :not_strict_decrease
    @test !verify_laurent_link_witness_candidate(column, R, nonwitness)

    nonidentity = merge(
        good_candidate,
        (; target_endpoint = good_candidate.source_endpoint),
    )
    @test !verify_laurent_link_witness_candidate(column, R, nonidentity)
end
