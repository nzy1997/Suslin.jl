using Test
using Suslin
using Oscar

if !(@isdefined case008_d14_laurent_post_descent_profile_report)
    include(joinpath(@__DIR__, "case008_d14_laurent_post_descent_profile.jl"))
end

const CASE008_D14_LINK_WITNESS_SOURCE_BOUNDARY = :case008_d14_post_descent
const CASE008_D14_LINK_WITNESS_BOUNDARY = :laurent_link_witness
const CASE008_D14_LINK_WITNESS_CONTEXT_STATUS = :link_witness_context
const CASE008_D14_LINK_WITNESS_REQUIRED_WITNESS_FIELDS = (
    :family,
    :pivot_index,
    :partner_index,
    :coefficient,
    :exponent,
    :ring_generators,
)

const CASE008_D14_LINK_WITNESS_REQUIRED_CONTEXT_FIELDS = (
    :case_id,
    :dimension,
    :ring_generators,
    :source_report_boundary,
    :source_boundary,
    :boundary,
    :source_report_status,
    :source_measure,
    :pivot_entry_index,
    :pivot_leading_exponent,
    :pivot_term_count,
    :candidate_partner_indices,
    :required_witness_fields,
    :status,
)

function _case008_d14_link_witness_context_has_required_fields(context)::Bool
    return all(
        field -> hasproperty(context, field),
        CASE008_D14_LINK_WITNESS_REQUIRED_CONTEXT_FIELDS,
    )
end

function _case008_d14_link_witness_candidate_partner_indices(
    dimension::Integer,
    pivot_index::Integer,
)
    1 <= pivot_index <= dimension ||
        throw(ArgumentError("pivot index $(pivot_index) is outside dimension $(dimension)"))
    return tuple((idx for idx in 1:dimension if idx != pivot_index)...)
end

function _case008_d14_link_witness_pivot_candidate(report)
    candidates = report.post_descent_leading_monomial_summary.candidates
    isempty(candidates) &&
        throw(ArgumentError("post-descent leading-monomial summary has no pivot candidates"))
    pivot = first(candidates)
    pivot.entry_index == report.after_measure.leading_entry_index ||
        throw(ArgumentError("pivot entry does not match source measure leading entry"))
    pivot.leading_exponent == report.after_measure.leading_exponent ||
        throw(ArgumentError("pivot leading exponent does not match source measure"))
    support_entry = report.post_descent_support_summary.per_entry[pivot.entry_index]
    support_entry.index == pivot.entry_index ||
        throw(ArgumentError("pivot support summary entry mismatch"))
    support_entry.term_count == pivot.term_count ||
        throw(ArgumentError("pivot term count does not match support summary"))
    return pivot
end

function case008_d14_laurent_link_witness_context(
    report = case008_d14_laurent_post_descent_profile_report(),
)
    validation_result =
        validate_case008_d14_laurent_post_descent_profile_report(report)
    validation_result == :ok ||
        throw(ArgumentError("post-descent source report must validate before link-witness context construction; got $(validation_result)"))

    pivot = _case008_d14_link_witness_pivot_candidate(report)
    return (;
        case_id = report.case_id,
        dimension = report.dimension,
        ring_generators = report.ring_generators,
        source_report_boundary = report.source_boundary,
        source_boundary = CASE008_D14_LINK_WITNESS_SOURCE_BOUNDARY,
        boundary = CASE008_D14_LINK_WITNESS_BOUNDARY,
        source_report_status = report.status,
        source_measure = report.after_measure,
        pivot_entry_index = pivot.entry_index,
        pivot_leading_exponent = pivot.leading_exponent,
        pivot_term_count = pivot.term_count,
        candidate_partner_indices =
            _case008_d14_link_witness_candidate_partner_indices(
                report.dimension,
                pivot.entry_index,
            ),
        required_witness_fields = CASE008_D14_LINK_WITNESS_REQUIRED_WITNESS_FIELDS,
        status = CASE008_D14_LINK_WITNESS_CONTEXT_STATUS,
    )
end

function validate_case008_d14_laurent_link_witness_context(
    context,
    report = case008_d14_laurent_post_descent_profile_report(),
)::Symbol
    validate_case008_d14_laurent_post_descent_profile_report(report) == :ok ||
        return :invalid_source_report
    _case008_d14_link_witness_context_has_required_fields(context) ||
        return :missing_context_fields

    context.status == CASE008_D14_LINK_WITNESS_CONTEXT_STATUS ||
        return :wrong_status
    context.case_id == "case_008" || return :wrong_case
    context.dimension == 14 || return :wrong_dimension
    context.ring_generators == ("u", "v") || return :wrong_ring_generators
    context.source_report_boundary == :case008_d14_original ||
        return :wrong_source_report_boundary
    context.source_boundary == CASE008_D14_LINK_WITNESS_SOURCE_BOUNDARY ||
        return :wrong_source_boundary
    context.boundary == CASE008_D14_LINK_WITNESS_BOUNDARY ||
        return :wrong_boundary
    context.source_report_status == :post_descent_profile_report ||
        return :wrong_source_report_status
    context.source_measure == report.after_measure ||
        return :stale_source_measure
    context.source_measure.whole_support_count == 7378 ||
        return :wrong_source_measure
    context.source_measure.max_entry_terms == 3734 ||
        return :wrong_source_measure
    context.source_measure.valuation_span == (97, 92) ||
        return :wrong_source_measure
    context.pivot_entry_index == 10 || return :wrong_pivot_entry_index
    context.pivot_leading_exponent == (49, -5) ||
        return :wrong_pivot_leading_exponent
    context.pivot_term_count == 3692 || return :wrong_pivot_term_count
    context.candidate_partner_indices ==
        (1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 14) ||
        return :wrong_candidate_partner_indices
    context.required_witness_fields isa Tuple ||
        return :wrong_required_witness_fields
    all(
        field -> field in context.required_witness_fields,
        CASE008_D14_LINK_WITNESS_REQUIRED_WITNESS_FIELDS,
    ) || return :missing_required_witness_field
    context.required_witness_fields ==
        CASE008_D14_LINK_WITNESS_REQUIRED_WITNESS_FIELDS ||
        return :wrong_required_witness_fields

    expected = case008_d14_laurent_link_witness_context(report)
    context == expected || return :stale_context
    return :ok
end

function _case008_d14_replace_tuple_entry(values::Tuple, idx::Int, value)
    return ntuple(j -> j == idx ? value : values[j], length(values))
end

function _case008_d14_context_without_field(context::NamedTuple, field::Symbol)
    kept = tuple((name for name in keys(context) if name != field)...)
    return NamedTuple{kept}(tuple((getproperty(context, name) for name in kept)...))
end

@testset "case_008 d=14 Laurent link-witness context" begin
    runtests = read(joinpath(@__DIR__, "..", "runtests.jl"), String)
    @test occursin(
        "\"expert/case008_d14_laurent_link_witness_context.jl\"",
        runtests,
    )

    report = case008_d14_laurent_post_descent_profile_report()
    @test validate_case008_d14_laurent_post_descent_profile_report(report) == :ok

    context = case008_d14_laurent_link_witness_context(report)
    @test context.case_id == "case_008"
    pivot = first(report.post_descent_leading_monomial_summary.candidates)
    @test context.case_id == "case_008"
    @test context.dimension == 14
    @test context.ring_generators == ("u", "v")
    @test context.source_report_boundary == :case008_d14_original
    @test context.source_boundary == :case008_d14_post_descent
    @test context.boundary == :laurent_link_witness
    @test context.source_report_status == :post_descent_profile_report
    @test context.source_measure == report.after_measure
    @test context.source_measure.whole_support_count == 7378
    @test context.source_measure.max_entry_terms == 3734
    @test context.source_measure.valuation_span == (97, 92)
    @test context.pivot_entry_index == 10
    @test context.pivot_entry_index == pivot.entry_index
    @test context.pivot_leading_exponent == (49, -5)
    @test context.pivot_leading_exponent == pivot.leading_exponent
    @test context.pivot_term_count == 3692
    @test context.pivot_term_count == pivot.term_count
    @test context.candidate_partner_indices ==
          (1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 14)
    @test !(context.pivot_entry_index in context.candidate_partner_indices)
    @test context.required_witness_fields ==
          (:family, :pivot_index, :partner_index, :coefficient, :exponent, :ring_generators)
    @test context.status == :link_witness_context
    @test validate_case008_d14_laurent_link_witness_context(context, report) == :ok

    stale_report = merge(
        report,
        (;
            after_measure = merge(
                report.after_measure,
                (; whole_support_count = report.after_measure.whole_support_count + 1),
            ),
        ),
    )
    @test validate_case008_d14_laurent_post_descent_profile_report(stale_report) ==
          :stale_after_measure
    @test_throws ArgumentError case008_d14_laurent_link_witness_context(stale_report)

    wrong_status_report = merge(report, (; status = :stale_post_descent_report))
    @test validate_case008_d14_laurent_post_descent_profile_report(
        wrong_status_report,
    ) == :wrong_status
    @test_throws ArgumentError case008_d14_laurent_link_witness_context(
        wrong_status_report,
    )

    wrong_relation_report = merge(report, (; measure_relation = :not_strict_decrease))
    @test validate_case008_d14_laurent_post_descent_profile_report(
        wrong_relation_report,
    ) == :wrong_measure_relation
    @test_throws ArgumentError case008_d14_laurent_link_witness_context(
        wrong_relation_report,
    )

    swapped_generators_report = merge(report, (; ring_generators = ("v", "u")))
    @test validate_case008_d14_laurent_post_descent_profile_report(
        swapped_generators_report,
    ) == :wrong_ring_generators
    @test_throws ArgumentError case008_d14_laurent_link_witness_context(
        swapped_generators_report,
    )

    tampered_pivot = merge(pivot, (; entry_index = 9))
    tampered_summary = merge(
        report.post_descent_leading_monomial_summary,
        (;
            candidates = _case008_d14_replace_tuple_entry(
                report.post_descent_leading_monomial_summary.candidates,
                1,
                tampered_pivot,
            ),
        ),
    )
    tampered_pivot_report = merge(
        report,
        (; post_descent_leading_monomial_summary = tampered_summary),
    )
    @test validate_case008_d14_laurent_post_descent_profile_report(
        tampered_pivot_report,
    ) == :wrong_leading_monomial_summary
    @test_throws ArgumentError case008_d14_laurent_link_witness_context(
        tampered_pivot_report,
    )

    missing_partner_field_context = merge(
        context,
        (;
            required_witness_fields =
                (:family, :pivot_index, :coefficient, :exponent, :ring_generators),
        ),
    )
    @test validate_case008_d14_laurent_link_witness_context(
        missing_partner_field_context,
        report,
    ) == :missing_required_witness_field

    missing_schema_context =
        _case008_d14_context_without_field(context, :required_witness_fields)
    @test validate_case008_d14_laurent_link_witness_context(
        missing_schema_context,
        report,
    ) == :missing_context_fields
end
