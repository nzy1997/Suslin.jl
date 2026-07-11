using Test
using Suslin
using Oscar

if !(@isdefined case008_d14_laurent_descent_step_certificate)
    include(joinpath(@__DIR__, "laurent_descent_step_certificate.jl"))
end

const CASE008_D14_POST_DESCENT_SOURCE_BOUNDARY = :case008_d14_original
const CASE008_D14_POST_DESCENT_STATUS = :post_descent_profile_report

const CASE008_D14_EXPECTED_POST_DESCENT_ENTRY_TERM_COUNTS = (
    3661,
    3574,
    3554,
    3561,
    3595,
    3734,
    3622,
    3454,
    3489,
    3692,
    3675,
    3600,
    3693,
    3495,
)

const CASE008_D14_EXPECTED_POST_DESCENT_SUPPORT_SUMMARY = (
    generator_order = ("u", "v"),
    entry_count = 14,
    per_entry = (
        (index = 1, support_count = 3661, term_count = 3661, support_bounds = (min_exponents = (-47, -47), max_exponents = (48, 44))),
        (index = 2, support_count = 3574, term_count = 3574, support_bounds = (min_exponents = (-46, -48), max_exponents = (48, 44))),
        (index = 3, support_count = 3554, term_count = 3554, support_bounds = (min_exponents = (-47, -48), max_exponents = (48, 44))),
        (index = 4, support_count = 3561, term_count = 3561, support_bounds = (min_exponents = (-48, -47), max_exponents = (46, 42))),
        (index = 5, support_count = 3595, term_count = 3595, support_bounds = (min_exponents = (-46, -47), max_exponents = (48, 44))),
        (index = 6, support_count = 3734, term_count = 3734, support_bounds = (min_exponents = (-48, -48), max_exponents = (48, 44))),
        (index = 7, support_count = 3622, term_count = 3622, support_bounds = (min_exponents = (-47, -48), max_exponents = (48, 43))),
        (index = 8, support_count = 3454, term_count = 3454, support_bounds = (min_exponents = (-47, -46), max_exponents = (46, 43))),
        (index = 9, support_count = 3489, term_count = 3489, support_bounds = (min_exponents = (-47, -46), max_exponents = (46, 42))),
        (index = 10, support_count = 3692, term_count = 3692, support_bounds = (min_exponents = (-47, -48), max_exponents = (49, 43))),
        (index = 11, support_count = 3675, term_count = 3675, support_bounds = (min_exponents = (-48, -48), max_exponents = (48, 44))),
        (index = 12, support_count = 3600, term_count = 3600, support_bounds = (min_exponents = (-48, -48), max_exponents = (47, 44))),
        (index = 13, support_count = 3693, term_count = 3693, support_bounds = (min_exponents = (-47, -48), max_exponents = (48, 44))),
        (index = 14, support_count = 3495, term_count = 3495, support_bounds = (min_exponents = (-47, -46), max_exponents = (46, 43))),
    ),
    per_entry_support_counts = (3661, 3574, 3554, 3561, 3595, 3734, 3622, 3454, 3489, 3692, 3675, 3600, 3693, 3495),
    whole_column_support_count = 7378,
    whole_column_bounds = (min_exponents = (-48, -48), max_exponents = (49, 44)),
)

const CASE008_D14_EXPECTED_POST_DESCENT_VALUATION_SUMMARY =
    (u = (min = -48, max = 49), v = (min = -48, max = 44))

const CASE008_D14_EXPECTED_POST_DESCENT_LEADING_MONOMIAL_SUMMARY = (
    generator_order = ("u", "v"),
    order = :lexicographic_descending,
    tie_breaker = :entry_index_ascending,
    candidate_count = 14,
    candidates = (
        (entry_index = 10, leading_exponent = (49, -5), term_count = 3692, support_bounds = (min_exponents = (-47, -48), max_exponents = (49, 43))),
        (entry_index = 1, leading_exponent = (48, -4), term_count = 3661, support_bounds = (min_exponents = (-47, -47), max_exponents = (48, 44))),
        (entry_index = 2, leading_exponent = (48, -4), term_count = 3574, support_bounds = (min_exponents = (-46, -48), max_exponents = (48, 44))),
        (entry_index = 13, leading_exponent = (48, -4), term_count = 3693, support_bounds = (min_exponents = (-47, -48), max_exponents = (48, 44))),
        (entry_index = 3, leading_exponent = (48, -5), term_count = 3554, support_bounds = (min_exponents = (-47, -48), max_exponents = (48, 44))),
        (entry_index = 5, leading_exponent = (48, -5), term_count = 3595, support_bounds = (min_exponents = (-46, -47), max_exponents = (48, 44))),
        (entry_index = 6, leading_exponent = (48, -5), term_count = 3734, support_bounds = (min_exponents = (-48, -48), max_exponents = (48, 44))),
        (entry_index = 7, leading_exponent = (48, -5), term_count = 3622, support_bounds = (min_exponents = (-47, -48), max_exponents = (48, 43))),
        (entry_index = 11, leading_exponent = (48, -5), term_count = 3675, support_bounds = (min_exponents = (-48, -48), max_exponents = (48, 44))),
        (entry_index = 12, leading_exponent = (47, -5), term_count = 3600, support_bounds = (min_exponents = (-48, -48), max_exponents = (47, 44))),
        (entry_index = 9, leading_exponent = (46, -4), term_count = 3489, support_bounds = (min_exponents = (-47, -46), max_exponents = (46, 42))),
        (entry_index = 14, leading_exponent = (46, -4), term_count = 3495, support_bounds = (min_exponents = (-47, -46), max_exponents = (46, 43))),
        (entry_index = 4, leading_exponent = (46, -5), term_count = 3561, support_bounds = (min_exponents = (-48, -47), max_exponents = (46, 42))),
        (entry_index = 8, leading_exponent = (46, -5), term_count = 3454, support_bounds = (min_exponents = (-47, -46), max_exponents = (46, 43))),
    ),
)

const CASE008_D14_REQUIRED_POST_DESCENT_REPORT_FIELDS = (
    :case_id,
    :dimension,
    :source_boundary,
    :ring_generators,
    :operation_family,
    :operation,
    :replay_status,
    :before_measure,
    :after_measure,
    :measure_relation,
    :post_descent_profile,
    :post_descent_support_summary,
    :post_descent_valuation_summary,
    :post_descent_leading_monomial_summary,
    :status,
)

const CASE008_D14_POST_DESCENT_FIXTURE_VALIDATION_CACHE = IdDict{Any, Symbol}()
const CASE008_D14_REPLAYED_POST_DESCENT_DATA_CACHE = IdDict{Any, Any}()
const CASE008_D14_POST_DESCENT_PROFILE_REPORT_CACHE = IdDict{Any, Any}()

function _case008_d14_post_descent_fixture_validation(fixture)::Symbol
    haskey(CASE008_D14_POST_DESCENT_FIXTURE_VALIDATION_CACHE, fixture) &&
        return CASE008_D14_POST_DESCENT_FIXTURE_VALIDATION_CACHE[fixture]
    validation = ToricBuilderCase008D14ColumnBoundary.validate_boundary_fixture(
        fixture,
    )
    CASE008_D14_POST_DESCENT_FIXTURE_VALIDATION_CACHE[fixture] = validation
    return validation
end

function _case008_d14_cached_laurent_post_descent_profile_report()
    isempty(CASE008_D14_POST_DESCENT_PROFILE_REPORT_CACHE) && return nothing
    return first(values(CASE008_D14_POST_DESCENT_PROFILE_REPORT_CACHE))
end

function _case008_d14_post_descent_default_fixture()
    cached = _case008_d14_cached_laurent_post_descent_profile_report()
    cached !== nothing && return cached.fixture
    return ToricBuilderCase008D14ColumnBoundary.boundary_fixture()
end

function _post_descent_report_has_required_fields(report)::Bool
    return all(
        field -> hasproperty(report, field),
        CASE008_D14_REQUIRED_POST_DESCENT_REPORT_FIELDS,
    )
end

function _case008_d14_replayed_post_descent_data(fixture)
    haskey(CASE008_D14_REPLAYED_POST_DESCENT_DATA_CACHE, fixture) &&
        return CASE008_D14_REPLAYED_POST_DESCENT_DATA_CACHE[fixture]
    fixture_validation = _case008_d14_post_descent_fixture_validation(fixture)
    fixture_validation == :ok ||
        throw(ArgumentError("invalid case_008 d14 fixture: $(fixture_validation)"))

    operation = CASE008_D14_RECORDED_DESCENT_OPERATION
    before_measure = Suslin._laurent_descent_measure_from_column(
        fixture.failing_column,
        fixture.ring;
        case_id = fixture.case_id,
    )
    after_column = Suslin._replay_laurent_elementary_entry_addition(
        fixture.failing_column,
        fixture.ring,
        operation,
    )
    after_measure = Suslin._laurent_descent_measure_from_column(
        after_column,
        fixture.ring;
        case_id = fixture.case_id,
    )
    post_descent_profile = laurent_descent_step_profile(
        after_column,
        fixture.ring;
        case_id = fixture.case_id,
    )
    replay = (;
        operation,
        before_measure,
        after_column,
        after_measure,
        post_descent_profile,
    )
    CASE008_D14_REPLAYED_POST_DESCENT_DATA_CACHE[fixture] = replay
    return replay
end

function case008_d14_laurent_post_descent_profile_report()
    cached = _case008_d14_cached_laurent_post_descent_profile_report()
    cached !== nothing && return cached.report
    return case008_d14_laurent_post_descent_profile_report(
        _case008_d14_post_descent_default_fixture(),
    )
end

function case008_d14_laurent_post_descent_profile_report(fixture)
    haskey(CASE008_D14_POST_DESCENT_PROFILE_REPORT_CACHE, fixture) &&
        return CASE008_D14_POST_DESCENT_PROFILE_REPORT_CACHE[fixture].report
    replay = _case008_d14_replayed_post_descent_data(fixture)
    relation = strictly_decreases_laurent_measure(
        replay.before_measure,
        replay.after_measure,
    ) ? :strict_decrease : :not_strict_decrease
    report = (;
        case_id = fixture.case_id,
        dimension = length(fixture.failing_column),
        source_boundary = CASE008_D14_POST_DESCENT_SOURCE_BOUNDARY,
        ring_generators = _ring_generator_names(fixture.ring),
        operation_family = replay.operation.family,
        operation = replay.operation,
        replay_status = :ok,
        before_measure = replay.before_measure,
        after_measure = replay.after_measure,
        measure_relation = relation,
        post_descent_profile = replay.post_descent_profile,
        post_descent_support_summary =
            replay.post_descent_profile.newton_support_summary,
        post_descent_valuation_summary = replay.post_descent_profile.valuation_ranges,
        post_descent_leading_monomial_summary =
            replay.post_descent_profile.leading_monomial_candidates,
        status = CASE008_D14_POST_DESCENT_STATUS,
    )
    CASE008_D14_POST_DESCENT_PROFILE_REPORT_CACHE[fixture] = (;
        fixture,
        report,
    )
    return report
end

function validate_case008_d14_laurent_post_descent_profile_report(report)::Symbol
    return validate_case008_d14_laurent_post_descent_profile_report(
        report,
        _case008_d14_post_descent_default_fixture(),
    )
end

function validate_case008_d14_laurent_post_descent_profile_report(
    report,
    fixture,
)::Symbol
    fixture_validation = _case008_d14_post_descent_fixture_validation(fixture)
    fixture_validation == :ok || return :invalid_fixture
    _post_descent_report_has_required_fields(report) ||
        return :missing_report_fields
    report.status == CASE008_D14_POST_DESCENT_STATUS || return :wrong_status
    report.case_id == "case_008" || return :wrong_case
    report.dimension == 14 || return :wrong_dimension
    report.source_boundary == CASE008_D14_POST_DESCENT_SOURCE_BOUNDARY ||
        return :wrong_source_boundary
    report.ring_generators == ("u", "v") || return :wrong_ring_generators
    report.operation_family == :entry_addition || return :wrong_operation_family
    report.replay_status == :ok || return :wrong_replay_status
    report.measure_relation == :strict_decrease || return :wrong_measure_relation

    expected = case008_d14_laurent_post_descent_profile_report(fixture)
    report.operation == expected.operation || return :stale_operation
    report.before_measure == expected.before_measure ||
        return :stale_before_measure
    report.after_measure == expected.after_measure || return :stale_after_measure
    strictly_decreases_laurent_measure(report.before_measure, report.after_measure) ||
        return :not_strict_decrease
    report.post_descent_profile == expected.post_descent_profile ||
        return :stale_post_descent_profile
    report.post_descent_support_summary == expected.post_descent_support_summary ||
        return :wrong_support_summary
    report.post_descent_valuation_summary == expected.post_descent_valuation_summary ||
        return :wrong_valuation_summary
    report.post_descent_leading_monomial_summary ==
        expected.post_descent_leading_monomial_summary ||
        return :wrong_leading_monomial_summary
    report == expected || return :stale_report
    return :ok
end

@testset "case_008 d=14 Laurent post-descent profile report" begin
    fixture = ToricBuilderCase008D14ColumnBoundary.boundary_fixture()
    @test _case008_d14_post_descent_fixture_validation(fixture) == :ok

    report = case008_d14_laurent_post_descent_profile_report(fixture)
    @test report.case_id == "case_008"
    @test report.dimension == 14
    @test report.source_boundary == :case008_d14_original
    @test report.ring_generators == ("u", "v")
    @test report.operation_family == :entry_addition
    @test report.operation == CASE008_D14_RECORDED_DESCENT_OPERATION
    @test report.replay_status == :ok
    @test report.measure_relation == :strict_decrease
    @test report.status == :post_descent_profile_report
    @test strictly_decreases_laurent_measure(
        report.before_measure,
        report.after_measure,
    )
    @test report.after_measure.whole_support_count == 7378
    @test report.after_measure.max_entry_terms == 3734
    @test report.after_measure.valuation_span == (97, 92)
    @test report.after_measure.leading_exponent == (49, -5)
    @test report.after_measure.leading_entry_index == 10
    @test report.post_descent_profile.entry_term_counts ==
          CASE008_D14_EXPECTED_POST_DESCENT_ENTRY_TERM_COUNTS
    @test report.post_descent_support_summary ==
          CASE008_D14_EXPECTED_POST_DESCENT_SUPPORT_SUMMARY
    @test report.post_descent_valuation_summary ==
          CASE008_D14_EXPECTED_POST_DESCENT_VALUATION_SUMMARY
    @test report.post_descent_leading_monomial_summary ==
          CASE008_D14_EXPECTED_POST_DESCENT_LEADING_MONOMIAL_SUMMARY
    @test report.post_descent_support_summary.whole_column_support_count == 7378
    @test report.post_descent_support_summary.whole_column_bounds ==
          (; min_exponents = (-48, -48), max_exponents = (49, 44))
    @test first(report.post_descent_leading_monomial_summary.candidates) ==
          (;
              entry_index = 10,
              leading_exponent = (49, -5),
              term_count = 3692,
              support_bounds =
                  (; min_exponents = (-47, -48), max_exponents = (49, 43)),
          )
    @test last(report.post_descent_leading_monomial_summary.candidates) ==
          (;
              entry_index = 8,
              leading_exponent = (46, -5),
              term_count = 3454,
              support_bounds =
                  (; min_exponents = (-47, -46), max_exponents = (46, 43)),
          )
    replay = _case008_d14_replayed_post_descent_data(fixture)
    @test report.post_descent_profile ==
          laurent_descent_step_profile(
              replay.after_column,
              fixture.ring;
              case_id = fixture.case_id,
          )
    @test report.post_descent_support_summary ==
          report.post_descent_profile.newton_support_summary
    @test report.post_descent_valuation_summary ==
          report.post_descent_profile.valuation_ranges
    @test report.post_descent_leading_monomial_summary ==
          report.post_descent_profile.leading_monomial_candidates
    @test validate_case008_d14_laurent_post_descent_profile_report(
        report,
        fixture,
    ) == :ok

    tampered_operation = merge(report.operation, (; exponent = (0, 0)))
    @test validate_case008_d14_laurent_post_descent_profile_report(
        merge(report, (; operation = tampered_operation)),
        fixture,
    ) == :stale_operation

    wrong_operation_generators =
        merge(report.operation, (; ring_generators = ("v", "u")))
    @test validate_case008_d14_laurent_post_descent_profile_report(
        merge(report, (; operation = wrong_operation_generators)),
        fixture,
    ) == :stale_operation

    @test validate_case008_d14_laurent_post_descent_profile_report(
        merge(report, (; ring_generators = ("v", "u"))),
        fixture,
    ) == :wrong_ring_generators

    stale_after = merge(
        report,
        (;
            after_measure = merge(
                report.after_measure,
                (; whole_support_count = report.after_measure.whole_support_count + 1),
            ),
        ),
    )
    @test validate_case008_d14_laurent_post_descent_profile_report(
        stale_after,
        fixture,
    ) == :stale_after_measure

    stale_profile = merge(
        report,
        (;
            post_descent_profile = merge(
                report.post_descent_profile,
                (; nonzero_entries = report.post_descent_profile.nonzero_entries + 1),
            ),
        ),
    )
    @test validate_case008_d14_laurent_post_descent_profile_report(
        stale_profile,
        fixture,
    ) == :stale_post_descent_profile

    wrong_support = merge(
        report,
        (;
            post_descent_support_summary = merge(
                report.post_descent_support_summary,
                (;
                    whole_column_support_count =
                        report.post_descent_support_summary.whole_column_support_count + 1,
                ),
            ),
        ),
    )
    @test validate_case008_d14_laurent_post_descent_profile_report(
        wrong_support,
        fixture,
    ) == :wrong_support_summary

    wrong_valuation = merge(
        report,
        (;
            post_descent_valuation_summary =
                (u = (; min = -48, max = 48), v = (; min = -48, max = 44)),
        ),
    )
    @test validate_case008_d14_laurent_post_descent_profile_report(
        wrong_valuation,
        fixture,
    ) == :wrong_valuation_summary

    wrong_leading = merge(
        report,
        (;
            post_descent_leading_monomial_summary = merge(
                report.post_descent_leading_monomial_summary,
                (; order = :lexicographic_ascending),
            ),
        ),
    )
    @test validate_case008_d14_laurent_post_descent_profile_report(
        wrong_leading,
        fixture,
    ) == :wrong_leading_monomial_summary
end
