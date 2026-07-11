using Test

if !(@isdefined case008_d14_laurent_descent_profile)
    include(joinpath(@__DIR__, "case008_d14_laurent_descent_profile.jl"))
end

const CASE008_D14_MEASURE_COMPONENTS = Suslin._LAURENT_DESCENT_MEASURE_COMPONENTS

const CASE008_D14_REQUIRED_MEASURE_FIELDS = (
    :status,
    :order,
    :components,
    CASE008_D14_MEASURE_COMPONENTS...,
)

function _case008_d14_validated_measure_profile(profile, fixture)
    validation = validate_laurent_descent_profile(profile, fixture)
    validation == :ok ||
        throw(ArgumentError("invalid case_008 d14 Laurent descent profile: $(validation)"))
    return profile
end

function case008_d14_laurent_descent_measure(
    profile;
    fixture = ToricBuilderCase008D14ColumnBoundary.boundary_fixture(),
)
    checked = _case008_d14_validated_measure_profile(profile, fixture)
    return Suslin._laurent_descent_measure_from_column(
        fixture.failing_column,
        fixture.ring;
        case_id = checked.case_id,
    )
end

strictly_decreases_laurent_measure(before, after) =
    Suslin._strictly_decreases_laurent_measure(before, after)

@testset "case_008 d=14 Laurent descent measure contract" begin
    fixture = ToricBuilderCase008D14ColumnBoundary.boundary_fixture()
    profile = case008_d14_laurent_descent_profile(fixture)

    measure = case008_d14_laurent_descent_measure(profile; fixture)
    @test measure.case_id == "case_008"
    @test measure.dimension == 14
    @test measure.ring_generators == ("u", "v")
    @test measure.status == :measure_contract
    @test measure.order == :lexicographic_minimize
    @test measure.components == CASE008_D14_MEASURE_COMPONENTS
    @test measure.whole_support_count == 7387
    @test measure.max_entry_terms == 3734
    @test measure.valuation_span == (97, 93)
    @test measure.leading_exponent == (49, -5)
    @test measure.leading_entry_index == 10

    smaller = merge(
        measure,
        (; whole_support_count = measure.whole_support_count - 1),
    )
    equal = merge(measure, (;))
    larger_terms = merge(measure, (; max_entry_terms = measure.max_entry_terms + 1))
    larger_span = merge(measure, (; valuation_span = (98, 93)))
    smaller_terms = merge(measure, (; max_entry_terms = measure.max_entry_terms - 1))
    smaller_span = merge(measure, (; valuation_span = (97, 92)))
    smaller_leading = merge(measure, (; leading_exponent = (49, -6)))
    smaller_tie_breaker =
        merge(measure, (; leading_entry_index = measure.leading_entry_index - 1))
    reordered_components = (
        :leading_entry_index,
        :whole_support_count,
        :max_entry_terms,
        :valuation_span,
        :leading_exponent,
    )
    reordered_before = merge(measure, (; components = reordered_components))
    reordered_after = merge(
        measure,
        (;
            components = reordered_components,
            whole_support_count = measure.whole_support_count + 1,
            leading_entry_index = measure.leading_entry_index - 1,
        ),
    )

    @test strictly_decreases_laurent_measure(measure, smaller)
    @test strictly_decreases_laurent_measure(measure, smaller_terms)
    @test strictly_decreases_laurent_measure(measure, smaller_span)
    @test strictly_decreases_laurent_measure(measure, smaller_leading)
    @test strictly_decreases_laurent_measure(measure, smaller_tie_breaker)
    @test !strictly_decreases_laurent_measure(measure, equal)
    @test !strictly_decreases_laurent_measure(measure, larger_terms)
    @test !strictly_decreases_laurent_measure(measure, larger_span)
    @test !strictly_decreases_laurent_measure(reordered_before, reordered_after)

    supported_status = merge(profile, (; status = :supported))
    @test_throws ArgumentError case008_d14_laurent_descent_measure(
        supported_status;
        fixture,
    )

    swapped_generators = merge(profile, (; ring_generators = ("v", "u")))
    @test_throws ArgumentError case008_d14_laurent_descent_measure(
        swapped_generators;
        fixture,
    )

    stale_support = merge(
        profile,
        (;
            newton_support_summary = merge(
                profile.newton_support_summary,
                (;
                    whole_column_support_count =
                        profile.newton_support_summary.whole_column_support_count - 1,
                ),
            ),
        ),
    )
    @test_throws ArgumentError case008_d14_laurent_descent_measure(
        stale_support;
        fixture,
    )

    tampered_leading = merge(
        profile,
        (;
            leading_monomial_candidates = merge(
                profile.leading_monomial_candidates,
                (;
                    candidates = (
                        merge(
                            first(profile.leading_monomial_candidates.candidates),
                            (; leading_exponent = (48, -5)),
                        ),
                        Base.tail(profile.leading_monomial_candidates.candidates)...,
                    ),
                ),
            ),
        ),
    )
    @test_throws ArgumentError case008_d14_laurent_descent_measure(
        tampered_leading;
        fixture,
    )
end
