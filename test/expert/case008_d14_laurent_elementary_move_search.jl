using Test
using Suslin
using Oscar

if !(@isdefined case008_d14_laurent_descent_measure)
    include(joinpath(@__DIR__, "case008_d14_laurent_descent_measure_contract.jl"))
end

const CASE008_D14_ELEMENTARY_OPERATION_FAMILIES = (:entry_addition,)
const CASE008_D14_ELEMENTARY_EXPONENT_RADIUS = 1
const CASE008_D14_ELEMENTARY_COEFFICIENT_FAMILY = (1,)

const CASE008_D14_ELEMENTARY_OPERATION_FIELDS = (
    :family,
    :target_index,
    :source_index,
    :coefficient,
    :exponent,
)

const CASE008_D14_ELEMENTARY_CANDIDATE_FIELDS = (
    :operation,
    :before_measure,
    :after_measure,
)

function _checked_entry_index(index, n::Int, name::AbstractString)::Int
    index isa Integer ||
        throw(ArgumentError("$(name) must be an integer entry index"))
    checked = Int(index)
    1 <= checked <= n ||
        throw(ArgumentError("$(name) must be between 1 and $(n)"))
    return checked
end

function _checked_exponent_pair(exponent)::Tuple{Int, Int}
    (exponent isa Tuple || exponent isa AbstractVector) ||
        throw(ArgumentError("exponent must be a length-2 tuple or vector"))
    length(exponent) == 2 ||
        throw(ArgumentError("exponent must have exactly two entries"))
    exponent[1] isa Integer ||
        throw(ArgumentError("first exponent must be an integer"))
    exponent[2] isa Integer ||
        throw(ArgumentError("second exponent must be an integer"))
    return (Int(exponent[1]), Int(exponent[2]))
end

function _has_required_fields(value, fields)::Bool
    return all(field -> hasproperty(value, field), fields)
end

function replay_laurent_elementary_entry_addition(column, R, operation)
    operation.family == :entry_addition ||
        throw(ArgumentError("unsupported operation family $(repr(operation.family))"))
    n = length(column)
    target = _checked_entry_index(operation.target_index, n, "target_index")
    source = _checked_entry_index(operation.source_index, n, "source_index")
    target != source ||
        throw(ArgumentError("target_index and source_index must differ"))
    exponent = _checked_exponent_pair(operation.exponent)
    coefficient = R(operation.coefficient)
    generators = gens(R)
    monomial = coefficient * generators[1]^exponent[1] * generators[2]^exponent[2]
    transformed = copy(column)
    transformed[target] = transformed[target] + monomial * column[source]
    return transformed
end

function _case008_d14_entry_support_set(entry)::Set{Tuple{Int, Int}}
    support = Set{Tuple{Int, Int}}()
    for exponent in _entry_support(entry)
        push!(support, _checked_exponent_pair(exponent))
    end
    return support
end

function _case008_d14_column_support_sets(column)::Vector{Set{Tuple{Int, Int}}}
    return [_case008_d14_entry_support_set(entry) for entry in column]
end

function _case008_d14_shifted_support(
    support::Set{Tuple{Int, Int}},
    exponent::Tuple{Int, Int},
)::Set{Tuple{Int, Int}}
    return Set{Tuple{Int, Int}}(
        (term[1] + exponent[1], term[2] + exponent[2])
        for term in support
    )
end

function _case008_d14_symmetric_difference(
    left::Set{Tuple{Int, Int}},
    right::Set{Tuple{Int, Int}},
)::Set{Tuple{Int, Int}}
    result = copy(left)
    for term in right
        if term in result
            delete!(result, term)
        else
            push!(result, term)
        end
    end
    return result
end

function _case008_d14_update_entry_addition_supports(
    supports::Vector{Set{Tuple{Int, Int}}},
    R,
    operation,
)::Vector{Set{Tuple{Int, Int}}}
    operation.family == :entry_addition ||
        throw(ArgumentError("unsupported operation family $(repr(operation.family))"))
    n = length(supports)
    target = _checked_entry_index(operation.target_index, n, "target_index")
    source = _checked_entry_index(operation.source_index, n, "source_index")
    target != source ||
        throw(ArgumentError("target_index and source_index must differ"))
    exponent = _checked_exponent_pair(operation.exponent)

    updated = copy(supports)
    iszero(R(operation.coefficient)) && return updated
    shifted = _case008_d14_shifted_support(supports[source], exponent)
    updated[target] =
        _case008_d14_symmetric_difference(supports[target], shifted)
    return updated
end

function _case008_d14_support_bounds(support::Set{Tuple{Int, Int}})
    isempty(support) && return nothing
    return (;
        min_exponents = (
            minimum(term[1] for term in support),
            minimum(term[2] for term in support),
        ),
        max_exponents = (
            maximum(term[1] for term in support),
            maximum(term[2] for term in support),
        ),
    )
end

function _case008_d14_whole_support(
    supports::Vector{Set{Tuple{Int, Int}}},
)::Set{Tuple{Int, Int}}
    whole_support = Set{Tuple{Int, Int}}()
    for support in supports
        union!(whole_support, support)
    end
    return whole_support
end

function _case008_d14_valuation_span(whole_support::Set{Tuple{Int, Int}})
    bounds = _case008_d14_support_bounds(whole_support)
    bounds === nothing && return (0, 0)
    return (
        bounds.max_exponents[1] - bounds.min_exponents[1],
        bounds.max_exponents[2] - bounds.min_exponents[2],
    )
end

function _case008_d14_leading_monomial_candidates_from_supports(
    supports::Vector{Set{Tuple{Int, Int}}},
    generator_names,
)
    candidates = NamedTuple[]
    for (idx, support) in enumerate(supports)
        isempty(support) && continue
        push!(
            candidates,
            (;
                entry_index = idx,
                leading_exponent = maximum(support),
                term_count = length(support),
                support_bounds = _case008_d14_support_bounds(support),
            ),
        )
    end
    ordered = sort(candidates; lt = _candidate_sort_lt)
    return (;
        generator_order = generator_names,
        order = :lexicographic_descending,
        tie_breaker = :entry_index_ascending,
        candidate_count = length(ordered),
        candidates = Tuple(ordered),
    )
end

function _case008_d14_measure_from_supports(
    supports::Vector{Set{Tuple{Int, Int}}},
    generator_names;
    case_id = "case_008",
)
    whole_support = _case008_d14_whole_support(supports)
    leading = _case008_d14_leading_monomial_candidates_from_supports(
        supports,
        generator_names,
    )
    isempty(leading.candidates) &&
        throw(ArgumentError("cannot measure a column with no nonzero entries"))
    leading_first = first(leading.candidates)
    return (;
        case_id,
        dimension = length(supports),
        ring_generators = Tuple(generator_names),
        status = :measure_contract,
        order = :lexicographic_minimize,
        components = CASE008_D14_MEASURE_COMPONENTS,
        whole_support_count = length(whole_support),
        max_entry_terms = maximum(length, supports; init = 0),
        valuation_span = _case008_d14_valuation_span(whole_support),
        leading_exponent = leading_first.leading_exponent,
        leading_entry_index = leading_first.entry_index,
    )
end

function _case008_d14_measure_from_column(
    column,
    R;
    case_id = "case_008",
)
    supports = _case008_d14_column_support_sets(column)
    return _case008_d14_measure_from_supports(
        supports,
        _ring_generator_names(R);
        case_id,
    )
end

function _case008_d14_candidate_case_id(candidate)
    hasproperty(candidate.before_measure, :case_id) ||
        throw(ArgumentError("candidate before_measure must include case_id"))
    return candidate.before_measure.case_id
end

function verify_laurent_elementary_move_candidate(
    original_column,
    R,
    candidate,
)::Bool
    try
        _has_required_fields(candidate, CASE008_D14_ELEMENTARY_CANDIDATE_FIELDS) ||
            return false
        operation = candidate.operation
        _has_required_fields(operation, CASE008_D14_ELEMENTARY_OPERATION_FIELDS) ||
            return false
        if hasproperty(operation, :ring_generators)
            operation.ring_generators == _ring_generator_names(R) || return false
        end

        case_id = _case008_d14_candidate_case_id(candidate)
        before_measure = _case008_d14_measure_from_column(
            original_column,
            R;
            case_id,
        )
        transformed = replay_laurent_elementary_entry_addition(
            original_column,
            R,
            operation,
        )
        after_measure = _case008_d14_measure_from_column(
            transformed,
            R;
            case_id,
        )

        candidate.before_measure == before_measure || return false
        candidate.after_measure == after_measure || return false
        return strictly_decreases_laurent_measure(before_measure, after_measure)
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end

function _case008_d14_elementary_operation(
    target::Int,
    source::Int,
    exponent::Tuple{Int, Int},
    coefficient,
    generator_names,
)
    return (;
        family = :entry_addition,
        target_index = target,
        source_index = source,
        coefficient,
        exponent,
        ring_generators = Tuple(generator_names),
    )
end

function _case008_d14_search_bounds(n::Int)
    radius = CASE008_D14_ELEMENTARY_EXPONENT_RADIUS
    ordered_entry_pairs = Tuple(
        (; target_index = target, source_index = source)
        for target in 1:n
        for source in 1:n
        if target != source
    )
    exponents = Tuple(
        (a, b)
        for a in -radius:radius
        for b in -radius:radius
    )
    return (;
        target_indices = Tuple(1:n),
        source_indices = Tuple(1:n),
        ordered_entry_pairs,
        ordered_pair_constraint = :target_index_not_equal_source_index,
        ordered_pair_count = length(ordered_entry_pairs),
        exponent_radius = radius,
        exponent_vectors = exponents,
        coefficient_family = CASE008_D14_ELEMENTARY_COEFFICIENT_FAMILY,
        coefficient_ring = :GF2,
        scan_order = (
            :target_index_ascending,
            :source_index_ascending,
            :exponent_a_ascending,
            :exponent_b_ascending,
            :coefficient_family_order,
        ),
    )
end

function case008_d14_laurent_elementary_move_search_report(
    fixture = ToricBuilderCase008D14ColumnBoundary.boundary_fixture(),
)::NamedTuple
    validation = ToricBuilderCase008D14ColumnBoundary.validate_boundary_fixture(fixture)
    validation == :ok ||
        throw(ArgumentError("invalid case_008 d14 fixture: $(validation)"))

    column = fixture.failing_column
    R = fixture.ring
    generator_names = _ring_generator_names(R)
    profile = case008_d14_laurent_descent_profile(fixture)
    baseline_measure = case008_d14_laurent_descent_measure(profile; fixture)
    supports = _case008_d14_column_support_sets(column)
    support_measure = _case008_d14_measure_from_supports(
        supports,
        generator_names;
        case_id = fixture.case_id,
    )
    support_measure == baseline_measure ||
        throw(ArgumentError("support-set baseline measure does not match contract"))

    candidates = NamedTuple[]
    checked_operation_count = 0
    bounds = _case008_d14_search_bounds(length(column))

    for pair in bounds.ordered_entry_pairs
        target = pair.target_index
        source = pair.source_index
        for exponent in bounds.exponent_vectors
            for coefficient in bounds.coefficient_family
                checked_operation_count += 1
                operation = _case008_d14_elementary_operation(
                    target,
                    source,
                    exponent,
                    coefficient,
                    generator_names,
                )
                after_supports =
                    _case008_d14_update_entry_addition_supports(
                        supports,
                        R,
                        operation,
                    )
                after_measure = _case008_d14_measure_from_supports(
                    after_supports,
                    generator_names;
                    case_id = fixture.case_id,
                )
                strictly_decreases_laurent_measure(
                    baseline_measure,
                    after_measure,
                ) || continue

                candidate = (;
                    operation,
                    before_measure = baseline_measure,
                    after_measure,
                )
                verify_laurent_elementary_move_candidate(
                    column,
                    R,
                    candidate,
                ) || continue
                push!(candidates, candidate)
            end
        end
    end

    return (;
        case_id = fixture.case_id,
        dimension = length(column),
        input_measure = baseline_measure,
        operation_families = CASE008_D14_ELEMENTARY_OPERATION_FAMILIES,
        search_bounds = bounds,
        status = isempty(candidates) ? :exhausted : :candidate_found,
        candidate_count = length(candidates),
        checked_operation_count,
        replay_verified_count = length(candidates),
        candidates = Tuple(candidates),
    )
end

@testset "case_008 d=14 Laurent elementary move search" begin
    runtests = read(joinpath(@__DIR__, "..", "runtests.jl"), String)
    @test occursin(
        "\"expert/case008_d14_laurent_elementary_move_search.jl\"",
        runtests,
    )

    fixture = ToricBuilderCase008D14ColumnBoundary.boundary_fixture()
    baseline_profile = case008_d14_laurent_descent_profile(fixture)
    baseline_measure = case008_d14_laurent_descent_measure(
        baseline_profile;
        fixture,
    )

    report = case008_d14_laurent_elementary_move_search_report()
    @test report.case_id == "case_008"
    @test report.dimension == 14
    @test report.input_measure == baseline_measure
    @test report.operation_families == (:entry_addition,)
    @test report.search_bounds.exponent_radius == 1
    @test report.search_bounds.coefficient_family == (1,)
    @test report.search_bounds.ordered_pair_constraint ==
          :target_index_not_equal_source_index
    @test length(report.search_bounds.ordered_entry_pairs) == 14 * 13
    @test all(
        pair -> pair.target_index != pair.source_index,
        report.search_bounds.ordered_entry_pairs,
    )
    @test first(report.search_bounds.ordered_entry_pairs) ==
          (; target_index = 1, source_index = 2)
    @test last(report.search_bounds.ordered_entry_pairs) ==
          (; target_index = 14, source_index = 13)
    @test report.checked_operation_count == 1638
    @test report.status in (:candidate_found, :exhausted)
    @test report.candidate_count == length(report.candidates)
    @test report.replay_verified_count == report.candidate_count
    @test all(
        candidate -> verify_laurent_elementary_move_candidate(
            fixture.failing_column,
            fixture.ring,
            candidate,
        ),
        report.candidates,
    )

    if report.status == :exhausted
        @test report.candidate_count == 0
        @test report.replay_verified_count == 0
        @test isempty(report.candidates)
    else
        @test report.candidate_count > 0
        @test all(
            candidate -> candidate.before_measure == baseline_measure,
            report.candidates,
        )
        @test all(
            candidate -> strictly_decreases_laurent_measure(
                candidate.before_measure,
                candidate.after_measure,
            ),
            report.candidates,
        )
    end
end

@testset "Laurent elementary move candidate verification controls" begin
    R, (u, v) = Suslin.suslin_laurent_polynomial_ring(GF(2), ["u", "v"])
    column = [u + v, u + v + one(R)]

    before_measure = _case008_d14_measure_from_column(
        column,
        R;
        case_id = "synthetic",
    )
    decreasing_operation = (;
        family = :entry_addition,
        target_index = 2,
        source_index = 1,
        coefficient = 1,
        exponent = (0, 0),
        ring_generators = ("u", "v"),
    )
    decreasing_column = replay_laurent_elementary_entry_addition(
        column,
        R,
        decreasing_operation,
    )
    decreasing_after_measure = _case008_d14_measure_from_column(
        decreasing_column,
        R;
        case_id = "synthetic",
    )
    good_candidate = (;
        operation = decreasing_operation,
        before_measure,
        after_measure = decreasing_after_measure,
    )
    @test verify_laurent_elementary_move_candidate(column, R, good_candidate)

    nondecreasing_operation = merge(
        decreasing_operation,
        (; target_index = 1, source_index = 2),
    )
    nondecreasing_column = replay_laurent_elementary_entry_addition(
        column,
        R,
        nondecreasing_operation,
    )
    nondecreasing_candidate = (;
        operation = nondecreasing_operation,
        before_measure,
        after_measure = _case008_d14_measure_from_column(
            nondecreasing_column,
            R;
            case_id = "synthetic",
        ),
    )
    @test !verify_laurent_elementary_move_candidate(
        column,
        R,
        nondecreasing_candidate,
    )

    malformed_candidate = merge(
        good_candidate,
        (;
            operation = merge(decreasing_operation, (; source_index = 3)),
        ),
    )
    @test !verify_laurent_elementary_move_candidate(column, R, malformed_candidate)

    stale_after_candidate = merge(
        good_candidate,
        (;
            after_measure = merge(
                decreasing_after_measure,
                (; max_entry_terms = decreasing_after_measure.max_entry_terms - 1),
            ),
        ),
    )
    @test !verify_laurent_elementary_move_candidate(column, R, stale_after_candidate)
end
