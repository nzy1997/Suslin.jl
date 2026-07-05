# Issue 317 Case008 D14 Laurent Descent Profile Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an expert-only deterministic Laurent descent profile for the checked-in `case_008 d=14` boundary fixture.

**Architecture:** Keep all profile logic in one expert test file so it remains test-only and cannot be mistaken for production reducer support. The profile helper computes Newton support, valuation ranges, leading-monomial candidates, and term-count summaries directly from the #315 fixture, then validates a supplied profile by recomputing the full profile from the fixture. Register the expert file in `test/runtests.jl` for static expert-suite coverage.

**Tech Stack:** Julia, Oscar Laurent polynomial elements, Suslin fixture helpers, Test stdlib.

## Global Constraints

- Do not modify production reducer code under `src/`.
- Do not implement Laurent ECP.
- Do not choose the final Laurent descent measure.
- Do not add Laurent link witnesses, endpoint reductions, Laurent normality/conjugation replay, or recursive Laurent peel integration.
- Do not make diagonal monomial balancing or polynomialization the primary route.
- Compute the profile from `test/fixtures/toricbuilder_case008_d14_column_boundary.jl`.
- Do not require a local ToricBuilder checkout.
- The profile must record `case_id = "case_008"`.
- The profile must record `dimension = 14`.
- The profile must record `ring_generators = ("u", "v")`.
- The profile must record `nonzero_entries = 14`.
- The profile must record `max_entry_terms = 3734`.
- The profile must include valuation ranges, Newton support summary, leading monomial candidates, candidate measure families, and `status = :profile_only`.
- Candidate measure families must include `:newton_support`, `:valuation`, and `:leading_monomial`.
- Leading candidate ordering metadata must include generator order, order, and tie-breaker.
- Tampering with the fixture column, support summary, term-count summary, ring-generator metadata, or status must make validation fail.
- A profile with `status = :supported` must be rejected.
- Required focused verification: `julia --project=. -e 'include("test/expert/case008_d14_laurent_descent_profile.jl")'`.
- Required package verification: `julia --project=. -e 'using Pkg; Pkg.test()'`.

---

## File Structure

- Create: `test/expert/case008_d14_laurent_descent_profile.jl`
  - Owns all profile construction, validation, positive checks, and negative controls.
- Modify: `test/runtests.jl`
  - Registers the new expert file near the existing Laurent/case008 expert checks.

### Task 1: Add the Expert Descent Profile

**Files:**
- Create: `test/expert/case008_d14_laurent_descent_profile.jl`
- Modify: `test/runtests.jl`

**Interfaces:**
- Consumes: `ToricBuilderCase008D14ColumnBoundary.boundary_fixture() -> NamedTuple`.
- Consumes: `ToricBuilderCase008D14ColumnBoundary.validate_boundary_fixture(fixture)::Symbol`.
- Produces: `case008_d14_laurent_descent_profile(fixture = boundary_fixture()) -> NamedTuple`.
- Produces: `validate_laurent_descent_profile(profile, fixture = boundary_fixture())::Symbol`.

- [ ] **Step 1: Write the failing expert profile test**

Create `test/expert/case008_d14_laurent_descent_profile.jl` with this structure:

```julia
using Test
using Suslin
using Oscar

include(joinpath(@__DIR__, "..", "fixtures", "toricbuilder_case008_d14_column_boundary.jl"))

const CASE008_D14_PROFILE_MEASURE_FAMILIES =
    (:newton_support, :valuation, :leading_monomial)

const CASE008_D14_REQUIRED_PROFILE_FIELDS = (
    :case_id,
    :dimension,
    :ring_generators,
    :nonzero_entries,
    :max_entry_terms,
    :entry_term_counts,
    :valuation_ranges,
    :newton_support_summary,
    :leading_monomial_candidates,
    :candidate_measure_families,
    :status,
)

function _profile_has_required_fields(profile)::Bool
    return all(field -> hasproperty(profile, field), CASE008_D14_REQUIRED_PROFILE_FIELDS)
end

function _ring_generator_names(R)
    return Tuple(string.(gens(R)))
end

function _entry_support(entry)
    iszero(entry) && return ()
    support = Tuple{Vararg{Int}}[]
    for exponent in exponents(entry)
        push!(support, Tuple(Int(exponent[idx]) for idx in eachindex(exponent)))
    end
    return Tuple(sort(support))
end

function _support_bounds(support)
    isempty(support) && return nothing
    dimension = length(first(support))
    return (;
        min_exponents = ntuple(idx -> minimum(exponent[idx] for exponent in support), dimension),
        max_exponents = ntuple(idx -> maximum(exponent[idx] for exponent in support), dimension),
    )
end

function _entry_support_summary(index::Int, entry)
    support = _entry_support(entry)
    return (;
        index,
        support_count = length(support),
        term_count = length(support),
        support_bounds = _support_bounds(support),
    )
end

function _whole_column_support(column)
    support = Tuple{Vararg{Int}}[]
    for entry in column
        append!(support, collect(_entry_support(entry)))
    end
    return Tuple(sort(unique(support)))
end

function _newton_support_summary(column, generator_names)
    per_entry = Tuple(
        _entry_support_summary(idx, entry)
        for (idx, entry) in enumerate(column)
    )
    whole_support = _whole_column_support(column)
    return (;
        generator_order = generator_names,
        entry_count = length(column),
        per_entry,
        per_entry_support_counts = Tuple(summary.support_count for summary in per_entry),
        whole_column_support_count = length(whole_support),
        whole_column_bounds = _support_bounds(whole_support),
    )
end

function _valuation_ranges(whole_bounds, generator_names)
    names = Tuple(Symbol.(generator_names))
    values = whole_bounds === nothing ?
        ntuple(_ -> (; min = nothing, max = nothing), length(generator_names)) :
        ntuple(
            idx -> (;
                min = whole_bounds.min_exponents[idx],
                max = whole_bounds.max_exponents[idx],
            ),
            length(generator_names),
        )
    return NamedTuple{names}(values)
end

function _candidate_sort_lt(left, right)::Bool
    left.leading_exponent == right.leading_exponent &&
        return left.entry_index < right.entry_index
    return isless(right.leading_exponent, left.leading_exponent)
end

function _leading_monomial_candidates(column, generator_names)
    candidates = NamedTuple[]
    for (idx, entry) in enumerate(column)
        support = _entry_support(entry)
        isempty(support) && continue
        push!(
            candidates,
            (;
                entry_index = idx,
                leading_exponent = maximum(support),
                term_count = length(support),
                support_bounds = _support_bounds(support),
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

function _computed_laurent_descent_profile(fixture)
    validation = ToricBuilderCase008D14ColumnBoundary.validate_boundary_fixture(fixture)
    validation == :ok ||
        throw(ArgumentError("invalid case_008 d14 fixture: $(validation)"))

    column = fixture.failing_column
    generator_names = _ring_generator_names(fixture.ring)
    support_summary = _newton_support_summary(column, generator_names)
    term_counts = support_summary.per_entry_support_counts

    return (;
        case_id = fixture.case_id,
        dimension = length(column),
        ring_generators = generator_names,
        nonzero_entries = count(!iszero, column),
        max_entry_terms = maximum(term_counts; init = 0),
        entry_term_counts = term_counts,
        valuation_ranges = _valuation_ranges(
            support_summary.whole_column_bounds,
            generator_names,
        ),
        newton_support_summary = support_summary,
        leading_monomial_candidates = _leading_monomial_candidates(
            column,
            generator_names,
        ),
        candidate_measure_families = CASE008_D14_PROFILE_MEASURE_FAMILIES,
        status = :profile_only,
    )
end

function case008_d14_laurent_descent_profile(
    fixture = ToricBuilderCase008D14ColumnBoundary.boundary_fixture(),
)
    return _computed_laurent_descent_profile(fixture)
end

function validate_laurent_descent_profile(
    profile,
    fixture = ToricBuilderCase008D14ColumnBoundary.boundary_fixture(),
)::Symbol
    fixture_validation =
        ToricBuilderCase008D14ColumnBoundary.validate_boundary_fixture(fixture)
    fixture_validation == :ok || return :invalid_fixture
    _profile_has_required_fields(profile) || return :missing_profile_fields
    profile.status == :profile_only || return :wrong_status
    profile.case_id == "case_008" || return :wrong_case
    profile.dimension == 14 || return :wrong_dimension
    profile.ring_generators == ("u", "v") || return :wrong_ring_generators
    profile.candidate_measure_families == CASE008_D14_PROFILE_MEASURE_FAMILIES ||
        return :wrong_candidate_measure_families

    expected = _computed_laurent_descent_profile(fixture)
    profile.nonzero_entries == expected.nonzero_entries ||
        return :wrong_term_count_summary
    profile.max_entry_terms == expected.max_entry_terms ||
        return :wrong_term_count_summary
    profile.entry_term_counts == expected.entry_term_counts ||
        return :wrong_term_count_summary
    profile.valuation_ranges == expected.valuation_ranges ||
        return :wrong_valuation_summary
    profile.newton_support_summary == expected.newton_support_summary ||
        return :wrong_support_summary
    profile.leading_monomial_candidates == expected.leading_monomial_candidates ||
        return :wrong_leading_monomial_candidates
    profile == expected || return :stale_profile
    return :ok
end

@testset "case_008 d=14 Laurent descent profile" begin
    runtests = read(joinpath(@__DIR__, "..", "runtests.jl"), String)
    @test occursin(
        "\"expert/case008_d14_laurent_descent_profile.jl\"",
        runtests,
    )

    fixture = ToricBuilderCase008D14ColumnBoundary.boundary_fixture()
    @test ToricBuilderCase008D14ColumnBoundary.validate_boundary_fixture(fixture) == :ok

    profile = case008_d14_laurent_descent_profile(fixture)
    @test profile.case_id == "case_008"
    @test profile.dimension == 14
    @test profile.ring_generators == ("u", "v")
    @test profile.nonzero_entries == 14
    @test profile.max_entry_terms == 3734
    @test profile.candidate_measure_families ==
          (:newton_support, :valuation, :leading_monomial)
    @test (:newton_support in profile.candidate_measure_families)
    @test (:valuation in profile.candidate_measure_families)
    @test (:leading_monomial in profile.candidate_measure_families)
    @test profile.status == :profile_only
    @test length(profile.entry_term_counts) == 14
    @test maximum(profile.entry_term_counts) == 3734
    @test count(>(0), profile.entry_term_counts) == 14

    support = profile.newton_support_summary
    @test support.generator_order == ("u", "v")
    @test support.entry_count == 14
    @test length(support.per_entry) == 14
    @test support.per_entry_support_counts == profile.entry_term_counts
    @test support.whole_column_support_count <= sum(profile.entry_term_counts)
    @test support.whole_column_bounds !== nothing
    @test length(support.whole_column_bounds.min_exponents) == 2
    @test length(support.whole_column_bounds.max_exponents) == 2

    @test keys(profile.valuation_ranges) == (:u, :v)
    @test profile.valuation_ranges.u.min == support.whole_column_bounds.min_exponents[1]
    @test profile.valuation_ranges.u.max == support.whole_column_bounds.max_exponents[1]
    @test profile.valuation_ranges.v.min == support.whole_column_bounds.min_exponents[2]
    @test profile.valuation_ranges.v.max == support.whole_column_bounds.max_exponents[2]

    leading = profile.leading_monomial_candidates
    @test leading.generator_order == ("u", "v")
    @test leading.order == :lexicographic_descending
    @test leading.tie_breaker == :entry_index_ascending
    @test leading.candidate_count == 14
    @test length(leading.candidates) == 14
    @test all(candidate -> candidate.term_count > 0, leading.candidates)
    for idx in 1:(length(leading.candidates) - 1)
        @test _candidate_sort_lt(leading.candidates[idx], leading.candidates[idx + 1])
    end
    @test validate_laurent_descent_profile(profile, fixture) == :ok

    corrupted_fixture =
        ToricBuilderCase008D14ColumnBoundary.corrupted_column_negative_control(fixture)
    @test validate_laurent_descent_profile(profile, corrupted_fixture) ==
          :invalid_fixture

    wrong_ring_generators = merge(profile, (; ring_generators = ("x", "y")))
    @test validate_laurent_descent_profile(wrong_ring_generators, fixture) ==
          :wrong_ring_generators

    wrong_terms = merge(profile, (; max_entry_terms = profile.max_entry_terms + 1))
    @test validate_laurent_descent_profile(wrong_terms, fixture) ==
          :wrong_term_count_summary

    wrong_support = merge(
        profile,
        (;
            newton_support_summary = merge(
                profile.newton_support_summary,
                (;
                    whole_column_support_count =
                        profile.newton_support_summary.whole_column_support_count + 1,
                ),
            ),
        ),
    )
    @test validate_laurent_descent_profile(wrong_support, fixture) ==
          :wrong_support_summary

    wrong_leading = merge(
        profile,
        (;
            leading_monomial_candidates = merge(
                profile.leading_monomial_candidates,
                (; order = :lexicographic_ascending),
            ),
        ),
    )
    @test validate_laurent_descent_profile(wrong_leading, fixture) ==
          :wrong_leading_monomial_candidates

    wrong_measures = merge(
        profile,
        (; candidate_measure_families = (:newton_support, :valuation)),
    )
    @test validate_laurent_descent_profile(wrong_measures, fixture) ==
          :wrong_candidate_measure_families

    supported_status = merge(profile, (; status = :supported))
    @test validate_laurent_descent_profile(supported_status, fixture) ==
          :wrong_status
end
```

- [ ] **Step 2: Register the expert test**

In `test/runtests.jl`, add the file next to the Laurent-native boundary diagnostic and case008 Laurent expert checks:

```julia
        "expert/laurent_column_reduction_diagnostics.jl",
        "expert/laurent_native_ecp_boundary_diagnostics.jl",
        "expert/case008_d14_laurent_descent_profile.jl",
        "expert/case008_d21_laurent_column_reduction.jl",
```

- [ ] **Step 3: Run the focused verification**

Run:

```bash
julia --project=. -e 'include("test/expert/case008_d14_laurent_descent_profile.jl")'
```

Expected: exits 0 and the new testset passes. If the registration assertion fails, Step 2 was missed.

- [ ] **Step 4: Run expert and package verification**

Run:

```bash
julia --project=. test/runtests.jl expert
julia --project=. -e 'using Pkg; Pkg.test()'
git diff --check
```

Expected: all commands exit 0. `Pkg.test()` should remain the default public/internal package suite; the new expert file is covered by the explicit expert command.

- [ ] **Step 5: Commit the implementation**

Run:

```bash
git add test/expert/case008_d14_laurent_descent_profile.jl test/runtests.jl docs/superpowers/plans/2026-07-05-issue-317-case008-d14-laurent-descent-profile.md
git commit -m "test: record case008 d14 Laurent descent profile"
```

Expected: a new implementation commit is created on the worker branch.

## Plan Self-Review

- Spec coverage: the task creates the expert-only profile, records all requested stable fields, includes candidate ordering metadata, validates by recomputation, and covers every requested negative control.
- Placeholder scan: no placeholder work remains; each step has exact files and commands.
- Type consistency: the produced helper names and return types match the test and validation steps.
