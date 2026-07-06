# Issue 321 Laurent Descent Measure Contract Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Define an expert-only replayable Laurent descent measure contract for the checked-in `case_008 d=14` profile.

**Architecture:** Add a new expert test file that includes the #317 profile helper, validates profiles through `validate_laurent_descent_profile`, derives a stable named-tuple measure, and exposes a strict lexicographic-minimize comparator. Register the new expert file in `test/runtests.jl`; do not touch production reducer code.

**Tech Stack:** Julia, Test stdlib, existing Oscar-backed test fixture/profile helpers.

## Global Constraints

- Keep the contract expert-only and test-only.
- Do not modify production code under `src/`.
- Do not implement Laurent link witnesses, endpoint reductions, recursive peel integration, or production support for `case_008 d=14`.
- Do not make diagonal monomial balancing or polynomialization part of this primary measure contract.
- Build the measure from the validated profile returned by `case008_d14_laurent_descent_profile()`.
- The real profile must pass `validate_laurent_descent_profile(profile, fixture) == :ok` before a measure object is accepted.
- The measure must record `case_id == "case_008"`.
- The measure must record `dimension == 14`.
- The measure must record `ring_generators == ("u", "v")`.
- The measure must record `status == :measure_contract`.
- The measure must record `order == :lexicographic_minimize`.
- The measure must record `components == (:whole_support_count, :max_entry_terms, :valuation_span, :leading_exponent, :leading_entry_index)`.
- The measure baseline must record `whole_support_count == 7387`.
- The measure baseline must record `max_entry_terms == 3734`.
- The measure baseline must record `valuation_span == (97, 93)`.
- The measure baseline must record `leading_exponent == (49, -5)`.
- The measure baseline must record `leading_entry_index == 10`.
- `strictly_decreases_laurent_measure(before, after)::Bool` must return true only when `after` is strictly smaller under the declared lexicographic-minimize component order.
- Tuple-valued components are compared lexicographically in generator order.
- `leading_entry_index` is only the final tie-breaker.
- A profile with `status != :profile_only`, swapped ring generators, stale support summary, or tampered leading-monomial metadata must be rejected before a measure object is accepted.
- Required focused verification: `julia --project=. -e 'include("test/expert/case008_d14_laurent_descent_measure_contract.jl")'`.
- Required package verification: `julia --project=. -e 'using Pkg; Pkg.test()'`.

---

## File Structure

- Create: `test/expert/case008_d14_laurent_descent_measure_contract.jl`
  - Owns measure construction, comparator helpers, positive checks, synthetic comparator controls, and profile rejection controls.
- Modify: `test/runtests.jl`
  - Registers the new expert file immediately after `expert/case008_d14_laurent_descent_profile.jl`.

### Task 1: Write the Failing Expert Contract Tests

**Files:**
- Create: `test/expert/case008_d14_laurent_descent_measure_contract.jl`

**Interfaces:**
- Consumes: `case008_d14_laurent_descent_profile(fixture = boundary_fixture()) -> NamedTuple`.
- Consumes: `validate_laurent_descent_profile(profile, fixture = boundary_fixture())::Symbol`.
- Expects later task to produce: `case008_d14_laurent_descent_measure(profile; fixture = boundary_fixture()) -> NamedTuple`.
- Expects later task to produce: `strictly_decreases_laurent_measure(before, after)::Bool`.

- [ ] **Step 1: Add the failing test file**

Create `test/expert/case008_d14_laurent_descent_measure_contract.jl` with these tests and no implementation for `case008_d14_laurent_descent_measure` or `strictly_decreases_laurent_measure` yet:

```julia
using Test

if !(@isdefined case008_d14_laurent_descent_profile)
    include(joinpath(@__DIR__, "case008_d14_laurent_descent_profile.jl"))
end

const CASE008_D14_MEASURE_COMPONENTS = (
    :whole_support_count,
    :max_entry_terms,
    :valuation_span,
    :leading_exponent,
    :leading_entry_index,
)

@testset "case_008 d=14 Laurent descent measure contract" begin
    runtests = read(joinpath(@__DIR__, "..", "runtests.jl"), String)
    @test occursin(
        "\"expert/case008_d14_laurent_descent_measure_contract.jl\"",
        runtests,
    )

    fixture = ToricBuilderCase008D14ColumnBoundary.boundary_fixture()
    profile = case008_d14_laurent_descent_profile(fixture)
    @test validate_laurent_descent_profile(profile, fixture) == :ok

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

    @test strictly_decreases_laurent_measure(measure, smaller)
    @test !strictly_decreases_laurent_measure(measure, equal)
    @test !strictly_decreases_laurent_measure(measure, larger_terms)
    @test !strictly_decreases_laurent_measure(measure, larger_span)

    supported_status = merge(profile, (; status = :supported))
    @test validate_laurent_descent_profile(supported_status, fixture) == :wrong_status
    @test_throws ArgumentError case008_d14_laurent_descent_measure(
        supported_status;
        fixture,
    )

    swapped_generators = merge(profile, (; ring_generators = ("v", "u")))
    @test validate_laurent_descent_profile(swapped_generators, fixture) ==
          :wrong_ring_generators
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
    @test validate_laurent_descent_profile(stale_support, fixture) ==
          :wrong_support_summary
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
    @test validate_laurent_descent_profile(tampered_leading, fixture) ==
          :wrong_leading_monomial_candidates
    @test_throws ArgumentError case008_d14_laurent_descent_measure(
        tampered_leading;
        fixture,
    )
end
```

- [ ] **Step 2: Run the focused test and verify RED**

Run:

```bash
julia --project=. -e 'include("test/expert/case008_d14_laurent_descent_measure_contract.jl")'
```

Expected: nonzero exit. The failure must mention the missing runner registration and/or `UndefVarError: case008_d14_laurent_descent_measure not defined`, proving the new contract does not yet exist.

### Task 2: Implement the Measure Contract and Register It

**Files:**
- Modify: `test/expert/case008_d14_laurent_descent_measure_contract.jl`
- Modify: `test/runtests.jl`

**Interfaces:**
- Produces: `case008_d14_laurent_descent_measure(profile; fixture = boundary_fixture()) -> NamedTuple`.
- Produces: `strictly_decreases_laurent_measure(before, after)::Bool`.

- [ ] **Step 1: Add the minimal implementation above the testset**

Insert this code below `CASE008_D14_MEASURE_COMPONENTS`:

```julia
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

function _case008_d14_valuation_span(profile)
    names = Tuple(Symbol.(profile.ring_generators))
    return ntuple(
        idx -> begin
            range = getproperty(profile.valuation_ranges, names[idx])
            range.max - range.min
        end,
        length(names),
    )
end

function case008_d14_laurent_descent_measure(
    profile;
    fixture = ToricBuilderCase008D14ColumnBoundary.boundary_fixture(),
)
    checked = _case008_d14_validated_measure_profile(profile, fixture)
    leading = first(checked.leading_monomial_candidates.candidates)
    return (;
        case_id = checked.case_id,
        dimension = checked.dimension,
        ring_generators = checked.ring_generators,
        status = :measure_contract,
        order = :lexicographic_minimize,
        components = CASE008_D14_MEASURE_COMPONENTS,
        whole_support_count =
            checked.newton_support_summary.whole_column_support_count,
        max_entry_terms = checked.max_entry_terms,
        valuation_span = _case008_d14_valuation_span(checked),
        leading_exponent = leading.leading_exponent,
        leading_entry_index = leading.entry_index,
    )
end

function _has_laurent_measure_fields(measure)::Bool
    return all(field -> hasproperty(measure, field), CASE008_D14_REQUIRED_MEASURE_FIELDS)
end

function _laurent_measure_component_values(measure)
    return Tuple(getproperty(measure, component) for component in measure.components)
end

function strictly_decreases_laurent_measure(before, after)::Bool
    _has_laurent_measure_fields(before) || return false
    _has_laurent_measure_fields(after) || return false
    before.status == :measure_contract || return false
    after.status == :measure_contract || return false
    before.order == :lexicographic_minimize || return false
    after.order == before.order || return false
    after.components == before.components || return false
    return isless(
        _laurent_measure_component_values(after),
        _laurent_measure_component_values(before),
    )
end
```

- [ ] **Step 2: Register the new expert test**

In `test/runtests.jl`, insert:

```julia
"expert/case008_d14_laurent_descent_measure_contract.jl",
```

immediately after:

```julia
"expert/case008_d14_laurent_descent_profile.jl",
```

- [ ] **Step 3: Run focused verification and verify GREEN**

Run:

```bash
julia --project=. -e 'include("test/expert/case008_d14_laurent_descent_measure_contract.jl")'
```

Expected: exit 0.

- [ ] **Step 4: Run expert group coverage**

Run:

```bash
julia --project=. test/runtests.jl expert
```

Expected: exit 0.

- [ ] **Step 5: Run package verification**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: exit 0.

- [ ] **Step 6: Run whitespace check**

Run:

```bash
git diff --check
```

Expected: exit 0.

## Plan Self-Review

- Spec coverage: the plan covers the exact measure fields, baseline values,
  validation gate, strict comparator, negative controls, runner registration,
  and required verification commands from the design and issue body.
- Placeholder scan: no deferred implementation placeholders remain.
- Type consistency: the same function names and component names are used in
  tests, implementation, and verification steps.

## Automatic Decisions

- Execution approach: Subagent-Driven selected because it is the recommended
  writing-plans option and Agent Desk is non-interactive.
- If subagent execution cannot write usable changes in this sandbox, execute
  the same task steps inline and record that fallback in the decision log.
