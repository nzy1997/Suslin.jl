# Issue 327 Laurent Descent Internal Helpers Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Promote Laurent descent measure, strict-decrease, replay, and certificate-validation mechanics from expert tests into case-agnostic internal helpers.

**Architecture:** Add private helpers to `src/algorithm/column_reduction.jl` near the Laurent diagnostic helpers so later diagnostics can call them without public exports. Internal tests exercise the production helpers directly; expert tests keep fixture-specific profiles/search reports while delegating shared replay and measure logic to `Suslin._...` helpers.

**Tech Stack:** Julia, Oscar Laurent polynomial rings, `Test`, existing `Suslin` internals.

## Global Constraints

- Do not expose a public API or add exports.
- Production helpers must stay case-agnostic: no `case_008`, `d=14`, or baseline constants in `src`.
- Do not change `diagnose_unimodular_column_reduction` behavior except as needed to compile shared helpers.
- Validator must recompute before and after measures from the input column and replayed operation instead of trusting supplied summaries.
- The internal test must recompute the checked-in `case_008 d=14` baseline values: whole support count `7387`, maximum entry terms `3734`, valuation span `(97, 93)`, leading exponent `(49, -5)`, leading entry index `10`.
- The internal test must replay `target_index = 1`, `source_index = 2`, `coefficient = 1`, `exponent = (-1, 1)` and prove strict decrease.
- Negative controls must reject swapped ring generators, malformed source/target indices, target equal to source, stale supplied measures, a zero-coefficient non-decreasing operation, and an operation whose after-measure is not recomputed from replay.

---

### Task 1: Promote Laurent Descent Helpers And Tests

**Files:**
- Modify: `src/algorithm/column_reduction.jl`
- Create: `test/internal/laurent_descent_measure_helpers.jl`
- Modify: `test/runtests.jl`
- Modify: `test/expert/case008_d14_laurent_descent_measure_contract.jl`
- Modify: `test/expert/case008_d14_laurent_elementary_move_search.jl`
- Modify: `test/expert/laurent_descent_step_certificate.jl`

**Interfaces:**
- Produces: `Suslin._laurent_descent_measure_from_column(column, R; case_id = nothing)::NamedTuple`
- Produces: `Suslin._strictly_decreases_laurent_measure(before, after)::Bool`
- Produces: `Suslin._replay_laurent_elementary_entry_addition(column, R, operation)::Vector`
- Produces: `Suslin._validate_laurent_descent_step_certificate(cert, column, R)::Symbol`
- Consumes: Existing `_is_laurent_polynomial_ring`, `_coerce_into_ring`, `gens`, `exponents`, and Laurent polynomial arithmetic.

- [ ] **Step 1: Write the failing internal test and register it**

Create `test/internal/laurent_descent_measure_helpers.jl` with tests that call only production helpers for shared mechanics:

```julia
using Test
using Suslin
using Oscar

include(joinpath(@__DIR__, "..", "fixtures", "toricbuilder_case008_d14_column_boundary.jl"))

const INTERNAL_D14_OPERATION = (;
    family = :entry_addition,
    target_index = 1,
    source_index = 2,
    coefficient = 1,
    exponent = (-1, 1),
    ring_generators = ("u", "v"),
)

function _internal_descent_certificate(column, R, operation; case_id = "case_008")
    before = Suslin._laurent_descent_measure_from_column(column, R; case_id)
    after_column = Suslin._replay_laurent_elementary_entry_addition(column, R, operation)
    after = Suslin._laurent_descent_measure_from_column(after_column, R; case_id)
    return (;
        case_id,
        dimension = length(column),
        ring_generators = Tuple(string.(gens(R))),
        operation,
        before_measure = before,
        after_measure = after,
        status = :descent_step_certificate,
        replay_status = :ok,
        measure_relation = Suslin._strictly_decreases_laurent_measure(before, after) ?
            :strict_decrease : :not_strict_decrease,
    )
end

@testset "internal Laurent descent measure helpers" begin
    runtests = read(joinpath(@__DIR__, "..", "runtests.jl"), String)
    @test occursin("\"internal/laurent_descent_measure_helpers.jl\"", runtests)

    fixture = ToricBuilderCase008D14ColumnBoundary.boundary_fixture()
    @test ToricBuilderCase008D14ColumnBoundary.validate_boundary_fixture(fixture) == :ok

    measure = Suslin._laurent_descent_measure_from_column(
        fixture.failing_column,
        fixture.ring;
        case_id = fixture.case_id,
    )
    @test measure.case_id == "case_008"
    @test measure.dimension == 14
    @test measure.ring_generators == ("u", "v")
    @test measure.status == :measure_contract
    @test measure.order == :lexicographic_minimize
    @test measure.components == (
        :whole_support_count,
        :max_entry_terms,
        :valuation_span,
        :leading_exponent,
        :leading_entry_index,
    )
    @test measure.whole_support_count == 7387
    @test measure.max_entry_terms == 3734
    @test measure.valuation_span == (97, 93)
    @test measure.leading_exponent == (49, -5)
    @test measure.leading_entry_index == 10

    after_column = Suslin._replay_laurent_elementary_entry_addition(
        fixture.failing_column,
        fixture.ring,
        INTERNAL_D14_OPERATION,
    )
    after_measure = Suslin._laurent_descent_measure_from_column(
        after_column,
        fixture.ring;
        case_id = fixture.case_id,
    )
    @test Suslin._strictly_decreases_laurent_measure(measure, after_measure)

    cert = _internal_descent_certificate(
        fixture.failing_column,
        fixture.ring,
        INTERNAL_D14_OPERATION;
        case_id = fixture.case_id,
    )
    @test Suslin._validate_laurent_descent_step_certificate(
        cert,
        fixture.failing_column,
        fixture.ring,
    ) == :ok

    swapped_cert = merge(cert, (; ring_generators = ("v", "u")))
    @test Suslin._validate_laurent_descent_step_certificate(
        swapped_cert,
        fixture.failing_column,
        fixture.ring,
    ) == :wrong_ring_generators

    swapped_operation = merge(INTERNAL_D14_OPERATION, (; ring_generators = ("v", "u")))
    swapped_operation_cert = merge(cert, (; operation = swapped_operation))
    @test Suslin._validate_laurent_descent_step_certificate(
        swapped_operation_cert,
        fixture.failing_column,
        fixture.ring,
    ) == :wrong_ring_generators

    bad_source = merge(INTERNAL_D14_OPERATION, (; source_index = 0))
    @test_throws ArgumentError Suslin._replay_laurent_elementary_entry_addition(
        fixture.failing_column,
        fixture.ring,
        bad_source,
    )
    bad_source_cert = merge(cert, (; operation = bad_source))
    @test Suslin._validate_laurent_descent_step_certificate(
        bad_source_cert,
        fixture.failing_column,
        fixture.ring,
    ) == :malformed_operation

    equal_indices = merge(INTERNAL_D14_OPERATION, (; source_index = 1))
    @test_throws ArgumentError Suslin._replay_laurent_elementary_entry_addition(
        fixture.failing_column,
        fixture.ring,
        equal_indices,
    )
    equal_indices_cert = merge(cert, (; operation = equal_indices))
    @test Suslin._validate_laurent_descent_step_certificate(
        equal_indices_cert,
        fixture.failing_column,
        fixture.ring,
    ) == :malformed_operation

    stale_before = merge(
        cert,
        (;
            before_measure = merge(
                cert.before_measure,
                (; whole_support_count = cert.before_measure.whole_support_count - 1),
            ),
        ),
    )
    @test Suslin._validate_laurent_descent_step_certificate(
        stale_before,
        fixture.failing_column,
        fixture.ring,
    ) == :stale_before_measure

    zero_operation = merge(INTERNAL_D14_OPERATION, (; coefficient = 0))
    zero_cert = _internal_descent_certificate(
        fixture.failing_column,
        fixture.ring,
        zero_operation;
        case_id = fixture.case_id,
    )
    zero_claim = merge(zero_cert, (; measure_relation = :strict_decrease))
    @test Suslin._validate_laurent_descent_step_certificate(
        zero_claim,
        fixture.failing_column,
        fixture.ring,
    ) == :not_strict_decrease

    stale_after_operation = merge(INTERNAL_D14_OPERATION, (; exponent = (0, 0)))
    stale_after_cert = merge(cert, (; operation = stale_after_operation))
    @test Suslin._validate_laurent_descent_step_certificate(
        stale_after_cert,
        fixture.failing_column,
        fixture.ring,
    ) == :stale_after_measure
end
```

Also add `"internal/laurent_descent_measure_helpers.jl"` to the `"internal"` group in `test/runtests.jl` immediately after `"internal/toricbuilder_case008_d14_column_boundary.jl"`.

- [ ] **Step 2: Run the red test and confirm it fails for missing internals**

Run:

```bash
julia --project=. -e 'include("test/internal/laurent_descent_measure_helpers.jl")'
```

Expected: FAIL with `UndefVarError` for `Suslin._laurent_descent_measure_from_column` or the first missing promoted helper.

- [ ] **Step 3: Implement internal helpers in `src/algorithm/column_reduction.jl`**

Add private constants and helper functions near the existing Laurent diagnostic helpers:

```julia
const _LAURENT_DESCENT_MEASURE_COMPONENTS = (
    :whole_support_count,
    :max_entry_terms,
    :valuation_span,
    :leading_exponent,
    :leading_entry_index,
)

const _LAURENT_DESCENT_OPERATION_FIELDS = (
    :family,
    :target_index,
    :source_index,
    :coefficient,
    :exponent,
)

const _LAURENT_DESCENT_CERTIFICATE_FIELDS = (
    :case_id,
    :dimension,
    :ring_generators,
    :operation,
    :before_measure,
    :after_measure,
    :status,
    :replay_status,
    :measure_relation,
)
```

Implement helper bodies that:

- require `R` to be a two-generator Laurent polynomial ring;
- extract sorted exponent tuples with `exponents(entry)`;
- compute whole support, max entry term count, valuation span from whole support bounds, leading exponent by lexicographic descending order, and entry-index ascending tie-breaker;
- include `case_id` in the measure only when the keyword is not `nothing`;
- compare measure tuples using the fixed component order;
- replay `target += coefficient * prod(gens(R)[i]^exponent[i]) * source`;
- convert all malformed operation data into `ArgumentError` in replay and into stable rejection symbols in validation;
- recompute both stored measures from the original column and replayed after-column before accepting a certificate.

- [ ] **Step 4: Run the internal test and fix implementation bugs**

Run:

```bash
julia --project=. -e 'include("test/internal/laurent_descent_measure_helpers.jl")'
```

Expected: exit 0.

- [ ] **Step 5: Update expert tests to reuse internal helpers where practical**

Make these narrow substitutions:

```julia
const CASE008_D14_MEASURE_COMPONENTS = Suslin._LAURENT_DESCENT_MEASURE_COMPONENTS

function case008_d14_laurent_descent_measure(profile; fixture = ...)
    checked = _case008_d14_validated_measure_profile(profile, fixture)
    return Suslin._laurent_descent_measure_from_column(
        fixture.failing_column,
        fixture.ring;
        case_id = checked.case_id,
    )
end

strictly_decreases_laurent_measure(before, after) =
    Suslin._strictly_decreases_laurent_measure(before, after)

replay_laurent_elementary_entry_addition(column, R, operation) =
    Suslin._replay_laurent_elementary_entry_addition(column, R, operation)

_case008_d14_measure_from_column(column, R; case_id = "case_008") =
    Suslin._laurent_descent_measure_from_column(column, R; case_id)
```

In `test/expert/laurent_descent_step_certificate.jl`, delegate certificate validation to `Suslin._validate_laurent_descent_step_certificate` and keep any profile-specific checks only if they remain necessary for the existing expert shell.

- [ ] **Step 6: Run required targeted verification**

Run:

```bash
julia --project=. -e 'include("test/internal/laurent_descent_measure_helpers.jl")'
julia --project=. -e 'include("test/expert/laurent_descent_step_certificate.jl")'
```

Expected: both commands exit 0.

- [ ] **Step 7: Run broader verification**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
git diff --check
```

Expected: both commands exit 0. Julia may print pre-existing world-age warnings from fixture catalog loading; do not treat warnings as failures unless the process exits nonzero.

- [ ] **Step 8: Commit**

Run:

```bash
git add src/algorithm/column_reduction.jl test/internal/laurent_descent_measure_helpers.jl test/runtests.jl test/expert/case008_d14_laurent_descent_measure_contract.jl test/expert/case008_d14_laurent_elementary_move_search.jl test/expert/laurent_descent_step_certificate.jl docs/superpowers/plans/2026-07-07-issue-327-laurent-descent-internal-helpers.md
git commit -m "Promote Laurent descent helpers"
```

Expected: commit succeeds with only source, test, and plan files staged.

## Self Review

Spec coverage: the plan covers internal measure shape, strict-decrease predicate, replay helper, validator recomputation, internal test registration, expert reuse, and all issue verification commands. Placeholder scan: no placeholder terms are used. Type consistency: helper names and field names match the design and issue body.
