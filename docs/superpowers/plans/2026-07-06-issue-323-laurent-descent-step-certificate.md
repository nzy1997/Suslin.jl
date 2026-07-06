# Issue 323 Laurent Descent-Step Certificate Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an expert-only replayable Laurent descent-step certificate shell.

**Architecture:** Create one expert test/helper file that reuses the #322 operation replay and #321 measure contract, constructs certificate named tuples from replayed columns, and validates certificates by recomputing before/after profiles and measures. Register the expert test; do not touch production reducer code.

**Tech Stack:** Julia, Test stdlib, Suslin/Oscar test fixtures, existing d14 Laurent profile, measure, and search-report helpers.

## Global Constraints

- Keep the certificate shell expert-only and test-only.
- Do not modify production code under `src/`.
- Do not implement Laurent link witnesses, endpoint reductions, normality/conjugation replay, determinant normalization, or recursive peel integration.
- Certificate status must be `:descent_step_certificate`.
- Certificate replay status must be `:ok`.
- Certificate measure relation must be `:strict_decrease`.
- Stable certificate fields must include `case_id`, `dimension`, `ring_generators`, `operation`, `before_measure`, `after_measure`, `status`, `replay_status`, and `measure_relation`.
- The certificate may also include `before_profile` and `after_profile` so stale profile payloads can be rejected.
- The operation schema must match #322: `family = :entry_addition`, `target_index`, `source_index`, `coefficient`, `exponent = (a, b)`, and `ring_generators`.
- Operation replay semantics must be `target_entry <- target_entry + coefficient * u^a * v^b * source_entry`.
- The validator must be `validate_laurent_descent_step_certificate(cert, column, R)::Symbol`.
- The validator must return `:ok` only after replaying the operation and recomputing the after-measure from the replayed column.
- The validator must reject a tampered operation, stale after-profile, equal-measure step, wrong-ring certificate, malformed operation schema, and wrong certificate status.
- Include at least one known-good synthetic Laurent certificate.
- Include a real `case_008 d=14` certificate constructor because the merged #322 search report found replay-verified candidates.
- Register `expert/laurent_descent_step_certificate.jl` in `test/runtests.jl` immediately after `expert/case008_d14_laurent_elementary_move_search.jl`.
- Required focused verification: `julia --project=. -e 'include("test/expert/laurent_descent_step_certificate.jl")'`.
- Required package verification: `julia --project=. -e 'using Pkg; Pkg.test()'`.

---

## File Structure

- Create: `test/expert/laurent_descent_step_certificate.jl`
  - Owns the certificate profile builder, certificate constructor, validator, synthetic controls, and `case008_d14_laurent_descent_step_certificate()`.
- Modify: `test/runtests.jl`
  - Registers the certificate expert test after the #322 bounded move search test.

### Task 1: Replayable Laurent Descent-Step Certificate Shell

**Files:**
- Create: `test/expert/laurent_descent_step_certificate.jl`
- Modify: `test/runtests.jl`

**Interfaces:**
- Consumes: `case008_d14_laurent_descent_profile(fixture) -> NamedTuple`.
- Consumes: `strictly_decreases_laurent_measure(before, after)::Bool`.
- Consumes: `replay_laurent_elementary_entry_addition(column, R, operation) -> Vector`.
- Consumes: `_case008_d14_measure_from_column(column, R; case_id) -> NamedTuple`.
- Produces: `laurent_descent_step_profile(column, R; case_id) -> NamedTuple`.
- Produces: `laurent_descent_step_certificate(column, R, before_profile, operation) -> NamedTuple`.
- Produces: `case008_d14_laurent_descent_step_certificate(fixture = ToricBuilderCase008D14ColumnBoundary.boundary_fixture()) -> NamedTuple`.
- Produces: `validate_laurent_descent_step_certificate(cert, column, R)::Symbol`.

- [ ] **Step 1: Write the failing certificate tests and runner registration**

Create `test/expert/laurent_descent_step_certificate.jl` with this initial test-focused shape:

```julia
using Test
using Suslin
using Oscar

if !(@isdefined case008_d14_laurent_elementary_move_search_report)
    include(joinpath(@__DIR__, "case008_d14_laurent_elementary_move_search.jl"))
end

@testset "Laurent descent-step certificate shell" begin
    runtests = read(joinpath(@__DIR__, "..", "runtests.jl"), String)
    @test occursin(
        "\"expert/laurent_descent_step_certificate.jl\"",
        runtests,
    )

    R, (u, v) = Suslin.suslin_laurent_polynomial_ring(GF(2), ["u", "v"])
    column = [u + v, u + v + one(R)]
    before_profile = laurent_descent_step_profile(
        column,
        R;
        case_id = "synthetic",
    )
    operation = (;
        family = :entry_addition,
        target_index = 2,
        source_index = 1,
        coefficient = 1,
        exponent = (0, 0),
        ring_generators = ("u", "v"),
    )

    cert = laurent_descent_step_certificate(
        column,
        R,
        before_profile,
        operation,
    )

    @test cert.case_id == "synthetic"
    @test cert.dimension == 2
    @test cert.ring_generators == ("u", "v")
    @test cert.operation == operation
    @test cert.status == :descent_step_certificate
    @test cert.replay_status == :ok
    @test cert.measure_relation == :strict_decrease
    @test strictly_decreases_laurent_measure(
        cert.before_measure,
        cert.after_measure,
    )
    @test validate_laurent_descent_step_certificate(cert, column, R) == :ok

    tampered_operation = merge(operation, (; exponent = (1, 0)))
    tampered_cert = merge(cert, (; operation = tampered_operation))
    @test validate_laurent_descent_step_certificate(
        tampered_cert,
        column,
        R,
    ) == :stale_after_profile

    stale_after_profile = merge(
        cert.after_profile,
        (; nonzero_entries = cert.after_profile.nonzero_entries + 1),
    )
    stale_after_cert = merge(cert, (; after_profile = stale_after_profile))
    @test validate_laurent_descent_step_certificate(
        stale_after_cert,
        column,
        R,
    ) == :stale_after_profile

    equal_operation = merge(operation, (; coefficient = 0))
    equal_cert = _laurent_descent_step_certificate_from_replay(
        column,
        R,
        before_profile,
        equal_operation;
        require_strict = false,
    )
    equal_claim = merge(equal_cert, (; measure_relation = :strict_decrease))
    @test validate_laurent_descent_step_certificate(
        equal_claim,
        column,
        R,
    ) == :not_strict_decrease

    wrong_ring_cert = merge(cert, (; ring_generators = ("v", "u")))
    @test validate_laurent_descent_step_certificate(
        wrong_ring_cert,
        column,
        R,
    ) == :wrong_ring_generators

    malformed_operation = (;
        family = :entry_addition,
        target_index = 2,
        source_index = 1,
        coefficient = 1,
        ring_generators = ("u", "v"),
    )
    malformed_cert = merge(cert, (; operation = malformed_operation))
    @test validate_laurent_descent_step_certificate(
        malformed_cert,
        column,
        R,
    ) == :malformed_operation

    wrong_status = merge(cert, (; status = :profile_only))
    @test validate_laurent_descent_step_certificate(
        wrong_status,
        column,
        R,
    ) == :wrong_status
end

@testset "case_008 d=14 Laurent descent-step certificate" begin
    fixture = ToricBuilderCase008D14ColumnBoundary.boundary_fixture()
    cert = case008_d14_laurent_descent_step_certificate(fixture)
    @test cert.case_id == "case_008"
    @test cert.dimension == 14
    @test cert.ring_generators == ("u", "v")
    @test cert.operation.family == :entry_addition
    @test cert.operation.target_index == 1
    @test cert.operation.source_index == 2
    @test cert.operation.coefficient == 1
    @test cert.operation.exponent == (-1, 1)
    @test cert.operation.ring_generators == ("u", "v")
    @test cert.status == :descent_step_certificate
    @test cert.replay_status == :ok
    @test cert.measure_relation == :strict_decrease
    @test validate_laurent_descent_step_certificate(
        cert,
        fixture.failing_column,
        fixture.ring,
    ) == :ok
end
```

Modify `test/runtests.jl` by adding:

```julia
"expert/laurent_descent_step_certificate.jl",
```

immediately after:

```julia
"expert/case008_d14_laurent_elementary_move_search.jl",
```

Run:

```bash
julia --project=. -e 'include("test/expert/laurent_descent_step_certificate.jl")'
```

Expected: fail before implementation with `UndefVarError: laurent_descent_step_profile not defined` or an equivalent missing certificate helper.

- [ ] **Step 2: Implement the certificate helpers**

Add these constants and helper functions above the testsets in `test/expert/laurent_descent_step_certificate.jl`:

```julia
const LAURENT_DESCENT_STEP_CERTIFICATE_FIELDS = (
    :case_id,
    :dimension,
    :ring_generators,
    :operation,
    :before_profile,
    :after_profile,
    :before_measure,
    :after_measure,
    :status,
    :replay_status,
    :measure_relation,
)

const LAURENT_DESCENT_STEP_OPERATION_FIELDS = (
    :family,
    :target_index,
    :source_index,
    :coefficient,
    :exponent,
    :ring_generators,
)

function laurent_descent_step_profile(column, R; case_id)
    generator_names = _ring_generator_names(R)
    support_summary = _newton_support_summary(column, generator_names)
    term_counts = support_summary.per_entry_support_counts
    return (;
        case_id,
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

function _laurent_descent_step_operation_status(operation, n::Int, R)::Symbol
    _has_required_fields(operation, LAURENT_DESCENT_STEP_OPERATION_FIELDS) ||
        return :malformed_operation
    operation.family == :entry_addition || return :malformed_operation
    operation.ring_generators == _ring_generator_names(R) ||
        return :wrong_ring_generators
    try
        target = _checked_entry_index(operation.target_index, n, "target_index")
        source = _checked_entry_index(operation.source_index, n, "source_index")
        target != source || return :malformed_operation
        _checked_exponent_pair(operation.exponent)
        R(operation.coefficient)
    catch err
        err isa InterruptException && rethrow()
        return :malformed_operation
    end
    return :ok
end

function _laurent_descent_step_certificate_from_replay(
    column,
    R,
    before_profile,
    operation;
    require_strict::Bool = true,
)
    hasproperty(before_profile, :case_id) ||
        throw(ArgumentError("before_profile must include case_id"))
    operation_status =
        _laurent_descent_step_operation_status(operation, length(column), R)
    operation_status == :ok ||
        throw(ArgumentError("invalid Laurent descent operation: $(operation_status)"))

    expected_before_profile = laurent_descent_step_profile(
        column,
        R;
        case_id = before_profile.case_id,
    )
    before_profile == expected_before_profile ||
        throw(ArgumentError("before_profile is stale for the input column"))
    before_measure = _case008_d14_measure_from_column(
        column,
        R;
        case_id = before_profile.case_id,
    )
    after_column = replay_laurent_elementary_entry_addition(
        column,
        R,
        operation,
    )
    after_profile = laurent_descent_step_profile(
        after_column,
        R;
        case_id = before_profile.case_id,
    )
    after_measure = _case008_d14_measure_from_column(
        after_column,
        R;
        case_id = before_profile.case_id,
    )
    relation = strictly_decreases_laurent_measure(
        before_measure,
        after_measure,
    ) ? :strict_decrease : :not_strict_decrease
    require_strict && relation == :not_strict_decrease &&
        throw(ArgumentError("operation does not strictly decrease the measure"))
    return (;
        case_id = before_profile.case_id,
        dimension = length(column),
        ring_generators = _ring_generator_names(R),
        operation,
        before_profile,
        after_profile,
        before_measure,
        after_measure,
        status = :descent_step_certificate,
        replay_status = :ok,
        measure_relation = relation,
    )
end

function laurent_descent_step_certificate(column, R, before_profile, operation)
    return _laurent_descent_step_certificate_from_replay(
        column,
        R,
        before_profile,
        operation;
        require_strict = true,
    )
end

const CASE008_D14_RECORDED_DESCENT_OPERATION = (;
    family = :entry_addition,
    target_index = 1,
    source_index = 2,
    coefficient = 1,
    exponent = (-1, 1),
    ring_generators = ("u", "v"),
)

function case008_d14_laurent_descent_step_certificate(
    fixture = ToricBuilderCase008D14ColumnBoundary.boundary_fixture(),
)
    profile = case008_d14_laurent_descent_profile(fixture)
    return laurent_descent_step_certificate(
        fixture.failing_column,
        fixture.ring,
        profile,
        CASE008_D14_RECORDED_DESCENT_OPERATION,
    )
end
```

- [ ] **Step 3: Implement the validator**

Add this validator above the testsets:

```julia
function validate_laurent_descent_step_certificate(cert, column, R)::Symbol
    try
        _has_required_fields(cert, LAURENT_DESCENT_STEP_CERTIFICATE_FIELDS) ||
            return :missing_certificate_fields
        cert.status == :descent_step_certificate || return :wrong_status
        cert.replay_status == :ok || return :wrong_replay_status
        cert.measure_relation == :strict_decrease ||
            return :wrong_measure_relation
        cert.dimension == length(column) || return :wrong_dimension

        ring_generators = _ring_generator_names(R)
        cert.ring_generators == ring_generators ||
            return :wrong_ring_generators

        operation_status = _laurent_descent_step_operation_status(
            cert.operation,
            length(column),
            R,
        )
        operation_status == :ok || return operation_status

        hasproperty(cert.before_profile, :case_id) ||
            return :missing_profile_fields
        cert.case_id == cert.before_profile.case_id || return :wrong_case

        expected_before_profile = laurent_descent_step_profile(
            column,
            R;
            case_id = cert.case_id,
        )
        cert.before_profile == expected_before_profile ||
            return :stale_before_profile

        expected_before_measure = _case008_d14_measure_from_column(
            column,
            R;
            case_id = cert.case_id,
        )
        cert.before_measure == expected_before_measure ||
            return :stale_before_measure

        after_column = replay_laurent_elementary_entry_addition(
            column,
            R,
            cert.operation,
        )
        expected_after_profile = laurent_descent_step_profile(
            after_column,
            R;
            case_id = cert.case_id,
        )
        cert.after_profile == expected_after_profile ||
            return :stale_after_profile

        expected_after_measure = _case008_d14_measure_from_column(
            after_column,
            R;
            case_id = cert.case_id,
        )
        cert.after_measure == expected_after_measure ||
            return :stale_after_measure

        strictly_decreases_laurent_measure(
            expected_before_measure,
            expected_after_measure,
        ) || return :not_strict_decrease
        return :ok
    catch err
        err isa InterruptException && rethrow()
        return :operation_replay_failed
    end
end
```

- [ ] **Step 4: Run focused red/green verification**

Run:

```bash
julia --project=. -e 'include("test/expert/laurent_descent_step_certificate.jl")'
```

Expected: pass. The output must include:

- `Laurent descent-step certificate shell`;
- `case_008 d=14 Laurent descent-step certificate`;
- all negative controls passing.

- [ ] **Step 5: Run broader verification and commit**

Run:

```bash
git diff --check
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: both commands exit 0. Do not commit `Manifest.toml` if it is generated only as a local verification artifact and remains ignored by git.

Commit:

```bash
git add docs/superpowers/plans/2026-07-06-issue-323-laurent-descent-step-certificate.md test/expert/laurent_descent_step_certificate.jl test/runtests.jl
git commit -m "test: add Laurent descent-step certificate"
```

## Self-Review

- Spec coverage: task covers stable certificate fields, operation schema, ring-generator metadata, replay-based validation, stale after-profile rejection, equal-measure rejection, malformed operation rejection, wrong-ring rejection, wrong-status rejection, synthetic good case, d14 candidate case, runner registration, focused verification, and package verification.
- Placeholder scan: no placeholder markers, deferred implementation, or unresolved decision remains.
- Type consistency: produced helper names match the task interfaces and validator return values used by tests.
