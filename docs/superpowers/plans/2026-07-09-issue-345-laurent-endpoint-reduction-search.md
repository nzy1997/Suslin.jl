# Issue 345 Laurent Endpoint-Reduction Search Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an expert-only bounded `case_008 d=14` Laurent endpoint-reduction search report that consumes the replay-derived issue 341 context and records either replay-verified candidates or an exhausted bounded search.

**Architecture:** Keep all new logic in a single expert test file. Reuse the issue 341 endpoint context and its cached replay source, recompute endpoint metadata from replayed source/target columns, and register the focused expert test in `test/runtests.jl`.

**Tech Stack:** Julia, `Test`, `Suslin`, `Oscar`, existing expert-only replay helpers.

## Global Constraints

- Do not promote helpers into `src/` or expose production endpoint-reduction diagnostics.
- The report must rebuild from `case008_d14_laurent_endpoint_reduction_context()` and the replay-certified d14 Laurent link-witness certificate path.
- Fixed fields must include `case_id == "case_008"`, `dimension == 14`, `ring_generators == ("u", "v")`, `source_boundary == :case008_d14_link_witness_certificate`, `context_status == :endpoint_reduction_context`, `boundary == :laurent_endpoint_reduction`, `pivot_index == 10`, `partner_index == 1`, and `witness_exponent == (1, -1)`.
- Bounded search parameters are endpoint indices `(10,)`, source indices `(1,)`, operation families `(:paired_laurent_endpoint_entry_addition,)`, exponent bounds `((-1, -1), (1, 1))`, all nine exponent vectors in that box, and coefficient families `(1,)`.
- If `status == :candidate_found`, require `replay_verified_count == candidate_count > 0` and `next_boundary == :laurent_endpoint_reduction_certificate`.
- If `status == :exhausted`, require `candidate_count == 0`, `replay_verified_count == 0`, and `next_boundary == :laurent_endpoint_reduction_search_expansion`.
- Negative controls must reject swapped ring generators, stale source endpoint metadata, stale target endpoint metadata, malformed endpoint operations, wrong endpoint indices, missing required endpoint-reduction fields, and candidate metadata that was copied instead of recomputed from replay.

---

## File Structure

- Create `test/expert/case008_d14_laurent_endpoint_reduction_search.jl`.
  This file owns the test-only search report, endpoint-operation replay helpers,
  report validation, candidate validation, positive contract tests, and negative
  controls.
- Modify `test/runtests.jl`.
  Add the new expert file immediately after
  `expert/case008_d14_laurent_endpoint_reduction_context.jl`.

### Task 1: Expert Search Report

**Files:**
- Create: `test/expert/case008_d14_laurent_endpoint_reduction_search.jl`
- Modify: `test/runtests.jl`
- Test: `test/expert/case008_d14_laurent_endpoint_reduction_search.jl`

**Interfaces:**
- Consumes: `case008_d14_laurent_endpoint_reduction_context()`, `_case008_d14_endpoint_reduction_replay_source()`, and `validate_case008_d14_laurent_endpoint_reduction_context(context, fixture)` from `test/expert/case008_d14_laurent_endpoint_reduction_context.jl`.
- Produces: `case008_d14_laurent_endpoint_reduction_search_report(fixture = ToricBuilderCase008D14ColumnBoundary.boundary_fixture())::NamedTuple`, `validate_case008_d14_laurent_endpoint_reduction_search_report(report, fixture = ToricBuilderCase008D14ColumnBoundary.boundary_fixture())::Symbol`, and `validate_case008_d14_laurent_endpoint_reduction_candidate(context, candidate, replay_source; require_strict = true)::Symbol`.

- [ ] **Step 1: Write the failing expert test and registration**

Create `test/expert/case008_d14_laurent_endpoint_reduction_search.jl` with the imports, constants, helper names, and testsets below. The functions may initially be stubs returning empty NamedTuples or `:invalid_report`, but the tests must reference the final public test-only names so the focused command fails for missing behavior, not syntax.

```julia
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
    :status,
)

@testset "case_008 d=14 Laurent endpoint reduction search report" begin
    runtests = read(joinpath(@__DIR__, "..", "runtests.jl"), String)
    @test occursin(
        "\"expert/case008_d14_laurent_endpoint_reduction_search.jl\"",
        runtests,
    )

    report = case008_d14_laurent_endpoint_reduction_search_report()
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
    @test validate_case008_d14_laurent_endpoint_reduction_search_report(report) == :ok

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
    context = case008_d14_laurent_endpoint_reduction_context(
        replay_source.fixture;
        replay_source,
    )
    report = case008_d14_laurent_endpoint_reduction_search_report(
        replay_source.fixture,
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
    ) == :stale_target_endpoint
    @test validate_case008_d14_laurent_endpoint_reduction_search_report(
        merge(
            report,
            (;
                required_endpoint_reduction_fields =
                    (:family, :endpoint_index, :ring_generators),
            ),
        ),
    ) == :missing_required_endpoint_reduction_field
    @test validate_case008_d14_laurent_endpoint_reduction_search_report(
        _case008_d14_endpoint_reduction_without_field(
            report,
            :required_endpoint_reduction_fields,
        ),
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
    ) == :wrong_endpoint_index

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
    ) == :stale_candidate_target_endpoint
end
```

In `test/runtests.jl`, add:

```julia
"expert/case008_d14_laurent_endpoint_reduction_search.jl",
```

immediately after:

```julia
"expert/case008_d14_laurent_endpoint_reduction_context.jl",
```

- [ ] **Step 2: Run the focused command and confirm RED**

Run:

```bash
julia --project=. -e 'include("test/expert/case008_d14_laurent_endpoint_reduction_search.jl")'
```

Expected: FAIL because `case008_d14_laurent_endpoint_reduction_search_report`
and helper functions are not implemented yet. Syntax and include errors must be
fixed before proceeding; the observed failure should be missing behavior.

- [ ] **Step 3: Implement report and validation helpers**

Implement these helpers in the same expert file, above the testsets:

```julia
function _case008_d14_endpoint_reduction_has_fields(value, fields)::Bool
    return all(field -> hasproperty(value, field), fields)
end

function _case008_d14_endpoint_reduction_without_field(value::NamedTuple, field::Symbol)
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

function _case008_d14_endpoint_reduction_operation_status(endpoint_operation, n::Int, R)::Symbol
    _case008_d14_endpoint_reduction_has_fields(
        endpoint_operation,
        CASE008_D14_ENDPOINT_REDUCTION_REQUIRED_FIELDS,
    ) || return :malformed_endpoint_operation
    endpoint_operation.family == :paired_laurent_endpoint_entry_addition ||
        return :malformed_endpoint_operation
    endpoint_operation.ring_generators == ("u", "v") ||
        return :wrong_ring_generators
    endpoint_operation.endpoint_index in CASE008_D14_ENDPOINT_REDUCTION_ENDPOINT_INDICES ||
        return :wrong_endpoint_index
    _case008_d14_endpoint_reduction_has_fields(
        endpoint_operation.operation,
        (:family, :target_index, :source_index, :coefficient, :exponent, :ring_generators),
    ) || return :malformed_endpoint_operation
    endpoint_operation.operation.target_index == endpoint_operation.endpoint_index ||
        return :wrong_endpoint_index
    endpoint_operation.operation.ring_generators == endpoint_operation.ring_generators ||
        return :wrong_ring_generators
    Suslin._laurent_descent_operation_status(endpoint_operation.operation, n, R) == :ok ||
        return :malformed_endpoint_operation
    return :ok
end

function _case008_d14_endpoint_metadata_from_column(column, R, endpoint_index::Int; case_id)
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
    endpoint_index = Int(endpoint_operation.endpoint_index)
    source_before = _case008_d14_endpoint_metadata_from_column(
        source_column,
        R,
        endpoint_index;
        case_id = context.case_id,
    )
    target_before = _case008_d14_endpoint_metadata_from_column(
        target_column,
        R,
        endpoint_index;
        case_id = context.case_id,
    )
    replayed_source_column = Suslin._replay_laurent_elementary_entry_addition(
        source_column,
        R,
        endpoint_operation.operation,
    )
    replayed_target_column = Suslin._replay_laurent_elementary_entry_addition(
        target_column,
        R,
        endpoint_operation.operation,
    )
    source_endpoint = _case008_d14_endpoint_metadata_from_column(
        replayed_source_column,
        R,
        endpoint_index;
        case_id = context.case_id,
    )
    target_endpoint = _case008_d14_endpoint_metadata_from_column(
        replayed_target_column,
        R,
        endpoint_index;
        case_id = context.case_id,
    )
    source_measure_relation =
        _case008_d14_endpoint_relation(source_before, source_endpoint)
    target_measure_relation =
        _case008_d14_endpoint_relation(target_before, target_endpoint)
    candidate_status =
        source_measure_relation == :strict_decrease &&
        target_measure_relation == :strict_decrease ?
        :strict_endpoint_decrease :
        :not_endpoint_reduction
    require_strict && candidate_status != :strict_endpoint_decrease &&
        throw(ArgumentError("endpoint operation does not strictly decrease both endpoint measures"))
    return (;
        endpoint_operation,
        source_endpoint,
        target_endpoint,
        source_measure_relation,
        target_measure_relation,
        replay_status = :ok,
        status = candidate_status,
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
            return operation_status == :ok ? :wrong_candidate_replay_status : operation_status
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
        candidate.status == expected.status || return :wrong_candidate_status
        require_strict && expected.status != :strict_endpoint_decrease &&
            return :not_endpoint_reduction
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
        throw(ArgumentError("endpoint-reduction context must validate; got $(context_validation)"))
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
                        nothing
                    end
                    candidate === nothing && continue
                    validation =
                        validate_case008_d14_laurent_endpoint_reduction_candidate(
                            context,
                            candidate,
                            replay_source,
                        )
                    validation == :ok || continue
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
        report.source_indices == CASE008_D14_ENDPOINT_REDUCTION_SOURCE_INDICES ||
            return :wrong_source_indices
        report.operation_families ==
            CASE008_D14_ENDPOINT_REDUCTION_SEARCH_OPERATION_FAMILIES ||
            return :wrong_operation_families
        report.operation_semantics ==
            CASE008_D14_ENDPOINT_REDUCTION_OPERATION_SEMANTICS ||
            return :wrong_operation_semantics
        report.exponent_bounds == CASE008_D14_ENDPOINT_REDUCTION_EXPONENT_BOUNDS ||
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
```

- [ ] **Step 4: Run focused verification**

Run:

```bash
julia --project=. -e 'include("test/expert/case008_d14_laurent_endpoint_reduction_search.jl")'
```

Expected: PASS, with the report either `:candidate_found` satisfying the
candidate branch or `:exhausted` satisfying the exhausted branch.

- [ ] **Step 5: Run package verification and commit**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: PASS for the package test suite.

Then run:

```bash
git diff --check
git status --short
```

Expected: no whitespace errors; only the new expert test, `test/runtests.jl`,
and this plan file are modified or untracked.

Commit:

```bash
git add docs/superpowers/plans/2026-07-09-issue-345-laurent-endpoint-reduction-search.md test/expert/case008_d14_laurent_endpoint_reduction_search.jl test/runtests.jl
git commit -m "Add issue 345 Laurent endpoint search"
```
