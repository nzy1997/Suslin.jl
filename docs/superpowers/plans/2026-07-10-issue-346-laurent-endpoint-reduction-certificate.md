# Issue 346 Laurent Endpoint-Reduction Certificate Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an expert-only replayable Laurent endpoint-reduction certificate shell with synthetic replay coverage and conditional `case_008 d=14` status consumption.

**Architecture:** Keep all new helpers in `test/expert/laurent_endpoint_reduction_certificate.jl`. The certificate builder replays one endpoint operation on the supplied column, recomputes source and target endpoint metadata, and the validator recomputes the same data from `cert.operation` rather than trusting endpoint fields. The d14 summary consumes #345 and only reports `:endpoint_reduction_certificate` when the first replay-verified candidate validates through this shell.

**Tech Stack:** Julia, Suslin internal Laurent replay helpers, Oscar, `Test`.

## Global Constraints

- Do not promote helpers into internal code in this issue.
- Do not update production diagnostics.
- Do not implement Laurent normality/conjugation replay, recursive peel integration, determinant normalization, or full `case_008` support.
- Do not count a synthetic-only certificate as real d14 endpoint-reduction progress.
- The focused command must be `julia --project=. -e 'include("test/expert/laurent_endpoint_reduction_certificate.jl")'`.
- The full package command required by Agent Desk is `julia --project=. -e 'using Pkg; Pkg.test()'`.

---

### Task 1: Synthetic Certificate Shell

**Files:**
- Create: `test/expert/laurent_endpoint_reduction_certificate.jl`

**Interfaces:**
- Consumes: `Suslin._laurent_descent_measure_from_column`, `Suslin._laurent_link_endpoint_metadata`, `Suslin._replay_laurent_elementary_entry_addition`, `Suslin._laurent_descent_operation_status`, `Suslin._strictly_decreases_laurent_measure`.
- Produces:
  - `laurent_endpoint_reduction_certificate(context, endpoint_operation, column, R; source_endpoint = nothing, require_strict = true)::NamedTuple`
  - `validate_laurent_endpoint_reduction_certificate(cert, column, R)::Symbol`

- [ ] **Step 1: Write the failing synthetic tests**

Create `test/expert/laurent_endpoint_reduction_certificate.jl` with imports, field constants, helper field checks, and a `@testset "Laurent endpoint-reduction certificate shell"` that expects this API:

```julia
using Test
using Suslin
using Oscar

const LAURENT_ENDPOINT_REDUCTION_CERTIFICATE_STATUS =
    :endpoint_reduction_certificate
const LAURENT_ENDPOINT_REDUCTION_CERTIFICATE_NEXT_BOUNDARY =
    :laurent_normality_replay
const LAURENT_ENDPOINT_REDUCTION_OPERATION_FIELDS = (
    :family,
    :endpoint_index,
    :operation,
    :ring_generators,
)
const LAURENT_ENDPOINT_REDUCTION_CERTIFICATE_FIELDS = (
    :case_id,
    :dimension,
    :ring_generators,
    :context_status,
    :operation,
    :source_endpoint,
    :target_endpoint,
    :replay_status,
    :identity_status,
    :next_boundary,
    :status,
)

function _laurent_endpoint_reduction_has_fields(value, fields)::Bool
    return all(field -> hasproperty(value, field), fields)
end

function _laurent_endpoint_reduction_without_field(value::NamedTuple, field::Symbol)
    kept = tuple((name for name in keys(value) if name != field)...)
    return NamedTuple{kept}(tuple((getproperty(value, name) for name in kept)...))
end

@testset "Laurent endpoint-reduction certificate shell" begin
    R, (u, v) = Suslin.suslin_laurent_polynomial_ring(GF(2), ["u", "v"])
    column = [one(R), u + one(R)]
    context = (;
        case_id = "synthetic",
        dimension = 2,
        ring_generators = ("u", "v"),
        status = :endpoint_reduction_context,
    )
    operation = (;
        family = :paired_laurent_endpoint_entry_addition,
        endpoint_index = 2,
        operation = (;
            family = :entry_addition,
            target_index = 2,
            source_index = 1,
            coefficient = 1,
            exponent = (1, 0),
            ring_generators = ("u", "v"),
        ),
        ring_generators = ("u", "v"),
    )

    source_endpoint = _laurent_endpoint_metadata_from_column(
        column,
        R,
        2;
        case_id = "synthetic",
    )
    cert = laurent_endpoint_reduction_certificate(
        context,
        operation,
        column,
        R;
        source_endpoint,
    )

    @test cert.case_id == "synthetic"
    @test cert.dimension == 2
    @test cert.ring_generators == ("u", "v")
    @test cert.context_status == :endpoint_reduction_context
    @test cert.operation == operation
    @test cert.source_endpoint.leading_exponent == (1, 0)
    @test cert.target_endpoint.leading_exponent == (0, 0)
    @test cert.replay_status == :ok
    @test cert.identity_status == :verified
    @test cert.next_boundary == :laurent_normality_replay
    @test cert.status == :endpoint_reduction_certificate
    @test validate_laurent_endpoint_reduction_certificate(cert, column, R) == :ok
end
```

- [ ] **Step 2: Run the synthetic test to verify RED**

Run: `julia --project=. -e 'include("test/expert/laurent_endpoint_reduction_certificate.jl")'`

Expected: FAIL because `_laurent_endpoint_metadata_from_column` or `laurent_endpoint_reduction_certificate` is not defined.

- [ ] **Step 3: Implement the minimal certificate shell**

Add these helper functions above the testset in `test/expert/laurent_endpoint_reduction_certificate.jl`:

```julia
function _laurent_endpoint_metadata_from_column(column, R, endpoint_index::Int; case_id)
    measure = Suslin._laurent_descent_measure_from_column(column, R; case_id)
    return Suslin._laurent_link_endpoint_metadata(
        column[endpoint_index],
        R,
        endpoint_index,
        measure;
        case_id,
    )
end

function _laurent_endpoint_reduction_operation_status(endpoint_operation, n::Int, R)::Symbol
    _laurent_endpoint_reduction_has_fields(
        endpoint_operation,
        LAURENT_ENDPOINT_REDUCTION_OPERATION_FIELDS,
    ) || return :malformed_endpoint_operation
    endpoint_operation.family in (
        :laurent_endpoint_entry_addition,
        :paired_laurent_endpoint_entry_addition,
    ) || return :malformed_endpoint_operation
    endpoint_operation.ring_generators == Suslin._laurent_descent_ring_generators(R) ||
        return :wrong_ring_generators
    _laurent_endpoint_reduction_has_fields(
        endpoint_operation.operation,
        (:family, :target_index, :source_index, :coefficient, :exponent, :ring_generators),
    ) || return :malformed_endpoint_operation
    endpoint_operation.operation.target_index == endpoint_operation.endpoint_index ||
        return :wrong_endpoint_index
    endpoint_operation.operation.ring_generators == endpoint_operation.ring_generators ||
        return :wrong_ring_generators
    Suslin._laurent_descent_operation_status(endpoint_operation.operation, n, R) ==
        :ok || return :malformed_endpoint_operation
    1 <= Int(endpoint_operation.endpoint_index) <= n || return :wrong_endpoint_index
    return :ok
end

function _laurent_endpoint_reduction_replay(context, endpoint_operation, column, R; require_strict::Bool)
    status = _laurent_endpoint_reduction_operation_status(
        endpoint_operation,
        length(column),
        R,
    )
    status == :ok ||
        throw(ArgumentError("invalid Laurent endpoint operation: $(status)"))
    endpoint_index = Int(endpoint_operation.endpoint_index)
    source_endpoint = _laurent_endpoint_metadata_from_column(
        column,
        R,
        endpoint_index;
        case_id = context.case_id,
    )
    replayed_column = Suslin._replay_laurent_elementary_entry_addition(
        column,
        R,
        endpoint_operation.operation,
    )
    target_endpoint = _laurent_endpoint_metadata_from_column(
        replayed_column,
        R,
        endpoint_index;
        case_id = context.case_id,
    )
    relation = Suslin._strictly_decreases_laurent_measure(
        source_endpoint.column_measure,
        target_endpoint.column_measure,
    ) ? :strict_decrease : :not_strict_decrease
    require_strict && relation != :strict_decrease &&
        throw(ArgumentError("endpoint operation does not strictly decrease the endpoint measure"))
    return (; endpoint_operation, source_endpoint, target_endpoint, relation)
end

function laurent_endpoint_reduction_certificate(
    context,
    endpoint_operation,
    column,
    R;
    source_endpoint = nothing,
    require_strict::Bool = true,
)
    _laurent_endpoint_reduction_has_fields(
        context,
        (:case_id, :dimension, :ring_generators, :status),
    ) || throw(ArgumentError("context must include case_id, dimension, ring_generators, and status"))
    context.status == :endpoint_reduction_context ||
        throw(ArgumentError("context must have status :endpoint_reduction_context"))
    context.dimension == length(column) ||
        throw(ArgumentError("context dimension does not match source column"))
    context.ring_generators == Suslin._laurent_descent_ring_generators(R) ||
        throw(ArgumentError("context ring_generators do not match the ring"))
    replay = _laurent_endpoint_reduction_replay(
        context,
        endpoint_operation,
        column,
        R;
        require_strict,
    )
    source_endpoint === nothing || source_endpoint == replay.source_endpoint ||
        throw(ArgumentError("source endpoint is stale for the input column"))
    return (;
        case_id = context.case_id,
        dimension = context.dimension,
        ring_generators = context.ring_generators,
        context_status = context.status,
        operation = replay.endpoint_operation,
        source_endpoint = replay.source_endpoint,
        target_endpoint = replay.target_endpoint,
        endpoint_measure_relation = replay.relation,
        replay_status = :ok,
        identity_status = :verified,
        next_boundary = LAURENT_ENDPOINT_REDUCTION_CERTIFICATE_NEXT_BOUNDARY,
        status = LAURENT_ENDPOINT_REDUCTION_CERTIFICATE_STATUS,
    )
end
```

- [ ] **Step 4: Implement the validator and negative controls**

Add `validate_laurent_endpoint_reduction_certificate(cert, column, R)::Symbol` below the builder and extend the synthetic testset with controls for tampered operation, swapped ring generators, stale source metadata, stale target metadata, malformed fields, wrong endpoint index, and non-strict replay:

```julia
function validate_laurent_endpoint_reduction_certificate(cert, column, R)::Symbol
    try
        _laurent_endpoint_reduction_has_fields(
            cert,
            LAURENT_ENDPOINT_REDUCTION_CERTIFICATE_FIELDS,
        ) || return :missing_certificate_fields
        cert.status == LAURENT_ENDPOINT_REDUCTION_CERTIFICATE_STATUS ||
            return :wrong_status
        cert.replay_status == :ok || return :wrong_replay_status
        cert.identity_status == :verified || return :wrong_identity_status
        cert.next_boundary == LAURENT_ENDPOINT_REDUCTION_CERTIFICATE_NEXT_BOUNDARY ||
            return :wrong_next_boundary
        cert.dimension == length(column) || return :wrong_dimension
        cert.ring_generators == Suslin._laurent_descent_ring_generators(R) ||
            return :wrong_ring_generators
        cert.context_status == :endpoint_reduction_context ||
            return :wrong_context_status

        context = (;
            case_id = cert.case_id,
            dimension = cert.dimension,
            ring_generators = cert.ring_generators,
            status = cert.context_status,
        )
        operation_status = _laurent_endpoint_reduction_operation_status(
            cert.operation,
            length(column),
            R,
        )
        operation_status == :ok || return operation_status
        replay = _laurent_endpoint_reduction_replay(
            context,
            cert.operation,
            column,
            R;
            require_strict = false,
        )
        replay.source_endpoint == cert.source_endpoint ||
            return :stale_source_endpoint
        replay.target_endpoint == cert.target_endpoint ||
            return :stale_target_endpoint
        hasproperty(cert, :endpoint_measure_relation) ||
            return :missing_certificate_fields
        cert.endpoint_measure_relation == replay.relation ||
            return :wrong_endpoint_measure_relation
        replay.relation == :strict_decrease || return :not_endpoint_reduction
        expected = laurent_endpoint_reduction_certificate(
            context,
            cert.operation,
            column,
            R;
            source_endpoint = cert.source_endpoint,
        )
        cert == expected || return :stale_certificate
        return :ok
    catch err
        err isa InterruptException && rethrow()
        return :operation_replay_failed
    end
end
```

- [ ] **Step 5: Run the synthetic test to verify GREEN**

Run: `julia --project=. -e 'include("test/expert/laurent_endpoint_reduction_certificate.jl")'`

Expected: PASS for the synthetic certificate and all negative controls.

- [ ] **Step 6: Commit Task 1**

Run:

```bash
git add test/expert/laurent_endpoint_reduction_certificate.jl
git commit -m "Add Laurent endpoint certificate shell"
```

### Task 2: d14 Search Report Consumption and Test Registration

**Files:**
- Modify: `test/expert/laurent_endpoint_reduction_certificate.jl`
- Modify: `test/runtests.jl`

**Interfaces:**
- Consumes: `case008_d14_laurent_endpoint_reduction_search_report(fixture)`, `validate_case008_d14_laurent_endpoint_reduction_search_report(report, fixture)`, `_case008_d14_endpoint_reduction_replay_source()`, `_case008_d14_endpoint_reduction_columns(replay_source)`, and `case008_d14_laurent_endpoint_reduction_context(fixture; replay_source)`.
- Produces: `case008_d14_laurent_endpoint_reduction_certificate_summary(fixture)::NamedTuple`.

- [ ] **Step 1: Write the failing d14 summary tests**

At the top of `test/expert/laurent_endpoint_reduction_certificate.jl`, include the #345 search report when missing:

```julia
if !(@isdefined case008_d14_laurent_endpoint_reduction_search_report)
    include(joinpath(@__DIR__, "case008_d14_laurent_endpoint_reduction_search.jl"))
end
```

Add a `@testset "case_008 d=14 Laurent endpoint-reduction certificate summary"` that calls `case008_d14_laurent_endpoint_reduction_certificate_summary()` and asserts:

```julia
summary = case008_d14_laurent_endpoint_reduction_certificate_summary()
@test summary.case_id == "case_008"
@test summary.dimension == 14
@test summary.search_status in (:candidate_found, :exhausted)
@test summary.d14_status in (
    :endpoint_reduction_certificate,
    :endpoint_reduction_search_expansion,
)

if summary.search_status == :candidate_found
    @test summary.d14_status == :endpoint_reduction_certificate
    @test summary.search_next_boundary == :laurent_endpoint_reduction_certificate
    @test summary.candidate_count > 0
    @test summary.replay_verified_count == summary.candidate_count
    @test summary.certificate.case_id == "case_008"
    @test summary.certificate.status == :endpoint_reduction_certificate
    @test summary.certificate.next_boundary == :laurent_normality_replay
    @test validate_laurent_endpoint_reduction_certificate(
        summary.certificate,
        summary.source_column,
        summary.ring,
    ) == :ok
else
    @test summary.d14_status == :endpoint_reduction_search_expansion
    @test summary.search_next_boundary ==
          :laurent_endpoint_reduction_search_expansion
    @test summary.candidate_count == 0
    @test summary.replay_verified_count == 0
    @test !hasproperty(summary, :certificate)
end
```

- [ ] **Step 2: Run the focused test to verify RED**

Run: `julia --project=. -e 'include("test/expert/laurent_endpoint_reduction_certificate.jl")'`

Expected: FAIL because `case008_d14_laurent_endpoint_reduction_certificate_summary` is not defined.

- [ ] **Step 3: Implement d14 summary consumption**

Add this helper before the d14 testset:

```julia
const CASE008_D14_ENDPOINT_REDUCTION_CERTIFICATE_SUMMARY_CACHE = Ref{Any}(nothing)

function case008_d14_laurent_endpoint_reduction_certificate_summary(
    fixture = ToricBuilderCase008D14ColumnBoundary.boundary_fixture(),
)
    cached = CASE008_D14_ENDPOINT_REDUCTION_CERTIFICATE_SUMMARY_CACHE[]
    cached !== nothing && cached.fixture == fixture && return cached.summary

    report = case008_d14_laurent_endpoint_reduction_search_report(fixture)
    report_validation =
        validate_case008_d14_laurent_endpoint_reduction_search_report(
            report,
            fixture,
        )
    report_validation == :ok ||
        throw(ArgumentError("d14 endpoint-reduction search report must validate; got $(report_validation)"))

    replay_source = _case008_d14_endpoint_reduction_replay_source(fixture)
    columns = _case008_d14_endpoint_reduction_columns(replay_source)
    common = (;
        case_id = report.case_id,
        dimension = report.dimension,
        search_status = report.status,
        search_next_boundary = report.next_boundary,
        candidate_count = report.candidate_count,
        replay_verified_count = report.replay_verified_count,
        source_column = columns.source_column,
        ring = columns.ring,
    )

    if report.status == :candidate_found
        context = case008_d14_laurent_endpoint_reduction_context(
            fixture;
            replay_source,
        )
        candidate = first(report.candidates)
        cert = laurent_endpoint_reduction_certificate(
            context,
            candidate.endpoint_operation,
            columns.source_column,
            columns.ring;
            source_endpoint = report.source_endpoint,
        )
        validate_laurent_endpoint_reduction_certificate(
            cert,
            columns.source_column,
            columns.ring,
        ) == :ok ||
            throw(ArgumentError("first d14 endpoint-reduction candidate did not validate through the certificate shell"))
        cert.target_endpoint == candidate.source_endpoint ||
            throw(ArgumentError("d14 certificate target endpoint must match the replayed source-side candidate endpoint"))
        summary = merge(
            common,
            (;
                d14_status = :endpoint_reduction_certificate,
                certificate = cert,
            ),
        )
        CASE008_D14_ENDPOINT_REDUCTION_CERTIFICATE_SUMMARY_CACHE[] = (;
            fixture,
            summary,
        )
        return summary
    end

    summary = merge(
        common,
        (; d14_status = :endpoint_reduction_search_expansion),
    )
    CASE008_D14_ENDPOINT_REDUCTION_CERTIFICATE_SUMMARY_CACHE[] = (;
        fixture,
        summary,
    )
    return summary
end
```

- [ ] **Step 4: Register the expert test**

Add `"expert/laurent_endpoint_reduction_certificate.jl"` to the expert list in
`test/runtests.jl` immediately after
`"expert/case008_d14_laurent_endpoint_reduction_search.jl"`.

- [ ] **Step 5: Run focused verification**

Run: `julia --project=. -e 'include("test/expert/laurent_endpoint_reduction_certificate.jl")'`

Expected: PASS. The d14 branch must assert certificate progress only for `:candidate_found`; exhausted reports remain at search expansion.

- [ ] **Step 6: Run full package verification**

Run: `julia --project=. -e 'using Pkg; Pkg.test()'`

Expected: PASS.

- [ ] **Step 7: Commit Task 2**

Run:

```bash
git add test/expert/laurent_endpoint_reduction_certificate.jl test/runtests.jl docs/superpowers/plans/2026-07-10-issue-346-laurent-endpoint-reduction-certificate.md
git commit -m "Add d14 Laurent endpoint certificate summary"
```

## Self-Review

Spec coverage: Task 1 covers the synthetic certificate object, replay validator,
and negative controls. Task 2 covers conditional #345 consumption, d14 summary
status, test registration, focused verification, and full package verification.

Placeholder scan: no unresolved placeholder language or incomplete steps are
intended in this plan.

Type consistency: all produced helper names match the issue interface and the
test-only naming pattern used by the #323 and #335 certificate shells.

## Execution Choice

Plan complete and saved to
`docs/superpowers/plans/2026-07-10-issue-346-laurent-endpoint-reduction-certificate.md`.

1. Subagent-Driven (recommended) - dispatch a fresh subagent per task and review between tasks.
2. Inline Execution - execute tasks in this session using executing-plans with checkpoints.

Non-interactive Agent Desk decision: choose Subagent-Driven because it is marked
recommended and the user requested the standing answer policy.
