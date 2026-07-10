# Issue 347 Laurent Endpoint-Reduction Internals Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Promote generic Laurent endpoint-reduction replay, candidate, and certificate helpers from expert tests into internal algorithm code.

**Architecture:** Add underscore-prefixed helpers next to the existing Laurent descent and link-witness helpers in `src/algorithm/column_reduction.jl`. Keep `case_008 d=14` bounds and fixture constants in expert tests, and make the expert search/certificate files delegate generic replay and validation to the new internals.

**Tech Stack:** Julia, Oscar, existing Suslin internal Laurent descent/link-witness helpers, `Test`.

## Global Constraints

- Do not expose a public API.
- Do not update `diagnose_unimodular_column_reduction` except as needed to compile shared helpers.
- Keep production helper code case-agnostic; no `case_008 d=14` fixture constants in `src/algorithm/column_reduction.jl`.
- Validate endpoint-reduction certificates by replay and recomputation, not by trusting supplied endpoint fields.
- Run focused verification:
  `julia --project=. -e 'include("test/internal/laurent_endpoint_reduction_helpers.jl")'`
- Run focused verification:
  `julia --project=. -e 'include("test/expert/laurent_endpoint_reduction_certificate.jl")'`
- Run final package verification:
  `julia --project=. -e 'using Pkg; Pkg.test()'`

---

## File Structure

- Modify `src/algorithm/column_reduction.jl`: add generic endpoint-reduction constants, status, replay, candidate, certificate, and validator helpers near the Laurent link-witness helpers.
- Create `test/internal/laurent_endpoint_reduction_helpers.jl`: focused synthetic and conditional d14 coverage for the promoted internals.
- Modify `test/runtests.jl`: include the new internal test after `internal/laurent_link_witness_helpers.jl`.
- Modify `test/expert/case008_d14_laurent_endpoint_reduction_search.jl`: keep d14-specific search bounds but call the promoted internals for generic status, replay, and candidate verification.
- Modify `test/expert/laurent_endpoint_reduction_certificate.jl`: keep test-facing helper names as wrappers around promoted internals.

### Task 1: Add Internal RED Test

**Files:**
- Create: `test/internal/laurent_endpoint_reduction_helpers.jl`
- Modify: `test/runtests.jl`

**Interfaces:**
- Consumes: planned internal names from the issue: `_laurent_endpoint_reduction_status`, `_replay_laurent_endpoint_reduction`, `_laurent_endpoint_reduction_candidate_from_replay`, `_verify_laurent_endpoint_reduction_candidate`, `_laurent_endpoint_reduction_certificate_from_replay`, `_validate_laurent_endpoint_reduction_certificate`.
- Produces: a failing test command that proves those helpers are not yet promoted.

- [ ] **Step 1: Write the failing internal test**

Create `test/internal/laurent_endpoint_reduction_helpers.jl` with these test sections:

```julia
using Test
using Suslin
using Oscar

if !(@isdefined case008_d14_laurent_endpoint_reduction_search_report)
    include(joinpath(@__DIR__, "..", "expert", "case008_d14_laurent_endpoint_reduction_search.jl"))
end

function _internal_endpoint_without_field(value::NamedTuple, field::Symbol)
    kept = tuple((name for name in keys(value) if name != field)...)
    return NamedTuple{kept}(tuple((getproperty(value, name) for name in kept)...))
end

@testset "internal Laurent endpoint-reduction helpers" begin
    runtests = read(joinpath(@__DIR__, "..", "runtests.jl"), String)
    @test occursin("\"internal/laurent_endpoint_reduction_helpers.jl\"", runtests)

    R, (u, v) = Suslin.suslin_laurent_polynomial_ring(GF(2), ["u", "v"])
    source_column = [one(R), u + one(R)]
    target_column = [one(R), u + v]
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

    @test Suslin._laurent_endpoint_reduction_status(operation, length(source_column), R) == :ok
    replay = Suslin._replay_laurent_endpoint_reduction(
        source_column,
        R,
        operation;
        case_id = context.case_id,
    )
    @test replay.source_endpoint.leading_exponent == (1, 0)
    @test replay.target_endpoint.leading_exponent == (0, 0)
    @test replay.relation == :strict_decrease

    candidate = Suslin._laurent_endpoint_reduction_candidate_from_replay(
        source_column,
        target_column,
        R,
        operation;
        case_id = context.case_id,
    )
    @test candidate.source_endpoint.leading_exponent == (0, 0)
    @test candidate.target_endpoint.leading_exponent == (0, 1)
    @test candidate.source_measure_relation == :strict_decrease
    @test candidate.target_measure_relation == :strict_decrease
    @test candidate.identity_status == :verified
    @test candidate.status == :strict_endpoint_decrease
    @test Suslin._verify_laurent_endpoint_reduction_candidate(
        source_column,
        target_column,
        R,
        candidate,
    )

    cert = Suslin._laurent_endpoint_reduction_certificate_from_replay(
        context,
        operation,
        source_column,
        R;
        source_endpoint = replay.source_endpoint,
    )
    @test cert.case_id == "synthetic"
    @test cert.target_endpoint.leading_exponent == (0, 0)
    @test cert.endpoint_measure_relation == :strict_decrease
    @test cert.next_boundary == :laurent_normality_replay
    @test cert.status == :endpoint_reduction_certificate
    @test Suslin._validate_laurent_endpoint_reduction_certificate(
        cert,
        source_column,
        R,
    ) == :ok

    missing_field = _internal_endpoint_without_field(operation, :family)
    @test Suslin._laurent_endpoint_reduction_status(
        missing_field,
        length(source_column),
        R,
    ) == :malformed_endpoint_operation

    wrong_generators = merge(operation, (; ring_generators = ("v", "u")))
    @test Suslin._laurent_endpoint_reduction_status(
        wrong_generators,
        length(source_column),
        R,
    ) == :wrong_ring_generators

    wrong_endpoint_index = merge(operation, (; endpoint_index = 1))
    @test Suslin._laurent_endpoint_reduction_status(
        wrong_endpoint_index,
        length(source_column),
        R,
    ) == :wrong_endpoint_index

    malformed_nested = merge(
        operation,
        (; operation = _internal_endpoint_without_field(operation.operation, :coefficient)),
    )
    @test Suslin._laurent_endpoint_reduction_status(
        malformed_nested,
        length(source_column),
        R,
    ) == :malformed_endpoint_operation

    stale_source = merge(
        cert,
        (;
            source_endpoint = merge(
                cert.source_endpoint,
                (; term_count = cert.source_endpoint.term_count + 1),
            ),
        ),
    )
    @test Suslin._validate_laurent_endpoint_reduction_certificate(
        stale_source,
        source_column,
        R,
    ) == :stale_source_endpoint

    stale_target = merge(
        cert,
        (;
            target_endpoint = merge(
                cert.target_endpoint,
                (; term_count = cert.target_endpoint.term_count + 1),
            ),
        ),
    )
    @test Suslin._validate_laurent_endpoint_reduction_certificate(
        stale_target,
        source_column,
        R,
    ) == :stale_target_endpoint

    tampered_operation = merge(
        cert,
        (;
            operation = merge(
                operation,
                (; operation = merge(operation.operation, (; coefficient = 0))),
            ),
        ),
    )
    @test Suslin._validate_laurent_endpoint_reduction_certificate(
        tampered_operation,
        source_column,
        R,
    ) == :stale_target_endpoint

    nonstrict_operation = merge(
        operation,
        (; operation = merge(operation.operation, (; exponent = (0, -1)))),
    )
    nonstrict_cert = Suslin._laurent_endpoint_reduction_certificate_from_replay(
        context,
        nonstrict_operation,
        source_column,
        R;
        source_endpoint = replay.source_endpoint,
        require_strict = false,
    )
    @test Suslin._validate_laurent_endpoint_reduction_certificate(
        nonstrict_cert,
        source_column,
        R,
    ) == :not_endpoint_reduction
    nonstrict_candidate = Suslin._laurent_endpoint_reduction_candidate_from_replay(
        source_column,
        target_column,
        R,
        nonstrict_operation;
        case_id = context.case_id,
        require_strict = false,
    )
    @test !Suslin._verify_laurent_endpoint_reduction_candidate(
        source_column,
        target_column,
        R,
        nonstrict_candidate,
    )
end

@testset "internal d14 Laurent endpoint-reduction candidate" begin
    report = case008_d14_laurent_endpoint_reduction_search_report()
    replay_source = _case008_d14_endpoint_reduction_replay_source()
    columns = _case008_d14_endpoint_reduction_columns(replay_source)

    if report.status == :candidate_found
        candidate = first(report.candidates)
        @test Suslin._verify_laurent_endpoint_reduction_candidate(
            columns.source_column,
            columns.target_column,
            columns.ring,
            candidate,
        )
        expected = Suslin._laurent_endpoint_reduction_candidate_from_replay(
            columns.source_column,
            columns.target_column,
            columns.ring,
            candidate.endpoint_operation;
            case_id = report.case_id,
        )
        @test candidate == expected

        context = case008_d14_laurent_endpoint_reduction_context(
            replay_source.fixture;
            replay_source,
        )
        cert = Suslin._laurent_endpoint_reduction_certificate_from_replay(
            context,
            candidate.endpoint_operation,
            columns.source_column,
            columns.ring;
            source_endpoint = report.source_endpoint,
        )
        @test cert.target_endpoint == candidate.source_endpoint
        @test Suslin._validate_laurent_endpoint_reduction_certificate(
            cert,
            columns.source_column,
            columns.ring,
        ) == :ok
    else
        @test report.status == :exhausted
        @test report.candidate_count == 0
        @test report.replay_verified_count == 0
        @test report.next_boundary == :laurent_endpoint_reduction_search_expansion
    end
end
```

Add `"internal/laurent_endpoint_reduction_helpers.jl",` immediately after
`"internal/laurent_link_witness_helpers.jl",` in `test/runtests.jl`.

- [ ] **Step 2: Run the RED test**

Run:

```bash
julia --project=. -e 'include("test/internal/laurent_endpoint_reduction_helpers.jl")'
```

Expected: nonzero exit with `UndefVarError` for
`_laurent_endpoint_reduction_status` or the next missing promoted helper.

### Task 2: Promote Generic Internal Helpers

**Files:**
- Modify: `src/algorithm/column_reduction.jl`

**Interfaces:**
- Consumes: `_laurent_descent_has_fields`, `_laurent_descent_ring_generators`, `_laurent_descent_checked_entry_index`, `_laurent_descent_operation_status`, `_replay_laurent_elementary_entry_addition`, `_laurent_link_endpoint_metadata`, `_strictly_decreases_laurent_measure`.
- Produces: the six internal endpoint-reduction helper functions listed in Task 1.

- [ ] **Step 1: Add constants near the Laurent link-witness constants**

Add endpoint-reduction field/status constants after
`_LAURENT_LINK_WITNESS_CERTIFICATE_FIELDS`:

```julia
const _LAURENT_ENDPOINT_REDUCTION_OPERATION_FIELDS = (
    :family,
    :endpoint_index,
    :operation,
    :ring_generators,
)
const _LAURENT_ENDPOINT_REDUCTION_OPERATION_FAMILIES = (
    :laurent_endpoint_entry_addition,
    :paired_laurent_endpoint_entry_addition,
)
const _LAURENT_ENDPOINT_REDUCTION_CANDIDATE_FIELDS = (
    :endpoint_operation,
    :source_endpoint,
    :target_endpoint,
    :source_measure_relation,
    :target_measure_relation,
    :replay_status,
    :identity_status,
    :status,
)
const _LAURENT_ENDPOINT_REDUCTION_CERTIFICATE_STATUS =
    :endpoint_reduction_certificate
const _LAURENT_ENDPOINT_REDUCTION_CERTIFICATE_NEXT_BOUNDARY =
    :laurent_normality_replay
const _LAURENT_ENDPOINT_REDUCTION_CERTIFICATE_FIELDS = (
    :case_id,
    :dimension,
    :ring_generators,
    :context_status,
    :operation,
    :source_endpoint,
    :target_endpoint,
    :endpoint_measure_relation,
    :replay_status,
    :identity_status,
    :next_boundary,
    :status,
)
```

- [ ] **Step 2: Add helper implementations after `_validate_laurent_link_witness_certificate`**

Insert these internal functions:

```julia
function _laurent_endpoint_metadata_from_column(column, R, endpoint_index::Int; case_id)
    measure = _laurent_descent_measure_from_column(column, R; case_id)
    return _laurent_link_endpoint_metadata(
        column[endpoint_index],
        R,
        endpoint_index,
        measure;
        case_id,
    )
end

function _laurent_endpoint_reduction_status(endpoint_operation, n::Int, R)::Symbol
    _require_two_generator_laurent_ring(R)
    _laurent_descent_has_fields(
        endpoint_operation,
        _LAURENT_ENDPOINT_REDUCTION_OPERATION_FIELDS,
    ) || return :malformed_endpoint_operation
    endpoint_operation.family in _LAURENT_ENDPOINT_REDUCTION_OPERATION_FAMILIES ||
        return :malformed_endpoint_operation
    endpoint_operation.ring_generators == _laurent_descent_ring_generators(R) ||
        return :wrong_ring_generators

    endpoint_index = try
        _laurent_descent_checked_entry_index(
            endpoint_operation.endpoint_index,
            n,
            "endpoint_index",
        )
    catch err
        err isa InterruptException && rethrow()
        return :wrong_endpoint_index
    end

    _laurent_descent_has_fields(
        endpoint_operation.operation,
        _LAURENT_DESCENT_OPERATION_FIELDS,
    ) || return :malformed_endpoint_operation
    endpoint_operation.operation.target_index == endpoint_index ||
        return :wrong_endpoint_index
    endpoint_operation.operation.ring_generators ==
        endpoint_operation.ring_generators ||
        return :wrong_ring_generators

    operation_status = _laurent_descent_operation_status(
        endpoint_operation.operation,
        n,
        R,
    )
    operation_status == :ok || return operation_status == :wrong_ring_generators ?
                                  :wrong_ring_generators :
                                  :malformed_endpoint_operation
    return :ok
end

function _replay_laurent_endpoint_reduction(
    column,
    R,
    endpoint_operation;
    case_id,
    require_strict::Bool = true,
)
    status = _laurent_endpoint_reduction_status(
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
        case_id,
    )
    replayed_column = _replay_laurent_elementary_entry_addition(
        column,
        R,
        endpoint_operation.operation,
    )
    target_endpoint = _laurent_endpoint_metadata_from_column(
        replayed_column,
        R,
        endpoint_index;
        case_id,
    )
    relation = _strictly_decreases_laurent_measure(
        source_endpoint.column_measure,
        target_endpoint.column_measure,
    ) ? :strict_decrease : :not_strict_decrease
    require_strict && relation != :strict_decrease &&
        throw(
            ArgumentError(
                "endpoint operation does not strictly decrease the endpoint measure",
            ),
        )
    return (;
        endpoint_operation,
        replayed_column,
        source_endpoint,
        target_endpoint,
        relation,
    )
end

function _laurent_endpoint_reduction_candidate_from_replay(
    source_column,
    target_column,
    R,
    endpoint_operation;
    case_id,
    require_strict::Bool = true,
)
    length(source_column) == length(target_column) ||
        throw(ArgumentError("source and target columns must have the same length"))
    source_replay = _replay_laurent_endpoint_reduction(
        source_column,
        R,
        endpoint_operation;
        case_id,
        require_strict = false,
    )
    target_replay = _replay_laurent_endpoint_reduction(
        target_column,
        R,
        endpoint_operation;
        case_id,
        require_strict = false,
    )
    candidate_status =
        source_replay.relation == :strict_decrease &&
        target_replay.relation == :strict_decrease ?
        :strict_endpoint_decrease : :not_endpoint_reduction
    require_strict && candidate_status != :strict_endpoint_decrease &&
        throw(
            ArgumentError(
                "endpoint operation does not strictly decrease both endpoint measures",
            ),
        )
    return (;
        endpoint_operation,
        source_endpoint = source_replay.target_endpoint,
        target_endpoint = target_replay.target_endpoint,
        source_measure_relation = source_replay.relation,
        target_measure_relation = target_replay.relation,
        replay_status = :ok,
        identity_status = :verified,
        status = candidate_status,
    )
end

function _verify_laurent_endpoint_reduction_candidate(
    source_column,
    target_column,
    R,
    candidate,
)::Bool
    try
        _laurent_descent_has_fields(
            candidate,
            _LAURENT_ENDPOINT_REDUCTION_CANDIDATE_FIELDS,
        ) || return false
        candidate.replay_status == :ok || return false
        candidate.identity_status == :verified || return false
        candidate.source_measure_relation == :strict_decrease || return false
        candidate.target_measure_relation == :strict_decrease || return false
        candidate.status == :strict_endpoint_decrease || return false
        _laurent_endpoint_reduction_status(
            candidate.endpoint_operation,
            length(source_column),
            R,
        ) == :ok || return false
        expected = _laurent_endpoint_reduction_candidate_from_replay(
            source_column,
            target_column,
            R,
            candidate.endpoint_operation;
            case_id = candidate.source_endpoint.case_id,
        )
        return candidate == expected
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end

function _laurent_endpoint_reduction_certificate_from_replay(
    context,
    endpoint_operation,
    column,
    R;
    source_endpoint = nothing,
    require_strict::Bool = true,
)
    _laurent_descent_has_fields(
        context,
        (:case_id, :dimension, :ring_generators, :status),
    ) || throw(ArgumentError("context must include case_id, dimension, ring_generators, and status"))
    context.status == :endpoint_reduction_context ||
        throw(ArgumentError("context must have status :endpoint_reduction_context"))
    context.dimension == length(column) ||
        throw(ArgumentError("context dimension does not match source column"))
    context.ring_generators == _laurent_descent_ring_generators(R) ||
        throw(ArgumentError("context ring_generators do not match the ring"))

    replay = _replay_laurent_endpoint_reduction(
        column,
        R,
        endpoint_operation;
        case_id = context.case_id,
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
        next_boundary = _LAURENT_ENDPOINT_REDUCTION_CERTIFICATE_NEXT_BOUNDARY,
        status = _LAURENT_ENDPOINT_REDUCTION_CERTIFICATE_STATUS,
    )
end

function _validate_laurent_endpoint_reduction_certificate(cert, column, R)::Symbol
    try
        _require_two_generator_laurent_ring(R)
        _laurent_descent_has_fields(
            cert,
            _LAURENT_ENDPOINT_REDUCTION_CERTIFICATE_FIELDS,
        ) || return :missing_certificate_fields
        cert.status == _LAURENT_ENDPOINT_REDUCTION_CERTIFICATE_STATUS ||
            return :wrong_status
        cert.replay_status == :ok || return :wrong_replay_status
        cert.identity_status == :verified || return :wrong_identity_status
        cert.next_boundary == _LAURENT_ENDPOINT_REDUCTION_CERTIFICATE_NEXT_BOUNDARY ||
            return :wrong_next_boundary
        cert.dimension == length(column) || return :wrong_dimension
        cert.ring_generators == _laurent_descent_ring_generators(R) ||
            return :wrong_ring_generators
        cert.context_status == :endpoint_reduction_context ||
            return :wrong_context_status

        operation_status = _laurent_endpoint_reduction_status(
            cert.operation,
            length(column),
            R,
        )
        operation_status == :ok || return operation_status
        replay = _replay_laurent_endpoint_reduction(
            column,
            R,
            cert.operation;
            case_id = cert.case_id,
            require_strict = false,
        )
        replay.source_endpoint == cert.source_endpoint ||
            return :stale_source_endpoint
        replay.target_endpoint == cert.target_endpoint ||
            return :stale_target_endpoint
        cert.endpoint_measure_relation == replay.relation ||
            return :wrong_endpoint_measure_relation
        replay.relation == :strict_decrease || return :not_endpoint_reduction

        context = (;
            case_id = cert.case_id,
            dimension = cert.dimension,
            ring_generators = cert.ring_generators,
            status = cert.context_status,
        )
        expected = _laurent_endpoint_reduction_certificate_from_replay(
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

- [ ] **Step 3: Run the internal test**

Run:

```bash
julia --project=. -e 'include("test/internal/laurent_endpoint_reduction_helpers.jl")'
```

Expected: exit 0 after the new helper implementation.

### Task 3: Reuse Internals from Expert Search and Certificate Tests

**Files:**
- Modify: `test/expert/case008_d14_laurent_endpoint_reduction_search.jl`
- Modify: `test/expert/laurent_endpoint_reduction_certificate.jl`

**Interfaces:**
- Consumes: the new internal helpers from Task 2.
- Produces: expert test files that keep only d14-specific bounds and test-facing wrappers.

- [ ] **Step 1: Update the search candidate field list**

In `CASE008_D14_ENDPOINT_REDUCTION_CANDIDATE_FIELDS`, add
`:identity_status` between `:replay_status` and `:status`.

- [ ] **Step 2: Make `_case008_d14_endpoint_reduction_operation_status` delegate generic shape checks**

At the top of `_case008_d14_endpoint_reduction_operation_status`, call:

```julia
generic_status = Suslin._laurent_endpoint_reduction_status(
    endpoint_operation,
    n,
    R,
)
generic_status == :ok || return generic_status
```

Then keep only d14-specific checks for family, endpoint index set, source index
set, exponent bounds/vectors, and coefficient family. Remove duplicate generic
required-field, ring-generator, target/endpoint, nested-shape, and elementary
operation status checks from that function.

- [ ] **Step 3: Make the d14 candidate constructor delegate replay**

Replace the body of `_case008_d14_endpoint_reduction_candidate_from_replay` with
case-specific status validation followed by:

```julia
return Suslin._laurent_endpoint_reduction_candidate_from_replay(
    source_column,
    target_column,
    R,
    endpoint_operation;
    case_id = context.case_id,
    require_strict,
)
```

Keep the existing non-strict error string by ensuring the internal helper throws
`"endpoint operation does not strictly decrease both endpoint measures"`.

- [ ] **Step 4: Make the d14 candidate validator use internal verification**

After comparing the candidate to the expected replay result, add:

```julia
require_strict && !Suslin._verify_laurent_endpoint_reduction_candidate(
    columns.source_column,
    columns.target_column,
    columns.ring,
    candidate,
) && return :identity_replay_failed
```

Keep d14-specific report checks unchanged.

- [ ] **Step 5: Replace certificate test-only helpers with wrappers**

In `test/expert/laurent_endpoint_reduction_certificate.jl`, keep
`_laurent_endpoint_reduction_without_field` for tests. Replace the local helper
implementations with wrappers:

```julia
function _laurent_endpoint_metadata_from_column(column, R, endpoint_index::Int; case_id)
    return Suslin._laurent_endpoint_metadata_from_column(
        column,
        R,
        endpoint_index;
        case_id,
    )
end

function _laurent_endpoint_reduction_operation_status(endpoint_operation, n::Int, R)::Symbol
    return Suslin._laurent_endpoint_reduction_status(endpoint_operation, n, R)
end

function _laurent_endpoint_reduction_replay(
    context,
    endpoint_operation,
    column,
    R;
    require_strict::Bool,
)
    return Suslin._replay_laurent_endpoint_reduction(
        column,
        R,
        endpoint_operation;
        case_id = context.case_id,
        require_strict,
    )
end

function laurent_endpoint_reduction_certificate(
    context,
    endpoint_operation,
    column,
    R;
    source_endpoint = nothing,
    require_strict::Bool = true,
)
    return Suslin._laurent_endpoint_reduction_certificate_from_replay(
        context,
        endpoint_operation,
        column,
        R;
        source_endpoint,
        require_strict,
    )
end

function validate_laurent_endpoint_reduction_certificate(cert, column, R)::Symbol
    return Suslin._validate_laurent_endpoint_reduction_certificate(cert, column, R)
end
```

Remove now-duplicated local constants and local replay/validation bodies that
are no longer used by tests.

- [ ] **Step 6: Run focused expert verification**

Run:

```bash
julia --project=. -e 'include("test/expert/laurent_endpoint_reduction_certificate.jl")'
```

Expected: exit 0.

### Task 4: Final Verification and Commit

**Files:**
- Modify: all files changed by Tasks 1-3.

**Interfaces:**
- Consumes: completed code and tests.
- Produces: one implementation commit ready for PR review.

- [ ] **Step 1: Run focused internal verification**

Run:

```bash
julia --project=. -e 'include("test/internal/laurent_endpoint_reduction_helpers.jl")'
```

Expected: exit 0.

- [ ] **Step 2: Run focused expert verification**

Run:

```bash
julia --project=. -e 'include("test/expert/laurent_endpoint_reduction_certificate.jl")'
```

Expected: exit 0.

- [ ] **Step 3: Run full package verification**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: exit 0.

- [ ] **Step 4: Review changed files**

Run:

```bash
git diff --stat HEAD
git diff --check
```

Expected: no whitespace errors; changed files are limited to the issue-347 docs,
internal helpers, endpoint expert tests, and `test/runtests.jl`.

- [ ] **Step 5: Commit implementation**

Run:

```bash
git add src/algorithm/column_reduction.jl test/internal/laurent_endpoint_reduction_helpers.jl test/expert/case008_d14_laurent_endpoint_reduction_search.jl test/expert/laurent_endpoint_reduction_certificate.jl test/runtests.jl docs/superpowers/plans/2026-07-10-issue-347-laurent-endpoint-reduction-internals.md
git commit -m "Promote Laurent endpoint reduction internals"
```
