# Issue 341 Laurent Endpoint-Reduction Context Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a replay-derived expert-only Laurent endpoint-reduction context for the certified `case_008 d=14` link-witness certificate.

**Architecture:** The context lives in one expert test file and consumes the existing replay chain from `case008_d14_laurent_link_witness_certificate_summary`. Validation recomputes the certified source and target endpoint metadata through internal link-witness certificate helpers, so no diagnostic detail or trusted endpoint copy is used.

**Tech Stack:** Julia, Test stdlib, Oscar, Suslin internal expert helpers.

## Global Constraints

- Do not expose a public API.
- Do not read endpoint metadata from issue 340 production diagnostic detail.
- Do not search for or claim an endpoint reduction.
- The fixed witness fields are `pivot_index == 10`, `partner_index == 1`, and `witness_exponent == (1, -1)`.
- The unresolved endpoint boundary is `:laurent_endpoint_reduction`.

---

### Task 1: Expert Endpoint Context And Validator

**Files:**
- Create: `test/expert/case008_d14_laurent_endpoint_reduction_context.jl`
- Modify: `test/runtests.jl`

**Interfaces:**
- Consumes: `case008_d14_laurent_link_witness_certificate_summary(fixture)`, `_case008_d14_link_witness_source_data(fixture)`, and `validate_laurent_link_witness_certificate(cert, column, R)::Symbol`.
- Produces: `_case008_d14_endpoint_reduction_replay_source(fixture)::NamedTuple`, `case008_d14_laurent_endpoint_reduction_context(fixture = ToricBuilderCase008D14ColumnBoundary.boundary_fixture(); replay_source = _case008_d14_endpoint_reduction_replay_source(fixture))::NamedTuple`, and `validate_case008_d14_laurent_endpoint_reduction_context(context, fixture = ToricBuilderCase008D14ColumnBoundary.boundary_fixture(); replay_source = _case008_d14_endpoint_reduction_replay_source(fixture))::Symbol`.

- [ ] **Step 1: Write the expert test file with the intended API and assertions**

Create `test/expert/case008_d14_laurent_endpoint_reduction_context.jl` with these top-level pieces:

```julia
using Test
using Suslin
using Oscar

if !(@isdefined case008_d14_laurent_link_witness_certificate_summary)
    include(joinpath(@__DIR__, "laurent_link_witness_certificate.jl"))
end

const CASE008_D14_ENDPOINT_REDUCTION_SOURCE_BOUNDARY =
    :case008_d14_link_witness_certificate
const CASE008_D14_ENDPOINT_REDUCTION_BOUNDARY =
    :laurent_endpoint_reduction
const CASE008_D14_ENDPOINT_REDUCTION_CONTEXT_STATUS =
    :endpoint_reduction_context
const CASE008_D14_ENDPOINT_REDUCTION_REQUIRED_FIELDS =
    (:family, :endpoint_index, :operation, :ring_generators)
const CASE008_D14_ENDPOINT_REDUCTION_CONTEXT_FIELDS = (
    :case_id,
    :dimension,
    :ring_generators,
    :source_boundary,
    :boundary,
    :link_witness_status,
    :witness_family,
    :pivot_index,
    :partner_index,
    :witness_exponent,
    :source_endpoint,
    :target_endpoint,
    :required_endpoint_reduction_fields,
    :status,
)
```

- [ ] **Step 2: Add small NamedTuple helpers**

Add helper functions for field checks, field omission, and tuple replacement:

```julia
function _case008_d14_endpoint_context_has_required_fields(context)::Bool
    return all(
        field -> hasproperty(context, field),
        CASE008_D14_ENDPOINT_REDUCTION_CONTEXT_FIELDS,
    )
end

function _case008_d14_endpoint_without_field(value::NamedTuple, field::Symbol)
    kept = tuple((name for name in keys(value) if name != field)...)
    return NamedTuple{kept}(tuple((getproperty(value, name) for name in kept)...))
end
```

Add a replay-source helper so the verification command reconstructs the certified certificate once per testset and does not repeat the expensive d14 search for every negative-control mutation:

```julia
function _case008_d14_endpoint_reduction_replay_source(
    fixture = ToricBuilderCase008D14ColumnBoundary.boundary_fixture(),
)
    summary = case008_d14_laurent_link_witness_certificate_summary(fixture)
    summary.d14_status == :link_witness_certificate ||
        throw(ArgumentError("case_008 d14 requires a certified Laurent link-witness certificate before endpoint context construction"))
    certificate_validation = validate_laurent_link_witness_certificate(
        summary.certificate,
        summary.source_column,
        summary.ring,
    )
    certificate_validation == :ok ||
        throw(ArgumentError("case_008 d14 link-witness certificate must validate; got $(certificate_validation)"))
    return (; summary, certificate_validation)
end
```

- [ ] **Step 3: Implement the context builder**

The builder must reject a d14 search path that lacks a certified link witness and must validate the certificate before returning endpoint metadata:

```julia
function case008_d14_laurent_endpoint_reduction_context(
    fixture = ToricBuilderCase008D14ColumnBoundary.boundary_fixture(),
    ;
    replay_source = _case008_d14_endpoint_reduction_replay_source(fixture),
)
    summary = replay_source.summary
    replay_source.certificate_validation == :ok ||
        throw(ArgumentError("case_008 d14 link-witness certificate must validate; got $(replay_source.certificate_validation)"))
    witness = summary.certificate.witness
    return (;
        case_id = summary.certificate.case_id,
        dimension = summary.certificate.dimension,
        ring_generators = summary.certificate.ring_generators,
        source_boundary = CASE008_D14_ENDPOINT_REDUCTION_SOURCE_BOUNDARY,
        boundary = CASE008_D14_ENDPOINT_REDUCTION_BOUNDARY,
        link_witness_status = summary.certificate.status,
        witness_family = witness.family,
        pivot_index = witness.pivot_index,
        partner_index = witness.partner_index,
        witness_exponent = witness.exponent,
        source_endpoint = summary.certificate.source_endpoint,
        target_endpoint = summary.certificate.target_endpoint,
        required_endpoint_reduction_fields =
            CASE008_D14_ENDPOINT_REDUCTION_REQUIRED_FIELDS,
        status = CASE008_D14_ENDPOINT_REDUCTION_CONTEXT_STATUS,
    )
end
```

- [ ] **Step 4: Implement the validator**

The validator recomputes the expected context from the replayed certificate and returns stable symbols for negative controls:

```julia
function validate_case008_d14_laurent_endpoint_reduction_context(
    context,
    fixture = ToricBuilderCase008D14ColumnBoundary.boundary_fixture(),
    ;
    replay_source = _case008_d14_endpoint_reduction_replay_source(fixture),
)::Symbol
    try
        _case008_d14_endpoint_context_has_required_fields(context) ||
            return :missing_context_fields
        summary = replay_source.summary
        summary.d14_status == :link_witness_certificate ||
            return :missing_link_witness_certificate
        replay_source.certificate_validation == :ok ||
            return :invalid_link_witness_certificate

        context.case_id == "case_008" || return :wrong_case
        context.dimension == 14 || return :wrong_dimension
        context.ring_generators == ("u", "v") || return :wrong_ring_generators
        context.source_boundary == CASE008_D14_ENDPOINT_REDUCTION_SOURCE_BOUNDARY ||
            return :wrong_source_boundary
        context.boundary == CASE008_D14_ENDPOINT_REDUCTION_BOUNDARY ||
            return :wrong_boundary
        context.link_witness_status == :link_witness_certificate ||
            return :wrong_link_witness_status
        context.witness_family == :two_entry_laurent_combination ||
            return :wrong_witness_family
        context.pivot_index == 10 || return :wrong_pivot_index
        context.partner_index == 1 || return :wrong_partner_index
        context.witness_exponent == (1, -1) || return :wrong_witness_exponent
        context.source_endpoint == summary.certificate.source_endpoint ||
            return :stale_source_endpoint
        context.target_endpoint == summary.certificate.target_endpoint ||
            return :stale_target_endpoint
        context.required_endpoint_reduction_fields isa Tuple ||
            return :wrong_required_endpoint_reduction_fields
        all(
            field -> field in context.required_endpoint_reduction_fields,
            CASE008_D14_ENDPOINT_REDUCTION_REQUIRED_FIELDS,
        ) || return :missing_required_endpoint_reduction_field
        context.required_endpoint_reduction_fields ==
            CASE008_D14_ENDPOINT_REDUCTION_REQUIRED_FIELDS ||
            return :wrong_required_endpoint_reduction_fields
        context.status == CASE008_D14_ENDPOINT_REDUCTION_CONTEXT_STATUS ||
            return :wrong_status

        expected = case008_d14_laurent_endpoint_reduction_context(
            fixture;
            replay_source,
        )
        context == expected || return :stale_context
        return :ok
    catch err
        err isa InterruptException && rethrow()
        return :invalid_context
    end
end
```

- [ ] **Step 5: Add positive and negative controls**

Use two `@testset`s. The positive test asserts every fixed issue field and endpoint recomputation. The negative test mutates the accepted context and asserts these validator symbols:

```julia
@test validate_case008_d14_laurent_endpoint_reduction_context(
    merge(context, (; link_witness_status = :stale_link_witness_certificate)),
    fixture,
) == :wrong_link_witness_status

@test validate_case008_d14_laurent_endpoint_reduction_context(
    merge(context, (; ring_generators = ("v", "u"))),
    fixture,
) == :wrong_ring_generators

@test validate_case008_d14_laurent_endpoint_reduction_context(
    merge(context, (; source_endpoint = merge(context.source_endpoint, (; term_count = context.source_endpoint.term_count + 1)))),
    fixture,
) == :stale_source_endpoint

@test validate_case008_d14_laurent_endpoint_reduction_context(
    merge(context, (; target_endpoint = merge(context.target_endpoint, (; term_count = context.target_endpoint.term_count + 1)))),
    fixture,
) == :stale_target_endpoint

@test validate_case008_d14_laurent_endpoint_reduction_context(
    merge(context, (; witness_exponent = (0, 0))),
    fixture,
) == :wrong_witness_exponent

@test validate_case008_d14_laurent_endpoint_reduction_context(
    merge(context, (; required_endpoint_reduction_fields = (:family, :endpoint_index, :operation))),
    fixture,
) == :missing_required_endpoint_reduction_field
```

- [ ] **Step 6: Register the expert test**

In `test/runtests.jl`, insert the new expert file immediately after `"expert/laurent_link_witness_certificate.jl"`:

```julia
"expert/case008_d14_laurent_endpoint_reduction_context.jl",
```

- [ ] **Step 7: Run targeted verification**

Run:

```bash
julia --project=. -e 'include("test/expert/case008_d14_laurent_endpoint_reduction_context.jl")'
```

Expected: command exits 0.

- [ ] **Step 8: Run full package verification**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: command exits 0 and reports public/internal tests passing.

- [ ] **Step 9: Check formatting hygiene**

Run:

```bash
git diff --check
```

Expected: no output and exit 0.

- [ ] **Step 10: Commit implementation**

Run:

```bash
git add docs/superpowers/plans/2026-07-09-issue-341-laurent-endpoint-reduction-context.md test/expert/case008_d14_laurent_endpoint_reduction_context.jl test/runtests.jl
git commit -m "Add issue 341 Laurent endpoint context"
```

Expected: commit succeeds with the plan, expert test, and suite registration.
