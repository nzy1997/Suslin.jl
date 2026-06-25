# Issue 133 Laurent Column Reduction Diagnostics Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an internal structured diagnostic for exact unimodular-column reduction failures, including stable unsupported Laurent failure codes.

**Architecture:** Keep the diagnostic beside `reduce_unimodular_column` and reuse existing stage helpers in reducer order. The reducer keeps its existing throwing behavior; tests inspect the unexported diagnostic through `Suslin.diagnose_unimodular_column_reduction`.

**Tech Stack:** Julia, Oscar, Suslin internal reducer helpers, Test stdlib.

## Global Constraints

- Do not export `diagnose_unimodular_column_reduction` from `src/Suslin.jl`.
- Do not add a new Laurent reduction algorithm.
- Do not change `reduce_unimodular_column` throwing behavior.
- Unsupported validated Laurent failures must use `failure_code == :unsupported_laurent_column_family`.
- Non-unimodular inputs must use a precondition failure, not `:unsupported_laurent_column_family`.
- The issue verification command must be `julia --project=. -e 'include("test/expert/laurent_column_reduction_diagnostics.jl")'`.
- The full verification command must be `julia --project=. -e 'using Pkg; Pkg.test()'`.

---

### Task 1: Add Internal Diagnostic And Expert Tests

**Files:**
- Modify: `src/algorithm/column_reduction.jl`
- Create: `test/expert/laurent_column_reduction_diagnostics.jl`
- Modify: `test/runtests.jl`

**Interfaces:**
- Consumes: `reduce_unimodular_column(v, R)`, `_validated_unimodular_column(v, R)`, `_reduce_laurent_unimodular_column_certificate(column, R)`, `_reduce_polynomial_unimodular_column_exact_certificate(column, R)`, `_reduce_after_monicity_normalization_certificate(column, R)`, `_reduce_via_supported_three_block_certificate(column, R)`, `ToricBuilderCase010ColumnBoundary.boundary_fixture()`.
- Produces: `diagnose_unimodular_column_reduction(v::AbstractVector, R)` returning a named tuple with fields `status`, `failure_code`, `ring_profile`, `column_length`, `attempted_stages`, and `message`.

- [ ] **Step 1: Write the failing expert test**

Create `test/expert/laurent_column_reduction_diagnostics.jl`:

```julia
using Test
using Suslin
using Oscar

include(joinpath(@__DIR__, "..", "fixtures", "toricbuilder_case010_column_boundary.jl"))

@testset "Laurent column reduction diagnostics" begin
    fixture = ToricBuilderCase010ColumnBoundary.boundary_fixture()

    unsupported = Suslin.diagnose_unimodular_column_reduction(
        fixture.failing_column,
        fixture.ring,
    )
    @test unsupported.status == :unsupported
    @test unsupported.failure_code == :unsupported_laurent_column_family
    @test unsupported.column_length == 5
    @test unsupported.ring_profile.kind == :laurent_polynomial
    @test unsupported.ring_profile.generators == ("u", "v")
    @test occursin("unsupported exact unimodular column reduction", unsupported.message)
    for stage in (:unit_entry, :witness_unit, :monicity_normalization, :three_entry_block)
        @test stage in unsupported.attempted_stages
    end
    @test :laurent_normalization in unsupported.attempted_stages

    d = nrows(fixture.normalized_matrix)
    supported_column = [fixture.normalized_matrix[row, d] for row in 1:d]
    supported = Suslin.diagnose_unimodular_column_reduction(supported_column, fixture.ring)
    @test supported.status == :supported
    @test supported.failure_code === nothing
    @test supported.column_length == d
    @test !isempty(supported.attempted_stages)
    Suslin.reduce_unimodular_column(supported_column, fixture.ring)

    negative = ToricBuilderCase010ColumnBoundary.non_unimodular_negative_control(fixture)
    precondition = Suslin.diagnose_unimodular_column_reduction(
        negative.failing_column,
        negative.ring,
    )
    @test precondition.status == :precondition_failed
    @test precondition.failure_code == :not_unimodular
    @test precondition.failure_code != :unsupported_laurent_column_family
    @test isempty(precondition.attempted_stages)
end
```

- [ ] **Step 2: Run the test to verify RED**

Run:

```bash
julia --project=. -e 'include("test/expert/laurent_column_reduction_diagnostics.jl")'
```

Expected: FAIL with `UndefVarError: diagnose_unimodular_column_reduction not defined`.

- [ ] **Step 3: Add the diagnostic implementation**

In `src/algorithm/column_reduction.jl`, insert these helpers after
`ecp_column_reduction_certificate` and before `_ecp_column_reduction_certificate_validated`:

```julia
function diagnose_unimodular_column_reduction(v::AbstractVector, R)
    ring_profile = _column_reduction_ring_profile(R)
    column_length = length(v)
    validation = _diagnose_unimodular_column_preconditions(v, R, ring_profile, column_length)
    validation.status == :ok || return validation.diagnostic
    return _diagnose_unimodular_column_reduction_validated(
        validation.column,
        R,
        ring_profile,
        column_length,
    )
end

function _column_reduction_diagnostic(
    status::Symbol,
    failure_code,
    ring_profile,
    column_length::Int,
    attempted_stages,
    message::AbstractString,
)
    return (;
        status,
        failure_code,
        ring_profile,
        column_length,
        attempted_stages = tuple(attempted_stages...),
        message = String(message),
    )
end

function _column_reduction_ring_profile(R)
    kind = try
        _is_laurent_polynomial_ring(R) ? :laurent_polynomial : :polynomial
    catch err
        err isa InterruptException && rethrow()
        :unknown
    end
    coefficient_ring = try
        string(base_ring(R))
    catch err
        err isa InterruptException && rethrow()
        ""
    end
    generator_names = try
        tuple(string.(gens(R))...)
    catch err
        err isa InterruptException && rethrow()
        ()
    end
    return (; kind, coefficient_ring, generators = generator_names)
end

function _diagnose_unimodular_column_preconditions(v::AbstractVector, R, ring_profile, column_length::Int)
    try
        Base.require_one_based_indexing(v)
    catch err
        err isa InterruptException && rethrow()
        return (;
            status = :precondition_failed,
            diagnostic = _column_reduction_diagnostic(
                :precondition_failed,
                :unsupported_indexing,
                ring_profile,
                column_length,
                Symbol[],
                "v must use one-based indexing",
            ),
        )
    end

    if column_length < 3
        return (;
            status = :precondition_failed,
            diagnostic = _column_reduction_diagnostic(
                :precondition_failed,
                :column_too_short,
                ring_profile,
                column_length,
                Symbol[],
                "v must have length at least 3",
            ),
        )
    end

    column = try
        [_coerce_into_ring(R, v[idx], "v[$idx]") for idx in 1:column_length]
    catch err
        err isa InterruptException && rethrow()
        return (;
            status = :precondition_failed,
            diagnostic = _column_reduction_diagnostic(
                :precondition_failed,
                :column_not_over_ring,
                ring_profile,
                column_length,
                Symbol[],
                _column_reduction_error_message(err),
            ),
        )
    end

    is_unimodular = try
        is_unimodular_column(column, R)
    catch err
        err isa InterruptException && rethrow()
        return (;
            status = :precondition_failed,
            diagnostic = _column_reduction_diagnostic(
                :precondition_failed,
                :unimodularity_check_failed,
                ring_profile,
                column_length,
                Symbol[],
                _column_reduction_error_message(err),
            ),
        )
    end

    is_unimodular || return (;
        status = :precondition_failed,
        diagnostic = _column_reduction_diagnostic(
            :precondition_failed,
            :not_unimodular,
            ring_profile,
            column_length,
            Symbol[],
            "v must be a unimodular column",
        ),
    )
    return (; status = :ok, column)
end

function _diagnose_unimodular_column_reduction_validated(column::AbstractVector, R, ring_profile, column_length::Int)
    attempted = Symbol[]
    result = _is_laurent_polynomial_ring(R) ?
        _diagnose_laurent_unimodular_column_reduction(column, R, attempted) :
        _diagnose_polynomial_unimodular_column_reduction(column, R, attempted)

    if result.supported
        return _column_reduction_diagnostic(
            :supported,
            nothing,
            ring_profile,
            column_length,
            attempted,
            "exact unimodular column reduction is supported by $(result.stage)",
        )
    end

    failure_code = _is_laurent_polynomial_ring(R) ?
        :unsupported_laurent_column_family :
        :unsupported_polynomial_column_family
    return _column_reduction_diagnostic(
        :unsupported,
        failure_code,
        ring_profile,
        column_length,
        attempted,
        _unsupported_unimodular_column_reduction_message(column, R),
    )
end

function _diagnose_laurent_unimodular_column_reduction(column::AbstractVector, R, attempted::Vector{Symbol})
    push!(attempted, :unit_entry)
    unit_idx = findfirst(is_unit, column)
    unit_idx !== nothing && return (; supported = true, stage = :unit_entry)

    push!(attempted, :laurent_normalization)
    normalization = normalize_laurent_object(column)
    poly_column = normalization.normalized_object
    P = normalization.metadata.polynomial_ring
    is_unimodular_column(poly_column, P) ||
        return (; supported = false, stage = nothing)
    return _diagnose_polynomial_unimodular_column_reduction(poly_column, P, attempted)
end

function _diagnose_polynomial_unimodular_column_reduction(column::AbstractVector, R, attempted::Vector{Symbol})
    small = _diagnose_exact_small_column_reduction(column, R, attempted)
    small.supported && return small

    if length(column) > 3
        push!(attempted, :three_entry_block)
        block = _reduce_via_supported_three_block_certificate(column, R)
        block !== nothing && return (; supported = true, stage = :three_entry_block)
    end

    return (; supported = false, stage = nothing)
end

function _diagnose_exact_small_column_reduction(column::AbstractVector, R, attempted::Vector{Symbol})
    supported = _diagnose_supported_unimodular_column_reduction(column, R, attempted)
    supported.supported && return supported

    if _has_at_least_two_generators(R)
        push!(attempted, :monicity_normalization)
        normalized = _reduce_after_monicity_normalization_certificate(column, R)
        normalized !== nothing && return (; supported = true, stage = :monicity_normalization)
    end

    return (; supported = false, stage = nothing)
end

function _diagnose_supported_unimodular_column_reduction(column::AbstractVector, R, attempted::Vector{Symbol})
    push!(attempted, :unit_entry)
    unit_idx = findfirst(is_unit, column)
    unit_idx !== nothing && return (; supported = true, stage = :unit_entry)

    push!(attempted, :witness_unit)
    witness = _unimodular_witness(column, R)
    witness_unit_idx = findfirst(is_unit, witness)
    witness_unit_idx !== nothing && return (; supported = true, stage = :witness_unit)

    return (; supported = false, stage = nothing)
end

function _column_reduction_error_message(err)
    return sprint(showerror, err)
end
```

- [ ] **Step 4: Add the reducer message helper**

Replace `_throw_unsupported_unimodular_column_reduction(column, R)` with a
shared message helper:

```julia
function _unsupported_unimodular_column_reduction_message(column::AbstractVector, R)
    if _is_laurent_polynomial_ring(R)
        return "unsupported exact unimodular column reduction for Laurent-normalized column of length $(length(column))"
    end
    return "unsupported exact unimodular column reduction for column of length $(length(column)) over $(R)"
end

function _throw_unsupported_unimodular_column_reduction(column::AbstractVector, R)
    throw(ArgumentError(_unsupported_unimodular_column_reduction_message(column, R)))
end
```

- [ ] **Step 5: Register the expert test**

In `test/runtests.jl`, insert `"expert/laurent_column_reduction_diagnostics.jl",`
immediately after `"expert/laurent_column_peel_issue38.jl",`.

- [ ] **Step 6: Run focused verification**

Run:

```bash
julia --project=. -e 'include("test/expert/laurent_column_reduction_diagnostics.jl")'
```

Expected: PASS with no failures.

- [ ] **Step 7: Run full verification**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: PASS for public and internal groups.

- [ ] **Step 8: Commit**

Run:

```bash
git add src/algorithm/column_reduction.jl test/expert/laurent_column_reduction_diagnostics.jl test/runtests.jl docs/superpowers/plans/2026-06-25-issue-133-laurent-column-reduction-diagnostics.md
git commit -m "feat: diagnose Laurent column reduction failures"
```

Expected: commit succeeds and contains the diagnostic, expert test, test registration, and this plan.

---

## Self-Review

- Spec coverage: Task 1 covers the internal diagnostic entry point, stable Laurent failure code, stage-attempt recording, supported Laurent negative control, non-unimodular precondition handling, focused verification, and full verification.
- Placeholder scan: no open-ended placeholders are used.
- Type consistency: the plan consistently uses `diagnose_unimodular_column_reduction(v, R)` and a named tuple with the required fields.
