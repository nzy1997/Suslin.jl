# Issue 134 Case 010 Laurent Column Reduction Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a certified Laurent one-step unit-creation stage that reduces the extracted `case_010` length-5 boundary column to `[0, 0, 0, 0, 1]`.

**Architecture:** Extend `src/algorithm/column_reduction.jl` with a narrow Laurent stage before the existing Laurent-normalization path. The stage creates a literal unit entry via one exact elementary row operation, then delegates to the existing unit-entry certificate and replay machinery. Tests cover the new `case_010` acceptance path, replay tampering rejection, and diagnostics after the former unsupported boundary becomes supported.

**Tech Stack:** Julia, Oscar, Suslin internal column-reduction certificates, Test stdlib.

## Global Constraints

- Repository has no `AGENTS.md`, `CLAUDE.md`, `GEMINI.md`, or `CONVENTIONS.md`; follow `README.md` test commands and existing Suslin style.
- Do not update the `case_010` report status.
- Do not claim support for arbitrary Laurent unimodular columns.
- Do not export new public APIs.
- Keep the new algorithm in `src/algorithm/column_reduction.jl`.
- The new stage must record enough metadata for `verify_ecp_column_reduction` replay to reject tampering.
- The issue verification command must be `julia --project=. -e 'include("test/expert/case010_laurent_column_reduction.jl")'`.
- The full verification command must be `julia --project=. -e 'using Pkg; Pkg.test()'`.

---

## File Structure

- Modify `src/algorithm/column_reduction.jl`: add the Laurent one-step unit-creation stage, hook it into Laurent reduction and diagnostics, and add exact replay support for the new stage kind.
- Modify `test/fixtures/toricbuilder_case010_column_boundary.jl`: stop treating the former unsupported reducer error as part of fixture validity, because this issue makes the boundary reducible.
- Modify `test/internal/toricbuilder_case010_column_boundary.jl`: update the boundary fixture test to assert the diagnostic is now supported by the new stage.
- Modify `test/expert/laurent_column_reduction_diagnostics.jl`: keep unsupported Laurent diagnostics covered with a separate unsupported column, and assert `case_010` is supported through `:laurent_unit_creation`.
- Create `test/expert/case010_laurent_column_reduction.jl`: focused issue acceptance test and negative controls.
- Modify `test/runtests.jl`: register the new expert test.

---

### Task 1: Certified Case 010 Laurent Unit-Creation Reduction

**Files:**
- Modify: `src/algorithm/column_reduction.jl`
- Modify: `test/fixtures/toricbuilder_case010_column_boundary.jl`
- Modify: `test/internal/toricbuilder_case010_column_boundary.jl`
- Modify: `test/expert/laurent_column_reduction_diagnostics.jl`
- Create: `test/expert/case010_laurent_column_reduction.jl`
- Modify: `test/runtests.jl`

**Interfaces:**
- Consumes: `ToricBuilderCase010ColumnBoundary.boundary_fixture()`, `reduce_unimodular_column(v, R)`, `ecp_column_reduction_certificate(v, R)`, `diagnose_unimodular_column_reduction(v, R)`, existing unit-entry certificate helpers, `verify_ecp_column_reduction(certificate)`.
- Produces: internal `_reduce_via_laurent_unit_creation_certificate(column, R)` returning `nothing` or `(; factors, stage)` with `stage.kind == :laurent_unit_creation`.

- [ ] **Step 1: Write the failing expert acceptance test**

Create `test/expert/case010_laurent_column_reduction.jl` with this content:

```julia
using Test
using Suslin
using Oscar

include(joinpath(@__DIR__, "..", "fixtures", "toricbuilder_case010_column_boundary.jl"))

function _case010_reduction_product(factors, R, n::Int)
    product = identity_matrix(R, n)
    for factor in factors
        product *= factor
    end
    return product
end

function _case010_apply_factors(factors, column, R)
    n = length(column)
    return _case010_reduction_product(factors, R, n) * matrix(R, n, 1, collect(column))
end

function _case010_target_column(R, n::Int)
    target = zero_matrix(R, n, 1)
    target[n, 1] = one(R)
    return target
end

function _case010_tamper_first_factor(factors, R, n::Int)
    tampered = copy(factors)
    tampered[1] = identity_matrix(R, n)
    return tampered
end

function _case010_tamper_certificate_first_factor(cert)
    tampered = _case010_tamper_first_factor(cert.factors, cert.ring, length(cert.original_column))
    return Suslin.ECPColumnReductionCertificate(
        cert.original_column,
        cert.ring,
        cert.stages,
        tampered,
        cert.final_column,
        cert.verification,
    )
end

@testset "case_010 Laurent column reduction" begin
    fixture = ToricBuilderCase010ColumnBoundary.boundary_fixture()
    R = fixture.ring
    column = fixture.failing_column
    target = _case010_target_column(R, length(column))

    factors = Suslin.reduce_unimodular_column(column, R)
    @test _case010_apply_factors(factors, column, R) == target

    tampered_factors = _case010_tamper_first_factor(factors, R, length(column))
    @test _case010_apply_factors(tampered_factors, column, R) != target

    certificate = Suslin.ecp_column_reduction_certificate(column, R)
    @test Suslin.verify_ecp_column_reduction(certificate)
    @test any(stage -> stage.kind == :laurent_unit_creation, certificate.stages)
    @test _case010_apply_factors(certificate.factors, certificate.original_column, certificate.ring) == target
    @test !Suslin.verify_ecp_column_reduction(_case010_tamper_certificate_first_factor(certificate))

    diagnostic = Suslin.diagnose_unimodular_column_reduction(column, R)
    @test diagnostic.status == :supported
    @test diagnostic.failure_code === nothing
    @test :laurent_unit_creation in diagnostic.attempted_stages
end
```

- [ ] **Step 2: Run the focused test to verify RED**

Run:

```bash
julia --project=. -e 'include("test/expert/case010_laurent_column_reduction.jl")'
```

Expected: FAIL before implementation. The failure should be the current unsupported exact Laurent column-reduction error for the fixture column, or `SystemError`/`LoadError` if the new test file has not yet been created correctly. Fix test-file syntax until the failure is the unsupported reducer behavior.

- [ ] **Step 3: Add the Laurent unit-creation stage implementation**

In `src/algorithm/column_reduction.jl`, insert these helpers after `_witness_unit_reduction_certificate_stage`:

```julia
function _try_laurent_divexact(numerator, denominator)
    iszero(denominator) && return nothing
    try
        return divexact(numerator, denominator)
    catch err
        err isa InterruptException && rethrow()
        (err isa ErrorException || err isa ArgumentError || err isa MethodError) && return nothing
        rethrow()
    end
end

function _laurent_unit_creation_candidate(column::AbstractVector, R)
    _is_laurent_polynomial_ring(R) || return nothing
    target_unit = one(R)

    for pivot_idx in eachindex(column), source_idx in eachindex(column)
        pivot_idx == source_idx && continue
        coeff = _try_laurent_divexact(target_unit - column[pivot_idx], column[source_idx])
        coeff === nothing && continue
        coeff == zero(R) && continue
        column[pivot_idx] + coeff * column[source_idx] == target_unit || continue
        return (;
            pivot_index = pivot_idx,
            source_index = source_idx,
            target_unit,
            creation_coefficient = coeff,
        )
    end

    return nothing
end

function _laurent_unit_creation_factors(n::Int, pivot_idx::Int, source_idx::Int, coeff, R)
    return [elementary_matrix(n, pivot_idx, source_idx, coeff, R)]
end

function _reduce_via_laurent_unit_creation_certificate(column::AbstractVector, R)
    candidate = _laurent_unit_creation_candidate(column, R)
    candidate === nothing && return nothing
    return _laurent_unit_creation_certificate_stage(
        column,
        R,
        candidate.pivot_index,
        candidate.source_index,
        candidate.target_unit,
        candidate.creation_coefficient,
    )
end

function _laurent_unit_creation_certificate_stage(
    column::AbstractVector,
    R,
    pivot_idx::Int,
    source_idx::Int,
    target_unit,
    creation_coefficient,
)
    n = length(column)
    creation_factors = _laurent_unit_creation_factors(n, pivot_idx, source_idx, creation_coefficient, R)
    created_column_matrix = _apply_reduction_factors(creation_factors, column, R)
    created_column = collect(_ecp_matrix_column_to_tuple(created_column_matrix))
    created_column[pivot_idx] == target_unit || error("internal Laurent unit-creation stage did not create the target unit")

    unit_stage = _unit_entry_reduction_certificate_stage(created_column, pivot_idx, R).stage
    factors = _checked_reduction_factors(
        vcat(unit_stage.factors, creation_factors),
        column,
        R,
        "Laurent unit-creation reduction",
    )
    stage = (;
        kind = :laurent_unit_creation,
        input_column = _ecp_column_tuple(column),
        pivot_index = pivot_idx,
        source_index = source_idx,
        target_unit,
        creation_coefficient,
        creation_factors,
        created_column = tuple(created_column...),
        unit_stage,
        factors,
        output_column = _apply_reduction_factors(factors, column, R),
    )
    return (; factors, stage)
end
```

Then modify `_reduce_laurent_unimodular_column_certificate(column, R)` so the beginning becomes:

```julia
function _reduce_laurent_unimodular_column_certificate(column::AbstractVector, R)
    unit_idx = findfirst(is_unit, column)
    unit_idx !== nothing && return _unit_entry_reduction_certificate_stage(column, unit_idx, R)

    unit_creation = _reduce_via_laurent_unit_creation_certificate(column, R)
    unit_creation !== nothing && return unit_creation

    normalization = normalize_laurent_object(column)
```

Modify `_diagnose_laurent_unimodular_column_reduction(column, R, attempted)` so it records and checks the new stage before Laurent normalization:

```julia
function _diagnose_laurent_unimodular_column_reduction(column::AbstractVector, R, attempted::Vector{Symbol})
    push!(attempted, :unit_entry)
    unit_idx = findfirst(is_unit, column)
    unit_idx !== nothing && return (; supported = true, stage = :unit_entry)

    push!(attempted, :laurent_unit_creation)
    unit_creation = _reduce_via_laurent_unit_creation_certificate(column, R)
    unit_creation !== nothing && return (; supported = true, stage = :laurent_unit_creation)

    push!(attempted, :laurent_normalization)
    normalization = normalize_laurent_object(column)
    poly_column = normalization.normalized_object
    P = normalization.metadata.polynomial_ring
    return _diagnose_polynomial_unimodular_column_reduction(poly_column, P, attempted)
end
```

Finally, add this branch to `_ecp_replay_stage(stage, input_column, R)` after the `:witness_unit` branch and before `:monicity_normalization`:

```julia
    elseif stage.kind == :laurent_unit_creation
        invalid_replay = (; ok = false, factors = Any[], output_column = matrix(R, length(input_column), 1, collect(input_column)))
        required_keys = (
            :kind,
            :input_column,
            :pivot_index,
            :source_index,
            :target_unit,
            :creation_coefficient,
            :creation_factors,
            :created_column,
            :unit_stage,
            :factors,
            :output_column,
        )
        _ecp_stage_keys_ok(stage, required_keys) || return invalid_replay
        _is_laurent_polynomial_ring(R) || return invalid_replay

        pivot_idx = stage.pivot_index
        source_idx = stage.source_index
        pivot_idx isa Integer && source_idx isa Integer || return invalid_replay
        1 <= pivot_idx <= length(input_column) || return invalid_replay
        1 <= source_idx <= length(input_column) || return invalid_replay
        pivot_idx == source_idx && return invalid_replay
        stage.target_unit == one(R) || return invalid_replay

        coeff = _try_laurent_divexact(stage.target_unit - input_column[pivot_idx], input_column[source_idx])
        coeff === nothing && return invalid_replay
        creation_factors = _laurent_unit_creation_factors(length(input_column), pivot_idx, source_idx, coeff, R)
        created_column = collect(_ecp_matrix_column_to_tuple(_apply_reduction_factors(creation_factors, input_column, R)))
        unit_replay = _ecp_replay_stage(stage.unit_stage, created_column, R)
        expected_factors = vcat(unit_replay.factors, creation_factors)
        expected_output = _apply_reduction_factors(expected_factors, input_column, R)
        ok = stage.input_column == _ecp_column_tuple(input_column) &&
            stage.creation_coefficient == coeff &&
            created_column[pivot_idx] == stage.target_unit &&
            stage.created_column == tuple(created_column...) &&
            unit_replay.ok &&
            _ecp_factor_sequences_equal(stage.creation_factors, creation_factors) &&
            _ecp_factor_sequences_equal(stage.factors, expected_factors) &&
            stage.output_column == expected_output
        return (; ok, factors = expected_factors, output_column = expected_output)
```

- [ ] **Step 4: Update the boundary fixture validation**

In `test/fixtures/toricbuilder_case010_column_boundary.jl`, remove the `_column_reduction_diagnostic_status` helper and replace the end of `validate_boundary_fixture(fixture)::Symbol`:

```julia
    _column_entries_are_over_ring(fixture.failing_column, R) || return :wrong_ring
    Suslin.is_unimodular_column(fixture.failing_column, R) || return :not_unimodular

    diagnostic_status = _column_reduction_diagnostic_status(
        fixture.failing_column,
        R,
        fixture.expected_diagnostic,
    )
    diagnostic_status == :expected_diagnostic || return diagnostic_status

    return :ok
end
```

with:

```julia
    _column_entries_are_over_ring(fixture.failing_column, R) || return :wrong_ring
    Suslin.is_unimodular_column(fixture.failing_column, R) || return :not_unimodular

    return :ok
end
```

In `test/internal/toricbuilder_case010_column_boundary.jl`, replace the reducer-error assertion block:

```julia
    err = try
        Suslin.reduce_unimodular_column(fixture.failing_column, fixture.ring)
        nothing
    catch caught
        caught
    end
    @test err isa ArgumentError
    @test occursin(fixture.expected_diagnostic, sprint(showerror, err))
```

with:

```julia
    diagnostic = Suslin.diagnose_unimodular_column_reduction(
        fixture.failing_column,
        fixture.ring,
    )
    @test diagnostic.status == :supported
    @test :laurent_unit_creation in diagnostic.attempted_stages
```

- [ ] **Step 5: Update the diagnostics expert test**

In `test/expert/laurent_column_reduction_diagnostics.jl`, replace the initial unsupported `case_010` block with:

```julia
    case010 = Suslin.diagnose_unimodular_column_reduction(
        fixture.failing_column,
        fixture.ring,
    )
    @test case010.status == :supported
    @test case010.failure_code === nothing
    @test case010.column_length == 5
    @test case010.ring_profile.kind == :laurent_polynomial
    @test case010.ring_profile.generators == ("u", "v")
    @test :laurent_unit_creation in case010.attempted_stages

    unsupported_ring, (x, y) = Suslin.suslin_laurent_polynomial_ring(GF(2), ["x", "y"])
    unsupported_column = [x + y, x * y + x + one(unsupported_ring), x^2 + x * y + y^2 + one(unsupported_ring)]
    unsupported = Suslin.diagnose_unimodular_column_reduction(
        unsupported_column,
        unsupported_ring,
    )
    @test unsupported.status == :unsupported
    @test unsupported.failure_code == :unsupported_laurent_column_family
    @test unsupported.column_length == 3
    @test unsupported.ring_profile.kind == :laurent_polynomial
    @test unsupported.ring_profile.generators == ("x", "y")
    @test occursin("unsupported exact unimodular column reduction", unsupported.message)
    for stage in (:unit_entry, :laurent_unit_creation, :laurent_normalization, :witness_unit, :monicity_normalization)
        @test stage in unsupported.attempted_stages
    end
```

Also delete the old block that expected `Suslin.reduce_unimodular_column(fixture.failing_column, fixture.ring)` to throw for `case_010`.

- [ ] **Step 6: Register the new expert test**

In `test/runtests.jl`, insert:

```julia
        "expert/case010_laurent_column_reduction.jl",
```

immediately after:

```julia
        "expert/laurent_column_reduction_diagnostics.jl",
```

- [ ] **Step 7: Run focused verification to verify GREEN**

Run:

```bash
julia --project=. -e 'include("test/expert/case010_laurent_column_reduction.jl")'
```

Expected: PASS with the `case_010 Laurent column reduction` testset and no failures.

- [ ] **Step 8: Run affected expert/internal checks**

Run:

```bash
julia --project=. -e 'include("test/internal/toricbuilder_case010_column_boundary.jl")'
julia --project=. -e 'include("test/expert/laurent_column_reduction_diagnostics.jl")'
```

Expected: both commands PASS with no failures.

- [ ] **Step 9: Run package verification**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: PASS for the package test entry point.

- [ ] **Step 10: Commit**

Run:

```bash
git add src/algorithm/column_reduction.jl test/fixtures/toricbuilder_case010_column_boundary.jl test/internal/toricbuilder_case010_column_boundary.jl test/expert/laurent_column_reduction_diagnostics.jl test/expert/case010_laurent_column_reduction.jl test/runtests.jl docs/superpowers/plans/2026-06-25-issue-134-case010-laurent-column-reduction.md
git commit -m "feat: reduce case010 Laurent column"
```

Expected: commit succeeds and contains the implementation, focused acceptance test, affected-test updates, and this plan.

---

## Self-Review

- Spec coverage: Task 1 covers the certified stage, replay metadata, tamper rejection, issue acceptance command, diagnostics update, boundary fixture update, and full package verification.
- Incomplete-marker scan: no open-ended implementation markers are used.
- Type consistency: the plan consistently uses `_reduce_via_laurent_unit_creation_certificate(column, R)`, `:laurent_unit_creation`, and the existing `ECPColumnReductionCertificate` replay flow.
