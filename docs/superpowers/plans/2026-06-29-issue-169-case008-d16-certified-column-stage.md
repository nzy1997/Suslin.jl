# Issue 169 Case008 D16 Certified Column Stage Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Reduce the `case_008` `d=16` Laurent column through a certified row-preconditioned elementary stage.

**Architecture:** Add a Laurent-only row-preconditioning stage that composes one elementary row addition with an existing verified transformed-column reduction. Certificate replay recomputes the precondition factor and transformed certificate exactly.

**Tech Stack:** Julia, Oscar, existing Suslin column reduction and ECP certificate helpers.

## Global Constraints

- Public inputs remain `reduce_unimodular_column(fixture.failing_column, fixture.ring)` and `ecp_column_reduction_certificate(fixture.failing_column, fixture.ring)`.
- The stored stage name must be algebraic, not `:case008_special_case`.
- Produce replayable elementary factors and certificate verification.
- Diagnostics for the supported column must report `status == :supported`.
- The negative control must keep unsupported non-unimodular inputs on precondition failure, not the new stage.
- Do not claim full `case_008` certificate success unless the bounded report verifies it.
- Do not broaden to arbitrary Laurent columns without an exact predicate-backed test.

---

## File Structure

- `test/expert/case008_d16_laurent_column_reduction.jl`: new red acceptance test with replay, diagnostics, and tamper controls.
- `test/expert/laurent_column_reduction_diagnostics.jl`: add d=16 supported diagnostic coverage while preserving existing unsupported coverage.
- `test/fixtures/toricbuilder_case008_d16_column_boundary.jl`: update expected diagnostic metadata once the reducer is green.
- `test/internal/toricbuilder_case008_d16_column_boundary.jl`: update fixture validator expectations once supported.
- `test/runtests.jl`: add the new expert test to the expert group.
- `src/algorithm/column_reduction.jl`: add the certified stage and replay support.

---

### Task 1: Add Red D16 Acceptance Coverage

**Files:**
- Create: `test/expert/case008_d16_laurent_column_reduction.jl`
- Modify: `test/expert/laurent_column_reduction_diagnostics.jl`
- Modify: `test/runtests.jl`

**Interfaces:**
- Consumes: `ToricBuilderCase008D16ColumnBoundary.boundary_fixture()`, `Suslin.reduce_unimodular_column`, `Suslin.ecp_column_reduction_certificate`, `Suslin.verify_ecp_column_reduction`, `Suslin.diagnose_unimodular_column_reduction`.
- Produces: failing acceptance tests for the d=16 stage and diagnostics.

- [ ] **Step 1: Write the failing expert acceptance test**

Create `test/expert/case008_d16_laurent_column_reduction.jl` with helpers matching the d=21 test shape:

```julia
using Test
using Suslin
using Oscar

include(joinpath(@__DIR__, "..", "fixtures", "toricbuilder_case008_d16_column_boundary.jl"))

function _case008_d16_reduction_product(factors, R, n::Int)
    product = identity_matrix(R, n)
    for factor in factors
        product *= factor
    end
    return product
end

function _case008_d16_apply_factors(factors, column, R)
    n = length(column)
    return _case008_d16_reduction_product(factors, R, n) * matrix(R, n, 1, collect(column))
end

function _case008_d16_target_column(R, n::Int)
    target = zero_matrix(R, n, 1)
    target[n, 1] = one(R)
    return target
end

function _case008_d16_tamper_first_factor(factors, R, n::Int)
    tampered = copy(factors)
    tampered[1] = identity_matrix(R, n)
    return tampered
end

function _case008_d16_tamper_certificate_first_factor(cert)
    tampered = _case008_d16_tamper_first_factor(
        cert.factors,
        cert.ring,
        length(cert.original_column),
    )
    return Suslin.ECPColumnReductionCertificate(
        cert.original_column,
        cert.ring,
        cert.stages,
        tampered,
        cert.final_column,
        cert.verification,
    )
end

function _case008_d16_tamper_stage_coefficient(cert)
    stages = collect(cert.stages)
    stage_idx = findfirst(stage -> stage.kind == :laurent_elementary_row_preconditioning, stages)
    stage_idx === nothing && error("missing row-preconditioning stage")
    stage = stages[stage_idx]
    stages[stage_idx] = merge(stage, (; coefficient = zero(cert.ring)))
    return Suslin.ECPColumnReductionCertificate(
        cert.original_column,
        cert.ring,
        tuple(stages...),
        cert.factors,
        cert.final_column,
        cert.verification,
    )
end

@testset "case_008 d=16 Laurent column reduction" begin
    fixture = ToricBuilderCase008D16ColumnBoundary.boundary_fixture()
    R = fixture.ring
    column = fixture.failing_column
    target = _case008_d16_target_column(R, length(column))

    factors = Suslin.reduce_unimodular_column(column, R)
    @test _case008_d16_apply_factors(factors, column, R) == target

    tampered_factors = _case008_d16_tamper_first_factor(factors, R, length(column))
    @test _case008_d16_apply_factors(tampered_factors, column, R) != target

    certificate = Suslin.ecp_column_reduction_certificate(column, R)
    @test Suslin.verify_ecp_column_reduction(certificate)
    @test any(stage -> stage.kind == :laurent_elementary_row_preconditioning, certificate.stages)
    @test _case008_d16_apply_factors(
        certificate.factors,
        certificate.original_column,
        certificate.ring,
    ) == target
    @test !Suslin.verify_ecp_column_reduction(
        _case008_d16_tamper_certificate_first_factor(certificate),
    )
    @test !Suslin.verify_ecp_column_reduction(
        _case008_d16_tamper_stage_coefficient(certificate),
    )

    diagnostic = Suslin.diagnose_unimodular_column_reduction(column, R)
    @test diagnostic.status == :supported
    @test diagnostic.failure_code === nothing
    @test diagnostic.column_length == 16
    @test diagnostic.ring_profile.kind == :laurent_polynomial
    @test diagnostic.ring_profile.generators == ("u", "v")
    @test :laurent_elementary_row_preconditioning in diagnostic.attempted_stages

    negative = ToricBuilderCase008D16ColumnBoundary.non_unimodular_negative_control(fixture)
    negative_diagnostic = Suslin.diagnose_unimodular_column_reduction(
        negative.failing_column,
        negative.ring,
    )
    @test negative_diagnostic.status == :precondition_failed
    @test negative_diagnostic.failure_code == :not_unimodular
    @test isempty(negative_diagnostic.attempted_stages)
    @test_throws ArgumentError Suslin.reduce_unimodular_column(
        negative.failing_column,
        negative.ring,
    )
end
```

- [ ] **Step 2: Add diagnostics coverage**

In `test/expert/laurent_column_reduction_diagnostics.jl`, include the d=16 fixture and assert the new stage detail:

```julia
include(joinpath(@__DIR__, "..", "fixtures", "toricbuilder_case008_d16_column_boundary.jl"))
```

Inside the Laurent diagnostics testset, after the d=21 assertions, add:

```julia
    case008_d16_fixture = ToricBuilderCase008D16ColumnBoundary.boundary_fixture()
    case008_d16 = Suslin.diagnose_unimodular_column_reduction(
        case008_d16_fixture.failing_column,
        case008_d16_fixture.ring,
    )
    @test case008_d16.status == :supported
    @test case008_d16.failure_code === nothing
    @test case008_d16.column_length == 16
    @test :laurent_elementary_row_preconditioning in case008_d16.attempted_stages
    case008_d16_preconditioned =
        _diagnostic_stage_detail(case008_d16, :laurent_elementary_row_preconditioning)
    @test case008_d16_preconditioned !== nothing
    @test case008_d16_preconditioned.outcome == :supported
    @test case008_d16_preconditioned.target_index isa Integer
    @test case008_d16_preconditioned.source_index isa Integer
    @test case008_d16_preconditioned.coefficient == one(case008_d16_fixture.ring)
```

- [ ] **Step 3: Register the expert test**

Add `"expert/case008_d16_laurent_column_reduction.jl"` to the expert list in
`test/runtests.jl`.

- [ ] **Step 4: Run RED**

Run:

```bash
julia --project=. -e 'include("test/expert/case008_d16_laurent_column_reduction.jl")'
```

Expected: FAIL because `reduce_unimodular_column` still throws unsupported exact unimodular column reduction for the d=16 column.

- [ ] **Step 5: Commit**

```bash
git add test/expert/case008_d16_laurent_column_reduction.jl test/expert/laurent_column_reduction_diagnostics.jl test/runtests.jl
git commit -m "test: cover case008 d16 laurent reduction"
```

---

### Task 2: Implement Certified Row-Preconditioned Laurent Stage

**Files:**
- Modify: `src/algorithm/column_reduction.jl`
- Modify: `test/fixtures/toricbuilder_case008_d16_column_boundary.jl`
- Modify: `test/internal/toricbuilder_case008_d16_column_boundary.jl`

**Interfaces:**
- Consumes: red tests from Task 1 and existing ECP certificate replay helpers.
- Produces: `:laurent_elementary_row_preconditioning` reducer stage whose factors reduce the d=16 column to `e_16`.

- [ ] **Step 1: Split the Laurent base reducer**

In `src/algorithm/column_reduction.jl`, keep the existing direct Laurent logic in a helper that does not recursively try the new stage:

```julia
function _reduce_laurent_unimodular_column_base_certificate(column::AbstractVector, R)
    unit_idx = findfirst(is_unit, column)
    unit_idx !== nothing && return _unit_entry_reduction_certificate_stage(column, unit_idx, R)

    unit_creation = _reduce_via_laurent_unit_creation_certificate(column, R)
    unit_creation !== nothing && return unit_creation

    witness_unit = _reduce_via_laurent_witness_unit_certificate(column, R)
    witness_unit !== nothing && return witness_unit

    normalization = normalize_laurent_object(column)
    poly_column = normalization.normalized_object
    P = normalization.metadata.polynomial_ring
    is_unimodular_column(poly_column, P) || return nothing

    poly_result = _reduce_polynomial_unimodular_column_exact_certificate(poly_column, P)
    poly_result === nothing && return nothing

    polynomial_certificate = _ecp_certificate_from_stage(poly_column, P, poly_result.stage)
    lifted_factors = [_lift_polynomial_reduction_factor(factor, R) for factor in polynomial_certificate.factors]
    shift = only(normalization.metadata.shift_monomials)
    inverse_shift = only(normalization.metadata.inverse_shift_monomials)
    normalization_factors = _unit_normalization_factors(length(column), inverse_shift, shift, R)
    factors = _checked_reduction_factors(
        vcat(normalization_factors, lifted_factors),
        column,
        R,
        "Laurent normalization reduction",
    )

    stage = (;
        kind = :laurent_normalization,
        input_column = _ecp_column_tuple(column),
        normalization,
        normalized_column = _ecp_column_tuple(poly_column),
        polynomial_certificate,
        lifted_factors,
        shift,
        inverse_shift,
        normalization_factors,
        factors,
        output_column = _apply_reduction_factors(factors, column, R),
    )
    return (; factors, stage)
end
```

Then make `_reduce_laurent_unimodular_column_certificate` call the base helper first, then the new stage:

```julia
function _reduce_laurent_unimodular_column_certificate(column::AbstractVector, R)
    base = _reduce_laurent_unimodular_column_base_certificate(column, R)
    base !== nothing && return base

    row_preconditioned = _reduce_via_laurent_elementary_row_preconditioning_certificate(column, R)
    row_preconditioned !== nothing && return row_preconditioned

    return nothing
end
```

- [ ] **Step 2: Add row-preconditioning candidate search and stage construction**

Add these helpers near the Laurent reducer functions:

```julia
_laurent_row_preconditioning_coefficients(R) = (one(R),)

function _laurent_row_preconditioning_candidate(column::AbstractVector, R)
    _is_laurent_polynomial_ring(R) || return nothing
    length(column) >= 6 || return nothing

    n = length(column)
    for coeff in _laurent_row_preconditioning_coefficients(R)
        coeff == zero(R) && continue
        for target_idx in 1:n, source_idx in 1:n
            target_idx == source_idx && continue
            precondition_factor = elementary_matrix(n, target_idx, source_idx, coeff, R)
            transformed_column = collect(_ecp_matrix_column_to_tuple(
                precondition_factor * matrix(R, n, 1, collect(column)),
            ))
            transformed_result =
                _reduce_laurent_unimodular_column_base_certificate(transformed_column, R)
            transformed_result === nothing && continue
            transformed_certificate =
                _ecp_certificate_from_stage(transformed_column, R, transformed_result.stage)
            return (;
                target_index = target_idx,
                source_index = source_idx,
                coefficient = coeff,
                precondition_factor,
                transformed_column,
                transformed_certificate,
            )
        end
    end

    return nothing
end

function _reduce_via_laurent_elementary_row_preconditioning_certificate(column::AbstractVector, R)
    candidate = _laurent_row_preconditioning_candidate(column, R)
    candidate === nothing && return nothing

    factors = _checked_reduction_factors(
        vcat(candidate.transformed_certificate.factors, [candidate.precondition_factor]),
        column,
        R,
        "Laurent elementary row-preconditioning reduction",
    )
    stage = (;
        kind = :laurent_elementary_row_preconditioning,
        input_column = _ecp_column_tuple(column),
        target_index = candidate.target_index,
        source_index = candidate.source_index,
        coefficient = candidate.coefficient,
        precondition_factor = candidate.precondition_factor,
        transformed_column = tuple(candidate.transformed_column...),
        transformed_certificate = candidate.transformed_certificate,
        factors,
        output_column = _apply_reduction_factors(factors, column, R),
    )
    return (; factors, stage)
end
```

- [ ] **Step 3: Add exact stage replay**

Add an `_ecp_replay_stage` branch for `stage.kind == :laurent_elementary_row_preconditioning` that:

```julia
    elseif stage.kind == :laurent_elementary_row_preconditioning
        invalid_replay = (; ok = false, factors = Any[], output_column = matrix(R, length(input_column), 1, collect(input_column)))
        _ecp_stage_keys_ok(
            stage,
            (:kind, :input_column, :target_index, :source_index, :coefficient, :precondition_factor, :transformed_column, :transformed_certificate, :factors, :output_column),
        ) || return invalid_replay
        _is_laurent_polynomial_ring(R) || return invalid_replay

        n = length(input_column)
        target_idx = stage.target_index
        source_idx = stage.source_index
        target_idx isa Integer && source_idx isa Integer || return invalid_replay
        1 <= target_idx <= n || return invalid_replay
        1 <= source_idx <= n || return invalid_replay
        target_idx == source_idx && return invalid_replay
        stage.coefficient in _laurent_row_preconditioning_coefficients(R) || return invalid_replay

        precondition_factor = elementary_matrix(n, target_idx, source_idx, stage.coefficient, R)
        transformed_matrix = precondition_factor * matrix(R, n, 1, collect(input_column))
        transformed_column = collect(_ecp_matrix_column_to_tuple(transformed_matrix))
        transformed_certificate_ok =
            verify_ecp_column_reduction(stage.transformed_certificate) &&
            stage.transformed_certificate.original_column == transformed_column &&
            stage.transformed_certificate.ring == R
        expected_factors = transformed_certificate_ok ?
            vcat(stage.transformed_certificate.factors, [precondition_factor]) :
            Any[]
        expected_output = transformed_certificate_ok ?
            _apply_reduction_factors(expected_factors, input_column, R) :
            matrix(R, n, 1, collect(input_column))
        ok = stage.input_column == _ecp_column_tuple(input_column) &&
            stage.precondition_factor == precondition_factor &&
            stage.transformed_column == tuple(transformed_column...) &&
            transformed_certificate_ok &&
            _ecp_factor_sequences_equal(stage.factors, expected_factors) &&
            stage.output_column == expected_output
        return (; ok, factors = expected_factors, output_column = expected_output)
```

- [ ] **Step 4: Add diagnostics for the new stage**

In `_diagnose_laurent_unimodular_column_reduction`, after the witness-unit
attempt and before normalization, call
`_laurent_row_preconditioning_candidate(column, R)`. On success, push a
`:laurent_elementary_row_preconditioning` detail with `outcome = :supported`,
indices, coefficient, and transformed stage kind, then return supported. On
failure, push `outcome = :no_row_preconditioning_candidate`.

- [ ] **Step 5: Update d=16 fixture metadata and internal test**

In `test/fixtures/toricbuilder_case008_d16_column_boundary.jl`, replace
`EXPECTED_DIAGNOSTIC.status` with `:supported`, `failure_code` with `nothing`,
and add fields for the row-preconditioning outcome. Update
`_diagnostic_stage_profile_matches` to check the new detail.

In `test/internal/toricbuilder_case008_d16_column_boundary.jl`, update assertions
that expected the old unsupported boundary to expect supported diagnostics and
the new algebraic stage.

- [ ] **Step 6: Run GREEN**

Run:

```bash
julia --project=. -e 'include("test/expert/case008_d16_laurent_column_reduction.jl")'
julia --project=. -e 'include("test/expert/laurent_column_reduction_diagnostics.jl")'
julia --project=. -e 'include("test/internal/toricbuilder_case008_d16_column_boundary.jl")'
```

Expected: all commands exit 0.

- [ ] **Step 7: Commit**

```bash
git add src/algorithm/column_reduction.jl test/fixtures/toricbuilder_case008_d16_column_boundary.jl test/internal/toricbuilder_case008_d16_column_boundary.jl
git commit -m "feat: reduce case008 d16 laurent column"
```

---

### Task 3: Run Required Verification And Report Boundary

**Files:**
- No source edits expected.
- Output: `/tmp/qblock-case008-after-d16.md`

**Interfaces:**
- Consumes: completed reducer from Tasks 1-2.
- Produces: verification evidence for the PR.

- [ ] **Step 1: Run targeted expert acceptance**

```bash
julia --project=. -e 'include("test/expert/case008_d16_laurent_column_reduction.jl")'
```

Expected: exit 0 and exact reduction to `e_16`.

- [ ] **Step 2: Run bounded report**

```bash
julia --project=. scripts/report_toricbuilder_cache_q_blocks.jl --exercise=case_008 --timeout-seconds=180 --output=/tmp/qblock-case008-after-d16.md
```

Expected: exit 0. Inspect `/tmp/qblock-case008-after-d16.md` and confirm it no longer stops at the d=16 `unsupported_laurent_column_family` boundary. It may reach a later structured boundary or pass the larger route.

- [ ] **Step 3: Run default package tests**

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: exit 0.

- [ ] **Step 4: Commit any verification-only metadata**

If no files changed, do not commit. If fixture metadata or docs need a small
correction found by verification, commit only those files.
