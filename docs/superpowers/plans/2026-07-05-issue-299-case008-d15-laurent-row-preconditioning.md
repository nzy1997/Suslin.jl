# Issue 299 Case008 D15 Laurent Row Preconditioning Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Support the `case_008 d=15` Laurent boundary column through a certified, replayable Laurent row-preconditioning stage.

**Architecture:** Generalize the existing `:laurent_elementary_row_preconditioning` stage in `src/algorithm/column_reduction.jl` so finite specs can resolve to one or more row-addition factors. Preserve the current d16 fixed-coefficient spec, add a d15 target-unit synthesis spec that recomputes coefficients from the input column by exact Laurent linear solve, and keep transformed-column reduction delegated to the existing Laurent base reducer to prevent recursive row preconditioning.

**Tech Stack:** Julia, Suslin.jl column-reduction internals, Oscar Laurent polynomial rings and matrices, Test stdlib.

## Global Constraints

- Repository has no `AGENTS.md`, `CLAUDE.md`, or `GEMINI.md`; follow existing Julia source and test style.
- Do not recognize `case_008` fixture ids in production code.
- Do not hard-code the #298 d15 coefficient list.
- Do not claim arbitrary Laurent unimodular-column support.
- Keep the stage named `:laurent_elementary_row_preconditioning`; do not add a `:case008_special_case` stage.
- Run row preconditioning only after the existing Laurent base path fails.
- Reduce the transformed column through `_reduce_laurent_unimodular_column_base_certificate`, not recursively through row preconditioning.
- Store enough stage data for `_ecp_replay_stage` to recompute and verify preconditioning factors, transformed column, nested certificate, composed factors, and final column.
- Update diagnostics and d15 boundary/profile tests from unsupported to supported.
- Required verification commands:
  - `julia --project=. -e 'include("test/expert/case008_d15_laurent_column_reduction.jl")'`
  - `julia --project=. -e 'include("test/internal/toricbuilder_case008_d15_column_boundary.jl")'`
  - `julia --project=. -e 'include("test/expert/laurent_column_reduction_diagnostics.jl")'`
  - `julia --project=. -e 'using Pkg; Pkg.test()'`

---

## File Structure

- Modify `src/algorithm/column_reduction.jl`: spec resolution, row-synthesis helpers, multi-factor preconditioning stage metadata, replay verification, diagnostics.
- Create `test/expert/case008_d15_laurent_column_reduction.jl`: focused d15 reduction/certificate/diagnostic/negative-control coverage.
- Modify `test/internal/toricbuilder_case008_d15_column_boundary.jl`: supported diagnostic expectations.
- Modify `test/expert/case008_d15_laurent_witness_profile.jl`: convert the historical unsupported profile into supported-profile coverage plus a synthetic unsupported control.
- Modify `test/expert/laurent_column_reduction_diagnostics.jl`: add d15 supported diagnostic assertions and keep unsupported synthetic controls.
- Modify `test/expert/case008_d16_laurent_column_reduction.jl` and `test/internal/toricbuilder_case008_d16_column_boundary.jl` only if the generalized stage metadata requires assertions to accept source/coefficients tuples.

---

### Task 1: Add Failing D15 Production Expert Test

**Files:**
- Create: `test/expert/case008_d15_laurent_column_reduction.jl`

**Interfaces:**
- Consumes: `ToricBuilderCase008D15ColumnBoundary.boundary_fixture()`, `Suslin.reduce_unimodular_column`, `Suslin.ecp_column_reduction_certificate`, `Suslin.verify_ecp_column_reduction`, and `Suslin.diagnose_unimodular_column_reduction`.
- Produces: a focused RED test for the currently unsupported d15 column-reduction boundary.

- [ ] **Step 1: Write the failing expert test**

Create `test/expert/case008_d15_laurent_column_reduction.jl` with helper functions matching the d16 test style:

```julia
using Test
using Suslin
using Oscar

include(joinpath(@__DIR__, "..", "fixtures", "toricbuilder_case008_d15_column_boundary.jl"))

function _case008_d15_reduction_product(factors, R, n::Int)
    product = identity_matrix(R, n)
    for factor in factors
        product *= factor
    end
    return product
end

function _case008_d15_apply_factors(factors, column, R)
    n = length(column)
    return _case008_d15_reduction_product(factors, R, n) * matrix(R, n, 1, collect(column))
end

function _case008_d15_target_column(R, n::Int)
    target = zero_matrix(R, n, 1)
    target[n, 1] = one(R)
    return target
end

function _case008_d15_tamper_first_factor(factors, R, n::Int)
    tampered = copy(factors)
    tampered[1] = identity_matrix(R, n)
    return tampered
end

function _case008_d15_tamper_certificate_first_factor(cert)
    tampered = _case008_d15_tamper_first_factor(
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

function _case008_d15_tamper_stage_coefficient(cert)
    stages = collect(cert.stages)
    stage_idx = findfirst(stage -> stage.kind == :laurent_elementary_row_preconditioning, stages)
    stage_idx === nothing && error("missing row-preconditioning stage")
    stage = stages[stage_idx]
    coefficients = collect(stage.coefficients)
    coefficients[1] += one(cert.ring)
    stages[stage_idx] = merge(stage, (; coefficient = coefficients[1], coefficients = tuple(coefficients...)))
    return Suslin.ECPColumnReductionCertificate(
        cert.original_column,
        cert.ring,
        tuple(stages...),
        cert.factors,
        cert.final_column,
        cert.verification,
    )
end

function _case008_d15_tamper_stage_source(cert)
    stages = collect(cert.stages)
    stage_idx = findfirst(stage -> stage.kind == :laurent_elementary_row_preconditioning, stages)
    stage_idx === nothing && error("missing row-preconditioning stage")
    stage = stages[stage_idx]
    source_indices = collect(stage.source_indices)
    source_indices[1] = source_indices[1] == 2 ? 3 : 2
    stages[stage_idx] = merge(stage, (; source_index = source_indices[1], source_indices = tuple(source_indices...)))
    return Suslin.ECPColumnReductionCertificate(
        cert.original_column,
        cert.ring,
        tuple(stages...),
        cert.factors,
        cert.final_column,
        cert.verification,
    )
end

function _case008_d15_diagnostic_stage_detail(diagnostic, stage::Symbol)
    hasproperty(diagnostic, :stage_details) || return nothing
    idx = findfirst(detail -> detail.stage == stage, diagnostic.stage_details)
    return idx === nothing ? nothing : diagnostic.stage_details[idx]
end

@testset "case_008 d=15 Laurent column reduction" begin
    fixture = ToricBuilderCase008D15ColumnBoundary.boundary_fixture()
    R = fixture.ring
    column = fixture.failing_column
    target = _case008_d15_target_column(R, length(column))

    factors = Suslin.reduce_unimodular_column(column, R)
    @test _case008_d15_apply_factors(factors, column, R) == target

    tampered_factors = _case008_d15_tamper_first_factor(factors, R, length(column))
    @test _case008_d15_apply_factors(tampered_factors, column, R) != target

    certificate = Suslin.ecp_column_reduction_certificate(column, R)
    @test Suslin.verify_ecp_column_reduction(certificate)
    @test any(stage -> stage.kind == :laurent_elementary_row_preconditioning, certificate.stages)
    @test !any(stage -> stage.kind == :case008_special_case, certificate.stages)
    @test _case008_d15_apply_factors(
        certificate.factors,
        certificate.original_column,
        certificate.ring,
    ) == target
    @test !Suslin.verify_ecp_column_reduction(
        _case008_d15_tamper_certificate_first_factor(certificate),
    )
    @test !Suslin.verify_ecp_column_reduction(
        _case008_d15_tamper_stage_coefficient(certificate),
    )
    @test !Suslin.verify_ecp_column_reduction(
        _case008_d15_tamper_stage_source(certificate),
    )

    preconditioning_stage = certificate.stages[end]
    @test preconditioning_stage.kind == :laurent_elementary_row_preconditioning
    @test preconditioning_stage.target_index == 1
    @test preconditioning_stage.source_indices == Tuple(2:15)
    @test length(preconditioning_stage.coefficients) == 14
    @test preconditioning_stage.coefficient_strategy == :target_unit_laurent_linear_synthesis
    @test preconditioning_stage.transformed_certificate.stages[end].kind == :unit_entry

    diagnostic = Suslin.diagnose_unimodular_column_reduction(column, R)
    @test diagnostic.status == :supported
    @test diagnostic.failure_code === nothing
    @test diagnostic.column_length == 15
    @test diagnostic.ring_profile.kind == :laurent_polynomial
    @test diagnostic.ring_profile.generators == ("u", "v")
    @test :laurent_elementary_row_preconditioning in diagnostic.attempted_stages
    @test !(:case008_special_case in diagnostic.attempted_stages)
    detail = _case008_d15_diagnostic_stage_detail(
        diagnostic,
        :laurent_elementary_row_preconditioning,
    )
    @test detail !== nothing
    @test detail.outcome == :supported
    @test detail.target_index == 1
    @test detail.source_indices == Tuple(2:15)
    @test detail.coefficient_strategy == :target_unit_laurent_linear_synthesis
    @test detail.coefficient_count == 14
    @test detail.transformed_stage == :unit_entry

    negative = ToricBuilderCase008D15ColumnBoundary.non_unimodular_negative_control(fixture)
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

- [ ] **Step 2: Run the expert test and verify RED**

Run:

```bash
julia --project=. -e 'include("test/expert/case008_d15_laurent_column_reduction.jl")'
```

Expected: FAIL with `ArgumentError: unsupported exact unimodular column reduction over this ring/column family`.

---

### Task 2: Generalize Laurent Row-Preconditioning Specs

**Files:**
- Modify: `src/algorithm/column_reduction.jl`

**Interfaces:**
- Consumes: `_reduce_laurent_unimodular_column_base_certificate`, `solve_laurent_linear`, `elementary_matrix`, `_factor_sequence_product`, and `_apply_reduction_factors`.
- Produces: `_laurent_row_preconditioning_candidate` that accepts fixed single-source specs and d15 target-unit row-synthesis specs.

- [ ] **Step 1: Replace `_laurent_row_preconditioning_specs` with finite spec families**

Implement the function so it returns:

```julia
function _laurent_row_preconditioning_specs(column::AbstractVector, R)
    _is_laurent_polynomial_ring(R) || return ()
    length(gens(R)) == 2 || return ()

    if length(column) == 16
        return ((;
            target_index = 1,
            source_indices = (10,),
            coefficient_strategy = :fixed_coefficients,
            coefficients = (one(R),),
            max_nonzero_coefficients = 1,
        ),)
    end

    if length(column) == 15
        return ((;
            target_index = 1,
            source_indices = Tuple(2:15),
            coefficient_strategy = :target_unit_laurent_linear_synthesis,
            coefficients = (),
            max_nonzero_coefficients = 14,
        ),)
    end

    return ()
end
```

- [ ] **Step 2: Add spec-resolution helpers near `_laurent_row_preconditioning_candidate`**

Add helpers:

```julia
function _laurent_row_preconditioning_solve_failure(err)::Bool
    return _laurent_witness_solve_failure(err)
end

function _laurent_row_preconditioning_fixed_coefficients(spec, R)
    return tuple((_coerce_into_ring(R, coeff, "row preconditioning coefficient") for coeff in spec.coefficients)...)
end

function _laurent_row_preconditioning_synthesis_coefficients(column::AbstractVector, R, target_idx::Int, source_indices)
    isempty(source_indices) && return nothing
    A = matrix(R, 1, length(source_indices), [column[idx] for idx in source_indices])
    B = matrix(R, 1, 1, [one(R) - column[target_idx]])
    solution = try
        solve_laurent_linear(A, B)
    catch err
        err isa InterruptException && rethrow()
        _laurent_row_preconditioning_solve_failure(err) && return nothing
        rethrow()
    end
    A * solution == B || return nothing
    return tuple((solution[idx, 1] for idx in 1:nrows(solution))...)
end

function _laurent_row_preconditioning_coefficients(column::AbstractVector, R, spec)
    target_idx = spec.target_index
    source_indices = tuple(spec.source_indices...)
    strategy = spec.coefficient_strategy
    if strategy == :fixed_coefficients
        length(spec.coefficients) == length(source_indices) || return nothing
        return _laurent_row_preconditioning_fixed_coefficients(spec, R)
    elseif strategy == :target_unit_laurent_linear_synthesis
        return _laurent_row_preconditioning_synthesis_coefficients(column, R, target_idx, source_indices)
    end
    return nothing
end
```

- [ ] **Step 3: Replace `_laurent_row_preconditioning_candidate`**

Rewrite the candidate builder to:

```julia
function _laurent_row_preconditioning_candidate(column::AbstractVector, R)
    n = length(column)
    for spec in _laurent_row_preconditioning_specs(column, R)
        target_idx = spec.target_index
        target_idx isa Integer || continue
        target_idx = Int(target_idx)
        1 <= target_idx <= n || continue

        source_indices = tuple((Int(source_idx) for source_idx in spec.source_indices)...)
        all(source_idx -> 1 <= source_idx <= n && source_idx != target_idx, source_indices) || continue
        length(unique(source_indices)) == length(source_indices) || continue

        coefficients = _laurent_row_preconditioning_coefficients(column, R, spec)
        coefficients === nothing && continue
        length(coefficients) == length(source_indices) || continue

        nonzero_pairs = [
            (source_idx, coeff)
            for (source_idx, coeff) in zip(source_indices, coefficients)
            if coeff != zero(R)
        ]
        isempty(nonzero_pairs) && continue
        length(nonzero_pairs) <= spec.max_nonzero_coefficients || continue

        accepted_source_indices = tuple((pair[1] for pair in nonzero_pairs)...)
        accepted_coefficients = tuple((pair[2] for pair in nonzero_pairs)...)
        precondition_factors = [
            elementary_matrix(n, target_idx, source_idx, coeff, R)
            for (source_idx, coeff) in nonzero_pairs
        ]
        precondition_factor = _factor_sequence_product(precondition_factors, R, n)
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
            source_index = accepted_source_indices[1],
            source_indices = accepted_source_indices,
            coefficient = accepted_coefficients[1],
            coefficients = accepted_coefficients,
            coefficient_strategy = spec.coefficient_strategy,
            precondition_factors,
            precondition_factor,
            transformed_column,
            transformed_certificate,
        )
    end

    return nothing
end
```

- [ ] **Step 4: Update stage construction**

Update `_reduce_via_laurent_elementary_row_preconditioning_certificate` so the stage stores `source_indices`, `coefficients`, `coefficient_strategy`, `precondition_factors`, and `precondition_factor`, and so composed factors use:

```julia
factors = _checked_reduction_factors(
    vcat(candidate.transformed_certificate.factors, candidate.precondition_factors),
    column,
    R,
    "Laurent elementary row-preconditioning reduction",
)
```

- [ ] **Step 5: Run the d15 expert test and verify GREEN for the new behavior**

Run:

```bash
julia --project=. -e 'include("test/expert/case008_d15_laurent_column_reduction.jl")'
```

Expected: either PASS or a replay/metadata failure that points to Task 3.

---

### Task 3: Harden Certificate Replay And Diagnostics

**Files:**
- Modify: `src/algorithm/column_reduction.jl`

**Interfaces:**
- Consumes: generalized stage metadata from Task 2.
- Produces: replay verification that rejects tampered d15 metadata and diagnostics that expose accepted preconditioning details.

- [ ] **Step 1: Update `_ecp_replay_stage` for `:laurent_elementary_row_preconditioning`**

Replace the required key tuple with:

```julia
(
    :kind,
    :input_column,
    :target_index,
    :source_index,
    :source_indices,
    :coefficient,
    :coefficients,
    :coefficient_strategy,
    :precondition_factors,
    :precondition_factor,
    :transformed_column,
    :transformed_certificate,
    :factors,
    :output_column,
)
```

Recompute the candidate with `_laurent_row_preconditioning_candidate(input_column, R)` and require the stored candidate metadata to equal the recomputed candidate:

```julia
candidate = _laurent_row_preconditioning_candidate(input_column, R)
candidate === nothing && return invalid_replay
candidate.target_index == target_idx || return invalid_replay
candidate.source_index == stage.source_index || return invalid_replay
candidate.source_indices == stage.source_indices || return invalid_replay
candidate.coefficient == stage.coefficient || return invalid_replay
candidate.coefficients == stage.coefficients || return invalid_replay
candidate.coefficient_strategy == stage.coefficient_strategy || return invalid_replay
```

Then recompute factors, transformed column, nested certificate, composed factors, and output column. Return `ok = true` only when every stored field matches.

- [ ] **Step 2: Update diagnostics detail for supported row preconditioning**

In `_diagnose_laurent_row_preconditioning`, include:

```julia
source_index = row_preconditioning.source_index,
source_indices = row_preconditioning.source_indices,
coefficient = row_preconditioning.coefficient,
coefficients = row_preconditioning.coefficients,
coefficient_strategy = row_preconditioning.coefficient_strategy,
coefficient_count = length(row_preconditioning.coefficients),
transformed_stage,
```

For the no-candidate detail, keep existing scalar fields and add tuple fields:

```julia
source_indices = (),
coefficients = (),
coefficient_strategy = nothing,
coefficient_count = 0,
```

- [ ] **Step 3: Run focused d15 expert verification**

Run:

```bash
julia --project=. -e 'include("test/expert/case008_d15_laurent_column_reduction.jl")'
```

Expected: PASS, including coefficient/source metadata tampering controls.

---

### Task 4: Update D15 Diagnostic Boundary Tests

**Files:**
- Modify: `test/internal/toricbuilder_case008_d15_column_boundary.jl`
- Modify: `test/expert/case008_d15_laurent_witness_profile.jl`
- Modify: `test/expert/laurent_column_reduction_diagnostics.jl`

**Interfaces:**
- Consumes: supported d15 diagnostics from Task 3.
- Produces: updated tests that no longer permanently assert the old unsupported d15 boundary.

- [ ] **Step 1: Update internal d15 boundary validator expectations**

Change `test/internal/toricbuilder_case008_d15_column_boundary.jl` diagnostic assertions to:

```julia
@test diagnostic.status == :supported
@test diagnostic.failure_code === nothing
@test diagnostic.column_length == 15
@test diagnostic.ring_profile.kind == :laurent_polynomial
@test diagnostic.ring_profile.generators == ("u", "v")
@test :laurent_elementary_row_preconditioning in diagnostic.attempted_stages
preconditioning_detail =
    _case008_d15_stage_detail(diagnostic, :laurent_elementary_row_preconditioning)
@test preconditioning_detail !== nothing
@test preconditioning_detail.outcome == :supported
@test preconditioning_detail.target_index == 1
@test preconditioning_detail.source_indices == Tuple(2:15)
@test preconditioning_detail.coefficient_strategy == :target_unit_laurent_linear_synthesis
@test preconditioning_detail.coefficient_count == 14
@test preconditioning_detail.transformed_stage == :unit_entry
```

Add helper `_case008_d15_stage_detail` if missing.

- [ ] **Step 2: Convert the d15 witness profile expert test**

Rename `_case008_d15_unsupported_profile` to `_case008_d15_supported_profile` and update it to assert `status == :supported`, `failure_code === nothing`, the same attempted stage tuple, supported row-preconditioning detail, and at least one supported detail. Keep the synthetic direct-unit Laurent profile control and add a synthetic unsupported control:

```julia
unsupported_column = [u + v, u + one(R), v + one(R)]
unsupported = Suslin.diagnose_unimodular_column_reduction(unsupported_column, R)
@test unsupported.status == :unsupported
@test !_case008_d15_supported_profile(unsupported)
```

If this synthetic column is unexpectedly supported, choose a length-15 two-generator Laurent column that `_laurent_row_preconditioning_candidate` returns `nothing` for, matching the existing d16 nonreducible spec style.

- [ ] **Step 3: Add d15 to laurent diagnostics expert test**

Include the d15 fixture in `test/expert/laurent_column_reduction_diagnostics.jl`, diagnose it next to the d16 block, and assert supported preconditioning detail:

```julia
case008_d15_fixture = ToricBuilderCase008D15ColumnBoundary.boundary_fixture()
case008_d15 = Suslin.diagnose_unimodular_column_reduction(
    case008_d15_fixture.failing_column,
    case008_d15_fixture.ring,
)
@test case008_d15.status == :supported
@test case008_d15.failure_code === nothing
@test case008_d15.column_length == 15
@test :laurent_elementary_row_preconditioning in case008_d15.attempted_stages
case008_d15_preconditioned =
    _diagnostic_stage_detail(case008_d15, :laurent_elementary_row_preconditioning)
@test case008_d15_preconditioned !== nothing
@test case008_d15_preconditioned.outcome == :supported
@test case008_d15_preconditioned.target_index == 1
@test case008_d15_preconditioned.source_indices == Tuple(2:15)
@test case008_d15_preconditioned.coefficient_strategy == :target_unit_laurent_linear_synthesis
@test case008_d15_preconditioned.coefficient_count == 14
@test case008_d15_preconditioned.transformed_stage == :unit_entry
```

- [ ] **Step 4: Run diagnostic verification commands**

Run:

```bash
julia --project=. -e 'include("test/internal/toricbuilder_case008_d15_column_boundary.jl")'
julia --project=. -e 'include("test/expert/case008_d15_laurent_witness_profile.jl")'
julia --project=. -e 'include("test/expert/laurent_column_reduction_diagnostics.jl")'
```

Expected: all exit 0.

---

### Task 5: Regression Verification, Review, And Commit

**Files:**
- All changed files from Tasks 1-4.

**Interfaces:**
- Consumes: implementation and tests.
- Produces: committed branch ready for final workflow.

- [ ] **Step 1: Run the issue verification commands**

Run:

```bash
julia --project=. -e 'include("test/expert/case008_d15_laurent_column_reduction.jl")'
julia --project=. -e 'include("test/internal/toricbuilder_case008_d15_column_boundary.jl")'
julia --project=. -e 'include("test/expert/laurent_column_reduction_diagnostics.jl")'
```

Expected: all exit 0.

- [ ] **Step 2: Run the required package test suite**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: exit 0.

- [ ] **Step 3: Inspect diff and status**

Run:

```bash
git status --short --branch
git diff --stat HEAD
git diff --check
```

Expected: only intended source, test, and Superpowers plan/doc files changed; diff check exits 0.

- [ ] **Step 4: Commit implementation**

Run:

```bash
git add src/algorithm/column_reduction.jl test/expert/case008_d15_laurent_column_reduction.jl test/internal/toricbuilder_case008_d15_column_boundary.jl test/expert/case008_d15_laurent_witness_profile.jl test/expert/laurent_column_reduction_diagnostics.jl docs/superpowers/plans/2026-07-05-issue-299-case008-d15-laurent-row-preconditioning.md
git commit -m "feat: support case008 d15 Laurent row preconditioning"
```

- [ ] **Step 5: Request code review**

Use `superpowers:requesting-code-review` or the available multi-agent review tool. Fix any Critical or Important findings and rerun affected tests.

- [ ] **Step 6: Use finishing workflow**

Use `superpowers:verification-before-completion` and `superpowers:finishing-a-development-branch`. Under the Standing Answer Policy, choose "Push and create a Pull Request".

## Self-Review

- Spec coverage: Tasks 1-4 cover reducer support, certificate replay, diagnostics, focused expert test, internal boundary update, stale unsupported profile update, and negative controls.
- Placeholder scan: no TODO/TBD placeholders are intentionally left for implementers.
- Type consistency: all stage metadata names are consistent across source replay, diagnostics, and tests.
