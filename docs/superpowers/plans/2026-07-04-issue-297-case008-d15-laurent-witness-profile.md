# Issue 297 Case008 D15 Laurent Witness Profile Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an expert-only profile test that records the current unsupported Laurent reducer boundary for the `case_008 d=15` fixture.

**Architecture:** Add one focused expert test file beside the existing `d=16` witness profile. The test uses `Suslin.diagnose_unimodular_column_reduction(column, R)` as the source of truth, keeps all helper logic local, and asserts both the unsupported fixture profile and a cheap supported direct-unit negative control.

**Tech Stack:** Julia, Test stdlib, Oscar, Suslin diagnostics, existing ToricBuilder fixture modules.

## Global Constraints

- Do not modify `src/`.
- Do not implement a reducer repair.
- Do not add this expert profile to the default `test/runtests.jl` suite.
- Use fixture module `test/fixtures/toricbuilder_case008_d15_column_boundary.jl`.
- Create expert test path `test/expert/case008_d15_laurent_witness_profile.jl`.
- The unsupported fixture must assert status `:unsupported`, failure code `:unsupported_laurent_column_family`, and column length `15`.
- The attempted stages must be exactly `(:unit_entry, :laurent_unit_creation, :laurent_witness_unit, :laurent_normalization, :laurent_elementary_row_preconditioning)`.
- The witness detail must report outcome `:witness_without_unit` and `witness_unit_index === nothing`.
- The normalization detail must report outcome `:normalized_not_unimodular`, `normalized_column_length == 15`, `normalized_status == :precondition_failed`, and `normalized_failure_code == :not_unimodular`.
- The row-preconditioning detail must report outcome `:no_row_preconditioning_candidate`.
- Required focused verification command: `julia --project=. -e 'include("test/expert/case008_d15_laurent_witness_profile.jl")'`.
- Required package verification command: `julia --project=. -e 'using Pkg; Pkg.test()'`.

---

### Task 1: Add Expert Diagnostic Profile Test

**Files:**
- Create: `test/expert/case008_d15_laurent_witness_profile.jl`

**Interfaces:**
- Consumes: `ToricBuilderCase008D15ColumnBoundary.boundary_fixture() -> NamedTuple`
- Consumes: `Suslin.diagnose_unimodular_column_reduction(column, R) -> diagnostic`
- Produces: local `_diagnostic_stage_detail(diagnostic, stage::Symbol)`
- Produces: local `_case008_d15_unsupported_profile(diagnostic)::Bool`

- [ ] **Step 1: Verify the expert test is absent**

Run:

```bash
test ! -e test/expert/case008_d15_laurent_witness_profile.jl
```

Expected: command exits `0`.

- [ ] **Step 2: Create the expert test file**

Create `test/expert/case008_d15_laurent_witness_profile.jl` with this content:

```julia
using Test
using Suslin
using Oscar

include(joinpath(@__DIR__, "..", "fixtures", "toricbuilder_case008_d15_column_boundary.jl"))

const CASE008_D15_EXPECTED_ATTEMPTED_STAGES = (
    :unit_entry,
    :laurent_unit_creation,
    :laurent_witness_unit,
    :laurent_normalization,
    :laurent_elementary_row_preconditioning,
)

function _diagnostic_stage_detail(diagnostic, stage::Symbol)
    idx = findfirst(detail -> detail.stage == stage, diagnostic.stage_details)
    return idx === nothing ? nothing : diagnostic.stage_details[idx]
end

function _case008_d15_unsupported_profile(diagnostic)::Bool
    diagnostic.status == :unsupported || return false
    diagnostic.failure_code == :unsupported_laurent_column_family || return false
    diagnostic.column_length == 15 || return false
    diagnostic.attempted_stages == CASE008_D15_EXPECTED_ATTEMPTED_STAGES || return false
    hasproperty(diagnostic, :stage_details) || return false
    length(diagnostic.stage_details) == length(CASE008_D15_EXPECTED_ATTEMPTED_STAGES) || return false
    any(detail -> detail.outcome == :supported, diagnostic.stage_details) && return false

    unit_entry = _diagnostic_stage_detail(diagnostic, :unit_entry)
    unit_entry !== nothing || return false
    unit_entry.outcome == :no_unit_entry || return false
    unit_entry.pivot_index === nothing || return false

    unit_creation = _diagnostic_stage_detail(diagnostic, :laurent_unit_creation)
    unit_creation !== nothing || return false
    unit_creation.outcome == :no_unit_creation_candidate || return false

    witness = _diagnostic_stage_detail(diagnostic, :laurent_witness_unit)
    witness !== nothing || return false
    witness.outcome == :witness_without_unit || return false
    witness.witness_unit_index === nothing || return false

    normalization = _diagnostic_stage_detail(diagnostic, :laurent_normalization)
    normalization !== nothing || return false
    normalization.outcome == :normalized_not_unimodular || return false
    normalization.normalized_column_length == 15 || return false
    normalization.normalized_status == :precondition_failed || return false
    normalization.normalized_failure_code == :not_unimodular || return false

    row_preconditioning =
        _diagnostic_stage_detail(diagnostic, :laurent_elementary_row_preconditioning)
    row_preconditioning !== nothing || return false
    row_preconditioning.outcome == :no_row_preconditioning_candidate || return false
    row_preconditioning.target_index === nothing || return false
    row_preconditioning.source_index === nothing || return false
    row_preconditioning.coefficient === nothing || return false
    row_preconditioning.transformed_stage === nothing || return false

    return true
end

@testset "case_008 d=15 Laurent witness profile" begin
    fixture = ToricBuilderCase008D15ColumnBoundary.boundary_fixture()
    @test ToricBuilderCase008D15ColumnBoundary.validate_boundary_fixture(fixture) == :ok

    diagnostic = Suslin.diagnose_unimodular_column_reduction(
        fixture.failing_column,
        fixture.ring,
    )

    @test diagnostic.status == :unsupported
    @test diagnostic.failure_code == :unsupported_laurent_column_family
    @test diagnostic.column_length == 15
    @test diagnostic.ring_profile.kind == :laurent_polynomial
    @test diagnostic.ring_profile.generators == ("u", "v")
    @test diagnostic.attempted_stages == CASE008_D15_EXPECTED_ATTEMPTED_STAGES
    @test hasproperty(diagnostic, :stage_details)
    @test diagnostic.stage_details isa Tuple
    @test all(detail -> detail isa NamedTuple, diagnostic.stage_details)
    @test length(diagnostic.stage_details) == length(CASE008_D15_EXPECTED_ATTEMPTED_STAGES)
    @test !any(detail -> detail.outcome == :supported, diagnostic.stage_details)

    unit_entry = _diagnostic_stage_detail(diagnostic, :unit_entry)
    @test unit_entry !== nothing
    @test unit_entry.outcome == :no_unit_entry
    @test unit_entry.pivot_index === nothing

    unit_creation = _diagnostic_stage_detail(diagnostic, :laurent_unit_creation)
    @test unit_creation !== nothing
    @test unit_creation.outcome == :no_unit_creation_candidate

    witness = _diagnostic_stage_detail(diagnostic, :laurent_witness_unit)
    @test witness !== nothing
    @test witness.outcome == :witness_without_unit
    @test witness.witness_unit_index === nothing

    normalization = _diagnostic_stage_detail(diagnostic, :laurent_normalization)
    @test normalization !== nothing
    @test normalization.outcome == :normalized_not_unimodular
    @test normalization.normalized_column_length == 15
    @test normalization.normalized_ring_kind == :polynomial
    @test normalization.normalized_status == :precondition_failed
    @test normalization.normalized_failure_code == :not_unimodular

    row_preconditioning =
        _diagnostic_stage_detail(diagnostic, :laurent_elementary_row_preconditioning)
    @test row_preconditioning !== nothing
    @test row_preconditioning.outcome == :no_row_preconditioning_candidate
    @test row_preconditioning.target_index === nothing
    @test row_preconditioning.source_index === nothing
    @test row_preconditioning.coefficient === nothing
    @test row_preconditioning.transformed_stage === nothing

    @test _case008_d15_unsupported_profile(diagnostic)
end

@testset "synthetic direct-unit Laurent profile control" begin
    R, (u, v) = Suslin.suslin_laurent_polynomial_ring(GF(2), ["u", "v"])
    direct_unit_column = [one(R), u + v, zero(R)]

    diagnostic = Suslin.diagnose_unimodular_column_reduction(direct_unit_column, R)

    @test diagnostic.status == :supported
    @test diagnostic.failure_code === nothing
    @test diagnostic.column_length == 3
    @test diagnostic.attempted_stages == (:unit_entry,)
    @test !_case008_d15_unsupported_profile(diagnostic)

    unit_entry = _diagnostic_stage_detail(diagnostic, :unit_entry)
    @test unit_entry !== nothing
    @test unit_entry.outcome == :supported
    @test unit_entry.pivot_index == 1
end
```

- [ ] **Step 3: Run the focused expert test**

Run:

```bash
julia --project=. -e 'include("test/expert/case008_d15_laurent_witness_profile.jl")'
```

Expected: command exits `0` and both testsets pass.

- [ ] **Step 4: Run the package test entry point**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: command exits `0`.

- [ ] **Step 5: Commit implementation**

Run:

```bash
git add test/expert/case008_d15_laurent_witness_profile.jl docs/superpowers/plans/2026-07-04-issue-297-case008-d15-laurent-witness-profile.md
git commit -m "test: profile case008 d15 Laurent witness boundary"
```

Expected: commit succeeds.

## Plan Self-Review

- Spec coverage: Task 1 creates the expert-only test, uses the #295 fixture,
  asserts all requested diagnostic fields, and includes the supported direct-unit
  negative control.
- Placeholder scan: no unfinished placeholders remain.
- Type consistency: helper names, diagnostic field names, and expected symbols
  match the probe output from `Suslin.diagnose_unimodular_column_reduction`.
