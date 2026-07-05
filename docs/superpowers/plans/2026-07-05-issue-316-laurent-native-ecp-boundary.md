# Issue 316 Laurent Native ECP Boundary Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Expose a structured terminal `:laurent_native_ecp_boundary` diagnostic for validated Laurent unimodular columns that exhaust existing Laurent reduction stages.

**Architecture:** Keep certificate reduction unchanged and modify only diagnostic reporting. The Laurent diagnostic path appends one terminal diagnostic-only stage after `:laurent_elementary_row_preconditioning` declines; tests verify that validated d14 reaches the new boundary, d15 remains supported, and precondition failures stop before stage attempts. The validated d14 test explicitly opts into large-support witness and normalized-delegation declines with structured diagnostic data instead of running unsupported solver/ideal-membership paths; default diagnostics still run existing supported stages.

**Tech Stack:** Julia, Oscar, Suslin internal reducer diagnostics, Test stdlib.

## Global Constraints

- New terminal stage symbol must be `:laurent_native_ecp_boundary`.
- Terminal stage detail must include `outcome = :staged_boundary`.
- Terminal stage detail must include `boundary = :laurent_native_ecp`.
- Terminal stage detail must include `requires_descent_measure = true`.
- Terminal stage detail must include `requires_link_witness = true`.
- Terminal stage detail must include `requires_endpoint_reduction = true`.
- Terminal stage detail must include `requires_laurent_normality_replay = true`.
- Terminal stage detail must include `requires_recursive_peel_integration = true`.
- Terminal stage detail must include `fallback_policy = :diagnostic_only`.
- Keep existing supported stages unchanged: `:unit_entry`, `:laurent_unit_creation`, `:laurent_witness_unit`, `:laurent_normalization`, delegated ordinary-polynomial stages where applicable, and `:laurent_elementary_row_preconditioning`.
- Do not route Laurent columns into `_reduce_via_general_ecp_pipeline_certificate`.
- Do not implement Laurent link witnesses, endpoint reductions, normality/conjugation replay, or recursive Laurent peel integration.
- Do not claim arbitrary Laurent `GL_n` support.
- Do not make diagonal monomial balancing or polynomialization the primary route.
- Add `assume_unimodular = true` only as an internal diagnostic keyword for already-validated fixtures; the default must still run unimodularity preconditions.
- Use an explicit `laurent_large_support_diagnostic_decline = true` opt-in with a conservative support term limit of 1000 so validated d14 records a large-support decline while default diagnostics keep existing supported routes.
- The supported `case_008 d=15` fixture must remain supported by `:laurent_elementary_row_preconditioning`.
- A non-unimodular Laurent column must fail with `:not_unimodular` before any Laurent-native ECP boundary stage is attempted.
- A large-support Laurent column with a unit witness must remain supported by `:laurent_witness_unit` by default and must not report the Laurent-native ECP boundary.
- The issue verification commands are `julia --project=. -e 'include("test/expert/laurent_native_ecp_boundary_diagnostics.jl")'` and `julia --project=. -e 'include("test/expert/laurent_column_reduction_diagnostics.jl")'`.
- The package verification command is `julia --project=. -e 'using Pkg; Pkg.test()'`.

---

## File Structure

- Modify: `src/algorithm/column_reduction.jl`
  - Owns the terminal Laurent-native ECP diagnostic stage helper and appends it after Laurent row preconditioning declines.
- Create: `test/expert/laurent_native_ecp_boundary_diagnostics.jl`
  - Owns the new d14 boundary contract, d15 supported negative control, and non-unimodular precondition negative control.
- Modify: `test/expert/laurent_column_reduction_diagnostics.jl`
  - Extends the existing unsupported Laurent diagnostic test to assert the new terminal stage and detail while keeping existing supported checks stable.
- Modify: `test/runtests.jl`
  - Registers the new expert test in the static expert group used by CI.

### Task 1: Add Laurent-Native ECP Boundary Diagnostics

**Files:**
- Modify: `src/algorithm/column_reduction.jl`
- Create: `test/expert/laurent_native_ecp_boundary_diagnostics.jl`
- Modify: `test/expert/laurent_column_reduction_diagnostics.jl`
- Modify: `test/runtests.jl`

**Interfaces:**
- Consumes: `Suslin.diagnose_unimodular_column_reduction(column, R) -> NamedTuple`.
- Consumes: `_column_reduction_stage_detail(stage::Symbol, R, outcome::Symbol; kwargs...) -> NamedTuple`.
- Consumes: `_diagnose_laurent_row_preconditioning(column, R, attempted, details)`.
- Produces: attempted diagnostic stage `:laurent_native_ecp_boundary`.
- Produces: `stage_details` entry with fields listed in Global Constraints.

- [ ] **Step 1: Write the failing Laurent-native boundary expert test**

Create `test/expert/laurent_native_ecp_boundary_diagnostics.jl`:

```julia
using Test
using Suslin
using Oscar

include(joinpath(@__DIR__, "..", "fixtures", "toricbuilder_case008_d14_column_boundary.jl"))
include(joinpath(@__DIR__, "..", "fixtures", "toricbuilder_case008_d15_column_boundary.jl"))

function _diagnostic_stage_detail(diagnostic, stage::Symbol)
    idx = findfirst(detail -> detail.stage == stage, diagnostic.stage_details)
    return idx === nothing ? nothing : diagnostic.stage_details[idx]
end

function _assert_laurent_native_ecp_boundary_detail(detail)
    @test detail !== nothing
    @test detail.outcome == :staged_boundary
    @test detail.boundary == :laurent_native_ecp
    @test detail.requires_descent_measure == true
    @test detail.requires_link_witness == true
    @test detail.requires_endpoint_reduction == true
    @test detail.requires_laurent_normality_replay == true
    @test detail.requires_recursive_peel_integration == true
    @test detail.fallback_policy == :diagnostic_only
end

@testset "Laurent native ECP boundary diagnostics" begin
    d14_fixture = ToricBuilderCase008D14ColumnBoundary.boundary_fixture()
    @test ToricBuilderCase008D14ColumnBoundary.validate_boundary_fixture(d14_fixture) == :ok

    d14 = Suslin.diagnose_unimodular_column_reduction(
        d14_fixture.failing_column,
        d14_fixture.ring,
        assume_unimodular = true,
        laurent_large_support_diagnostic_decline = true,
    )
    @test d14.status == :unsupported
    @test d14.failure_code == :unsupported_laurent_column_family
    @test d14.column_length == 14
    @test :laurent_elementary_row_preconditioning in d14.attempted_stages
    @test :laurent_native_ecp_boundary in d14.attempted_stages
    preconditioning_idx = findfirst(==(:laurent_elementary_row_preconditioning), d14.attempted_stages)
    boundary_idx = findfirst(==(:laurent_native_ecp_boundary), d14.attempted_stages)
    @test preconditioning_idx !== nothing
    @test boundary_idx !== nothing
    @test boundary_idx > preconditioning_idx
    @test length(d14.stage_details) == length(d14.attempted_stages)
    d14_witness = _diagnostic_stage_detail(d14, :laurent_witness_unit)
    @test d14_witness !== nothing
    @test d14_witness.outcome == :witness_support_too_large
    @test d14_witness.max_entry_term_count == 3734
    d14_normalization = _diagnostic_stage_detail(d14, :laurent_normalization)
    @test d14_normalization !== nothing
    @test d14_normalization.outcome == :delegation_declined_large_support
    @test d14_normalization.normalized_status == :declined
    @test d14_normalization.normalized_failure_code == :support_too_large
    @test d14_normalization.max_entry_term_count == 3734
    _assert_laurent_native_ecp_boundary_detail(
        _diagnostic_stage_detail(d14, :laurent_native_ecp_boundary),
    )

    d15_fixture = ToricBuilderCase008D15ColumnBoundary.boundary_fixture()
    d15 = Suslin.diagnose_unimodular_column_reduction(
        d15_fixture.failing_column,
        d15_fixture.ring,
    )
    @test d15.status == :supported
    @test d15.failure_code === nothing
    @test :laurent_elementary_row_preconditioning in d15.attempted_stages
    @test !(:laurent_native_ecp_boundary in d15.attempted_stages)
    @test _diagnostic_stage_detail(d15, :laurent_native_ecp_boundary) === nothing
    d15_preconditioned =
        _diagnostic_stage_detail(d15, :laurent_elementary_row_preconditioning)
    @test d15_preconditioned !== nothing
    @test d15_preconditioned.outcome == :supported

    R, (u, _) = Suslin.suslin_laurent_polynomial_ring(GF(2), ["u", "v"])
    large_support_entry = sum(u^k for k in 1:1001)
    large_support_witness_column = [
        large_support_entry,
        large_support_entry + one(R),
        zero(R),
    ]
    large_support_witness = Suslin.diagnose_unimodular_column_reduction(
        large_support_witness_column,
        R,
    )
    @test large_support_witness.status == :supported
    @test large_support_witness.failure_code === nothing
    @test :laurent_witness_unit in large_support_witness.attempted_stages
    @test !(:laurent_native_ecp_boundary in large_support_witness.attempted_stages)
    large_support_witness_detail =
        _diagnostic_stage_detail(large_support_witness, :laurent_witness_unit)
    @test large_support_witness_detail !== nothing
    @test large_support_witness_detail.outcome == :supported
    @test large_support_witness_detail.witness_unit_index !== nothing

    non_unimodular =
        ToricBuilderCase008D14ColumnBoundary.non_unimodular_negative_control(d14_fixture)
    non_unimodular_diagnostic = Suslin.diagnose_unimodular_column_reduction(
        non_unimodular.failing_column,
        non_unimodular.ring,
    )
    @test non_unimodular_diagnostic.status == :precondition_failed
    @test non_unimodular_diagnostic.failure_code == :not_unimodular
    @test isempty(non_unimodular_diagnostic.attempted_stages)
    @test isempty(non_unimodular_diagnostic.stage_details)
    @test !(:laurent_native_ecp_boundary in non_unimodular_diagnostic.attempted_stages)
end
```

Add a self-registration assertion near the top of the testset:

```julia
    runtests = read(joinpath(@__DIR__, "..", "runtests.jl"), String)
    @test occursin(
        "\"expert/laurent_native_ecp_boundary_diagnostics.jl\"",
        runtests,
    )
```

- [ ] **Step 2: Update the existing diagnostic test expectations**

In `test/expert/laurent_column_reduction_diagnostics.jl`, in the d15 and d16
supported fixture checks, add:

```julia
    @test !(:laurent_native_ecp_boundary in case008_d15.attempted_stages)
    @test _diagnostic_stage_detail(case008_d15, :laurent_native_ecp_boundary) === nothing
```

and:

```julia
    @test !(:laurent_native_ecp_boundary in case008_d16.attempted_stages)
    @test _diagnostic_stage_detail(case008_d16, :laurent_native_ecp_boundary) === nothing
```

In the existing `unsupported` Laurent example, change the stage loop to include
the new terminal stage:

```julia
    for stage in (:unit_entry, :laurent_unit_creation, :laurent_witness_unit, :laurent_normalization, :witness_unit, :monicity_normalization, :laurent_elementary_row_preconditioning, :laurent_native_ecp_boundary)
        @test stage in unsupported.attempted_stages
    end
```

Then add this detail assertion after the existing `@test !any(detail -> detail.outcome == :supported, unsupported.stage_details)`:

```julia
    unsupported_boundary = _diagnostic_stage_detail(unsupported, :laurent_native_ecp_boundary)
    @test unsupported_boundary !== nothing
    @test unsupported_boundary.outcome == :staged_boundary
    @test unsupported_boundary.boundary == :laurent_native_ecp
    @test unsupported_boundary.requires_descent_measure == true
    @test unsupported_boundary.requires_link_witness == true
    @test unsupported_boundary.requires_endpoint_reduction == true
    @test unsupported_boundary.requires_laurent_normality_replay == true
    @test unsupported_boundary.requires_recursive_peel_integration == true
    @test unsupported_boundary.fallback_policy == :diagnostic_only
```

In `test/runtests.jl`, add:

```julia
        "expert/laurent_native_ecp_boundary_diagnostics.jl",
```

immediately after `expert/laurent_column_reduction_diagnostics.jl`.

- [ ] **Step 3: Run focused tests to verify RED**

Run:

```bash
julia --project=. -e 'include("test/expert/laurent_native_ecp_boundary_diagnostics.jl")'
julia --project=. -e 'include("test/expert/laurent_column_reduction_diagnostics.jl")'
julia --project=. test/runtests.jl expert
```

Expected: the commands fail after loading dependencies because
`:laurent_native_ecp_boundary` is not yet produced. If the local environment
lacks `Oscar`, record that the RED run is environment-blocked and continue with
the code change; do not change production code before the test files are
written.

- [ ] **Step 4: Add the terminal boundary stage helper**

In `src/algorithm/column_reduction.jl`, immediately after
`_column_reduction_stage_detail`, add:

```julia
function _laurent_native_ecp_boundary_stage_detail(R)
    return _column_reduction_stage_detail(
        :laurent_native_ecp_boundary,
        R,
        :staged_boundary;
        boundary = :laurent_native_ecp,
        requires_descent_measure = true,
        requires_link_witness = true,
        requires_endpoint_reduction = true,
        requires_laurent_normality_replay = true,
        requires_recursive_peel_integration = true,
        fallback_policy = :diagnostic_only,
    )
end
```

Add a small diagnostic support-size helper near the same area. It should count
entry coefficients defensively, return `nothing` when the maximum entry support
is at or below 1000 terms, and otherwise return a named tuple with
`max_entry_term_count` and `support_term_limit`.

- [ ] **Step 5: Append the boundary after Laurent row preconditioning declines**

In `_diagnose_laurent_row_preconditioning`, after the existing
`:no_row_preconditioning_candidate` detail is pushed and before the final
unsupported return, add:

```julia
    push!(attempted, :laurent_native_ecp_boundary)
    push!(details, _laurent_native_ecp_boundary_stage_detail(R))
    return (; supported = false, stage = :laurent_native_ecp_boundary)
```

Remove the old final `return (; supported = false, stage = nothing)` from that
function so each unsupported path after row preconditioning has exactly one
terminal boundary stage detail.

Also thread `assume_unimodular::Bool = false` through
`diagnose_unimodular_column_reduction` into
`_diagnose_unimodular_column_preconditions` as
`check_unimodularity = !assume_unimodular`. When false, the precondition helper
must still enforce indexing, length, and ring coercion before returning the
validated column.

Thread `laurent_large_support_diagnostic_decline::Bool = false` through
`diagnose_unimodular_column_reduction` into the Laurent diagnostic path. Before
calling `_laurent_unimodular_witness`, use the support-size helper only when
that flag is true to record `:witness_support_too_large` with
`witness_unit_index = nothing` for the validated d14 input. After Laurent
normalization succeeds, apply the same opt-in support-size check to the
normalized ordinary column and record
`:delegation_declined_large_support` with
`normalized_status = :declined` and
`normalized_failure_code = :support_too_large` before continuing to row
preconditioning.

- [ ] **Step 6: Run focused tests to verify GREEN**

Run:

```bash
julia --project=. -e 'include("test/expert/laurent_native_ecp_boundary_diagnostics.jl")'
julia --project=. -e 'include("test/expert/laurent_column_reduction_diagnostics.jl")'
```

Expected: both commands exit 0 in a dependency-complete environment.

- [ ] **Step 7: Run package and diff verification**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
julia --project=. test/runtests.jl expert
git diff --check
```

Expected: both commands exit 0 in a dependency-complete environment. If Julia
dependency resolution is blocked by the sandbox, record the exact missing
dependency or network error and keep `git diff --check` as the local syntax
gate.

- [ ] **Step 8: Commit**

Local git commits may be blocked in Agent Desk because the managed worktree
index lives outside the writable sandbox. If local commit is blocked, include
these changed files in the GitHub API commit during the PR creation step:

```bash
docs/superpowers/specs/2026-07-05-issue-316-laurent-native-ecp-boundary-design.md
docs/superpowers/plans/2026-07-05-issue-316-laurent-native-ecp-boundary.md
src/algorithm/column_reduction.jl
test/expert/laurent_native_ecp_boundary_diagnostics.jl
test/expert/laurent_column_reduction_diagnostics.jl
test/runtests.jl
```

Commit message:

```bash
git commit -m "test: expose Laurent native ECP boundary diagnostics"
```

## Plan Self-Review

- Spec coverage: the plan covers terminal stage visibility, required detail
  fields, d14 positive boundary, d15 supported negative control,
  non-unimodular precondition negative control, no Laurent ECP implementation,
  no ordinary general ECP route for Laurent columns, the validated-fixture
  `assume_unimodular` shortcut, the opt-in large-support diagnostic decline,
  the default large-support witness-unit regression, and static expert runner
  registration.
- Placeholder scan: no deferred implementation placeholders remain.
- Type consistency: stage symbols and named-tuple field names match the issue
  body and the existing diagnostic helper style.

## Automatic Execution Choice

Plan complete and saved to
`docs/superpowers/plans/2026-07-05-issue-316-laurent-native-ecp-boundary.md`.
Using option 1, Subagent-Driven, because it is marked recommended by the
writing-plans skill and the Standing Answer Policy says to choose recommended
options automatically.
