# Issue 165 Stage-Level Laurent Column Reduction Diagnostics Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add backward-compatible `stage_details` diagnostics for Laurent column-reduction stages, including the `case_008` d=16 witness-without-unit and normalized-not-unimodular boundary.

**Architecture:** Keep `diagnose_unimodular_column_reduction` as a named-tuple diagnostic API and add a parallel `stage_details` tuple. Thread a details vector through the existing diagnostic helper paths, use the existing Laurent witness and normalization helpers, and stop normalized ordinary diagnosis at `:not_unimodular` when the normalized column fails the ordinary unimodularity precondition.

**Tech Stack:** Julia, Oscar, existing Suslin column-reduction helpers, Test stdlib, ToricBuilder fixture modules.

## Global Constraints

- Keep the change additive and backward-compatible: existing `status`, `failure_code`, `attempted_stages`, and `message` fields must remain valid.
- Do not export a new public diagnostic type.
- Do not change `reduce_unimodular_column` behavior.
- Do not add a new reducer stage.
- `stage_details` must be a tuple of named tuples.
- Precondition failures must not contain fake successful stage details.
- The d=16 fixture must include `stage = :laurent_witness_unit`, `outcome = :witness_without_unit`, and `witness_unit_index === nothing`.
- The d=16 fixture must include a `:laurent_normalization` detail that records a length-16 normalized ordinary polynomial column with `normalized_status == :precondition_failed` and `normalized_failure_code == :not_unimodular`.
- A supported d=21 Laurent witness-unit fixture must still report `status == :supported`.
- The required verification commands are:
  - `julia --project=. -e 'include("test/expert/laurent_column_reduction_diagnostics.jl")'`
  - `julia --project=. -e 'include("test/internal/toricbuilder_case008_d16_column_boundary.jl")'`
  - `julia --project=. test/runtests.jl`
  - `julia --project=. -e 'using Pkg; Pkg.test()'`

---

## File Structure

- `test/expert/laurent_column_reduction_diagnostics.jl`: extend focused diagnostic coverage for `stage_details`, negative controls, d=21 supported witness-unit details, and unsupported Laurent details.
- `test/internal/toricbuilder_case008_d16_column_boundary.jl`: pin d=16 stage details for the Laurent witness-unit and Laurent normalization boundary.
- `test/fixtures/toricbuilder_case008_d16_column_boundary.jl`: make the fixture validator reject diagnostics that do not report the expected d=16 stage-detail profile.
- `src/algorithm/column_reduction.jl`: add the diagnostic detail plumbing and normalized ordinary-column precondition classification.

---

### Task 1: Add Failing Stage-Details Tests

**Files:**
- Modify: `test/expert/laurent_column_reduction_diagnostics.jl`
- Modify: `test/internal/toricbuilder_case008_d16_column_boundary.jl`

**Interfaces:**
- Consumes: `Suslin.diagnose_unimodular_column_reduction(column, R)`, `ToricBuilderCase008D21ColumnBoundary.boundary_fixture()`, `ToricBuilderCase008D16ColumnBoundary.boundary_fixture()`.
- Produces: failing tests that require `diagnostic.stage_details` and the issue #165 d=16 detail outcomes.

- [ ] **Step 1: Add test helpers to the expert diagnostics test**

In `test/expert/laurent_column_reduction_diagnostics.jl`, after the fixture includes, add:

```julia
function _diagnostic_stage_detail(diagnostic, stage::Symbol)
    idx = findfirst(detail -> detail.stage == stage, diagnostic.stage_details)
    return idx === nothing ? nothing : diagnostic.stage_details[idx]
end
```

- [ ] **Step 2: Extend the expert diagnostics assertions**

In the `case010`, `case008_d21`, `unsupported`, `supported`, and `precondition`
assertion blocks of `test/expert/laurent_column_reduction_diagnostics.jl`, add
these assertions:

```julia
    @test hasproperty(case010, :stage_details)
    @test length(case010.stage_details) == length(case010.attempted_stages)
    case010_unit_creation = _diagnostic_stage_detail(case010, :laurent_unit_creation)
    @test case010_unit_creation !== nothing
    @test case010_unit_creation.outcome == :supported
    @test case010_unit_creation.pivot_index isa Integer
```

```julia
    @test hasproperty(case008_d21, :stage_details)
    @test length(case008_d21.stage_details) == length(case008_d21.attempted_stages)
    case008_d21_witness = _diagnostic_stage_detail(case008_d21, :laurent_witness_unit)
    @test case008_d21_witness !== nothing
    @test case008_d21_witness.outcome == :supported
    @test case008_d21_witness.witness_unit_index isa Integer
```

```julia
    @test hasproperty(unsupported, :stage_details)
    @test length(unsupported.stage_details) == length(unsupported.attempted_stages)
    @test !any(detail -> detail.outcome == :supported, unsupported.stage_details)
```

```julia
    @test hasproperty(supported, :stage_details)
    @test length(supported.stage_details) == length(supported.attempted_stages)
```

```julia
    @test hasproperty(precondition, :stage_details)
    @test isempty(precondition.stage_details)
```

- [ ] **Step 3: Add the d=16 detail helper to the internal fixture test**

In `test/internal/toricbuilder_case008_d16_column_boundary.jl`, after the
constant path definition, add:

```julia
function _case008_d16_stage_detail(diagnostic, stage::Symbol)
    idx = findfirst(detail -> detail.stage == stage, diagnostic.stage_details)
    return idx === nothing ? nothing : diagnostic.stage_details[idx]
end
```

- [ ] **Step 4: Pin the d=16 detail outcomes**

In `test/internal/toricbuilder_case008_d16_column_boundary.jl`, immediately
after the existing ring-profile diagnostic assertions, add:

```julia
    @test hasproperty(diagnostic, :stage_details)
    @test length(diagnostic.stage_details) == length(diagnostic.attempted_stages)

    witness_detail = _case008_d16_stage_detail(diagnostic, :laurent_witness_unit)
    @test witness_detail !== nothing
    @test witness_detail.outcome == :witness_without_unit
    @test witness_detail.witness_unit_index === nothing

    normalization_detail = _case008_d16_stage_detail(diagnostic, :laurent_normalization)
    @test normalization_detail !== nothing
    @test normalization_detail.outcome == :normalized_not_unimodular
    @test normalization_detail.normalized_column_length == 16
    @test normalization_detail.normalized_ring_kind == :polynomial
    @test normalization_detail.normalized_status == :precondition_failed
    @test normalization_detail.normalized_failure_code == :not_unimodular
```

- [ ] **Step 5: Run the focused tests and verify RED**

Run:

```bash
julia --project=. -e 'include("test/expert/laurent_column_reduction_diagnostics.jl")'
julia --project=. -e 'include("test/internal/toricbuilder_case008_d16_column_boundary.jl")'
```

Expected: both commands fail with errors or test failures caused by missing
`stage_details` on the diagnostic named tuple.

- [ ] **Step 6: Commit the failing tests**

Run:

```bash
git add test/expert/laurent_column_reduction_diagnostics.jl test/internal/toricbuilder_case008_d16_column_boundary.jl
git commit -m "test: cover stage-level laurent diagnostics"
```

Expected: a commit containing only test changes.

---

### Task 2: Add Diagnostic Stage Details

**Files:**
- Modify: `src/algorithm/column_reduction.jl`
- Test: `test/expert/laurent_column_reduction_diagnostics.jl`
- Test: `test/internal/toricbuilder_case008_d16_column_boundary.jl`

**Interfaces:**
- Consumes: `_laurent_unimodular_witness(column, R)`, `_reduce_via_laurent_unit_creation_certificate(column, R)`, `_witness_unit_reduction_certificate_stage(column, witness, pivot_idx, R)`, `normalize_laurent_object(column)`, `is_unimodular_column(column, R)`.
- Produces: `stage_details` in every diagnostic named tuple, with one detail per attempted stage.

- [ ] **Step 1: Extend `_column_reduction_diagnostic`**

In `src/algorithm/column_reduction.jl`, replace `_column_reduction_diagnostic`
with:

```julia
function _column_reduction_diagnostic(
    status::Symbol,
    failure_code,
    ring_profile,
    column_length::Int,
    attempted_stages,
    message::AbstractString,
    stage_details = (),
)
    return (;
        status,
        failure_code,
        ring_profile,
        column_length,
        attempted_stages = tuple(attempted_stages...),
        message = String(message),
        stage_details = tuple(stage_details...),
    )
end
```

- [ ] **Step 2: Add stage-detail helpers after `_column_reduction_ring_profile`**

Insert:

```julia
function _column_reduction_ring_kind(R)
    return _column_reduction_ring_profile(R).kind
end

function _column_reduction_stage_detail(stage::Symbol, R, outcome::Symbol; kwargs...)
    return (;
        stage,
        ring_kind = _column_reduction_ring_kind(R),
        outcome,
        kwargs...,
    )
end
```

- [ ] **Step 3: Thread `details` through validated diagnosis**

In `_diagnose_unimodular_column_reduction_validated`, create `details = Any[]`,
pass it into the Laurent or polynomial helper, and pass it into
`_column_reduction_diagnostic` for both supported and unsupported returns:

```julia
    attempted = Symbol[]
    details = Any[]
    result = _is_laurent_polynomial_ring(R) ?
        _diagnose_laurent_unimodular_column_reduction(column, R, attempted, details) :
        _diagnose_polynomial_unimodular_column_reduction(column, R, attempted, details)
```

Each `_column_reduction_diagnostic(...)` call in this function must receive
`details` as the final argument.

- [ ] **Step 4: Replace `_diagnose_laurent_unimodular_column_reduction`**

Replace the function with:

```julia
function _diagnose_laurent_unimodular_column_reduction(
    column::AbstractVector,
    R,
    attempted::Vector{Symbol},
    details::Vector,
)
    push!(attempted, :unit_entry)
    unit_idx = findfirst(is_unit, column)
    if unit_idx !== nothing
        push!(details, _column_reduction_stage_detail(:unit_entry, R, :supported; pivot_index = unit_idx))
        return (; supported = true, stage = :unit_entry)
    end
    push!(details, _column_reduction_stage_detail(:unit_entry, R, :no_unit_entry; pivot_index = nothing))

    push!(attempted, :laurent_unit_creation)
    unit_creation = _reduce_via_laurent_unit_creation_certificate(column, R)
    if unit_creation !== nothing
        push!(
            details,
            _column_reduction_stage_detail(
                :laurent_unit_creation,
                R,
                :supported;
                pivot_index = unit_creation.stage.pivot_index,
                source_index = unit_creation.stage.source_index,
            ),
        )
        return (; supported = true, stage = :laurent_unit_creation)
    end
    push!(details, _column_reduction_stage_detail(:laurent_unit_creation, R, :no_unit_creation_candidate))

    push!(attempted, :laurent_witness_unit)
    witness = _laurent_unimodular_witness(column, R)
    if witness === nothing
        push!(
            details,
            _column_reduction_stage_detail(
                :laurent_witness_unit,
                R,
                :witness_unavailable;
                witness_unit_index = nothing,
            ),
        )
    else
        witness_unit_idx = findfirst(is_unit, witness)
        if witness_unit_idx === nothing
            push!(
                details,
                _column_reduction_stage_detail(
                    :laurent_witness_unit,
                    R,
                    :witness_without_unit;
                    witness_unit_index = nothing,
                ),
            )
        else
            _witness_unit_reduction_certificate_stage(column, witness, witness_unit_idx, R)
            push!(
                details,
                _column_reduction_stage_detail(
                    :laurent_witness_unit,
                    R,
                    :supported;
                    witness_unit_index = witness_unit_idx,
                ),
            )
            return (; supported = true, stage = :laurent_witness_unit)
        end
    end

    push!(attempted, :laurent_normalization)
    normalization = normalize_laurent_object(column)
    poly_column = normalization.normalized_object
    P = normalization.metadata.polynomial_ring

    normalized_unimodular = try
        is_unimodular_column(poly_column, P)
    catch err
        err isa InterruptException && rethrow()
        push!(
            details,
            _column_reduction_stage_detail(
                :laurent_normalization,
                R,
                :normalized_unimodularity_check_failed;
                normalized_column_length = length(poly_column),
                normalized_ring_kind = _column_reduction_ring_kind(P),
                normalized_status = :precondition_failed,
                normalized_failure_code = :unimodularity_check_failed,
                normalized_message = _column_reduction_error_message(err),
            ),
        )
        return (; supported = false, stage = nothing)
    end

    if !normalized_unimodular
        push!(
            details,
            _column_reduction_stage_detail(
                :laurent_normalization,
                R,
                :normalized_not_unimodular;
                normalized_column_length = length(poly_column),
                normalized_ring_kind = _column_reduction_ring_kind(P),
                normalized_status = :precondition_failed,
                normalized_failure_code = :not_unimodular,
            ),
        )
        return (; supported = false, stage = nothing)
    end

    ordinary_attempted = Symbol[]
    ordinary_details = Any[]
    result = _diagnose_polynomial_unimodular_column_reduction(poly_column, P, ordinary_attempted, ordinary_details)
    normalized_status = result.supported ? :supported : :unsupported
    normalized_failure_code = result.supported ? nothing : :unsupported_polynomial_column_family
    push!(
        details,
        _column_reduction_stage_detail(
            :laurent_normalization,
            R,
            :delegated_to_polynomial;
            normalized_column_length = length(poly_column),
            normalized_ring_kind = _column_reduction_ring_kind(P),
            normalized_status,
            normalized_failure_code,
        ),
    )
    append!(attempted, ordinary_attempted)
    append!(details, ordinary_details)
    return result
end
```

- [ ] **Step 5: Replace ordinary diagnostic helper signatures and details**

Update `_diagnose_polynomial_unimodular_column_reduction`,
`_diagnose_exact_small_column_reduction`, and
`_diagnose_supported_unimodular_column_reduction` to accept `details::Vector`.
Use the same stage order as today and append:

```julia
_column_reduction_stage_detail(:unit_entry, R, :supported; pivot_index = unit_idx)
_column_reduction_stage_detail(:unit_entry, R, :no_unit_entry; pivot_index = nothing)
_column_reduction_stage_detail(:witness_unit, R, :witness_unavailable; witness_unit_index = nothing)
_column_reduction_stage_detail(:witness_unit, R, :witness_without_unit; witness_unit_index = nothing)
_column_reduction_stage_detail(:witness_unit, R, :supported; witness_unit_index = witness_unit_idx)
_column_reduction_stage_detail(:monicity_normalization, R, :supported; normalized_column_length = length(column))
_column_reduction_stage_detail(:monicity_normalization, R, :no_monicity_normalization)
_column_reduction_stage_detail(:three_entry_block, R, :supported; block_indices = block.stage.indices, pivot_index = block.stage.indices[end])
_column_reduction_stage_detail(:three_entry_block, R, :no_supported_three_block)
```

Keep all current return values unchanged except for the extra detail threading.

- [ ] **Step 6: Run focused tests and verify GREEN**

Run:

```bash
julia --project=. -e 'include("test/expert/laurent_column_reduction_diagnostics.jl")'
julia --project=. -e 'include("test/internal/toricbuilder_case008_d16_column_boundary.jl")'
```

Expected: both commands exit 0.

- [ ] **Step 7: Commit the implementation**

Run:

```bash
git add src/algorithm/column_reduction.jl
git commit -m "feat: add stage-level column diagnostics"
```

Expected: a commit containing only `src/algorithm/column_reduction.jl`.

---

### Task 3: Pin The D16 Fixture Validator Profile

**Files:**
- Modify: `test/fixtures/toricbuilder_case008_d16_column_boundary.jl`
- Test: `test/internal/toricbuilder_case008_d16_column_boundary.jl`

**Interfaces:**
- Consumes: `diagnostic.stage_details` from Task 2.
- Produces: fixture validation that rejects missing or wrong d=16 stage-detail profile.

- [ ] **Step 1: Extend the expected diagnostic tuple**

In `test/fixtures/toricbuilder_case008_d16_column_boundary.jl`, replace
`EXPECTED_DIAGNOSTIC` with:

```julia
const EXPECTED_DIAGNOSTIC = (;
    status = :unsupported,
    failure_code = :unsupported_laurent_column_family,
    laurent_witness_outcome = :witness_without_unit,
    laurent_witness_unit_index = nothing,
    laurent_normalization_outcome = :normalized_not_unimodular,
    normalized_column_length = 16,
    normalized_status = :precondition_failed,
    normalized_failure_code = :not_unimodular,
)
```

- [ ] **Step 2: Add fixture-local detail helpers**

After `_diagnostic_matches_expected`, add:

```julia
function _diagnostic_stage_detail(diagnostic, stage::Symbol)
    hasproperty(diagnostic, :stage_details) || return nothing
    idx = findfirst(detail -> detail.stage == stage, diagnostic.stage_details)
    return idx === nothing ? nothing : diagnostic.stage_details[idx]
end

function _diagnostic_stage_profile_matches(diagnostic, expected)::Bool
    witness = _diagnostic_stage_detail(diagnostic, :laurent_witness_unit)
    witness === nothing && return false
    witness.outcome == expected.laurent_witness_outcome || return false
    witness.witness_unit_index === expected.laurent_witness_unit_index || return false

    normalization = _diagnostic_stage_detail(diagnostic, :laurent_normalization)
    normalization === nothing && return false
    normalization.outcome == expected.laurent_normalization_outcome || return false
    normalization.normalized_column_length == expected.normalized_column_length || return false
    normalization.normalized_status == expected.normalized_status || return false
    normalization.normalized_failure_code == expected.normalized_failure_code || return false

    return true
end
```

- [ ] **Step 3: Strengthen `_diagnostic_matches_expected`**

Replace `_diagnostic_matches_expected` with:

```julia
function _diagnostic_matches_expected(column, R, expected)::Bool
    diagnostic = Suslin.diagnose_unimodular_column_reduction(column, R)
    return diagnostic.status == expected.status &&
        diagnostic.failure_code == expected.failure_code &&
        diagnostic.column_length == FIRST_FAILING_PEEL_DIMENSION &&
        _diagnostic_stage_profile_matches(diagnostic, expected)
end
```

- [ ] **Step 4: Run the d=16 focused test**

Run:

```bash
julia --project=. -e 'include("test/internal/toricbuilder_case008_d16_column_boundary.jl")'
```

Expected: exits 0.

- [ ] **Step 5: Commit the fixture validator update**

Run:

```bash
git add test/fixtures/toricbuilder_case008_d16_column_boundary.jl
git commit -m "test: pin case008 d16 diagnostic profile"
```

Expected: a commit containing only the fixture validator change.

---

## Final Verification

After all tasks pass their focused checks, run:

```bash
julia --project=. -e 'include("test/expert/laurent_column_reduction_diagnostics.jl")'
julia --project=. -e 'include("test/internal/toricbuilder_case008_d16_column_boundary.jl")'
julia --project=. test/runtests.jl
julia --project=. -e 'using Pkg; Pkg.test()'
git diff --check origin/main..HEAD
```

Expected: all Julia commands exit 0, and `git diff --check origin/main..HEAD`
prints no output.

## Plan Self-Review

- Spec coverage: Task 1 defines failing tests, Task 2 implements additive
  `stage_details`, and Task 3 pins the d=16 fixture profile.
- Reserved-token scan: no incomplete-work markers are intended in this plan.
- Type consistency: all new diagnostic detail fields use symbols, integers, or
  `nothing` as required by the issue examples.
