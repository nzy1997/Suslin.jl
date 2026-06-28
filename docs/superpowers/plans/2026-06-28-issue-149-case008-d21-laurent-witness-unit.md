# Issue 149 Case008 D21 Laurent Witness Unit Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Reduce the extracted `case_008` `d=21` Laurent unimodular column through a certified witness-unit stage.

**Architecture:** Add a narrow Laurent route inside the existing column-reduction pipeline. The route solves a Laurent Bezout row system with `solve_laurent_linear`, requires a unit witness coefficient, and delegates factor construction and replay to the existing `:witness_unit` certificate stage.

**Tech Stack:** Julia, Oscar, existing Suslin column-reduction helpers, Test stdlib, ToricBuilder fixture files.

## Global Constraints

- Do not claim arbitrary Laurent unimodular-column support.
- Do not require full `case_008` certificate success if the repair exposes a later boundary.
- Do not change the public `elementary_factorization` contract for original Laurent `GL_n` inputs.
- Do not export new public APIs.
- Keep the stored replay stage as `:witness_unit`; use `:laurent_witness_unit` only as a diagnostics attempted-stage marker.
- The repair must be described by an algebraic predicate, not by the `case_008` fixture id.

---

## File Structure

- `test/expert/case008_d21_laurent_column_reduction.jl`: new expert acceptance test for issue #149, including factor replay and certificate tampering negative controls.
- `src/algorithm/column_reduction.jl`: add the private Laurent witness helper and route; wire it into Laurent reduction and diagnostics.
- `test/fixtures/toricbuilder_case008_d21_column_boundary.jl`: update the fixture's expected diagnostic from unsupported to supported.
- `test/internal/toricbuilder_case008_d21_column_boundary.jl`: update boundary validator assertions so the already-extracted fixture remains valid after the route is supported.
- `test/expert/laurent_column_reduction_diagnostics.jl`: add a `case_008 d=21` supported diagnostic assertion while preserving the existing unsupported Laurent negative coverage.

---

### Task 1: Add The Issue 149 Expert Acceptance Test

**Files:**
- Create: `test/expert/case008_d21_laurent_column_reduction.jl`

**Interfaces:**
- Consumes: `ToricBuilderCase008D21ColumnBoundary.boundary_fixture()`, `Suslin.reduce_unimodular_column(column, R)`, `Suslin.ecp_column_reduction_certificate(column, R)`, `Suslin.verify_ecp_column_reduction(certificate)`, `Suslin.diagnose_unimodular_column_reduction(column, R)`, `Suslin.ECPColumnReductionCertificate`.
- Produces: a failing expert test that defines the issue #149 behavioral contract before implementation.

- [ ] **Step 1: Create the failing expert test**

Create `test/expert/case008_d21_laurent_column_reduction.jl` with this exact content:

```julia
using Test
using Suslin
using Oscar

include(joinpath(@__DIR__, "..", "fixtures", "toricbuilder_case008_d21_column_boundary.jl"))

function _case008_d21_reduction_product(factors, R, n::Int)
    product = identity_matrix(R, n)
    for factor in factors
        product *= factor
    end
    return product
end

function _case008_d21_apply_factors(factors, column, R)
    n = length(column)
    return _case008_d21_reduction_product(factors, R, n) * matrix(R, n, 1, collect(column))
end

function _case008_d21_target_column(R, n::Int)
    target = zero_matrix(R, n, 1)
    target[n, 1] = one(R)
    return target
end

function _case008_d21_tamper_first_factor(factors, R, n::Int)
    tampered = copy(factors)
    tampered[1] = identity_matrix(R, n)
    return tampered
end

function _case008_d21_tamper_certificate_first_factor(cert)
    tampered = _case008_d21_tamper_first_factor(
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

@testset "case_008 d=21 Laurent column reduction" begin
    fixture = ToricBuilderCase008D21ColumnBoundary.boundary_fixture()
    R = fixture.ring
    column = fixture.failing_column
    target = _case008_d21_target_column(R, length(column))

    factors = Suslin.reduce_unimodular_column(column, R)
    @test _case008_d21_apply_factors(factors, column, R) == target

    tampered_factors = _case008_d21_tamper_first_factor(factors, R, length(column))
    @test _case008_d21_apply_factors(tampered_factors, column, R) != target

    certificate = Suslin.ecp_column_reduction_certificate(column, R)
    @test Suslin.verify_ecp_column_reduction(certificate)
    @test any(stage -> stage.kind == :witness_unit, certificate.stages)
    @test _case008_d21_apply_factors(
        certificate.factors,
        certificate.original_column,
        certificate.ring,
    ) == target
    @test !Suslin.verify_ecp_column_reduction(
        _case008_d21_tamper_certificate_first_factor(certificate),
    )

    diagnostic = Suslin.diagnose_unimodular_column_reduction(column, R)
    @test diagnostic.status == :supported
    @test diagnostic.failure_code === nothing
    @test diagnostic.column_length == 21
    @test diagnostic.ring_profile.kind == :laurent_polynomial
    @test diagnostic.ring_profile.generators == ("u", "v")
    @test :laurent_witness_unit in diagnostic.attempted_stages
end
```

- [ ] **Step 2: Run the new test and verify it fails for the old boundary**

Run:

```bash
julia --project=. -e 'include("test/expert/case008_d21_laurent_column_reduction.jl")'
```

Expected: FAIL with an `ArgumentError` containing `unsupported exact unimodular column reduction for Laurent-normalized column of length 21`.

- [ ] **Step 3: Commit the failing acceptance test**

Run:

```bash
git add test/expert/case008_d21_laurent_column_reduction.jl
git commit -m "test: cover case008 d21 laurent reduction"
```

Expected: a commit containing only `test/expert/case008_d21_laurent_column_reduction.jl`.

---

### Task 2: Implement The Laurent Witness-Unit Route

**Files:**
- Modify: `src/algorithm/column_reduction.jl`
- Test: `test/expert/case008_d21_laurent_column_reduction.jl`

**Interfaces:**
- Consumes: `solve_laurent_linear(A, B)`, `_witness_unit_reduction_certificate_stage(column, witness, pivot_idx, R)`, `_reduce_via_laurent_unit_creation_certificate(column, R)`, `normalize_laurent_object(column)`.
- Produces: `_laurent_unimodular_witness(column::AbstractVector, R)`, `_reduce_via_laurent_witness_unit_certificate(column::AbstractVector, R)`, and diagnostics attempted stage `:laurent_witness_unit`.

- [ ] **Step 1: Add helper functions after `_reduce_via_witness_unit`**

In `src/algorithm/column_reduction.jl`, find:

```julia
function _reduce_via_witness_unit(column::AbstractVector, witness::AbstractVector, pivot_idx::Int, R)
    return _witness_unit_reduction_certificate_stage(column, witness, pivot_idx, R).factors
end
```

Immediately after it, insert:

```julia
function _laurent_witness_solve_failure(err)::Bool
    err isa ErrorException || return false
    message = sprint(showerror, err)
    return occursin("No exact solution exists for A * U = B", message) ||
        occursin("not liftable to the given generating system", message)
end

function _laurent_unimodular_witness(column::AbstractVector, R)
    _is_laurent_polynomial_ring(R) || return nothing
    row = matrix(R, 1, length(column), collect(column))
    rhs = matrix(R, 1, 1, [one(R)])

    solution = try
        solve_laurent_linear(row, rhs)
    catch err
        err isa InterruptException && rethrow()
        _laurent_witness_solve_failure(err) && return nothing
        rethrow()
    end

    row * solution == rhs || return nothing
    return [solution[idx, 1] for idx in 1:nrows(solution)]
end

function _reduce_via_laurent_witness_unit_certificate(column::AbstractVector, R)
    witness = _laurent_unimodular_witness(column, R)
    witness === nothing && return nothing

    witness_unit_idx = findfirst(is_unit, witness)
    witness_unit_idx === nothing && return nothing

    return _witness_unit_reduction_certificate_stage(column, witness, witness_unit_idx, R)
end
```

- [ ] **Step 2: Wire the route into Laurent certificate reduction**

In `src/algorithm/column_reduction.jl`, replace the start of `_reduce_laurent_unimodular_column_certificate` with:

```julia
function _reduce_laurent_unimodular_column_certificate(column::AbstractVector, R)
    unit_idx = findfirst(is_unit, column)
    unit_idx !== nothing && return _unit_entry_reduction_certificate_stage(column, unit_idx, R)

    unit_creation = _reduce_via_laurent_unit_creation_certificate(column, R)
    unit_creation !== nothing && return unit_creation

    witness_unit = _reduce_via_laurent_witness_unit_certificate(column, R)
    witness_unit !== nothing && return witness_unit

    normalization = normalize_laurent_object(column)
```

Leave the rest of the function unchanged from `poly_column = normalization.normalized_object` through the final `return (; factors, stage)`.

- [ ] **Step 3: Wire diagnostics into the Laurent diagnostic path**

In `src/algorithm/column_reduction.jl`, replace the start of `_diagnose_laurent_unimodular_column_reduction` with:

```julia
function _diagnose_laurent_unimodular_column_reduction(column::AbstractVector, R, attempted::Vector{Symbol})
    push!(attempted, :unit_entry)
    unit_idx = findfirst(is_unit, column)
    unit_idx !== nothing && return (; supported = true, stage = :unit_entry)

    push!(attempted, :laurent_unit_creation)
    unit_creation = _reduce_via_laurent_unit_creation_certificate(column, R)
    unit_creation !== nothing && return (; supported = true, stage = :laurent_unit_creation)

    push!(attempted, :laurent_witness_unit)
    witness_unit = _reduce_via_laurent_witness_unit_certificate(column, R)
    witness_unit !== nothing && return (; supported = true, stage = :laurent_witness_unit)

    push!(attempted, :laurent_normalization)
    normalization = normalize_laurent_object(column)
```

Leave the rest of the function unchanged from `poly_column = normalization.normalized_object` through the call to `_diagnose_polynomial_unimodular_column_reduction`.

- [ ] **Step 4: Run the issue expert test and verify it passes**

Run:

```bash
julia --project=. -e 'include("test/expert/case008_d21_laurent_column_reduction.jl")'
```

Expected: PASS with testset `case_008 d=21 Laurent column reduction`.

- [ ] **Step 5: Run the existing case010 regression test**

Run:

```bash
julia --project=. -e 'include("test/expert/case010_laurent_column_reduction.jl")'
```

Expected: PASS. This confirms the new route did not steal the existing `:laurent_unit_creation` route for `case_010`.

- [ ] **Step 6: Commit the reducer implementation**

Run:

```bash
git add src/algorithm/column_reduction.jl
git commit -m "feat: add laurent witness unit reduction"
```

Expected: a commit containing only `src/algorithm/column_reduction.jl`.

---

### Task 3: Update Boundary Expectations And Run Focused Verification

**Files:**
- Modify: `test/fixtures/toricbuilder_case008_d21_column_boundary.jl`
- Modify: `test/internal/toricbuilder_case008_d21_column_boundary.jl`
- Modify: `test/expert/laurent_column_reduction_diagnostics.jl`
- Test: `test/internal/toricbuilder_case008_d21_column_boundary.jl`
- Test: `test/expert/laurent_column_reduction_diagnostics.jl`

**Interfaces:**
- Consumes: supported diagnostics from Task 2.
- Produces: updated fixture validation and diagnostics coverage that match the new supported route while preserving non-unimodular and unsupported-column negative controls.

- [ ] **Step 1: Update the fixture expected diagnostic**

In `test/fixtures/toricbuilder_case008_d21_column_boundary.jl`, replace:

```julia
        expected_diagnostic = (;
            status = :unsupported,
            failure_code = :unsupported_laurent_column_family,
        ),
```

with:

```julia
        expected_diagnostic = (;
            status = :supported,
            failure_code = nothing,
        ),
```

- [ ] **Step 2: Update the internal boundary test diagnostic assertions**

In `test/internal/toricbuilder_case008_d21_column_boundary.jl`, replace:

```julia
    @test diagnostic.status == :unsupported
    @test diagnostic.failure_code == :unsupported_laurent_column_family
```

with:

```julia
    @test diagnostic.status == :supported
    @test diagnostic.failure_code === nothing
    @test :laurent_witness_unit in diagnostic.attempted_stages
```

Leave the existing `column_length`, `ring_profile`, validator, and negative-control assertions in place.

- [ ] **Step 3: Add case008 supported diagnostics coverage**

In `test/expert/laurent_column_reduction_diagnostics.jl`, add this include after the existing `case010` fixture include:

```julia
include(joinpath(@__DIR__, "..", "fixtures", "toricbuilder_case008_d21_column_boundary.jl"))
```

Then, inside the `@testset "Laurent column reduction diagnostics" begin` block immediately after the existing `case010` diagnostic assertions, insert:

```julia
    case008 = ToricBuilderCase008D21ColumnBoundary.boundary_fixture()
    case008_d21 = Suslin.diagnose_unimodular_column_reduction(
        case008.failing_column,
        case008.ring,
    )
    @test case008_d21.status == :supported
    @test case008_d21.failure_code === nothing
    @test case008_d21.column_length == 21
    @test case008_d21.ring_profile.kind == :laurent_polynomial
    @test case008_d21.ring_profile.generators == ("u", "v")
    @test :laurent_witness_unit in case008_d21.attempted_stages
```

Keep the existing `unsupported_column` block unchanged so `:unsupported_laurent_column_family` remains covered.

- [ ] **Step 4: Run focused fixture and diagnostics tests**

Run:

```bash
julia --project=. -e 'include("test/internal/toricbuilder_case008_d21_column_boundary.jl")'
julia --project=. -e 'include("test/expert/laurent_column_reduction_diagnostics.jl")'
```

Expected: both commands PASS.

- [ ] **Step 5: Run the bounded case008 report**

Run:

```bash
julia --project=. scripts/report_toricbuilder_cache_q_blocks.jl --exercise=case_008 --timeout-seconds=180 --output=/tmp/qblock-case008-after-d21.md
```

Expected: command exits successfully. Inspect `/tmp/qblock-case008-after-d21.md` and confirm it does not contain `unsupported_laurent_column_family` at `current d=21`. It may contain a later structured boundary or a Laurent `GL_n` certificate pass.

- [ ] **Step 6: Run the package's default test entry point**

Run:

```bash
julia --project=. test/runtests.jl
```

Expected: PASS for the default `public` and `internal` groups.

- [ ] **Step 7: Commit expectation updates and verification coverage**

Run:

```bash
git add test/fixtures/toricbuilder_case008_d21_column_boundary.jl test/internal/toricbuilder_case008_d21_column_boundary.jl test/expert/laurent_column_reduction_diagnostics.jl
git commit -m "test: update case008 d21 boundary support"
```

Expected: a commit containing only the fixture/internal/diagnostics expectation updates.

---

## Final Verification

After all tasks are complete, run:

```bash
julia --project=. -e 'include("test/expert/case008_d21_laurent_column_reduction.jl")'
julia --project=. -e 'include("test/expert/case010_laurent_column_reduction.jl")'
julia --project=. -e 'include("test/internal/toricbuilder_case008_d21_column_boundary.jl")'
julia --project=. -e 'include("test/expert/laurent_column_reduction_diagnostics.jl")'
julia --project=. test/runtests.jl
julia --project=. scripts/report_toricbuilder_cache_q_blocks.jl --exercise=case_008 --timeout-seconds=180 --output=/tmp/qblock-case008-after-d21.md
```

Expected:

- the four focused `include(...)` commands pass;
- the default test runner passes;
- the bounded `case_008` report no longer records the old `current d=21` `unsupported_laurent_column_family` boundary.
