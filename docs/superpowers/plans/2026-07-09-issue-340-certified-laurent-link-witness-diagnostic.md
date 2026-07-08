# Issue 340 Certified Laurent Link-Witness Diagnostic Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Expose a fixture-bound certified Laurent link-witness diagnostic stage for the validated `case_008 d=14` path and advance the terminal diagnostic boundary to endpoint reduction.

**Architecture:** Reuse the internal Laurent link-witness certificate replay and validation helpers from #339. Gate the new production diagnostic on the existing certified descent step, the exact post-descent d14 support fingerprint, and successful certificate validation. Leave unsupported and tampered Laurent columns on the existing Laurent-native ECP boundary without claiming link-witness progress.

**Tech Stack:** Julia, Suslin.jl, Oscar, existing `NamedTuple` diagnostic details, existing expert Julia test suite.

## Global Constraints

- The new diagnostic stage is conservative and fixture-bound to validated `case_008 d=14`.
- Do not implement Laurent endpoint reductions, endpoint-reduction contexts, Laurent normality/conjugation replay, determinant normalization, recursive peel integration, or full `case_008` success.
- The stage detail remains plain data derived from the internal certificate validator.
- Tampered columns, unsupported Laurent columns, non-unimodular columns, and already-supported d15 columns must not report `:laurent_link_witness_certificate`.
- The terminal `:laurent_native_ecp_boundary` still requires endpoint reduction, Laurent normality replay, and recursive peel integration.

---

### Task 1: Certified Link-Witness Diagnostic Stage

**Files:**
- Modify: `src/algorithm/column_reduction.jl`
- Modify: `test/expert/laurent_native_ecp_boundary_diagnostics.jl`

**Interfaces:**
- Consumes: `_laurent_descent_step_diagnostic_certificate(column, R)`, `_replay_laurent_elementary_entry_addition(column, R, operation)`, `_laurent_link_witness_certificate_from_replay(context, witness, column, R)`, `_validate_laurent_link_witness_certificate(cert, column, R)`, `_column_reduction_stage_detail(stage, R, outcome; kwargs...)`.
- Produces: `_laurent_link_witness_diagnostic_certificate(column::AbstractVector, R) -> Union{Nothing, NamedTuple}`.
- Produces: `_laurent_link_witness_certificate_stage_detail(R, certificate) -> NamedTuple`.
- Produces: `_laurent_native_ecp_boundary_stage_detail(R; certified_descent_step::Bool = false, certified_link_witness::Bool = false) -> NamedTuple`.

- [ ] **Step 1: Write the failing focused diagnostic test**

Patch `test/expert/laurent_native_ecp_boundary_diagnostics.jl`.

Replace `_assert_laurent_native_ecp_boundary_detail` with:

```julia
function _assert_laurent_native_ecp_boundary_detail(
    detail;
    requires_descent_measure::Bool = true,
    certified_descent_scope = nothing,
    requires_link_witness::Bool = true,
    next_boundary = certified_descent_scope === nothing ? nothing : :laurent_link_witness,
)
    @test detail !== nothing
    @test detail.outcome == :staged_boundary
    @test detail.boundary == :laurent_native_ecp
    @test detail.requires_descent_measure == requires_descent_measure
    @test detail.certified_descent_scope == certified_descent_scope
    @test detail.next_boundary == next_boundary
    @test detail.requires_link_witness == requires_link_witness
    @test detail.requires_endpoint_reduction == true
    @test detail.requires_laurent_normality_replay == true
    @test detail.requires_recursive_peel_integration == true
    @test detail.fallback_policy == :diagnostic_only
end
```

In the d14 positive assertion block, immediately after:

```julia
descent_idx = findfirst(==(:laurent_descent_step_certificate), d14.attempted_stages)
```

add:

```julia
link_witness_idx =
    findfirst(==(:laurent_link_witness_certificate), d14.attempted_stages)
```

Immediately after:

```julia
@test descent_idx !== nothing
```

add:

```julia
@test link_witness_idx !== nothing
```

Immediately after:

```julia
@test boundary_idx > descent_idx
```

add:

```julia
@test link_witness_idx > descent_idx
@test boundary_idx > link_witness_idx
```

Immediately after the existing descent-stage assertions ending with:

```julia
@test Suslin._strictly_decreases_laurent_measure(
    descent.before_measure,
    descent.after_measure,
)
```

add:

```julia
    link_witness =
        _diagnostic_stage_detail(d14, :laurent_link_witness_certificate)
    @test link_witness !== nothing
    @test link_witness.outcome == :certified_link_witness
    @test link_witness.witness_family == :two_entry_laurent_combination
    @test link_witness.pivot_index == 10
    @test link_witness.partner_index == 1
    @test link_witness.coefficient == 1
    @test link_witness.exponent == (1, -1)
    @test link_witness.replay_status == :ok
    @test link_witness.identity_status == :verified
    @test link_witness.certificate_status == :link_witness_certificate
    @test link_witness.context_status == :link_witness_context
    @test link_witness.source_endpoint.status == :link_witness_endpoint_metadata
    @test link_witness.target_endpoint.status == :link_witness_endpoint_metadata
    @test link_witness.source_endpoint.case_id == "case_008"
    @test link_witness.target_endpoint.case_id == "case_008"
    @test link_witness.source_endpoint.entry_index == 10
    @test link_witness.target_endpoint.entry_index == 10
    @test link_witness.next_boundary == :laurent_endpoint_reduction
```

Update the d14 terminal boundary assertion to:

```julia
    _assert_laurent_native_ecp_boundary_detail(
        _diagnostic_stage_detail(d14, :laurent_native_ecp_boundary);
        requires_descent_measure = false,
        certified_descent_scope = :single_certified_step,
        requires_link_witness = false,
        next_boundary = :laurent_endpoint_reduction,
    )
```

In the d15 negative assertions, immediately after:

```julia
@test !(:laurent_descent_step_certificate in d15.attempted_stages)
```

add:

```julia
@test !(:laurent_link_witness_certificate in d15.attempted_stages)
```

In the non-unimodular negative assertions, immediately after:

```julia
@test !(:laurent_descent_step_certificate in non_unimodular_diagnostic.attempted_stages)
```

add:

```julia
@test !(:laurent_link_witness_certificate in non_unimodular_diagnostic.attempted_stages)
```

In the tampered d14 negative assertions, immediately after:

```julia
@test !(:laurent_descent_step_certificate in tampered.attempted_stages)
```

add:

```julia
@test !(:laurent_link_witness_certificate in tampered.attempted_stages)
```

- [ ] **Step 2: Run the focused test and verify it fails for the missing stage**

Run:

```bash
julia --project=. -e 'include("test/expert/laurent_native_ecp_boundary_diagnostics.jl")'
```

Expected: FAIL because `link_witness_idx` is `nothing`, `:laurent_link_witness_certificate` is absent, and the boundary still reports `requires_link_witness == true`.

- [ ] **Step 3: Add the certified d14 link-witness diagnostic helper**

Patch `src/algorithm/column_reduction.jl`.

Immediately after `_LAURENT_D14_CERTIFIED_DESCENT_AFTER_FINGERPRINT`, add:

```julia
const _LAURENT_D14_CERTIFIED_LINK_WITNESS_SOURCE_FINGERPRINT =
    _LAURENT_D14_CERTIFIED_DESCENT_AFTER_FINGERPRINT
const _LAURENT_D14_CERTIFIED_LINK_WITNESS = (;
    family = :two_entry_laurent_combination,
    pivot_index = 10,
    partner_index = 1,
    coefficient = 1,
    exponent = (1, -1),
    ring_generators = ("u", "v"),
)
```

Immediately after `_laurent_descent_step_certificate_stage_detail`, add:

```julia
function _laurent_link_witness_diagnostic_certificate(column::AbstractVector, R)
    try
        witness = _LAURENT_D14_CERTIFIED_LINK_WITNESS
        _laurent_link_witness_status(witness, length(column), R) == :ok ||
            return nothing
        _laurent_descent_column_support_fingerprint(column) ==
            _LAURENT_D14_CERTIFIED_LINK_WITNESS_SOURCE_FINGERPRINT ||
            return nothing

        context = (;
            case_id = _LAURENT_D14_CERTIFIED_DESCENT_CASE_ID,
            dimension = length(column),
            ring_generators = _laurent_descent_ring_generators(R),
            status = :link_witness_context,
        )
        certificate = _laurent_link_witness_certificate_from_replay(
            context,
            witness,
            column,
            R,
        )
        _validate_laurent_link_witness_certificate(certificate, column, R) ==
            :ok || return nothing
        return certificate
    catch err
        err isa InterruptException && rethrow()
        return nothing
    end
end

function _laurent_link_witness_certificate_stage_detail(R, certificate)
    witness = certificate.witness
    return _column_reduction_stage_detail(
        :laurent_link_witness_certificate,
        R,
        :certified_link_witness;
        witness_family = witness.family,
        pivot_index = witness.pivot_index,
        partner_index = witness.partner_index,
        coefficient = witness.coefficient,
        exponent = witness.exponent,
        source_endpoint = certificate.source_endpoint,
        target_endpoint = certificate.target_endpoint,
        replay_status = certificate.replay_status,
        identity_status = certificate.identity_status,
        certificate_status = certificate.status,
        context_status = certificate.context_status,
        next_boundary = certificate.next_boundary,
    )
end
```

- [ ] **Step 4: Thread the certified link witness into the Laurent diagnostic flow**

Replace `_laurent_native_ecp_boundary_stage_detail` with:

```julia
function _laurent_native_ecp_boundary_stage_detail(
    R;
    certified_descent_step::Bool = false,
    certified_link_witness::Bool = false,
)
    next_boundary = if certified_link_witness
        :laurent_endpoint_reduction
    elseif certified_descent_step
        :laurent_link_witness
    else
        nothing
    end
    return _column_reduction_stage_detail(
        :laurent_native_ecp_boundary,
        R,
        :staged_boundary;
        boundary = :laurent_native_ecp,
        requires_descent_measure = !certified_descent_step,
        certified_descent_scope = certified_descent_step ? :single_certified_step : nothing,
        next_boundary,
        requires_link_witness = !certified_link_witness,
        requires_endpoint_reduction = true,
        requires_laurent_normality_replay = true,
        requires_recursive_peel_integration = true,
        fallback_policy = :diagnostic_only,
    )
end
```

In `_diagnose_laurent_unimodular_column_reduction`, replace:

```julia
    descent_certificate = _laurent_descent_step_diagnostic_certificate(column, R)
    certified_descent_step = descent_certificate !== nothing
    if certified_descent_step
        push!(attempted, :laurent_descent_step_certificate)
        push!(
            details,
            _laurent_descent_step_certificate_stage_detail(R, descent_certificate),
        )
    end
    push!(attempted, :laurent_native_ecp_boundary)
    push!(
        details,
        _laurent_native_ecp_boundary_stage_detail(R; certified_descent_step),
    )
```

with:

```julia
    descent_certificate = _laurent_descent_step_diagnostic_certificate(column, R)
    certified_descent_step = descent_certificate !== nothing
    link_witness_certificate = nothing
    if certified_descent_step
        push!(attempted, :laurent_descent_step_certificate)
        push!(
            details,
            _laurent_descent_step_certificate_stage_detail(R, descent_certificate),
        )

        post_descent_column = _replay_laurent_elementary_entry_addition(
            column,
            R,
            descent_certificate.operation,
        )
        link_witness_certificate =
            _laurent_link_witness_diagnostic_certificate(post_descent_column, R)
        if link_witness_certificate !== nothing
            push!(attempted, :laurent_link_witness_certificate)
            push!(
                details,
                _laurent_link_witness_certificate_stage_detail(
                    R,
                    link_witness_certificate,
                ),
            )
        end
    end
    certified_link_witness = link_witness_certificate !== nothing
    push!(attempted, :laurent_native_ecp_boundary)
    push!(
        details,
        _laurent_native_ecp_boundary_stage_detail(
            R;
            certified_descent_step,
            certified_link_witness,
        ),
    )
```

- [ ] **Step 5: Run the focused test and verify it passes**

Run:

```bash
julia --project=. -e 'include("test/expert/laurent_native_ecp_boundary_diagnostics.jl")'
```

Expected: PASS. The d14 testset should report the certified link-witness stage and the negative controls should not report it.

- [ ] **Step 6: Run full package verification**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: PASS. Existing Julia world-age warnings from Quillen fixture loading may appear, but there should be no test failures.

- [ ] **Step 7: Run whitespace verification**

Run:

```bash
git diff --check
```

Expected: no output and exit code 0.

- [ ] **Step 8: Commit implementation**

Run:

```bash
git add src/algorithm/column_reduction.jl test/expert/laurent_native_ecp_boundary_diagnostics.jl docs/superpowers/plans/2026-07-09-issue-340-certified-laurent-link-witness-diagnostic.md
git commit -m "Expose issue 340 Laurent link witness diagnostic"
```

Expected: a commit containing only the implementation, focused test updates, and this plan.
