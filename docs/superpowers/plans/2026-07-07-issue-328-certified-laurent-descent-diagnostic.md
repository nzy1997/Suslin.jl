# Issue 328 Certified Laurent Descent Diagnostic Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Expose a conservative certified Laurent descent-step diagnostic stage before the remaining Laurent-native ECP boundary.

**Architecture:** Extend `src/algorithm/column_reduction.jl` with private diagnostic helpers that build and validate the single recorded d14 entry-addition certificate through existing internal Laurent descent helpers. Update the existing expert diagnostic test to assert the new stage, the revised boundary fields, and negative controls.

**Tech Stack:** Julia, Oscar Laurent polynomial rings, `Test`, existing `Suslin` internals.

## Global Constraints

- Do not expose a public API or add exports.
- Do not claim full Laurent ECP support or a complete Laurent descent algorithm.
- Production code must not key on a fixture module or test-only helper.
- The certified diagnostic stage must appear only after the column matches the recorded support evidence and validates through `_validate_laurent_descent_step_certificate`.
- The certified operation is `family = :entry_addition`, `target_index = 1`, `source_index = 2`, `coefficient = 1`, `exponent = (-1, 1)`, and `ring_generators = ("u", "v")`.
- Certified stage detail must include `outcome = :certified_descent_step`, `descent_scope = :single_certified_step`, `measure_relation = :strict_decrease`, `replay_status = :ok`, and `next_boundary = :laurent_link_witness`.
- The certified d14 terminal boundary must set `requires_descent_measure = false` while keeping link witness, endpoint reduction, Laurent normality replay, and recursive peel integration requirements true.
- Generic unsupported Laurent columns must continue to reach the generic Laurent-native ECP boundary without a certified descent-step diagnostic.
- Out of scope: Laurent link witnesses, endpoint reductions, normality/conjugation replay, determinant normalization, recursive peel integration, and full `case_008` success.

---

### Task 1: Certified Descent Diagnostic Stage

**Files:**
- Modify: `src/algorithm/column_reduction.jl`
- Modify: `test/expert/laurent_native_ecp_boundary_diagnostics.jl`
- Modify: `test/runtests.jl` only if a new expert test file is created; otherwise leave registration unchanged.

**Interfaces:**
- Consumes: `Suslin._laurent_descent_measure_from_column(column, R; case_id = nothing)::NamedTuple`
- Consumes: `Suslin._replay_laurent_elementary_entry_addition(column, R, operation)::Vector`
- Consumes: `Suslin._strictly_decreases_laurent_measure(before, after)::Bool`
- Consumes: `Suslin._validate_laurent_descent_step_certificate(cert, column, R)::Symbol`
- Produces: diagnostic attempted stage `:laurent_descent_step_certificate` with plain-data detail fields from the issue.
- Produces: revised `_laurent_native_ecp_boundary_stage_detail(R; certified_descent_step = false)` boundary detail.

- [ ] **Step 1: Write the failing diagnostic tests**

Update `test/expert/laurent_native_ecp_boundary_diagnostics.jl` so the d14 block asserts:

```julia
descent_idx = findfirst(==(:laurent_descent_step_certificate), d14.attempted_stages)
@test descent_idx !== nothing
@test boundary_idx > descent_idx
descent = _diagnostic_stage_detail(d14, :laurent_descent_step_certificate)
@test descent.outcome == :certified_descent_step
@test descent.descent_scope == :single_certified_step
@test descent.operation_family == :entry_addition
@test descent.target_index == 1
@test descent.source_index == 2
@test descent.coefficient == 1
@test descent.exponent == (-1, 1)
@test descent.measure_relation == :strict_decrease
@test descent.replay_status == :ok
@test descent.next_boundary == :laurent_link_witness
@test Suslin._strictly_decreases_laurent_measure(
    descent.before_measure,
    descent.after_measure,
)
```

Update the d14 boundary assertions to expect `requires_descent_measure == false`
and `certified_descent_scope == :single_certified_step`. Keep the existing
`requires_link_witness`, `requires_endpoint_reduction`,
`requires_laurent_normality_replay`, and
`requires_recursive_peel_integration` assertions true.

Add negative assertions:

```julia
@test !(:laurent_descent_step_certificate in d15.attempted_stages)
@test !(:laurent_descent_step_certificate in non_unimodular_diagnostic.attempted_stages)

tampered_column = copy(d14_fixture.failing_column)
tampered_column[1] = tampered_column[1] + one(d14_fixture.ring)
tampered = Suslin.diagnose_unimodular_column_reduction(
    tampered_column,
    d14_fixture.ring;
    assume_unimodular = true,
    laurent_large_support_diagnostic_decline = true,
)
@test tampered.status == :unsupported
@test !(:laurent_descent_step_certificate in tampered.attempted_stages)
tampered_boundary = _diagnostic_stage_detail(tampered, :laurent_native_ecp_boundary)
@test tampered_boundary.requires_descent_measure == true
```

- [ ] **Step 2: Run the red test and confirm it fails for missing stage**

Run:

```bash
julia --project=. -e 'include("test/expert/laurent_native_ecp_boundary_diagnostics.jl")'
```

Expected: FAIL because `:laurent_descent_step_certificate` is absent and the
d14 boundary still reports `requires_descent_measure == true`.

- [ ] **Step 3: Implement the private diagnostic helpers**

In `src/algorithm/column_reduction.jl`, add a private recorded operation:

```julia
const _LAURENT_D14_CERTIFIED_DESCENT_OPERATION = (;
    family = :entry_addition,
    target_index = 1,
    source_index = 2,
    coefficient = 1,
    exponent = (-1, 1),
    ring_generators = ("u", "v"),
)
```

Add compact recorded support fingerprints for the certified before and after
columns. The before fingerprint is:

```julia
(;
    entry_support_counts = (3573, 3574, 3554, 3561, 3595, 3734, 3622, 3454, 3489, 3692, 3675, 3600, 3693, 3495),
    support_checksum = 0x9fba30357fde4b25,
)
```

The after fingerprint is:

```julia
(;
    entry_support_counts = (3661, 3574, 3554, 3561, 3595, 3734, 3622, 3454, 3489, 3692, 3675, 3600, 3693, 3495),
    support_checksum = 0x6230882975e8c766,
)
```

Add `_laurent_descent_column_support_fingerprint(column)` using sorted Laurent
supports and a stable `UInt64` FNV-style checksum over entry indices, support
lengths, and exponent pairs. Then add a helper that attempts certificate
construction by checking the recorded support evidence and recomputing before
and after measures from the supplied column:

```julia
function _laurent_descent_step_diagnostic_certificate(column::AbstractVector, R)
    try
        operation = _LAURENT_D14_CERTIFIED_DESCENT_OPERATION
        _laurent_descent_operation_status(operation, length(column), R) == :ok ||
            return nothing
        _laurent_descent_column_support_fingerprint(column) ==
            _LAURENT_D14_CERTIFIED_DESCENT_BEFORE_FINGERPRINT || return nothing
        before_measure = _laurent_descent_measure_from_column(column, R; case_id = "case_008")
        after_column = _replay_laurent_elementary_entry_addition(column, R, operation)
        _laurent_descent_column_support_fingerprint(after_column) ==
            _LAURENT_D14_CERTIFIED_DESCENT_AFTER_FINGERPRINT || return nothing
        after_measure = _laurent_descent_measure_from_column(after_column, R; case_id = "case_008")
        _strictly_decreases_laurent_measure(before_measure, after_measure) ||
            return nothing
        certificate = (;
            case_id = "case_008",
            dimension = length(column),
            ring_generators = _laurent_descent_ring_generators(R),
            operation,
            before_measure,
            after_measure,
            status = :descent_step_certificate,
            replay_status = :ok,
            measure_relation = :strict_decrease,
        )
        _validate_laurent_descent_step_certificate(certificate, column, R) == :ok ||
            return nothing
        return certificate
    catch err
        err isa InterruptException && rethrow()
        return nothing
    end
end
```

Add `_laurent_descent_step_certificate_stage_detail(R, certificate)` that returns
`_column_reduction_stage_detail(:laurent_descent_step_certificate, R,
:certified_descent_step; ...)` with the required fields copied from
`certificate.operation` and `certificate`.

Update `_laurent_native_ecp_boundary_stage_detail` to accept
`; certified_descent_step::Bool = false` and set:

```julia
requires_descent_measure = !certified_descent_step
certified_descent_scope = certified_descent_step ? :single_certified_step : nothing
next_boundary = certified_descent_step ? :laurent_link_witness : nothing
```

Keep every other existing boundary field unchanged.

- [ ] **Step 4: Insert the stage in the diagnostic flow**

In `_diagnose_laurent_row_preconditioning`, after pushing the
`:no_row_preconditioning_candidate` row-preconditioning detail and before
pushing `:laurent_native_ecp_boundary`, call the helper:

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
    _laurent_native_ecp_boundary_stage_detail(
        R;
        certified_descent_step,
    ),
)
```

- [ ] **Step 5: Run the focused verification**

Run:

```bash
julia --project=. -e 'include("test/expert/laurent_native_ecp_boundary_diagnostics.jl")'
```

Expected: PASS, including the d14 certified stage and negative controls.

- [ ] **Step 6: Run the package verification**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: PASS. Existing Julia 1.12 world-age warnings from Quillen fixture
catalog loading may appear and are not failures if the exit code is 0.

- [ ] **Step 7: Commit**

Run:

```bash
git add docs/superpowers/specs/2026-07-07-issue-328-certified-laurent-descent-diagnostic-design.md docs/superpowers/plans/2026-07-07-issue-328-certified-laurent-descent-diagnostic.md src/algorithm/column_reduction.jl test/expert/laurent_native_ecp_boundary_diagnostics.jl test/runtests.jl
git commit -m "Expose certified Laurent descent diagnostic"
```

## Self Review

The plan covers the issue interface, exact stage detail fields, boundary detail
change, required negative controls, focused verification command, and package
verification command. It contains no placeholders and keeps implementation scope
to one diagnostic stage.
