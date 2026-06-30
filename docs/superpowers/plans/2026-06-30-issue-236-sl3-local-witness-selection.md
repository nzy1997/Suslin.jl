# Issue 236 SL3 Local Witness Selection Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an internal verified witness-selection layer that turns checked #235 ordinary-polynomial `SL_3` contexts into local Murthy special-form witnesses, or staged diagnostics.

**Architecture:** Add one internal witness record and selector/verifier in `src/algorithm/factorization.jl`. Reuse `_sl3_local_target_entries` and add a tiny monicity-proof helper in `src/algorithm/sl3_local.jl`. Cover already-special-form, supplied variable-change replay, staged metadata shells, and negative controls in a focused expert test.

**Tech Stack:** Julia, Oscar polynomial rings and matrices, existing `SL3RealizationInputContext`, existing `SL3LocalMurthyInputContext` shape helpers, `Test`.

## Global Constraints

- Do not export new names from `src/Suslin.jl`.
- Do not route public `elementary_factorization` through this selector.
- Only accept verified #235 ordinary-polynomial determinant-one `SL_3` contexts.
- Accept already-special-form matrices only when `_sl3_local_target_entries` succeeds and extracted `p` is monic in the selected variable.
- Accept supplied variable-change or normality/conjugation metadata only when replay data names the context matrix, selected variable alignment, and a replayable local-form target or entries.
- Do not implement broad automatic coordinate-change search, Murthy certificates, Quillen local evidence, global patches, public route dispatch, ECP, or recursive `SL_n`.
- Missing evidence must produce a staged diagnostic containing `missing supported local-form witness`.
- Contradictory supplied evidence must throw `ArgumentError` or make the verifier return `false` before any route code can consume it.
- Focused verification command is `julia --project=. -e 'include("test/expert/park_woodburn_sl3_witness_selection.jl")'`.
- Package verification command is `julia --project=. -e 'using Pkg; Pkg.test()'`.

---

## File Structure

- Modify `src/algorithm/sl3_local.jl`: add `_sl3_local_monicity_witness(p, var_idx::Int, R)` beside `_is_monic_in_variable`.
- Modify `src/algorithm/factorization.jl`: add `SL3LocalFormWitnessSelection`, selector helpers, constructor, and verifier after the #235 `SL3RealizationInputContext` helpers.
- Create `test/expert/park_woodburn_sl3_witness_selection.jl`: red-green expert coverage for #236.
- Modify `test/runtests.jl`: register the new expert test near the #235 driver-context test.

### Task 1: Add Red Witness-Selection Tests

**Files:**
- Create: `test/expert/park_woodburn_sl3_witness_selection.jl`
- Modify: `test/runtests.jl`

**Interfaces:**
- Consumes future `Suslin._select_sl3_local_form_witness` and `Suslin._verify_sl3_local_form_witness_selection`.
- Produces focused tests for positive, staged, and negative witness-selection behavior.

- [ ] **Step 1: Write the failing test file**

Create `test/expert/park_woodburn_sl3_witness_selection.jl` with:

```julia
using Test
using Oscar
using Suslin

const PARK_WOODBURN_SL3_DRIVER_CATALOG_PATH =
    joinpath(@__DIR__, "..", "fixtures", "park_woodburn_sl3_driver_cases.jl")

function _selection_as_namedtuple(selection)
    names = propertynames(selection)
    return NamedTuple{names}(Tuple(getproperty(selection, name) for name in names))
end

function _corrupt_selection(selection, updates)
    return Suslin.SL3LocalFormWitnessSelection(
        values(merge(_selection_as_namedtuple(selection), updates))...,
    )
end

function _selection_corruption_is_rejected(selection, updates)
    try
        corrupted = _corrupt_selection(selection, updates)
        return !Suslin._verify_sl3_local_form_witness_selection(corrupted)
    catch err
        return err isa ArgumentError
    end
end

function _catalog_metadata(entry; expected_status = entry.expected_status)
    return (;
        fixture_id = entry.id,
        role = entry.role,
        expected_status,
        consumer_issue_ids = entry.consumer_issue_ids,
    )
end

function _special_form_matrix(R, p, q, r, s)
    return matrix(R, [
        p q zero(R);
        r s zero(R);
        zero(R) zero(R) one(R)
    ])
end

function _captured_error(f)
    try
        f()
        return nothing
    catch err
        return err
    end
end

if !isdefined(Main, :ParkWoodburnSL3DriverFixtureCatalog)
    include(PARK_WOODBURN_SL3_DRIVER_CATALOG_PATH)
end

@testset "Park-Woodburn SL3 local witness selection" begin
    entries = ParkWoodburnSL3DriverFixtureCatalog.cases_by_id()

    fast_entry = entries["sl3-driver-univariate-fast-local-qq"]
    fast_context = Suslin._sl3_realization_input_context(
        fast_entry.matrix;
        selected_variable = fast_entry.selected_variable,
        catalog_metadata = _catalog_metadata(fast_entry),
    )
    fast_selection = Suslin._select_sl3_local_form_witness(fast_context)
    @test fast_selection.support_status == :supported
    @test fast_selection.replay_status == :replayed
    @test fast_selection.witness_source == :already_special_form
    @test fast_selection.selected_variable == fast_entry.selected_variable.generator
    @test fast_selection.selected_variable_index == fast_entry.selected_variable.index
    @test fast_selection.entries ==
          (; p = fast_entry.matrix[1, 1], q = fast_entry.matrix[1, 2],
             r = fast_entry.matrix[2, 1], s = fast_entry.matrix[2, 2])
    @test fast_selection.monicity_witness.is_monic == true
    @test fast_selection.monicity_witness.variable == fast_entry.selected_variable.generator
    @test Suslin._verify_sl3_local_form_witness_selection(fast_selection)

    R, (X, y) = Oscar.polynomial_ring(QQ, ["X", "y"])
    p = X + y + 1
    q = one(R)
    r = y
    s = one(R)
    multivariate_target = _special_form_matrix(R, p, q, r, s)
    @test det(multivariate_target) == one(R)
    multivariate_context = Suslin._sl3_realization_input_context(
        multivariate_target;
        selected_variable = (; name = "X", generator = X, index = 1, status = :passes),
        catalog_metadata = (; fixture_id = "issue-236-multivariate-special-form"),
        local_form_witness = (; entries = (; p, q, r, s)),
    )
    multivariate_selection = Suslin._select_sl3_local_form_witness(multivariate_context)
    @test multivariate_selection.support_status == :supported
    @test multivariate_selection.entries == (; p, q, r, s)
    @test multivariate_selection.local_form_matrix == multivariate_target
    @test multivariate_selection.monicity_witness.degree == degree(p, 1)
    @test multivariate_selection.monicity_witness.leading_coefficient == one(R)
    @test Suslin._verify_sl3_local_form_witness_selection(multivariate_selection)

    staged_entry = entries["sl3-driver-det-one-no-witness-staged-qq"]
    staged_context = Suslin._sl3_realization_input_context(
        staged_entry.matrix;
        selected_variable = staged_entry.selected_variable,
        catalog_metadata = _catalog_metadata(staged_entry),
    )
    staged_selection = Suslin._select_sl3_local_form_witness(staged_context)
    @test staged_selection.support_status == :staged
    @test staged_selection.replay_status == :missing
    @test staged_selection.entries === nothing
    @test occursin(
        "missing supported local-form witness",
        staged_selection.staged_diagnostic.reason,
    )
    @test Suslin._verify_sl3_local_form_witness_selection(staged_selection)

    source_matrix = staged_entry.matrix
    variable_change_metadata = (;
        replay_id = "issue-236-variable-change-replay",
        source_matrix,
        selected_variable = staged_entry.selected_variable.generator,
        local_form_matrix = multivariate_target,
        replay_steps = ((; kind = :supplied_variable_change, name = :identity_catalog_replay),),
    )
    variable_context = Suslin._sl3_realization_input_context(
        source_matrix;
        selected_variable = staged_entry.selected_variable,
        catalog_metadata = (; fixture_id = "issue-236-supplied-variable-change"),
        variable_change_metadata,
    )
    variable_selection = Suslin._select_sl3_local_form_witness(variable_context)
    @test variable_selection.support_status == :supported
    @test variable_selection.replay_status == :replayed
    @test variable_selection.witness_source == :variable_change
    @test variable_selection.variable_change_status == :replayed
    @test variable_selection.entries == (; p, q, r, s)
    @test variable_selection.variable_change_metadata == variable_change_metadata
    @test Suslin._verify_sl3_local_form_witness_selection(variable_selection)

    normality_shell_context = Suslin._sl3_realization_input_context(
        source_matrix;
        selected_variable = staged_entry.selected_variable,
        catalog_metadata = (; fixture_id = "issue-236-normality-shell"),
        normality_conjugation_metadata = (;
            replay_id = "normality-shell-without-local-form",
            source_matrix,
            replay_steps = ((; kind = :normality_shell),),
        ),
    )
    normality_shell_selection =
        Suslin._select_sl3_local_form_witness(normality_shell_context)
    @test normality_shell_selection.support_status == :staged
    @test normality_shell_selection.normality_conjugation_status == :recorded
    @test :normality_conjugation in
          normality_shell_selection.staged_diagnostic.partial_evidence
    @test Suslin._verify_sl3_local_form_witness_selection(normality_shell_selection)

    @test_throws ArgumentError Suslin._select_sl3_local_form_witness(
        multivariate_context;
        selected_variable = X + one(R),
    )

    nonmonic_p = 2 * X + y + one(R)
    nonmonic_target = _special_form_matrix(R, nonmonic_p, q, r, s)
    nonmonic_context = Suslin._sl3_realization_input_context(
        nonmonic_target;
        selected_variable = X,
        catalog_metadata = (; fixture_id = "issue-236-nonmonic-local-form"),
    )
    nonmonic_error =
        _captured_error(() -> Suslin._select_sl3_local_form_witness(nonmonic_context))
    @test nonmonic_error isa ArgumentError
    @test occursin("local-form p is not monic", sprint(showerror, nonmonic_error))

    @test_throws ArgumentError Suslin._select_sl3_local_form_witness(
        multivariate_context;
        local_form_witness = (; entries = (; p = p + one(R), q, r, s)),
    )

    @test_throws ArgumentError Suslin._select_sl3_local_form_witness(
        variable_context;
        variable_change_metadata = merge(
            variable_change_metadata,
            (; source_matrix = identity_matrix(parent(source_matrix[1, 1]), 3),),
        ),
    )

    missing_payload_context = Suslin._sl3_realization_input_context(
        source_matrix;
        selected_variable = staged_entry.selected_variable,
        catalog_metadata = (; fixture_id = "issue-236-missing-variable-payload"),
        variable_change_metadata = (;
            replay_id = "missing-variable-change-payload",
            source_matrix,
            selected_variable = staged_entry.selected_variable.generator,
            local_form_matrix = multivariate_target,
        ),
    )
    missing_payload_selection =
        Suslin._select_sl3_local_form_witness(missing_payload_context)
    @test missing_payload_selection.support_status == :staged
    @test missing_payload_selection.variable_change_status == :recorded

    @test _selection_corruption_is_rejected(
        multivariate_selection,
        (; entries = (; p = p + one(R), q, r, s)),
    )
    @test _selection_corruption_is_rejected(
        multivariate_selection,
        (; monicity_witness = merge(
            multivariate_selection.monicity_witness,
            (; is_monic = false,),
        )),
    )
    @test _selection_corruption_is_rejected(
        variable_selection,
        (; variable_change_status = :missing,),
    )
end
```

- [ ] **Step 2: Register the expert test**

Add this file to the expert group in `test/runtests.jl` immediately after
`"expert/park_woodburn_sl3_driver_context.jl"`:

```julia
"expert/park_woodburn_sl3_witness_selection.jl",
```

- [ ] **Step 3: Run RED verification**

Run:

```bash
julia --project=. -e 'include("test/expert/park_woodburn_sl3_witness_selection.jl")'
```

Expected: FAIL with `UndefVarError: SL3LocalFormWitnessSelection not defined`
or `UndefVarError: _select_sl3_local_form_witness not defined`.

### Task 2: Add Monicity Proof Helper

**Files:**
- Modify: `src/algorithm/sl3_local.jl`
- Test: `test/expert/park_woodburn_sl3_witness_selection.jl`

**Interfaces:**
- Produces `_sl3_local_monicity_witness(p, var_idx::Int, R)`.
- Consumes existing `_sl3_local_coefficient_in_variable_degree` and `_is_monic_in_variable`.

- [ ] **Step 1: Add the helper near `_is_monic_in_variable`**

Insert before `_is_monic_in_variable`:

```julia
function _sl3_local_monicity_witness(p, var_idx::Int, R)
    variable_count = length(collect(gens(R)))
    if var_idx < 1 || var_idx > variable_count || iszero(p)
        return (;
            variable = var_idx >= 1 && var_idx <= variable_count ? collect(gens(R))[var_idx] : nothing,
            variable_index = var_idx,
            degree = -1,
            leading_coefficient = zero(R),
            is_monic = false,
        )
    end

    target_degree = degree(p, var_idx)
    leading = target_degree < 0 ?
        zero(R) :
        _sl3_local_coefficient_in_variable_degree(p, var_idx, target_degree, R)
    return (;
        variable = collect(gens(R))[var_idx],
        variable_index = var_idx,
        degree = target_degree,
        leading_coefficient = leading,
        is_monic = leading == one(R),
    )
end
```

- [ ] **Step 2: Run focused test to keep RED scoped to selector**

Run:

```bash
julia --project=. -e 'include("test/expert/park_woodburn_sl3_witness_selection.jl")'
```

Expected: still FAIL on missing selector names, not on monicity helper syntax.

### Task 3: Implement Local Witness Selection

**Files:**
- Modify: `src/algorithm/factorization.jl`
- Test: `test/expert/park_woodburn_sl3_witness_selection.jl`

**Interfaces:**
- Produces `SL3LocalFormWitnessSelection`, `_select_sl3_local_form_witness`, and `_verify_sl3_local_form_witness_selection`.
- Consumes `SL3RealizationInputContext`, `_verify_sl3_realization_input_context`, `_sl3_local_target_entries`, `_sl3_local_monicity_witness`, `_sl3_realization_input_context_selected_variable`, and `_sl3_realization_input_context_has_replay_payload`.

- [ ] **Step 1: Add the record type**

Insert after `SL3RealizationInputContext`:

```julia
struct SL3LocalFormWitnessSelection
    context::SL3RealizationInputContext
    selected_variable
    selected_variable_index
    selected_variable_name
    entries
    local_form_matrix
    monicity_witness::NamedTuple
    local_form_witness
    variable_change_metadata
    variable_change_status::Symbol
    normality_conjugation_metadata
    normality_conjugation_status::Symbol
    replay_status::Symbol
    support_status::Symbol
    witness_source::Symbol
    staged_diagnostic::NamedTuple
    verification
end
```

- [ ] **Step 2: Add helper functions**

Add helpers after `_sl3_realization_input_context`:

```julia
function _sl3_local_witness_hint_or_context(hint, stored)
    return hint === nothing ? stored : hint
end

function _sl3_local_witness_selected_variable(context, hint)
    R = context.base_ring
    selected_hint = _sl3_local_witness_hint_or_context(hint, context.selected_variable)
    selected, selected_index, status =
        _sl3_realization_input_context_selected_variable(R, selected_hint)
    status == :passes ||
        throw(ArgumentError("SL_3 local witness selection requires a selected variable"))
    if context.selected_variable !== nothing && selected != context.selected_variable
        throw(ArgumentError("selected variable hint does not match the SL_3 context"))
    end
    if context.selected_variable_index !== nothing &&
            selected_index != context.selected_variable_index
        throw(ArgumentError("selected variable index does not match the SL_3 context"))
    end
    selected_name = string(collect(gens(R))[selected_index])
    return selected, selected_index, selected_name
end

function _sl3_local_witness_entries_from_data(data)
    data === nothing && return nothing
    entries = _sl3_realization_input_context_extract(data, (:entries, :local_form_entries))
    entries !== nothing && return entries
    if all(field -> hasproperty(data, field), (:p, :q, :r, :s))
        return (; p = data.p, q = data.q, r = data.r, s = data.s)
    end
    return nothing
end

function _sl3_local_witness_matrix_from_entries(R, entries)
    entries === nothing && return nothing
    return matrix(R, [
        entries.p entries.q zero(R);
        entries.r entries.s zero(R);
        zero(R) zero(R) one(R)
    ])
end

function _sl3_local_witness_extract_matrix(data)
    return _sl3_realization_input_context_extract(
        data,
        (:local_form_matrix, :special_form_matrix, :transformed_matrix, :local_target),
    )
end

function _sl3_local_witness_source_matrix(data)
    return _sl3_realization_input_context_extract(
        data,
        (:source_matrix, :input_matrix, :context_matrix, :original_matrix),
    )
end
```

Also add `_sl3_local_witness_checked_entries`, `_sl3_local_witness_metadata_status`,
`_sl3_local_form_witness_selection_fields`, and verifier helpers with this
behavior:

- `_sl3_local_witness_checked_entries(R, matrix_value, selected_variable,
  selected_index, source_label)` calls `_sl3_local_target_entries(matrix_value)`.
  It throws if the shape is not special form, if `det(matrix_value) != one(R)`,
  or if `_sl3_local_monicity_witness(entries.p, selected_index, R).is_monic`
  is false.
- It returns `(entries, matrix_value, monicity_witness)`.
- `_sl3_local_witness_metadata_status(context, metadata, selected_variable,
  selected_index, source_label)` returns `(status, entries, matrix, monicity)`
  where status is `:missing`, `:recorded`, or `:replayed`.
- Metadata is `:replayed` only when `_sl3_realization_input_context_has_replay_payload(metadata)`
  is true, its source matrix equals `context.matrix`, its selected variable
  field equals `selected_variable`, and it has local-form matrix or local-form
  entries that pass `_sl3_local_witness_checked_entries`.
- Metadata with a replay payload but no local-form target is `:recorded`.
- Metadata with local-form target but no replay payload is `:recorded`.
- Source matrix mismatch, selected variable mismatch, or malformed local-form
  target throws `ArgumentError`.

- [ ] **Step 3: Add constructor and verifier**

Implement:

```julia
function _select_sl3_local_form_witness(
    context::SL3RealizationInputContext;
    selected_variable = nothing,
    local_form_witness = nothing,
    variable_change_metadata = nothing,
    normality_conjugation_metadata = nothing,
)
    fields = _sl3_local_form_witness_selection_fields(
        context;
        selected_variable,
        local_form_witness,
        variable_change_metadata,
        normality_conjugation_metadata,
    )
    unchecked = SL3LocalFormWitnessSelection(values(merge(fields, (; verification = nothing,)))...)
    verification = _sl3_local_form_witness_selection_core_verification(unchecked)
    checked = SL3LocalFormWitnessSelection(values(merge(fields, (; verification,)))...)
    _verify_sl3_local_form_witness_selection(checked) ||
        error("internal SL_3 local witness selection verification failed")
    return checked
end

function _verify_sl3_local_form_witness_selection(selection)::Bool
    try
        return _sl3_local_form_witness_selection_verification(selection).overall_ok
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end
```

The field builder must choose the first supported source in this order:
already-special-form context matrix, replayed variable-change metadata,
replayed normality/conjugation metadata. If none is supported, it returns a
staged record with `entries = nothing`, `local_form_matrix = nothing`,
`replay_status = :missing`, `support_status = :staged`,
`witness_source = :staged`, and diagnostic reason
`"missing supported local-form witness"`.

The verifier recomputes fields using the stored context and stored metadata.
For supported records it must compare entries, matrix, monicity witness,
source statuses, replay/support status, witness source, staged diagnostic, and
stored verification. For staged records it must verify the staged fields
recompute exactly.

- [ ] **Step 4: Run GREEN focused verification**

Run:

```bash
julia --project=. -e 'include("test/expert/park_woodburn_sl3_witness_selection.jl")'
```

Expected: PASS.

### Task 4: Package Verification And Review

**Files:**
- Verify: `src/algorithm/factorization.jl`
- Verify: `src/algorithm/sl3_local.jl`
- Verify: `test/expert/park_woodburn_sl3_witness_selection.jl`
- Verify: `test/runtests.jl`

**Interfaces:**
- Confirms issue #236 behavior and package integration.

- [ ] **Step 1: Run focused verification**

Run:

```bash
julia --project=. -e 'include("test/expert/park_woodburn_sl3_witness_selection.jl")'
```

Expected: PASS.

- [ ] **Step 2: Run package verification**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: PASS.

- [ ] **Step 3: Run diff hygiene**

Run:

```bash
git diff --check
```

Expected: no output and exit code 0.

- [ ] **Step 4: Inspect final diff**

Run:

```bash
git status --short
git diff --stat origin/main...HEAD
```

Expected: no generated artifacts; diff limited to issue #236 spec/plan, source
helpers, focused expert test, and test registration.

## Plan Self-Review

- The plan covers each issue verification case and negative control.
- The implementation tasks follow TDD: red test, minimal helper, selector, then package verification.
- The selected-variable and replay-evidence boundaries are explicit.
- The plan does not add public dispatch, Murthy certificates, Quillen patches, coordinate-change search, ECP, or recursive `SL_n`.
