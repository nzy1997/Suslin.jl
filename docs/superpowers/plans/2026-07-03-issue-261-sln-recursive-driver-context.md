# Issue 261 SLn Recursive Driver Context Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an internal checked `SL_n` recursive driver input context and verifier for the staged #186 ordinary-polynomial driver boundary.

**Architecture:** Follow the existing checked-context pattern used by `SL3RealizationInputContext`: compute fields from inputs, store a verification summary, and verify by recomputing instead of trusting stored flags. Place the internal type and helpers in `src/algorithm/polynomial_column_peel.jl`, with expert tests backed by the #260 `park_woodburn_sln_driver_cases.jl` catalog.

**Tech Stack:** Julia, Oscar, Suslin internal matrix/ring helpers, Test stdlib.

## Global Constraints

- Repository has no `AGENTS.md`; follow README test commands.
- Keep the context internal; do not export it and do not route public `elementary_factorization` through it.
- Modify `src/algorithm/polynomial_column_peel.jl`.
- Add `test/expert/park_woodburn_sln_driver_context.jl`.
- Register the expert test in `test/runtests.jl`.
- Build on `_validate_factorization_matrix`, `_factorization_ring_profile`, `_polynomial_exact_field_backed_ring`, and `_require_polynomial_sl_determinant`.
- Recompute every stored field in the verifier, including last column, determinant status, ring profile, generator metadata, route provenance, staged diagnostic, and staged reason code.
- Staged reason codes must include `:missing_ecp_evidence`, `:missing_final_sl3_route`, `:missing_variable_metadata`, `:unsupported_coefficient_ring`, and `:determinant_not_one`.
- Do not reduce the last column, produce elementary factors, call #184 `SL_3` routes, change public route selection, or broaden Laurent/ToricBuilder support.
- Run `julia --project=. -e 'include("test/expert/park_woodburn_sln_driver_context.jl")'`.
- Run `julia --project=. -e 'using Pkg; Pkg.test()'`.

---

## File Structure

- Modify `src/algorithm/polynomial_column_peel.jl`: define `SLnRecursiveDriverInputContext`, constructor, recomputation helpers, staged diagnostic helper, and verifier.
- Create `test/expert/park_woodburn_sln_driver_context.jl`: #260-backed context tests and tamper gates.
- Modify `test/runtests.jl`: include the expert test near the existing Park-Woodburn column-peel and SL3 driver tests.

### Task 1: Red Expert Context Contract

**Files:**
- Create: `test/expert/park_woodburn_sln_driver_context.jl`
- Modify: `test/runtests.jl`

**Interfaces:**
- Consumes future `Suslin._sln_recursive_driver_input_context(A; kwargs...)`.
- Consumes future `Suslin._verify_sln_recursive_driver_input_context(context)`.
- Produces an expert regression contract for supported, staged, legacy, unsupported, and tampered contexts.

- [ ] **Step 1: Write the failing context test**

Create `test/expert/park_woodburn_sln_driver_context.jl` with helper functions:

```julia
using Test
using Oscar
using Suslin

const PARK_WOODBURN_SLN_DRIVER_CATALOG_PATH =
    joinpath(@__DIR__, "..", "fixtures", "park_woodburn_sln_driver_cases.jl")

function _sln_ctx_as_namedtuple(context)
    names = propertynames(context)
    return NamedTuple{names}(Tuple(getproperty(context, name) for name in names))
end

function _sln_ctx_replace_field(context, field::Symbol, value)
    fields = fieldnames(typeof(context))
    idx = findfirst(==(field), fields)
    idx === nothing && error("unknown SLn context field $(field)")
    values = [getfield(context, name) for name in fields]
    values[idx] = value
    return Suslin.SLnRecursiveDriverInputContext(values...)
end

function _sln_context_from_entry(entry; variable_order = entry.ring.generators)
    return Suslin._sln_recursive_driver_input_context(
        entry.matrix;
        variable_order = variable_order,
        selected_variable = isempty(entry.ring.generators) ? nothing : entry.ring.generators[1],
        ecp_witness_metadata = entry.peel_steps[1].last_column_ecp,
        final_route_metadata = entry.final_route,
        route_provenance_metadata = entry.route_provenance,
        catalog_id = entry.id,
    )
end
```

Add tests that load `ParkWoodburnSLnDriverFixtureCatalog.cases_by_id()` and assert:

```julia
@testset "Park-Woodburn SLn recursive driver context" begin
    if !isdefined(Main, :ParkWoodburnSLnDriverFixtureCatalog)
        include(PARK_WOODBURN_SLN_DRIVER_CATALOG_PATH)
    end
    entries = ParkWoodburnSLnDriverFixtureCatalog.cases_by_id()
    catalog = ParkWoodburnSLnDriverFixtureCatalog.catalog()
    negative = Dict(entry.id => entry for entry in catalog.negative_controls)

    mainline = entries["sln-driver-sl4-gf2-ecp-mainline"]
    mainline_ctx = _sln_context_from_entry(mainline)
    @test mainline_ctx.dimension == 4
    @test mainline_ctx.support_classification == :supported
    @test mainline_ctx.staged_reason_code === nothing
    @test mainline_ctx.ecp_evidence_status == :replayed
    @test mainline_ctx.final_route_evidence_status == :replayed
    @test mainline_ctx.route_provenance_status == :recorded
    @test mainline_ctx.last_column == [mainline.matrix[row, 4] for row in 1:4]
    @test mainline_ctx.determinant_status == :one
    @test mainline_ctx.exact_field_status == :supported
    @test Suslin._verify_sln_recursive_driver_input_context(mainline_ctx)

    multistep = entries["sln-driver-sl5-gf2-two-step"]
    multistep_ctx = _sln_context_from_entry(multistep)
    @test multistep_ctx.dimension == 5
    @test multistep_ctx.last_column == [multistep.matrix[row, 5] for row in 1:5]
    @test multistep_ctx.support_classification == :supported
    @test Suslin._verify_sln_recursive_driver_input_context(multistep_ctx)

    legacy = entries["sln-driver-legacy-recursive-column-peel-qq"]
    legacy_ctx = _sln_context_from_entry(legacy)
    @test legacy_ctx.support_classification == :staged
    @test legacy_ctx.staged_reason_code == :missing_ecp_evidence
    @test legacy_ctx.ecp_evidence_status == :missing
    @test Suslin._verify_sln_recursive_driver_input_context(legacy_ctx)

    staged = entries["sln-driver-staged-missing-final-sl3-qq"]
    staged_ctx = _sln_context_from_entry(staged)
    @test staged_ctx.support_classification == :staged
    @test staged_ctx.staged_reason_code == :missing_final_sl3_route
    @test staged_ctx.ecp_evidence_status == :replayed
    @test staged_ctx.final_route_evidence_status == :missing
    @test Suslin._verify_sln_recursive_driver_input_context(staged_ctx)

    unsupported = negative["sln-driver-negative-unsupported-coefficient-ring"]
    unsupported_ctx = _sln_context_from_entry(unsupported)
    @test unsupported_ctx.support_classification == :staged
    @test unsupported_ctx.staged_reason_code == :unsupported_coefficient_ring
    @test unsupported_ctx.exact_field_status == :unsupported
    @test Suslin._verify_sln_recursive_driver_input_context(unsupported_ctx)

    missing_variable_ctx = Suslin._sln_recursive_driver_input_context(
        staged.matrix;
        variable_order = nothing,
        ecp_witness_metadata = staged.peel_steps[1].last_column_ecp,
        final_route_metadata = staged.final_route,
        route_provenance_metadata = staged.route_provenance,
        catalog_id = staged.id,
    )
    @test missing_variable_ctx.staged_reason_code == :missing_variable_metadata

    det_bad = negative["sln-driver-negative-det-not-one"]
    det_bad_ctx = Suslin._sln_recursive_driver_input_context(
        det_bad.matrix;
        variable_order = det_bad.ring.generators,
        ecp_witness_metadata = det_bad.peel_steps[1].last_column_ecp,
        final_route_metadata = det_bad.final_route,
        route_provenance_metadata = det_bad.route_provenance,
        catalog_id = det_bad.id,
    )
    @test det_bad_ctx.staged_reason_code == :determinant_not_one

    for (field, value) in (
        (:determinant_status, :not_one),
        (:ring_profile, :tampered),
        (:generators, reverse(mainline_ctx.generators)),
        (:last_column, reverse(mainline_ctx.last_column)),
        (:route_provenance_metadata, merge(mainline_ctx.route_provenance_metadata, (; source = "tampered"))),
        (:staged_reason_code, :missing_ecp_evidence),
        (:staged_diagnostic, merge(mainline_ctx.staged_diagnostic, (; message = "tampered"))),
        (:verification, merge(mainline_ctx.verification, (; determinant_status_ok = false))),
    )
        @test !Suslin._verify_sln_recursive_driver_input_context(
            _sln_ctx_replace_field(mainline_ctx, field, value),
        )
    end

    @test_throws ArgumentError Suslin._sln_recursive_driver_input_context(identity_matrix(mainline.base_ring, 2))
end
```

- [ ] **Step 2: Register the expert test**

Add this file to the `expert` group in `test/runtests.jl` after `expert/park_woodburn_polynomial_column_peel.jl`:

```julia
"expert/park_woodburn_sln_driver_context.jl",
```

- [ ] **Step 3: Run RED verification**

Run:

```bash
julia --project=. -e 'include("test/expert/park_woodburn_sln_driver_context.jl")'
```

Expected: failure with `UndefVarError` or missing field errors for `SLnRecursiveDriverInputContext` / `_sln_recursive_driver_input_context`.

- [ ] **Step 4: Commit the red test**

```bash
git add test/expert/park_woodburn_sln_driver_context.jl test/runtests.jl
git commit -m "test: add issue 261 SLn driver context contract"
```

### Task 2: Context Implementation

**Files:**
- Modify: `src/algorithm/polynomial_column_peel.jl`

**Interfaces:**
- Consumes: `ECPColumnReductionCertificate`, `verify_ecp_column_reduction`, `_validate_factorization_matrix`, `_factorization_ring_profile`, `_polynomial_exact_field_backed_ring`, `_require_polynomial_sl_determinant`.
- Produces: `SLnRecursiveDriverInputContext`, `_sln_recursive_driver_input_context`, `_verify_sln_recursive_driver_input_context`.

- [ ] **Step 1: Add the internal type and constants**

Insert before `struct PolynomialColumnPeelStep`:

```julia
struct SLnRecursiveDriverInputContext
    matrix
    base_ring
    coefficient_ring
    dimension::Int
    ring_profile::Symbol
    exact_field_status::Symbol
    determinant
    determinant_status::Symbol
    generators::Tuple
    generator_names::Tuple
    variable_order
    variable_order_status::Symbol
    selected_variable
    selected_variable_index
    selected_variable_status::Symbol
    last_column::Vector
    last_column_profile::NamedTuple
    initial_unimodularity_witness
    initial_unimodularity_witness_status::Symbol
    ecp_witness_metadata
    ecp_evidence_status::Symbol
    final_route_metadata
    final_route_evidence_status::Symbol
    route_provenance_metadata::NamedTuple
    route_provenance_status::Symbol
    catalog_id
    support_classification::Symbol
    staged_reason_code
    staged_diagnostic::NamedTuple
    verification
end

const _SLN_RECURSIVE_DRIVER_STAGED_REASON_CODES = Set([
    :missing_ecp_evidence,
    :missing_final_sl3_route,
    :missing_variable_metadata,
    :unsupported_coefficient_ring,
    :determinant_not_one,
])
```

- [ ] **Step 2: Add field recomputation helpers**

Add helpers `_sln_recursive_driver_selected_variable`, `_sln_recursive_driver_variable_order`,
`_sln_recursive_driver_determinant_status`, `_sln_recursive_driver_last_column_profile`,
`_sln_recursive_driver_unimodularity_witness`, `_sln_recursive_driver_ecp_status`,
`_sln_recursive_driver_final_route_status`, `_sln_recursive_driver_route_provenance`,
and `_sln_recursive_driver_staged_diagnostic`.

Required behavior:

```julia
_sln_recursive_driver_determinant_status(A) =
    try
        _require_polynomial_sl_determinant(A)
        (:one, det(A))
    catch err
        err isa InterruptException && rethrow()
        err isa ArgumentError || rethrow()
        (:not_one, det(A))
    end
```

`_sln_recursive_driver_ecp_status(metadata, column, R, n)` must return
`:replayed` only for a verified `ECPColumnReductionCertificate` whose ring,
original column, and final target column match. It returns `:missing` for
`nothing`, `status = :missing`, `status = :absent`, or certificate mismatch.
It returns `:recorded` for shell metadata that has identifiers but no verified
certificate.

`_sln_recursive_driver_final_route_status(metadata, R)` must return
`:replayed` only when metadata has `status = :replayed`, has a non-empty route
identifier or replay payload, and any provided matrix is a `3 x 3` matrix over
`R`. It returns `:missing` for `nothing`, `status = :missing`, or `status =
:absent`, and `:recorded` otherwise.

- [ ] **Step 3: Add constructor fields and verifier**

Implement:

```julia
function _sln_recursive_driver_input_context_fields(
    A;
    variable_order = :auto,
    selected_variable = nothing,
    ecp_witness_metadata = nothing,
    final_route_metadata = nothing,
    route_provenance_metadata = nothing,
    catalog_id = nothing,
)
    dimension = _validate_factorization_matrix(A)
    dimension >= 4 || throw(ArgumentError("SL_n recursive driver context requires size at least 4"))
    R = base_ring(A)
    ring_profile = _factorization_ring_profile(R)
    ring_profile == :polynomial ||
        throw(ArgumentError("SL_n recursive driver context requires an ordinary polynomial base ring"))
    coefficient = coefficient_ring(R)
    exact_field_status = _polynomial_exact_field_backed_ring(R) ? :supported : :unsupported
    generators = Tuple(collect(gens(R)))
    generator_names = Tuple(string(generator) for generator in generators)
    normalized_order, variable_order_status =
        _sln_recursive_driver_variable_order(R, variable_order)
    selected, selected_index, selected_status =
        _sln_recursive_driver_selected_variable(R, selected_variable, normalized_order)
    determinant_status, determinant_value = _sln_recursive_driver_determinant_status(A)
    last_column = [A[row, dimension] for row in 1:dimension]
    last_column_profile = _sln_recursive_driver_last_column_profile(last_column, R)
    witness, witness_status =
        _sln_recursive_driver_unimodularity_witness(last_column, R, exact_field_status, determinant_status)
    ecp_status = _sln_recursive_driver_ecp_status(ecp_witness_metadata, last_column, R, dimension)
    final_status = _sln_recursive_driver_final_route_status(final_route_metadata, R)
    provenance, provenance_status =
        _sln_recursive_driver_route_provenance(route_provenance_metadata, catalog_id)
    staged_diagnostic = _sln_recursive_driver_staged_diagnostic(
        exact_field_status,
        determinant_status,
        variable_order_status,
        selected_status,
        ecp_status,
        final_status,
        provenance_status,
    )
    return (; matrix = A, base_ring = R, coefficient_ring = coefficient,
        dimension, ring_profile, exact_field_status, determinant = determinant_value,
        determinant_status, generators, generator_names, variable_order = normalized_order,
        variable_order_status, selected_variable = selected,
        selected_variable_index = selected_index, selected_variable_status = selected_status,
        last_column, last_column_profile, initial_unimodularity_witness = witness,
        initial_unimodularity_witness_status = witness_status,
        ecp_witness_metadata, ecp_evidence_status = ecp_status,
        final_route_metadata, final_route_evidence_status = final_status,
        route_provenance_metadata = provenance, route_provenance_status = provenance_status,
        catalog_id, support_classification = staged_diagnostic.status,
        staged_reason_code = staged_diagnostic.reason_code,
        staged_diagnostic)
end
```

Then add `_sln_recursive_driver_input_context_core_verification(context)` that
recomputes the fields from `context.matrix` and stored hints and returns a
`NamedTuple` of `*_ok` booleans plus `overall_core_ok`. Add
`_sln_recursive_driver_input_context_verification(context)` to compare stored
verification with the recomputed core, and
`_verify_sln_recursive_driver_input_context(context)::Bool`.

- [ ] **Step 4: Add public-in-module constructor**

Add:

```julia
function _sln_recursive_driver_input_context(A; kwargs...)
    fields = _sln_recursive_driver_input_context_fields(A; kwargs...)
    raw = SLnRecursiveDriverInputContext(values(merge(fields, (; verification = nothing,)))...)
    verification = _sln_recursive_driver_input_context_core_verification(raw)
    checked = SLnRecursiveDriverInputContext(values(merge(fields, (; verification,)))...)
    _verify_sln_recursive_driver_input_context(checked) ||
        error("internal SL_n recursive driver input context verification failed")
    return checked
end
```

- [ ] **Step 5: Run GREEN verification**

Run:

```bash
julia --project=. -e 'include("test/expert/park_woodburn_sln_driver_context.jl")'
```

Expected: exit 0.

- [ ] **Step 6: Commit implementation**

```bash
git add src/algorithm/polynomial_column_peel.jl
git commit -m "feat: add checked SLn recursive driver context"
```

### Task 3: Full Verification, Review, and PR Prep

**Files:**
- No planned file edits unless verification exposes an issue.

**Interfaces:**
- Consumes finished Tasks 1 and 2.
- Produces a verified branch ready for PR creation.

- [ ] **Step 1: Run targeted issue verification**

Run:

```bash
julia --project=. -e 'include("test/expert/park_woodburn_sln_driver_context.jl")'
```

Expected: exit 0.

- [ ] **Step 2: Run package verification**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: exit 0.

- [ ] **Step 3: Run diff hygiene**

Run:

```bash
git diff --check
```

Expected: exit 0.

- [ ] **Step 4: Review branch diff**

Run:

```bash
git diff --stat origin/main...HEAD
git diff --name-only origin/main...HEAD
```

Expected: changed files are limited to the Superpowers design/plan docs, the
new expert test, `test/runtests.jl`, and `src/algorithm/polynomial_column_peel.jl`.

- [ ] **Step 5: Commit any verification fixes**

If Step 1 or Step 2 required changes, run the targeted verification again and
commit only the changed implementation/test files:

```bash
git add src/algorithm/polynomial_column_peel.jl test/expert/park_woodburn_sln_driver_context.jl test/runtests.jl
git commit -m "test: harden issue 261 SLn driver context"
```

## Self-Review

Spec coverage: Task 1 covers supported, legacy, staged, unsupported, missing
variable metadata, determinant-not-one diagnostics, and tamper gates. Task 2
implements the internal context and verifier without public routing. Task 3
runs the requested targeted and package verification gates.

Placeholder scan: no task depends on unspecified files, unstated commands, or
open-ended implementation steps.

Type consistency: the planned context fields, constructor name, verifier names,
and test helper names are consistent across all tasks.
