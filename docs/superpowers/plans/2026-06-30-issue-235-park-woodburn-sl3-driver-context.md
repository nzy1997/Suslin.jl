# Issue 235 Park-Woodburn SL3 Driver Context Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an internal checked `SL_3` realization input context for ordinary polynomial matrices without changing public factorization dispatch.

**Architecture:** Add a non-exported context record, constructor, and verifier in `src/algorithm/factorization.jl`. The constructor reuses existing factorization precondition helpers, stores normalized metadata, classifies evidence availability, and keeps determinant-one multivariate inputs staged unless replayable evidence metadata is present. A focused expert test exercises #234 catalog cases and verifier corruption checks.

**Tech Stack:** Julia, Oscar exact polynomial rings and matrices, existing Suslin factorization helpers, #234 Park-Woodburn `SL_3` driver fixture catalog, `Test`.

## Global Constraints

- Do not route public `elementary_factorization` through the new context.
- Do not export the new context type or helper names from `src/Suslin.jl`.
- Only accept `3 x 3` exact field-backed ordinary polynomial determinant-one matrices.
- Recompute every stored verifier field: size, ring profile, exact-field status, determinant status, selected-variable membership in `gens(R)`, evidence availability, support/staged status, and staged diagnostics.
- A determinant-one multivariate `SL_3` input with no replayable local-form, variable-change, normality/conjugation, or Quillen/Murthy metadata must remain staged.
- Metadata with fixture ids but no replay target is recorded as partial evidence and must not make the context supported.
- Construction may throw `ArgumentError` for invalid input or contradictory catalog metadata.
- Focused verification command is `julia --project=. -e 'include("test/expert/park_woodburn_sl3_driver_context.jl")'`.
- Package verification command is `julia --project=. -e 'using Pkg; Pkg.test()'`.

---

## File Structure

- Modify `src/algorithm/factorization.jl`: add `SL3RealizationInputContext`, `_sl3_realization_input_context`, `_verify_sl3_realization_input_context`, and private helpers near the existing polynomial factorization route helper types.
- Create `test/expert/park_woodburn_sl3_driver_context.jl`: focused expert tests for positive, staged, and corruption cases.
- Modify `test/runtests.jl`: register the new expert test file near other Park-Woodburn route tests.

### Task 1: Add Red Expert Context Contract

**Files:**
- Create: `test/expert/park_woodburn_sl3_driver_context.jl`
- Modify: `test/runtests.jl`

**Interfaces:**
- Consumes future `Suslin._sl3_realization_input_context` and `Suslin._verify_sl3_realization_input_context`.
- Produces focused tests that later implementation must satisfy.

- [ ] **Step 1: Write the failing expert test**

Create `test/expert/park_woodburn_sl3_driver_context.jl` with:

```julia
using Test
using Oscar
using Suslin

const PARK_WOODBURN_SL3_DRIVER_CATALOG_PATH =
    joinpath(@__DIR__, "..", "fixtures", "park_woodburn_sl3_driver_cases.jl")

function _context_as_namedtuple(context)
    names = propertynames(context)
    return NamedTuple{names}(Tuple(getproperty(context, name) for name in names))
end

function _corrupt_context(context, updates)
    return Suslin.SL3RealizationInputContext(
        values(merge(_context_as_namedtuple(context), updates))...,
    )
end

function _catalog_metadata(entry)
    return (;
        fixture_id = entry.id,
        role = entry.role,
        expected_status = entry.expected_status,
        consumer_issue_ids = entry.consumer_issue_ids,
    )
end
```

The testset must include:

```julia
include(PARK_WOODBURN_SL3_DRIVER_CATALOG_PATH)
catalog = ParkWoodburnSL3DriverFixtureCatalog.catalog()
entries = ParkWoodburnSL3DriverFixtureCatalog.cases_by_id()
negative = Dict(entry.id => entry for entry in catalog.negative_controls)
```

Add assertions for these cases:

- `sl3-driver-univariate-fast-local-qq`: construct with selected variable,
  catalog metadata, and local-form witness; expect `support_status == :supported`,
  `evidence_status == :replayable`, `local_form_status == :replayed`,
  determinant status `:one`, exact-field status `:supported`, and verifier `true`.
- `sl3-driver-legacy-quillen-patched-substitution-qq`: construct with selected
  variable, catalog metadata, and `quillen_murthy_metadata = entry.upstream_evidence`;
  expect `support_status == :staged`, `evidence_status == :partial`,
  `quillen_murthy_status == :recorded`, and a staged diagnostic mentioning
  `:quillen_murthy`.
- `sl3-driver-multivariate-monic-special-form-qq`: construct with selected
  variable, catalog metadata, and local-form witness; expect supported
  replayed local-form evidence and verifier `true`.
- `sl3-driver-det-one-no-witness-staged-qq`: construct with selected variable
  and catalog metadata only; expect staged missing evidence and verifier `true`.
- negative controls: determinant not one, unsupported coefficient ring,
  selected variable not a generator, and supported-without-witness catalog
  metadata throw `ArgumentError`.
- corruption controls: selected variable, determinant status, ring profile,
  evidence status, and staged diagnostic changes make
  `_verify_sl3_realization_input_context` return `false`.

- [ ] **Step 2: Register the expert test**

Add `"expert/park_woodburn_sl3_driver_context.jl"` to the expert list in
`test/runtests.jl` immediately after `"expert/park_woodburn_route_certificate.jl"`.

- [ ] **Step 3: Run RED verification**

Run:

```bash
julia --project=. -e 'include("test/expert/park_woodburn_sl3_driver_context.jl")'
```

Expected: FAIL with `UndefVarError: SL3RealizationInputContext not defined` or
`UndefVarError: _sl3_realization_input_context not defined`.

- [ ] **Step 4: Commit the red test contract**

```bash
git add test/expert/park_woodburn_sl3_driver_context.jl test/runtests.jl
git commit -m "test: add sl3 realization context contract"
```

### Task 2: Implement Checked SL3 Realization Context

**Files:**
- Modify: `src/algorithm/factorization.jl`
- Test: `test/expert/park_woodburn_sl3_driver_context.jl`

**Interfaces:**
- Produces `SL3RealizationInputContext`, `_sl3_realization_input_context`, and `_verify_sl3_realization_input_context`.
- Keeps all helpers internal to `Suslin` by not adding exports.

- [ ] **Step 1: Add the context record**

Add this struct near `PolynomialFactorizationRouteCertificate`:

```julia
struct SL3RealizationInputContext
    matrix
    base_ring
    coefficient_ring
    size::Int
    ring_profile::Symbol
    generators::Tuple
    generator_names::Tuple
    selected_variable
    selected_variable_index
    selected_variable_status::Symbol
    determinant
    determinant_status::Symbol
    exact_field_status::Symbol
    catalog_metadata::NamedTuple
    local_form_witness
    local_form_status::Symbol
    variable_change_metadata
    variable_change_status::Symbol
    normality_conjugation_metadata
    normality_conjugation_status::Symbol
    quillen_murthy_metadata
    quillen_murthy_status::Symbol
    evidence_status::Symbol
    support_status::Symbol
    staged_diagnostic::NamedTuple
    verification
end
```

- [ ] **Step 2: Add constructor and verifier**

Implement:

```julia
function _sl3_realization_input_context(A; selected_variable=nothing,
        catalog_metadata=(;), local_form_witness=nothing,
        variable_change_metadata=nothing, normality_conjugation_metadata=nothing,
        quillen_murthy_metadata=nothing)
    fields = _sl3_realization_input_context_fields(A; selected_variable,
        catalog_metadata, local_form_witness, variable_change_metadata,
        normality_conjugation_metadata, quillen_murthy_metadata)
    context = SL3RealizationInputContext(values(merge(fields, (; verification = nothing,)))...)
    verification = _sl3_realization_input_context_core_verification(context)
    checked = SL3RealizationInputContext(values(merge(fields, (; verification,)))...)
    _verify_sl3_realization_input_context(checked) ||
        error("internal SL_3 realization input context verification failed")
    return checked
end

function _verify_sl3_realization_input_context(context)::Bool
    try
        return _sl3_realization_input_context_verification(context).overall_ok
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end
```

The field builder must call `_validate_factorization_matrix(A)`,
`_factorization_ring_profile(base_ring(A))`, and
`_require_polynomial_sl_determinant(A)`. It must throw `ArgumentError` unless
the matrix is exactly size 3, the ring profile is `:polynomial`, the
coefficient ring is an exact `Field`, and any selected variable is a generator.

- [ ] **Step 3: Add evidence helper behavior**

Implement helpers with these exact classifications:

- selected variable present and valid: `selected_variable_status = :passes`;
- selected variable absent: `selected_variable_status = :missing`;
- exact field-backed coefficient ring: `exact_field_status = :supported`;
- determinant one: `determinant_status = :one`;
- local-form witness that replays a monic matrix entry: `local_form_status = :replayed`;
- existing univariate fast-local shape without a witness: `local_form_status = :fast_local`;
- metadata with a replay id and matching target matrix: `:replayed`;
- metadata with a replay id but no target matrix: `:recorded`;
- missing metadata: `:missing`;
- any replayed or fast-local evidence: `evidence_status = :replayable` and `support_status = :supported`;
- recorded-only evidence: `evidence_status = :partial` and `support_status = :staged`;
- no evidence: `evidence_status = :missing` and `support_status = :staged`.

The staged diagnostic must be a named tuple containing `status`,
`missing_evidence`, `partial_evidence`, `selected_variable_status`,
`determinant_status`, `exact_field_status`, and `message`.

- [ ] **Step 4: Run GREEN focused verification**

Run:

```bash
julia --project=. -e 'include("test/expert/park_woodburn_sl3_driver_context.jl")'
```

Expected: PASS.

- [ ] **Step 5: Commit the implementation**

```bash
git add src/algorithm/factorization.jl
git commit -m "feat: add checked sl3 realization input context"
```

### Task 3: Verify Integration And Package Gate

**Files:**
- Verify: `src/algorithm/factorization.jl`
- Verify: `test/expert/park_woodburn_sl3_driver_context.jl`
- Verify: `test/runtests.jl`

**Interfaces:**
- Confirms the new context stays internal and package tests are still green.

- [ ] **Step 1: Run focused expert verification**

Run:

```bash
julia --project=. -e 'include("test/expert/park_woodburn_sl3_driver_context.jl")'
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

Expected: no untracked generated artifacts; diff limited to the spec, plan,
factorization context, expert test, and test registration.

## Plan Self-Review

- The plan covers every issue requirement and verification command.
- The task ordering follows TDD: failing expert test before production code.
- No public route selection or exported API is changed.
- The supported boundary requires replayable evidence; recorded ids alone stay
  staged.
