# Issue 262 ECP-Backed Polynomial Peel-Step Certificates Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Record and verify ECP-backed replay evidence for each ordinary-polynomial column-peel step.

**Architecture:** Extend `PolynomialColumnPeelStep` in `src/algorithm/polynomial_column_peel.jl` with explicit ECP evidence, route provenance, right-clearing coefficients, determinant/descent metadata, and per-step verification. Keep the existing factor-only constructor as compatibility, but make the certificate verifier require independently verified ECP evidence through `verify_ecp_column_reduction`.

**Tech Stack:** Julia, Oscar matrices, Suslin internal ECP verifier, Test stdlib.

## Global Constraints

- Repository has no `AGENTS.md`; follow README test commands.
- Preserve the existing factor-only `PolynomialColumnPeelStep` constructor as compatibility.
- Mainline `_polynomial_column_peel_step` must record an `ECPColumnReductionCertificate` from `ecp_column_reduction_certificate`.
- Verify ECP evidence through `verify_ecp_column_reduction`; do not infer ECP support from factor multiplication alone.
- Right-clearing factors must be deterministic from the bottom row of `B*A`; verifier recomputes them instead of trusting stored coefficients.
- The verifier must reject corrupted ECP evidence, ECP route provenance, one left factor, one right-clearing coefficient, the recorded last column, and the next block even if `product` is restored.
- Do not assemble multiple recursive steps, route public factorization, implement the #184 `SL_3` route, broaden Laurent/ToricBuilder support, or optimize factor counts.
- Run `julia --project=. -e 'include("test/expert/park_woodburn_sln_peel_step.jl")'`.
- Run `julia --project=. -e 'include("test/expert/park_woodburn_polynomial_column_peel.jl")'`.
- Run `julia --project=. -e 'using Pkg; Pkg.test()'`.

---

## File Structure

- Modify `src/algorithm/polynomial_column_peel.jl`: extend `PolynomialColumnPeelStep`, add step metadata/provenance helpers, add a step verifier, and wire certificate verification through it.
- Create `test/expert/park_woodburn_sln_peel_step.jl`: focused #260 multivariate step test and negative controls.
- Modify `test/expert/park_woodburn_polynomial_column_peel.jl`: update helper constructors and add whole-certificate negative controls for new fields.
- Modify `test/runtests.jl`: register the new expert test.

### Task 1: Red Peel-Step Contract

**Files:**
- Create: `test/expert/park_woodburn_sln_peel_step.jl`
- Modify: `test/runtests.jl`

**Interfaces:**
- Consumes future `PolynomialColumnPeelStep` fields: `ecp_evidence`, `ecp_route_provenance`, `right_clearing_coefficients`, `block_embedding_indices`, `determinant_metadata`, `descent_metadata`, `verification`.
- Consumes future `Suslin._polynomial_column_peel_step_verification(step)`.
- Produces expert coverage for a #260 multivariate `SL_4` peel step.

- [ ] **Step 1: Write the failing expert test**

Create `test/expert/park_woodburn_sln_peel_step.jl` with helper functions that rebuild a step and certificate by field name:

```julia
using Test
using Suslin
using Oscar

const PARK_WOODBURN_SLN_PEEL_STEP_CATALOG_PATH =
    joinpath(@__DIR__, "..", "fixtures", "park_woodburn_sln_driver_cases.jl")

function _sln_step_product(factors, R, n::Int)
    product = identity_matrix(R, n)
    for factor in factors
        product *= factor
    end
    return product
end

function _sln_step_target_column(R, n::Int)
    target = zero_matrix(R, n, 1)
    target[n, 1] = one(R)
    return target
end

function _sln_replace_step(step; kwargs...)
    values = Dict{Symbol, Any}(name => getfield(step, name) for name in fieldnames(typeof(step)))
    for (name, value) in kwargs
        values[name] = value
    end
    return Suslin.PolynomialColumnPeelStep((values[name] for name in fieldnames(typeof(step)))...)
end
```

The test body should load `ParkWoodburnSLnDriverFixtureCatalog`, build
`step = Suslin._polynomial_column_peel_step(entry.matrix)` for
`entry = entries["sln-driver-sl4-gf2-ecp-mainline"]`, and assert:

```julia
@test step.dimension == 4
@test step.ecp_evidence isa Suslin.ECPColumnReductionCertificate
@test step.left_certificate == step.ecp_evidence
@test Suslin.verify_ecp_column_reduction(step.ecp_evidence)
@test step.ecp_route_provenance.verifier == :verify_ecp_column_reduction
@test step.ecp_route_provenance.status == :verified
@test step.ecp_route_provenance.route in (:general_ecp_pipeline, :embedded_three_block, :witness_unit, :unit_entry, :monicity_normalization, :laurent_normalization, :laurent_elementary_row_preconditioning, :unknown)
@test _sln_step_product(step.left_factors, base_ring(step.input_matrix), step.dimension) *
      matrix(base_ring(step.input_matrix), step.dimension, 1, step.last_column) ==
      _sln_step_target_column(base_ring(step.input_matrix), step.dimension)
@test step.right_clearing_coefficients ==
      tuple((step.after_left_matrix[step.dimension, col] for col in 1:(step.dimension - 1))...)
@test step.peeled_matrix == block_embedding(step.next_block, step.dimension, collect(1:(step.dimension - 1)))
@test det(step.next_block) == one(base_ring(step.input_matrix))
@test step.determinant_metadata.next_block_determinant == one(base_ring(step.input_matrix))
@test step.descent_metadata.input_dimension == 4
@test step.descent_metadata.next_dimension == 3
@test Suslin._polynomial_column_peel_step_verification(step).overall_ok
```

Add negative controls that call `_sln_replace_step` to corrupt
`ecp_evidence`, `ecp_route_provenance`, `left_factors`,
`right_clearing_coefficients`, `last_column`, and `next_block`, and assert the
step verifier returns `false`.

- [ ] **Step 2: Register the expert test**

Add this file to the `expert` list in `test/runtests.jl` immediately before
`expert/park_woodburn_polynomial_column_peel.jl`:

```julia
"expert/park_woodburn_sln_peel_step.jl",
```

- [ ] **Step 3: Run RED verification**

Run:

```bash
julia --project=. -e 'include("test/expert/park_woodburn_sln_peel_step.jl")'
```

Expected: failure because the new step fields and verifier do not exist.

### Task 2: Step Metadata and Verifier

**Files:**
- Modify: `src/algorithm/polynomial_column_peel.jl`

**Interfaces:**
- Consumes: `ECPColumnReductionCertificate`, `verify_ecp_column_reduction`, `_factor_product`, `_expected_column_peel_right_factors`, `_factor_sequences_equal`.
- Produces: richer `PolynomialColumnPeelStep` records and `_polynomial_column_peel_step_verification(step)`.

- [ ] **Step 1: Extend the step struct and constructors**

Add fields after `left_certificate` and before `after_left_matrix`:

```julia
ecp_evidence
ecp_route_provenance::NamedTuple
```

Add fields after `right_factors` and before `peeled_matrix`:

```julia
right_clearing_coefficients::Tuple
```

Add fields before `end`:

```julia
block_embedding_indices::Vector{Int}
determinant_metadata::NamedTuple
descent_metadata::NamedTuple
verification
```

Update the legacy constructor to call a helper that derives deterministic
metadata from the supplied matrices and uses `nothing` for missing ECP evidence.

- [ ] **Step 2: Add ECP provenance and metadata helpers**

Add helpers near `_polynomial_column_peel_step`:

```julia
function _polynomial_column_peel_ecp_route_provenance(evidence)
    if evidence isa ECPColumnReductionCertificate && verify_ecp_column_reduction(evidence)
        terminal = evidence.stages[end]
        route = hasproperty(terminal, :route_metadata) && hasproperty(terminal.route_metadata, :route) ?
            terminal.route_metadata.route :
            (hasproperty(terminal, :kind) ? terminal.kind : :unknown)
        return (;
            source = :ecp_column_reduction_certificate,
            verifier = :verify_ecp_column_reduction,
            status = :verified,
            route,
            stage_kinds = tuple((hasproperty(stage, :kind) ? stage.kind : :unknown for stage in evidence.stages)...),
            factor_count = length(evidence.factors),
        )
    end
    return (;
        source = :missing_ecp_certificate,
        verifier = :verify_ecp_column_reduction,
        status = :missing,
        route = :unknown,
        stage_kinds = (),
        factor_count = 0,
    )
end

function _polynomial_column_peel_right_clearing_coefficients(after_left, d::Int)
    return tuple((after_left[d, col] for col in 1:(d - 1))...)
end
```

Add determinant/descent helpers that recompute from `input_matrix`,
`peeled_matrix`, and `next_block`.

- [ ] **Step 3: Build rich records in `_polynomial_column_peel_step`**

After creating `left_certificate`, set:

```julia
ecp_evidence = left_certificate
ecp_route_provenance = _polynomial_column_peel_ecp_route_provenance(ecp_evidence)
right_clearing_coefficients = _polynomial_column_peel_right_clearing_coefficients(after_left, d)
block_embedding_indices = collect(1:(d - 1))
determinant_metadata = _polynomial_column_peel_determinant_metadata(current, peeled, next_block)
descent_metadata = _polynomial_column_peel_descent_metadata(d)
```

Construct a provisional step with `verification = nothing`, compute
`verification = _polynomial_column_peel_step_core_verification(provisional)`,
require `verification.overall_core_ok`, and return a checked step with the
stored verification.

- [ ] **Step 4: Verify all step fields by recomputation**

Replace `_is_valid_polynomial_column_peel_step_data` internals with a wrapper
around a new `_polynomial_column_peel_step_core_verification(step)` for full
records, while preserving the existing data-only signature for tests.

The new core verifier must return a `NamedTuple` containing:

```julia
overall_core_ok
shape_ok
last_column_ok
ecp_evidence_ok
ecp_route_provenance_ok
left_factors_ok
after_left_ok
right_clearing_coefficients_ok
right_factors_ok
peeled_matrix_ok
block_embedding_ok
next_block_ok
determinant_metadata_ok
descent_metadata_ok
```

`_polynomial_column_peel_step_verification(step)` must merge
`stored_verification_ok = step.verification == core` and
`overall_ok = core.overall_core_ok && stored_verification_ok`.

- [ ] **Step 5: Run GREEN verification for the new step test**

Run:

```bash
julia --project=. -e 'include("test/expert/park_woodburn_sln_peel_step.jl")'
```

Expected: exit 0.

### Task 3: Whole-Certificate Hardening

**Files:**
- Modify: `src/algorithm/polynomial_column_peel.jl`
- Modify: `test/expert/park_woodburn_polynomial_column_peel.jl`

**Interfaces:**
- Consumes: `_polynomial_column_peel_step_verification(step)`.
- Produces: whole-certificate verification that rejects corrupted rich step metadata.

- [ ] **Step 1: Update whole-certificate verification**

In `_polynomial_column_peel_core_verification`, compute `steps_ok` by calling
`_polynomial_column_peel_step_verification(step).overall_ok` for each step.
Keep `left_certificates_ok` as a named compatibility component, but implement it
through the same ECP evidence verifier.

- [ ] **Step 2: Update test helper constructors**

Update `_pw_poly_replace_step` in
`test/expert/park_woodburn_polynomial_column_peel.jl` to accept and pass the
new fields by default:

```julia
ecp_evidence = step.ecp_evidence
ecp_route_provenance = step.ecp_route_provenance
right_clearing_coefficients = step.right_clearing_coefficients
block_embedding_indices = step.block_embedding_indices
determinant_metadata = step.determinant_metadata
descent_metadata = step.descent_metadata
verification = step.verification
```

- [ ] **Step 3: Add negative controls**

Add helpers to corrupt `ecp_evidence`, `ecp_route_provenance`, and
`right_clearing_coefficients`. In the existing testset, assert:

```julia
@test !Suslin._verify_polynomial_column_peel_certificate(_pw_poly_corrupt_ecp_evidence(recursive_cert))
@test !Suslin._verify_polynomial_column_peel_certificate(_pw_poly_corrupt_ecp_route_provenance(recursive_cert))
@test !Suslin._verify_polynomial_column_peel_certificate(_pw_poly_corrupt_right_clearing_coefficient(recursive_cert))
```

Each helper should return a certificate with `product = cert.product` restored,
matching the issue's negative-control requirement.

- [ ] **Step 4: Run polynomial peel verification**

Run:

```bash
julia --project=. -e 'include("test/expert/park_woodburn_polynomial_column_peel.jl")'
```

Expected: exit 0.

### Task 4: Full Verification and Commit

**Files:**
- Verify all modified files.

**Interfaces:**
- Produces a reviewable commit and PR-ready branch.

- [ ] **Step 1: Run requested expert commands**

Run:

```bash
julia --project=. -e 'include("test/expert/park_woodburn_sln_peel_step.jl")'
julia --project=. -e 'include("test/expert/park_woodburn_polynomial_column_peel.jl")'
```

Expected: both exit 0.

- [ ] **Step 2: Run package tests**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: exit 0.

- [ ] **Step 3: Review the diff**

Run:

```bash
git diff --check
git status --short
git diff -- src/algorithm/polynomial_column_peel.jl test/expert/park_woodburn_sln_peel_step.jl test/expert/park_woodburn_polynomial_column_peel.jl test/runtests.jl
```

Expected: no whitespace errors and only issue #262 files changed.

- [ ] **Step 4: Commit**

Run:

```bash
git add docs/superpowers/specs/2026-07-03-issue-262-ecp-backed-polynomial-peel-design.md docs/superpowers/plans/2026-07-03-issue-262-ecp-backed-polynomial-peel.md src/algorithm/polynomial_column_peel.jl test/expert/park_woodburn_sln_peel_step.jl test/expert/park_woodburn_polynomial_column_peel.jl test/runtests.jl
git commit -m "Implement #262: record ECP-backed polynomial peel steps"
```

Expected: commit created on the worker branch.

## Self-Review

This plan covers all issue requirements: ECP evidence storage, verifier-backed
ECP replay, route provenance, left/right factor replay, deterministic clearing
coefficients, block embedding, determinant/descent metadata, new #260
multivariate step coverage, negative controls, and the requested verification
commands. It preserves compatibility by retaining the legacy constructor while
making mainline certificate verification require verified ECP evidence.
