# Issue 84 ECP Column Certificates Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add replayable expert certificates for the existing ECP column reducer families while keeping `reduce_unimodular_column(v, R)` factor-returning.

**Architecture:** Extend `src/algorithm/column_reduction.jl` with a thin certificate/result path used by the legacy reducer. Stage constructors reuse the existing unit, witness-unit, monicity-normalization, embedded-block, and Laurent-normalization logic; `verify_ecp_column_reduction` replays every recorded field and checks the final product sends the original column to `e_n`.

**Tech Stack:** Julia, Oscar polynomial/Laurent rings, Suslin elementary matrix helpers, Test stdlib.

## Global Constraints

- Preserve existing public behavior: `reduce_unimodular_column(v, R)` returns only the factor sequence.
- Do not broaden the accepted set of columns.
- Do not optimize factor count.
- Extend `src/algorithm/column_reduction.jl`; do not add a parallel reducer.
- Keep the certificate layer thin and replay-driven.
- Every recorded certificate or stage field must participate in exact verification.
- New expert names are intentionally not exported from `src/Suslin.jl`.
- Certificate entry point is `ecp_column_reduction_certificate(v, R)`.
- Verifier entry point is `verify_ecp_column_reduction(cert)::Bool`.
- Certificate type is `ECPColumnReductionCertificate`.
- Focused verification command is `julia --project=. -e 'include("test/expert/ecp_column_certificate.jl")'`.
- Required package verification command is `julia --project=. -e 'using Pkg; Pkg.test()'`.
- Do not commit a generated `Manifest.toml` unless the repository already tracked one before this issue.

---

## File Structure

- Modify `src/algorithm/column_reduction.jl`: add `ECPColumnReductionCertificate`, certificate constructors, stage builders, replay verifier, and legacy reducer delegation.
- Create `test/expert/ecp_column_certificate.jl`: focused expert tests for successful reducer families, legacy API preservation, and tamper rejection.
- Modify `test/runtests.jl`: register the focused expert test in the expert group.

---

### Task 1: Expert Certificate Acceptance Tests

**Files:**
- Create: `test/expert/ecp_column_certificate.jl`
- Modify: `test/runtests.jl`

**Interfaces:**
- Consumes: existing `test/fixtures/ecp_column_cases.jl`, `Suslin.reduce_unimodular_column`, `Suslin.is_unimodular_column`, and helper matrix APIs.
- Produces: RED coverage for `Suslin.ECPColumnReductionCertificate`, `Suslin.ecp_column_reduction_certificate(v, R)`, and `Suslin.verify_ecp_column_reduction(cert)::Bool`.

- [ ] **Step 1: Write the failing expert test**

Create `test/expert/ecp_column_certificate.jl` with this content:

```julia
using Test
using Suslin
using Oscar

include(joinpath(@__DIR__, "..", "fixtures", "ecp_column_cases.jl"))

function _ecp_cert_column(entry)
    return [getproperty(entry.entries, name) for name in entry.column_order]
end

function _ecp_cert_target_column(R, n::Int)
    target = zero_matrix(R, n, 1)
    target[n, 1] = one(R)
    return target
end

function _ecp_cert_factor_product(factors, R, n::Int)
    product = identity_matrix(R, n)
    for factor in factors
        product *= factor
    end
    return product
end

function _ecp_cert_apply_factors(factors, column, R)
    return _ecp_cert_factor_product(factors, R, length(column)) * matrix(R, length(column), 1, collect(column))
end

function _assert_ecp_certificate_replays(cert)
    R = cert.ring
    @test cert isa Suslin.ECPColumnReductionCertificate
    @test Suslin.verify_ecp_column_reduction(cert)
    @test _ecp_cert_apply_factors(cert.factors, cert.original_column, R) == _ecp_cert_target_column(R, length(cert.original_column))
    @test cert.final_column == _ecp_cert_target_column(R, length(cert.original_column))
    legacy_factors = Suslin.reduce_unimodular_column(cert.original_column, R)
    @test legacy_factors isa Vector
    @test _ecp_cert_apply_factors(legacy_factors, cert.original_column, R) == _ecp_cert_target_column(R, length(cert.original_column))
end

function _tamper_first_factor(cert)
    factors = copy(cert.factors)
    factors[1] = identity_matrix(cert.ring, length(cert.original_column))
    return Suslin.ECPColumnReductionCertificate(
        cert.original_column,
        cert.ring,
        cert.stages,
        factors,
        cert.final_column,
        cert.verification,
    )
end

function _tamper_witness(cert)
    stages = collect(cert.stages)
    witness_idx = findfirst(stage -> stage.kind == :witness_unit, stages)
    witness_idx === nothing && error("certificate has no witness stage")
    stage = stages[witness_idx]
    witness = collect(stage.witness)
    witness[1] += one(cert.ring)
    stages[witness_idx] = merge(stage, (; witness = tuple(witness...)))
    return Suslin.ECPColumnReductionCertificate(
        cert.original_column,
        cert.ring,
        tuple(stages...),
        cert.factors,
        cert.final_column,
        cert.verification,
    )
end

function _tamper_monicity_inverse(cert)
    stages = collect(cert.stages)
    monic_idx = findfirst(stage -> stage.kind == :monicity_normalization, stages)
    monic_idx === nothing && error("certificate has no monicity stage")
    stage = stages[monic_idx]
    inverse_values = collect(stage.inverse_values)
    inverse_values[stage.variable_index] = inverse_values[stage.variable_index] + one(cert.ring)
    stages[monic_idx] = merge(stage, (; inverse_values = tuple(inverse_values...)))
    return Suslin.ECPColumnReductionCertificate(
        cert.original_column,
        cert.ring,
        tuple(stages...),
        cert.factors,
        cert.final_column,
        cert.verification,
    )
end

function _tamper_embedded_indices(cert)
    stages = collect(cert.stages)
    embedded_idx = findfirst(stage -> stage.kind == :embedded_three_block, stages)
    embedded_idx === nothing && error("certificate has no embedded block stage")
    stage = stages[embedded_idx]
    stages[embedded_idx] = merge(stage, (; indices = reverse(stage.indices)))
    return Suslin.ECPColumnReductionCertificate(
        cert.original_column,
        cert.ring,
        tuple(stages...),
        cert.factors,
        cert.final_column,
        cert.verification,
    )
end

function _tamper_laurent_shift(cert)
    stages = collect(cert.stages)
    laurent_idx = findfirst(stage -> stage.kind == :laurent_normalization, stages)
    laurent_idx === nothing && error("certificate has no Laurent normalization stage")
    stage = stages[laurent_idx]
    stages[laurent_idx] = merge(stage, (; shift = one(cert.ring),))
    return Suslin.ECPColumnReductionCertificate(
        cert.original_column,
        cert.ring,
        tuple(stages...),
        cert.factors,
        cert.final_column,
        cert.verification,
    )
end

@testset "ECP column reduction certificates" begin
    cases = ECPColumnFixtureCatalog.cases_by_id()

    unit_cert = Suslin.ecp_column_reduction_certificate(_ecp_cert_column(cases["ecp-unit-entry-gf2"]), cases["ecp-unit-entry-gf2"].ring.object)
    @test any(stage -> stage.kind == :unit_entry, unit_cert.stages)
    _assert_ecp_certificate_replays(unit_cert)

    witness_cert = Suslin.ecp_column_reduction_certificate(_ecp_cert_column(cases["ecp-witness-unit-gf2"]), cases["ecp-witness-unit-gf2"].ring.object)
    @test any(stage -> stage.kind == :witness_unit, witness_cert.stages)
    _assert_ecp_certificate_replays(witness_cert)

    monic_cert = Suslin.ecp_column_reduction_certificate(_ecp_cert_column(cases["ecp-variable-change-monic-gf2"]), cases["ecp-variable-change-monic-gf2"].ring.object)
    @test any(stage -> stage.kind == :monicity_normalization, monic_cert.stages)
    _assert_ecp_certificate_replays(monic_cert)

    embedded_cert = Suslin.ecp_column_reduction_certificate(_ecp_cert_column(cases["ecp-longer-embedded-block-gf2"]), cases["ecp-longer-embedded-block-gf2"].ring.object)
    @test any(stage -> stage.kind == :embedded_three_block, embedded_cert.stages)
    _assert_ecp_certificate_replays(embedded_cert)

    R, (x, y) = Suslin.suslin_laurent_polynomial_ring(GF(2), ["x", "y"])
    laurent_column = [
        x^-1 + x^-2 * y^2,
        x^-1 * y + x^-1 + x^-2,
        one(R) + x^-1 * y + x^-2 * y + x^-2,
        x^-1 + x^-2 * y,
        x^-1 * y + x^-2 * y^2,
        x^-2 * y + x^-1 * y^2,
    ]
    laurent_cert = Suslin.ecp_column_reduction_certificate(laurent_column, R)
    @test any(stage -> stage.kind == :laurent_normalization, laurent_cert.stages)
    _assert_ecp_certificate_replays(laurent_cert)

    @test !Suslin.verify_ecp_column_reduction(_tamper_first_factor(unit_cert))
    @test !Suslin.verify_ecp_column_reduction(_tamper_witness(witness_cert))
    @test !Suslin.verify_ecp_column_reduction(_tamper_monicity_inverse(monic_cert))
    @test !Suslin.verify_ecp_column_reduction(_tamper_embedded_indices(embedded_cert))
    @test !Suslin.verify_ecp_column_reduction(_tamper_laurent_shift(laurent_cert))
    @test !Suslin.verify_ecp_column_reduction((; original_column = unit_cert.original_column))

    unsupported = _ecp_cert_column(cases["ecp-unsupported-unimodular-gf2"])
    @test_throws ArgumentError Suslin.ecp_column_reduction_certificate(unsupported, cases["ecp-unsupported-unimodular-gf2"].ring.object)

    non_unimodular = _ecp_cert_column(cases["ecp-non-unimodular-gf2"])
    @test_throws ArgumentError Suslin.ecp_column_reduction_certificate(non_unimodular, cases["ecp-non-unimodular-gf2"].ring.object)
end
```

- [ ] **Step 2: Register the expert file**

Modify `test/runtests.jl` and add:

```julia
"expert/ecp_column_certificate.jl",
```

after:

```julia
"expert/unimodular_reduction_exact.jl",
```

- [ ] **Step 3: Run the focused test to verify RED**

Run:

```bash
julia --project=. -e 'include("test/expert/ecp_column_certificate.jl")'
```

Expected: FAIL with `UndefVarError` or equivalent because `ECPColumnReductionCertificate`, `ecp_column_reduction_certificate`, and `verify_ecp_column_reduction` do not exist.

- [ ] **Step 4: Commit the RED test**

```bash
git add test/expert/ecp_column_certificate.jl test/runtests.jl docs/superpowers/plans/2026-06-22-issue-84-ecp-column-certificates.md
git commit -m "test: cover ECP column certificates"
```

---

### Task 2: Reducer Certificate Path And Replay Verifier

**Files:**
- Modify: `src/algorithm/column_reduction.jl`

**Interfaces:**
- Consumes: existing `_validated_unimodular_column`, `_reduce_unit_witness_column`, `_unimodular_witness`, `_reduce_via_witness_unit`, `_reduce_after_monicity_normalization`, `_embedded_three_block_reduction`, `_reduce_laurent_unimodular_column`, `_factor_sequence_product`, `_apply_reduction_factors`, `_target_reduced_column`, `_substitute_matrix_entries`, `_unit_normalization_factors`, `_lift_polynomial_reduction_factor`, `normalize_laurent_object`, `verify_laurent_normalization`, and `block_embedding`.
- Produces: non-exported `ECPColumnReductionCertificate`.
- Produces: non-exported `ecp_column_reduction_certificate(v, R)`.
- Produces: non-exported `verify_ecp_column_reduction(cert)::Bool`.
- Produces: stage builders returning `(; factors, stage)` for every existing successful reducer family.

- [ ] **Step 1: Add the certificate type and public-preserving entry points**

Near the top of `src/algorithm/column_reduction.jl`, add:

```julia
struct ECPColumnReductionCertificate
    original_column
    ring
    stages
    factors::Vector
    final_column
    verification
end

function reduce_unimodular_column(v::AbstractVector, R)
    return ecp_column_reduction_certificate(v, R).factors
end

function ecp_column_reduction_certificate(v::AbstractVector, R)
    column = _validated_unimodular_column(v, R)
    result = _is_laurent_polynomial_ring(R) ?
        _reduce_laurent_unimodular_column_certificate(column, R) :
        _reduce_polynomial_unimodular_column_exact_certificate(column, R)
    result !== nothing || _throw_unsupported_unimodular_column_reduction(column, R)

    factors = _checked_reduction_factors(result.factors, column, R, "certificate reducer")
    stages = ((; kind = :validation, input_length = length(column), is_unimodular = true), result.stage)
    final_column = _apply_reduction_factors(factors, column, R)
    provisional = ECPColumnReductionCertificate(column, R, stages, factors, final_column, nothing)
    verification = _ecp_column_reduction_replay_summary(provisional)
    verification.overall_ok || error("internal ECP column reduction certificate verification failed")
    certificate = ECPColumnReductionCertificate(column, R, stages, factors, final_column, verification)
    verify_ecp_column_reduction(certificate) || error("internal ECP column reduction certificate storage verification failed")
    return certificate
end
```

Delete the old first `reduce_unimodular_column` body instead of keeping two methods with identical signatures.

- [ ] **Step 2: Add stage result helpers**

Add helper functions that build stage records and return `(; factors, stage)`:

```julia
_unit_entry_reduction_certificate_stage(column, pivot_idx::Int, R)
_witness_unit_reduction_certificate_stage(column, witness::AbstractVector, pivot_idx::Int, R)
_reduce_supported_unimodular_column_certificate(column, R)
_reduce_exact_small_column_certificate(column, R)
_reduce_polynomial_unimodular_column_exact_certificate(column, R)
_reduce_after_monicity_normalization_certificate(column, R)
_reduce_via_supported_three_block_certificate(column, R)
_embedded_three_block_reduction_certificate_stage(column, R, indices, subresult)
_reduce_laurent_unimodular_column_certificate(column, R)
_ecp_certificate_from_stage(column, R, stage)
```

These helpers must use the existing factor formulas. For example, the unit and
witness stages must compute:

```julia
factors = _reduce_unit_witness_column(column, pivot_idx, R)
stage = (;
    kind = :unit_entry,
    input_column = tuple(column...),
    pivot_index = pivot_idx,
    pivot_value = column[pivot_idx],
    pivot_inverse = inv(column[pivot_idx]),
    factors,
    output_column = _apply_reduction_factors(factors, column, R),
)
```

and:

```julia
unit_creation_factors = _witness_unit_creation_factors(column, witness, pivot_idx, R)
created_column = vec(collect(_apply_reduction_factors(unit_creation_factors, column, R)))
unit_stage = _unit_entry_reduction_certificate_stage(created_column, pivot_idx, R).stage
factors = vcat(unit_stage.factors, unit_creation_factors)
stage = (;
    kind = :witness_unit,
    input_column = tuple(column...),
    witness = tuple(witness...),
    pivot_index = pivot_idx,
    witness_unit = witness[pivot_idx],
    witness_unit_inverse = inv(witness[pivot_idx]),
    unit_creation_factors,
    created_column = tuple(created_column...),
    unit_stage,
    factors,
    output_column = _apply_reduction_factors(factors, column, R),
)
```

Add `_witness_unit_creation_factors(column, witness, pivot_idx, R)` and refactor
`_reduce_via_witness_unit` to call it so factor formulas cannot drift.

- [ ] **Step 3: Refactor existing reducer helpers to delegate to certificate helpers**

Change the existing factor-returning helpers so they unwrap certificate results:

```julia
function _reduce_supported_unimodular_column(column::AbstractVector, R)
    result = _reduce_supported_unimodular_column_certificate(column, R)
    return result === nothing ? nothing : result.factors
end

function _reduce_polynomial_unimodular_column_exact(column::AbstractVector, R)
    result = _reduce_polynomial_unimodular_column_exact_certificate(column, R)
    return result === nothing ? nothing : result.factors
end

function _reduce_exact_small_column(column::AbstractVector, R)
    result = _reduce_exact_small_column_certificate(column, R)
    return result === nothing ? nothing : result.factors
end

function _reduce_via_supported_three_block(column::AbstractVector, R)
    result = _reduce_via_supported_three_block_certificate(column, R)
    return result === nothing ? nothing : result.factors
end

function _reduce_after_monicity_normalization(column::AbstractVector, R)
    result = _reduce_after_monicity_normalization_certificate(column, R)
    return result === nothing ? nothing : result.factors
end

function _reduce_laurent_unimodular_column(column::AbstractVector, R)
    result = _reduce_laurent_unimodular_column_certificate(column, R)
    return result === nothing ? nothing : result.factors
end
```

- [ ] **Step 4: Add replay verifier helpers**

Add these exact helper names:

```julia
verify_ecp_column_reduction(certificate)::Bool
_ecp_column_reduction_replay_summary(certificate)
_ecp_replay_stages(certificate)
_ecp_replay_stage(stage, input_column, R)
_ecp_factor_sequences_equal(left, right)
_ecp_stage_keys_ok(stage, expected)
_ecp_column_tuple(column)
_ecp_matrix_column_to_tuple(column_matrix)
```

`verify_ecp_column_reduction` must be:

```julia
function verify_ecp_column_reduction(certificate)::Bool
    try
        replay = _ecp_column_reduction_replay_summary(certificate)
        return replay.overall_ok && certificate.verification == replay
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end
```

`_ecp_replay_stage` must dispatch on `stage.kind` and recompute the expected
factors for `:unit_entry`, `:witness_unit`, `:monicity_normalization`,
`:embedded_three_block`, and `:laurent_normalization`. Each branch must check
exact stage keys, exact input column, exact witness or substitution equations,
exact stored factors, and exact stage output. It returns
`(; ok, factors, output_column)`.

- [ ] **Step 5: Run focused test to verify GREEN**

Run:

```bash
julia --project=. -e 'include("test/expert/ecp_column_certificate.jl")'
```

Expected: PASS.

- [ ] **Step 6: Run regression commands**

Run:

```bash
julia --project=. -e 'include("test/expert/unimodular_reduction_exact.jl")'
julia --project=. -e 'include("test/internal/ecp_column_fixtures.jl")'
```

Expected: both commands exit 0.

- [ ] **Step 7: Commit implementation**

```bash
git add src/algorithm/column_reduction.jl
git commit -m "feat: add ECP column reduction certificates"
```

---

## Self-Review

- Spec coverage: Task 1 covers all required successful families, legacy behavior, and negative controls. Task 2 implements the certificate object, expert constructor, replay verifier, reducer delegation, and exact stage replay.
- Placeholder scan: no TBD, TODO, fill-in, or placeholder steps remain.
- Type consistency: the plan consistently uses `ECPColumnReductionCertificate`, `ecp_column_reduction_certificate(v, R)`, and `verify_ecp_column_reduction(cert)::Bool`.
