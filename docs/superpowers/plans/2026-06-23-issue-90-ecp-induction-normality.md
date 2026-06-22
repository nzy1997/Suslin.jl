# Issue 90 ECP Induction Normality Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a replayable ECP induction/normality certificate that turns a verified #88 link step plus a verified lower-variable reduction into elementary factors reducing the original column to `e_n`.

**Architecture:** Extend `src/algorithm/column_reduction.jl` with a narrow internal certificate layered on the existing link-step and column-reduction replay APIs. Normality support is fixture-backed and explicit: supplied witness data records a non-identity embedded `SL_2` elementary contribution, the conjugator derived from the lower-variable reduction, and replayed elementary factors from `realize_conjugate_elementary`.

**Tech Stack:** Julia, Oscar polynomial rings, existing `ECPLinkStepCertificate`, existing `ECPColumnReductionCertificate`, `realize_conjugate_elementary` from `src/algorithm/normality.jl`, Julia `Test`.

## Global Constraints

- Repository has no `AGENTS.md`.
- Base branch is `main`; worker branch is `agent/issue-90-implement-ecp-induction-and-normality-replay-for-run-1`.
- Extend `src/algorithm/column_reduction.jl`; do not add a parallel reducer.
- Add `test/expert/ecp_induction_normality.jl` and register it in `test/runtests.jl`.
- Preserve `reduce_unimodular_column(v, R)` public behavior.
- Keep new names non-exported expert/internal APIs accessed as `Suslin.<name>` in expert tests.
- Support ordinary polynomial rings and verified #88 link-step certificates only.
- Keep normality witness data explicit with `source = :supplied_normality_witness`.
- Stage-fail missing normality witness data, identity `SL_2` entries, wrong conjugators, unsupported lower reductions, and failed final replay with `ArgumentError`.
- Do not implement Quillen patching, the full Park-Woodburn matrix driver, Laurent determinant correction, factor-count optimization, or public reducer routing.
- Verification command required by issue #90: `julia --project=. -e 'include("test/expert/ecp_induction_normality.jl")'`.
- Verification command required by Agent Desk: `julia --project=. -e 'using Pkg; Pkg.test()'`.

---

## File Structure

- Modify `src/algorithm/column_reduction.jl`: add `ECPInductionNormalityCertificate`, constructor, verifier, lower-reduction normalization helpers, normality-witness replay helpers, final-factor replay helpers, and elementary-factor checks.
- Create `test/expert/ecp_induction_normality.jl`: fixture-backed positive and negative coverage for the supported induction/normality family.
- Modify `test/runtests.jl`: register `expert/ecp_induction_normality.jl` in the expert group.

### Task 1: Add Induction/Normality Expert Tests

**Files:**
- Create: `test/expert/ecp_induction_normality.jl`

**Interfaces:**
- Consumes: `test/fixtures/ecp_column_cases.jl`, #88-style witness helpers, existing `Suslin.ecp_link_step_certificate`, existing `Suslin.ecp_column_reduction_certificate`, and wished-for `Suslin.ecp_induction_normality_certificate` / `Suslin.verify_ecp_induction_normality_certificate`.
- Produces: failing tests for two fixture-backed columns, explicit non-identity normality data, final factor replay, and negative tamper controls.

- [ ] **Step 1: Write the failing test file**

Create `test/expert/ecp_induction_normality.jl` with this structure:

```julia
using Test
using Oscar
using Suslin

const ECP_COLUMN_CATALOG_PATH = joinpath(@__DIR__, "..", "fixtures", "ecp_column_cases.jl")
include(ECP_COLUMN_CATALOG_PATH)

function _case_by_id(id::AbstractString)
    return ECPColumnFixtureCatalog.cases_by_id()[id]
end

function _column(entry)
    return [getproperty(entry.entries, name) for name in entry.column_order]
end

function _gf2_link_witness(entry)
    R = entry.ring.object
    x, y = entry.ring.generators
    v = _column(entry)
    G = y * v[2] + v[3]
    return (;
        source = :supplied_link_witness,
        residue_probes = ((; id = :gf2_fixture_probe, kind = :deterministic_fixture, maximal_ideal_generators = (y,)),),
        tail_reductions = ((; probe_id = :gf2_fixture_probe, G, lifted_tail_coefficients = (y, one(R)), tilde_G = G),),
        resultants = (one(R),),
        bezout_coefficients = ((; f = x, h = one(R)),),
        coverage_multipliers = (one(R),),
        path_points = (zero(R), x),
    )
end

function _qq_link_witness(entry)
    R = entry.ring.object
    x, y = entry.ring.generators
    return (;
        source = :supplied_link_witness,
        residue_probes = (
            (; id = :qq_y_probe, kind = :deterministic_fixture, maximal_ideal_generators = (y,)),
            (; id = :qq_x_probe, kind = :deterministic_fixture, maximal_ideal_generators = (x,)),
        ),
        tail_reductions = (
            (; probe_id = :qq_y_probe, G = y, lifted_tail_coefficients = (zero(R), one(R)), tilde_G = y),
            (; probe_id = :qq_x_probe, G = x, lifted_tail_coefficients = (one(R), zero(R)), tilde_G = x),
        ),
        resultants = (y^2, y + one(R)),
        bezout_coefficients = ((; f = zero(R), h = y), (; f = one(R), h = -x)),
        coverage_multipliers = (one(R), one(R) - y),
        path_points = (zero(R), y^2 * x, x),
    )
end

function _factor_product(factors, R, n::Int)
    product = identity_matrix(R, n)
    for factor in factors
        product *= factor
    end
    return product
end

function _apply_factors(factors, column, R)
    return _factor_product(factors, R, length(column)) * matrix(R, length(column), 1, collect(column))
end

function _normality_witness(lower_certificate, entry)
    R = lower_certificate.ring
    n = length(lower_certificate.original_column)
    lower_product = _factor_product(lower_certificate.factors, R, n)
    return (;
        source = :supplied_normality_witness,
        conjugator = inv(lower_product),
        sl2_indices = (n, 1),
        sl2_entry = entry,
    )
end

function _replace_record_field(record, field::Symbol, value)
    fields = fieldnames(typeof(record))
    values = [getfield(record, name) for name in fields]
    idx = findfirst(==(field), fields)
    idx === nothing && error("unknown record field: $(field)")
    values[idx] = value
    return typeof(record)(values...)
end

function _replace_rewrite_field(rewrite, field::Symbol, value)
    haskey(rewrite, field) || error("unknown rewrite field: $(field)")
    return merge(rewrite, NamedTuple{(field,)}((value,)))
end

function _fixture_certificate(id::AbstractString, witness_builder, normality_entry_builder)
    entry = _case_by_id(id)
    R = entry.ring.object
    column = _column(entry)
    link = Suslin.ecp_link_step_certificate(
        column,
        R;
        variable_order = entry.ring.generators,
        selected_variable = entry.ring.generators[1],
        supplied_link_witness = witness_builder(entry),
    )
    lower = Suslin.ecp_column_reduction_certificate(collect(link.lower_variable_column), R)
    witness = _normality_witness(lower, normality_entry_builder(R, entry.ring.generators))
    certificate = Suslin.ecp_induction_normality_certificate(
        column,
        R;
        link_step = link,
        lower_reduction = lower,
        normality_witness = witness,
    )
    return (; entry, R, column, link, lower, witness, certificate)
end

@testset "ECP induction and normality replay" begin
    gf2 = _fixture_certificate(
        "ecp-variable-change-monic-gf2",
        _gf2_link_witness,
        (R, generators) -> generators[2] + one(R),
    )
    @test Suslin.verify_ecp_link_step_certificate(gf2.link)
    @test Suslin.verify_ecp_column_reduction(gf2.lower)
    @test gf2.certificate.verification.overall_ok
    @test gf2.certificate.normality_rewrite.sl2_block != identity_matrix(gf2.R, 2)
    @test gf2.certificate.normality_rewrite.fixed_lower_column_ok
    @test _apply_factors(gf2.certificate.final_factors, gf2.column, gf2.R) ==
          Suslin._target_reduced_column(gf2.R, length(gf2.column))
    @test Suslin.verify_ecp_induction_normality_certificate(gf2.certificate)

    qq = _fixture_certificate(
        "ecp-monic-first-entry-qq",
        _qq_link_witness,
        (R, generators) -> generators[1] + generators[2],
    )
    @test Suslin.verify_ecp_link_step_certificate(qq.link)
    @test Suslin.verify_ecp_column_reduction(qq.lower)
    @test qq.certificate.verification.overall_ok
    @test qq.certificate.normality_rewrite.sl2_block != identity_matrix(qq.R, 2)
    @test qq.certificate.normality_rewrite.rewrite_product == _factor_product(
        qq.certificate.normality_rewrite.rewrite_factors,
        qq.R,
        length(qq.column),
    )
    @test _apply_factors(qq.certificate.final_factors, qq.column, qq.R) ==
          Suslin._target_reduced_column(qq.R, length(qq.column))
    @test Suslin.verify_ecp_induction_normality_certificate(qq.certificate)

    tampered_lifted = copy(qq.certificate.lifted_lower_variable_factors)
    tampered_lifted[1] = identity_matrix(qq.R, length(qq.column))
    tampered_lifted_certificate = _replace_record_field(
        qq.certificate,
        :lifted_lower_variable_factors,
        tampered_lifted,
    )
    @test !Suslin.verify_ecp_induction_normality_certificate(tampered_lifted_certificate)

    tampered_rewrite_factors = copy(qq.certificate.normality_rewrite.rewrite_factors)
    tampered_rewrite_factors[1] = identity_matrix(qq.R, length(qq.column))
    tampered_rewrite = _replace_rewrite_field(
        qq.certificate.normality_rewrite,
        :rewrite_factors,
        tampered_rewrite_factors,
    )
    tampered_rewrite_certificate = _replace_record_field(
        qq.certificate,
        :normality_rewrite,
        tampered_rewrite,
    )
    @test !Suslin.verify_ecp_induction_normality_certificate(tampered_rewrite_certificate)

    identity_witness = merge(qq.witness, (; sl2_entry = zero(qq.R)))
    @test_throws ArgumentError Suslin.ecp_induction_normality_certificate(
        qq.column,
        qq.R;
        link_step = qq.link,
        lower_reduction = qq.lower,
        normality_witness = identity_witness,
    )
end
```

- [ ] **Step 2: Run the test to verify RED**

Run:

```bash
julia --project=. -e 'include("test/expert/ecp_induction_normality.jl")'
```

Expected: FAIL with `UndefVarError: ecp_induction_normality_certificate not defined`.

### Task 2: Implement The Certificate And Replay

**Files:**
- Modify: `src/algorithm/column_reduction.jl`

**Interfaces:**
- Consumes: `ECPLinkStepCertificate`, `verify_ecp_link_step_certificate`, `ECPColumnReductionCertificate`, `verify_ecp_column_reduction`, `_factor_sequence_product`, `_apply_reduction_factors`, `_target_reduced_column`, `_lift_polynomial_reduction_factor`, `_ecp_inverse_factor_sequence`, `_ecp_factor_sequences_equal`, `block_embedding`, `elementary_matrix`, `realize_conjugate_elementary`.
- Produces: `ECPInductionNormalityCertificate`, `ecp_induction_normality_certificate`, `verify_ecp_induction_normality_certificate`, and internal replay helpers.

- [ ] **Step 1: Add the record type**

Add after `ECPLinkStepCertificate`:

```julia
struct ECPInductionNormalityCertificate
    original_column
    ring
    link_step::ECPLinkStepCertificate
    lower_variable_column
    lower_reduction_certificate
    lower_variable_factors::Vector
    lifted_lower_variable_factors::Vector
    normality_witness
    normality_rewrite
    final_factors::Vector
    final_column
    verification
end
```

- [ ] **Step 2: Add the constructor**

Add after `ecp_link_step_certificate`:

```julia
function ecp_induction_normality_certificate(
    v::AbstractVector,
    R;
    link_step = nothing,
    lower_reduction = nothing,
    normality_witness = nothing,
)
    _is_laurent_polynomial_ring(R) &&
        throw(ArgumentError("ECP induction/normality currently supports ordinary polynomial columns only"))
    link_step === nothing &&
        throw(ArgumentError("ECP induction/normality requires a verified link-step certificate"))
    verify_ecp_link_step_certificate(link_step) ||
        throw(ArgumentError("ECP induction/normality requires a verified link-step certificate"))
    _same_base_ring(link_step.ring, R) ||
        throw(ArgumentError("ECP induction/normality input ring must match the link-step ring"))

    column = _validated_unimodular_column(v, R)
    tuple(column...) == link_step.original_column ||
        throw(ArgumentError("ECP induction/normality input column must match the link-step column"))

    lower_column = collect(link_step.lower_variable_column)
    lower_certificate, lower_factors = _ecp_verified_lower_reduction(lower_reduction, lower_column, R)
    lifted_lower_factors = [_lift_polynomial_reduction_factor(factor, R) for factor in lower_factors]
    normality_rewrite = _ecp_induction_normality_rewrite(
        normality_witness,
        lower_column,
        lifted_lower_factors,
        R,
    )
    final_factors = vcat(lifted_lower_factors, normality_rewrite.rewrite_factors, link_step.reduction_factors)
    final_column = _apply_reduction_factors(final_factors, column, R)
    provisional = ECPInductionNormalityCertificate(
        tuple(column...),
        R,
        link_step,
        tuple(lower_column...),
        lower_certificate,
        lower_factors,
        lifted_lower_factors,
        normality_witness,
        normality_rewrite,
        final_factors,
        final_column,
        nothing,
    )
    verification = _ecp_induction_normality_replay_summary(provisional)
    verification.overall_ok ||
        throw(ArgumentError("constructed ECP induction/normality certificate failed exact replay verification"))
    certificate = ECPInductionNormalityCertificate(
        provisional.original_column,
        provisional.ring,
        provisional.link_step,
        provisional.lower_variable_column,
        provisional.lower_reduction_certificate,
        provisional.lower_variable_factors,
        provisional.lifted_lower_variable_factors,
        provisional.normality_witness,
        provisional.normality_rewrite,
        provisional.final_factors,
        provisional.final_column,
        verification,
    )
    verify_ecp_induction_normality_certificate(certificate) ||
        throw(ArgumentError("stored ECP induction/normality certificate failed exact replay verification"))
    return certificate
end
```

- [ ] **Step 3: Add lower-reduction helpers**

Add:

```julia
function _ecp_verified_lower_reduction(lower_reduction, lower_column, R)
    if lower_reduction === nothing
        certificate = ecp_column_reduction_certificate(lower_column, R)
        return certificate, certificate.factors
    end

    if lower_reduction isa ECPColumnReductionCertificate
        verify_ecp_column_reduction(lower_reduction) ||
            throw(ArgumentError("lower-variable reduction certificate does not verify"))
        _same_base_ring(lower_reduction.ring, R) ||
            throw(ArgumentError("lower-variable reduction ring must match the original ring"))
        lower_reduction.original_column == lower_column ||
            throw(ArgumentError("lower-variable reduction column must match v(0)"))
        return lower_reduction, lower_reduction.factors
    end

    factors = collect(lower_reduction)
    _apply_reduction_factors(factors, lower_column, R) == _target_reduced_column(R, length(lower_column)) ||
        throw(ArgumentError("lower-variable factor sequence does not reduce v(0) to e_n"))
    return nothing, factors
end
```

- [ ] **Step 4: Add normality replay helpers**

Add:

```julia
function _ecp_induction_normality_rewrite(normality_witness, lower_column, lifted_lower_factors, R)
    normality_witness === nothing &&
        throw(ArgumentError("ECP induction/normality requires explicit normality witness data"))
    _ecp_normality_witness_keys_ok(normality_witness) ||
        throw(ArgumentError("normality witness must contain source, conjugator, sl2_indices, and sl2_entry"))
    normality_witness.source == :supplied_normality_witness ||
        throw(ArgumentError("normality witness must use source = :supplied_normality_witness"))

    n = length(lower_column)
    lower_product = _factor_sequence_product(lifted_lower_factors, R, n)
    conjugator = normality_witness.conjugator
    nrows(conjugator) == n && ncols(conjugator) == n && _same_base_ring(base_ring(conjugator), R) ||
        throw(ArgumentError("normality witness conjugator must be an n by n matrix over R"))
    conjugator * lower_product == identity_matrix(R, n) ||
        throw(ArgumentError("normality witness conjugator must invert the lifted lower-variable reduction"))

    sl2_indices = tuple(normality_witness.sl2_indices...)
    length(sl2_indices) == 2 || throw(ArgumentError("normality witness sl2_indices must contain two indices"))
    fixed_index = Int(sl2_indices[1])
    moving_index = Int(sl2_indices[2])
    fixed_index == n || throw(ArgumentError("normality witness must use e_n as the fixed SL_2 coordinate"))
    1 <= moving_index <= n && moving_index != fixed_index ||
        throw(ArgumentError("normality witness moving index must be distinct and in range"))
    entry = _coerce_into_ring(R, normality_witness.sl2_entry, "normality witness sl2_entry")
    entry != zero(R) || throw(ArgumentError("normality witness must record a non-identity SL_2 contribution"))

    sl2_block = elementary_matrix(2, 1, 2, entry, R)
    sl2_embedding = block_embedding(sl2_block, n, sl2_indices)
    rewrite_factors = try
        realize_conjugate_elementary(conjugator, fixed_index, moving_index, entry)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("normality witness could not be rewritten into elementary factors"))
    end
    rewrite_product = _factor_sequence_product(rewrite_factors, R, n)
    expected_rewrite_product = conjugator * sl2_embedding * lower_product
    lower_matrix = matrix(R, n, 1, collect(lower_column))
    fixed_lower_column_ok = rewrite_product * lower_matrix == lower_matrix
    rewrite_product_ok = rewrite_product == expected_rewrite_product
    return (;
        source = :supplied_normality_witness,
        conjugator,
        lower_product,
        sl2_indices,
        fixed_index,
        moving_index,
        sl2_entry = entry,
        sl2_block,
        sl2_embedding,
        rewrite_factors,
        rewrite_product,
        expected_rewrite_product,
        rewrite_product_ok,
        fixed_lower_column_ok,
        overall_ok = rewrite_product_ok && fixed_lower_column_ok,
    )
end

function _ecp_normality_witness_keys_ok(normality_witness)
    return propertynames(normality_witness) == (:source, :conjugator, :sl2_indices, :sl2_entry)
end
```

- [ ] **Step 5: Add verifier and replay summary**

Add:

```julia
function verify_ecp_induction_normality_certificate(certificate)::Bool
    try
        replay = _ecp_induction_normality_replay_summary(certificate)
        return replay.overall_ok && certificate.verification == replay
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end

function _ecp_induction_normality_replay_summary(certificate)
    R = certificate.ring
    n = length(certificate.original_column)
    link_step_ok = verify_ecp_link_step_certificate(certificate.link_step)
    input_ok = link_step_ok && certificate.original_column == certificate.link_step.original_column
    lower_variable_column_ok = link_step_ok &&
        certificate.lower_variable_column == certificate.link_step.lower_variable_column

    lower_certificate, lower_factors = try
        _ecp_verified_lower_reduction(
            certificate.lower_reduction_certificate === nothing ?
                certificate.lower_variable_factors :
                certificate.lower_reduction_certificate,
            collect(certificate.lower_variable_column),
            R,
        )
    catch err
        err isa InterruptException && rethrow()
        nothing, Any[]
    end
    lower_reduction_ok = _ecp_factor_sequences_equal(certificate.lower_variable_factors, lower_factors)
    lifted_lower_factors = [_lift_polynomial_reduction_factor(factor, R) for factor in lower_factors]
    lifted_lower_factors_ok = _ecp_factor_sequences_equal(certificate.lifted_lower_variable_factors, lifted_lower_factors)

    normality_rewrite = try
        _ecp_induction_normality_rewrite(
            certificate.normality_witness,
            collect(certificate.lower_variable_column),
            certificate.lifted_lower_variable_factors,
            R,
        )
    catch err
        err isa InterruptException && rethrow()
        nothing
    end
    normality_rewrite_ok = normality_rewrite !== nothing &&
        certificate.normality_rewrite == normality_rewrite &&
        normality_rewrite.overall_ok

    expected_final_factors = normality_rewrite === nothing ?
        Any[] :
        vcat(certificate.lifted_lower_variable_factors, normality_rewrite.rewrite_factors, certificate.link_step.reduction_factors)
    final_factors_ok = _ecp_factor_sequences_equal(certificate.final_factors, expected_final_factors)
    final_column = final_factors_ok ?
        _apply_reduction_factors(certificate.final_factors, collect(certificate.original_column), R) :
        zero_matrix(R, n, 1)
    final_column_ok = certificate.final_column == final_column
    final_reduction_ok = final_column == _target_reduced_column(R, n)
    final_factors_elementary_ok = all(factor -> _ecp_is_elementary_factor(factor, R, n), certificate.final_factors)
    overall_ok = link_step_ok &&
        input_ok &&
        lower_variable_column_ok &&
        lower_reduction_ok &&
        lifted_lower_factors_ok &&
        normality_rewrite_ok &&
        final_factors_ok &&
        final_column_ok &&
        final_reduction_ok &&
        final_factors_elementary_ok
    return (;
        overall_ok,
        link_step_ok,
        input_ok,
        lower_variable_column_ok,
        lower_reduction_ok,
        lifted_lower_factors_ok,
        normality_rewrite_ok,
        final_factors_ok,
        final_column_ok,
        final_reduction_ok,
        final_factors_elementary_ok,
    )
end

function _ecp_is_elementary_factor(factor, R, n::Int)
    nrows(factor) == n || return false
    ncols(factor) == n || return false
    _same_base_ring(base_ring(factor), R) || return false
    identity = identity_matrix(R, n)
    positions = [(row, col) for row in 1:n, col in 1:n if factor[row, col] != identity[row, col]]
    length(positions) == 1 || return false
    row, col = only(positions)
    return row != col
end
```

- [ ] **Step 6: Run the issue test to verify GREEN**

Run:

```bash
julia --project=. -e 'include("test/expert/ecp_induction_normality.jl")'
```

Expected: PASS.

### Task 3: Register Tests And Run Final Verification

**Files:**
- Modify: `test/runtests.jl`

**Interfaces:**
- Consumes: `test/expert/ecp_induction_normality.jl`.
- Produces: expert test registration and verified branch.

- [ ] **Step 1: Register the expert test**

Add `expert/ecp_induction_normality.jl` after `expert/ecp_link_step.jl` in the expert group:

```julia
"expert/ecp_link_witnesses.jl",
"expert/ecp_link_step.jl",
"expert/ecp_induction_normality.jl",
```

- [ ] **Step 2: Run issue verification**

Run:

```bash
julia --project=. -e 'include("test/expert/ecp_induction_normality.jl")'
```

Expected: PASS.

- [ ] **Step 3: Run expert harness**

Run:

```bash
julia --project=. test/runtests.jl expert
```

Expected: PASS.

- [ ] **Step 4: Run package verification**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: PASS.

- [ ] **Step 5: Check whitespace**

Run:

```bash
git diff --check
```

Expected: no output and exit 0.

## Plan Self-Review

- Spec coverage: The plan covers link-step consumption, lower-variable reduction verification, lifting, explicit normality/conjugation replay, non-identity `SL_2` data, final factor replay, negative controls, and required verification commands.
- Placeholder scan: No placeholder markers remain; every task names exact files, functions, and commands.
- Type consistency: The produced names are `ECPInductionNormalityCertificate`, `ecp_induction_normality_certificate`, and `verify_ecp_induction_normality_certificate`; tests use the same names through `Suslin.<name>`.
