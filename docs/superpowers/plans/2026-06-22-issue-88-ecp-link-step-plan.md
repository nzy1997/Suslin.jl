# Issue 88 ECP Link Step Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a replayable Park-Woodburn ECP link-step certificate that maps the validated path from `v(0)` to `v(X)` and records inverse factors reducing `v(X)` to `v(0)`.

**Architecture:** Extend the existing internal ECP reducer module with a narrow link-step certificate layered on top of #87 link witnesses and #84 column-reduction replay certificates. The first supported family replays every link identity, records verified identity `SL_2` blocks, and transports endpoint path columns with replayed elementary factors derived from endpoint certificates.

**Tech Stack:** Julia, Oscar polynomial rings, existing `src/algorithm/column_reduction.jl` helpers, `elementary_matrix` / `block_embedding` matrix conventions, Julia `Test`.

## Global Constraints

- Repository has no `AGENTS.md`.
- Base branch is `main`; worker branch is `agent/issue-88-implement-the-park-woodburn-ecp-link-step-operat-run-1`.
- Extend `src/algorithm/column_reduction.jl`; do not add a parallel reducer.
- Add `test/expert/ecp_link_step.jl` and register it in `test/runtests.jl`.
- Preserve `reduce_unimodular_column(v, R)` public behavior.
- Keep new names non-exported expert/internal APIs accessed as `Suslin.<name>` in expert tests.
- Support ordinary polynomial rings and validated #87 `:supplied_link_witness` records only.
- Stage-fail unsupported path columns or unsupported non-identity `SL_2` families with `ArgumentError`.
- Stage-fail verified supplied witnesses whose residue-probe signature is not one of the explicitly supported #87 fixture families.
- Do not perform lower-variable induction or normality rewrite.
- Do not implement Quillen patching, factor-count optimization, random search, or public reducer routing.
- Verification command required by issue #88: `julia --project=. -e 'include("test/expert/ecp_link_step.jl")'`.
- Verification command required by Agent Desk: `julia --project=. -e 'using Pkg; Pkg.test()'`.

---

## File Structure

- Modify `src/algorithm/column_reduction.jl`: add `ECPLinkStepCertificate`, constructor, replay verifier, path evaluation helpers, segment transport helpers, inverse factor helpers, and segment replay helpers.
- Create `test/expert/ecp_link_step.jl`: fixture-backed positive and negative coverage for the supported link-step family.
- Modify `test/runtests.jl`: register `expert/ecp_link_step.jl` in the expert group.

### Task 1: Add Link-Step Expert Tests

**Files:**
- Create: `test/expert/ecp_link_step.jl`

**Interfaces:**
- Consumes: `test/fixtures/ecp_column_cases.jl`, witness helper shapes from `test/expert/ecp_link_witnesses.jl`, and wished-for `Suslin.ecp_link_step_certificate` / `Suslin.verify_ecp_link_step_certificate`.
- Produces: failing tests for path columns, segment maps, composed maps, lower-variable obligation, and negative tamper controls.

- [ ] **Step 1: Write the failing test**

Create `test/expert/ecp_link_step.jl` with fixture helper functions copied from the #87 expert witness test and add this assertion pattern:

```julia
@testset "ECP link step certificate replays path transport" begin
    gf2_entry = _case_by_id("ecp-variable-change-monic-gf2")
    gf2_column = _column(gf2_entry)
    gf2_record = Suslin.ecp_link_step_certificate(
        gf2_column,
        gf2_entry.ring.object;
        variable_order = gf2_entry.ring.generators,
        selected_variable = gf2_entry.ring.generators[1],
        supplied_link_witness = _gf2_link_witness(gf2_entry),
    )
    @test gf2_record.verification.overall_ok == true
    @test length(gf2_record.path_columns) == 2
    @test length(gf2_record.segments) == 1
    @test gf2_record.lower_variable_column == gf2_record.path_columns[1]
    @test gf2_record.transformed_column == tuple(gf2_column...)
    @test Suslin.verify_ecp_link_step_certificate(gf2_record)

    qq_entry = _case_by_id("ecp-monic-first-entry-qq")
    qq_column = _column(qq_entry)
    qq_record = Suslin.ecp_link_step_certificate(
        qq_column,
        qq_entry.ring.object;
        variable_order = qq_entry.ring.generators,
        selected_variable = qq_entry.ring.generators[1],
        supplied_link_witness = _qq_link_witness(qq_entry),
    )
    @test length(qq_record.path_columns) == 3
    @test length(qq_record.segments) == 2
    @test any(segment -> segment.delta != zero(qq_entry.ring.object), qq_record.segments)
    @test qq_record.verification.composed_forward_ok == true
    @test qq_record.verification.composed_reduction_ok == true
    @test Suslin.verify_ecp_link_step_certificate(qq_record)
end
```

- [ ] **Step 2: Verify RED**

Run:

```bash
julia --project=. -e 'include("test/expert/ecp_link_step.jl")'
```

Expected: FAIL with `UndefVarError: ecp_link_step_certificate not defined`.

- [ ] **Step 3: Add negative controls to the same test file**

Add helpers to reconstruct records by replacing a segment field and assert:

```julia
@test !Suslin.verify_ecp_link_step_certificate(_tamper_segment(qq_record, 1, :sl2_block, identity_matrix(R, 2) + elementary_matrix(2, 1, 2, one(R), R)))
@test !Suslin.verify_ecp_link_step_certificate(_tamper_segment(qq_record, 1, :elementary_factors, qq_record.segments[1].elementary_factors[1:(end - 1)]))
@test_throws ArgumentError Suslin.ecp_link_step_certificate(qq_column, R; variable_order = qq_entry.ring.generators, selected_variable = x, supplied_link_witness = _mutate_witness(qq_witness; path_points = (zero(R), x * y^2 + one(R), x)))
```

### Task 2: Implement Link-Step Certificate And Replay

**Files:**
- Modify: `src/algorithm/column_reduction.jl`

**Interfaces:**
- Consumes: `ECPLinkWitnessRecord`, `ecp_link_witness`, `verify_ecp_link_witness`, `_ecp_link_witness_replay_summary`, `_ecp_certificate_from_stage`, `ecp_column_reduction_certificate`, `_factor_sequence_product`, `_apply_reduction_factors`, `_ecp_factor_sequences_equal`, `_substitute_matrix_entries`, `block_embedding`, and `identity_matrix`.
- Produces: `ECPLinkStepCertificate`, `ecp_link_step_certificate`, `verify_ecp_link_step_certificate`, and internal replay helpers.

- [ ] **Step 1: Add the record type**

Add near `ECPLinkWitnessRecord`:

```julia
struct ECPLinkStepCertificate
    original_column
    ring
    link_witness::ECPLinkWitnessRecord
    path_points
    path_columns
    segments
    lower_variable_column
    transformed_column
    forward_factors::Vector
    reduction_factors::Vector
    verification
end
```

- [ ] **Step 2: Add constructor skeleton**

Implement:

```julia
function ecp_link_step_certificate(
    v::AbstractVector,
    R;
    link_witness = nothing,
    variable_order = tuple(gens(R)...),
    selected_variable = nothing,
    selected_monic_index::Integer = 1,
    supplied_link_witness = nothing,
)
    _is_laurent_polynomial_ring(R) &&
        throw(ArgumentError("ECP link steps currently support ordinary polynomial columns only"))
    witness = link_witness === nothing ?
        ecp_link_witness(v, R; variable_order, selected_variable, selected_monic_index, supplied_link_witness) :
        link_witness
    verify_ecp_link_witness(witness) ||
        throw(ArgumentError("ECP link step requires a verified Park-Woodburn link witness"))
    column = _validated_unimodular_column(v, R)
    tuple(column...) == witness.original_column ||
        throw(ArgumentError("ECP link step input column must match the link witness column"))
    # Build path columns, segments, composed factors, and stored verification.
end
```

- [ ] **Step 3: Add path and segment helpers**

Add helpers with these names and behaviors:

```julia
function _ecp_evaluate_at_selected_variable(values, selected_variable_index::Int, point, R)
    substitutions = collect(gens(R))
    substitutions[selected_variable_index] = point
    return tuple((_coerce_into_ring(R, evaluate(value, substitutions), "path column entry") for value in values)...)
end

function _ecp_link_step_segment(witness::ECPLinkWitnessRecord, idx::Int, path_columns)
    # Verify supported fixture family, delta, tail, Bezout, divisibility,
    # endpoint reductions, SL2 identity, and transport factors.
end
```

Segment records must include `index`, `from_path_point`, `to_path_point`,
`delta`, `from_column`, `to_column`, `sl2_block`, `sl2_embedding`,
`elementary_factors`, `forward_factors`, `inverse_factors`, `support_family`,
`endpoint_transport_matrix`, `from_certificate`, `to_certificate`,
`link_identity`, and `verification`.

- [ ] **Step 4: Add factor inversion and composition helpers**

Add:

```julia
function _ecp_inverse_elementary_factor(factor)
    R = base_ring(factor)
    n = nrows(factor)
    identity = identity_matrix(R, n)
    positions = [(row, col) for row in 1:n, col in 1:n if factor[row, col] != identity[row, col]]
    length(positions) == 1 || throw(ArgumentError("expected an elementary factor with one off-diagonal entry"))
    row, col = only(positions)
    row != col || throw(ArgumentError("expected an off-diagonal elementary factor"))
    return elementary_matrix(n, row, col, -factor[row, col], R)
end

function _ecp_inverse_factor_sequence(factors)
    return [_ecp_inverse_elementary_factor(factor) for factor in reverse(factors)]
end
```

Use these helpers to set `segment.forward_factors = vcat(_ecp_inverse_factor_sequence(to_cert.factors), from_cert.factors)` and `segment.inverse_factors = _ecp_inverse_factor_sequence(segment.forward_factors)`.
The constructor must first check the witness probe IDs and only allow the
known #87 fixture signatures `(:gf2_fixture_probe,)` and
`(:qq_y_probe, :qq_x_probe)`, recording
`:supplied_fixture_identity_sl2_endpoint_transport`; all other verified
supplied witnesses must throw `ArgumentError`.

- [ ] **Step 5: Add replay verifier**

Implement `verify_ecp_link_step_certificate(certificate)::Bool` and `_ecp_link_step_replay_summary(certificate)`. The replay summary must recompute path columns and every segment from the stored link witness, compare stored fields to recomputed fields, and verify:

```julia
segment_forward_ok = _apply_reduction_factors(segment.forward_factors, collect(segment.from_column), R) == matrix(R, n, 1, collect(segment.to_column))
segment_inverse_ok = _apply_reduction_factors(segment.inverse_factors, collect(segment.to_column), R) == matrix(R, n, 1, collect(segment.from_column))
composed_forward_ok = _apply_reduction_factors(certificate.forward_factors, collect(certificate.lower_variable_column), R) == matrix(R, n, 1, collect(certificate.transformed_column))
composed_reduction_ok = _apply_reduction_factors(certificate.reduction_factors, collect(certificate.transformed_column), R) == matrix(R, n, 1, collect(certificate.lower_variable_column))
```

- [ ] **Step 6: Verify GREEN**

Run:

```bash
julia --project=. -e 'include("test/expert/ecp_link_step.jl")'
```

Expected: PASS.

### Task 3: Register Expert Test And Final Checks

**Files:**
- Modify: `test/runtests.jl`

**Interfaces:**
- Consumes: `test/expert/ecp_link_step.jl`.
- Produces: expert test registration and verified branch.

- [ ] **Step 1: Register the expert test**

Add `"expert/ecp_link_step.jl"` to the expert group next to `"expert/ecp_link_witnesses.jl"`.

- [ ] **Step 2: Run issue verification**

Run:

```bash
julia --project=. -e 'include("test/expert/ecp_link_step.jl")'
```

Expected: PASS.

- [ ] **Step 3: Run expert harness**

Run:

```bash
julia --project=. test/runtests.jl expert
```

Expected: PASS.

- [ ] **Step 4: Run full package verification**

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

- Spec coverage: The plan covers path columns, per-segment `SL_2` and elementary contributions, replay data, lower-variable obligation, composed verification, narrow supported-family staging, and negative controls.
- Placeholder scan: No placeholder markers remain; skeleton steps name exact functions and expected behavior.
- Type consistency: The produced names are `ECPLinkStepCertificate`, `ecp_link_step_certificate`, and `verify_ecp_link_step_certificate`; tests use the same names through `Suslin.<name>`.
