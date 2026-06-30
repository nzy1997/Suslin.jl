# Issue 219 Murthy Adapter Consumption Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Consume verified #211 Murthy Quillen-adapter records as supplied local evidence and assemble a #218 evidence-backed global patch.

**Architecture:** Add an internal consumer certificate beside the existing Quillen evidence records in `src/algorithm/quillen_induction.jl`. The consumer verifies Murthy adapter replay and context alignment, converts ordinary adapters to #214 local sequence certificates, then delegates global assembly to `assemble_quillen_patch_from_local_evidence`.

**Tech Stack:** Julia, Oscar exact ordinary polynomial rings and matrices, existing Suslin Murthy local replay helpers, #214 local sequence certificates, #218 supplied-evidence patch assembly, `Test`.

## Global Constraints

- Input is an ordinary-polynomial `SL_3` local-form matrix `A`, selected variable `X`, and a nonempty collection of verified `MurthyQuillenLocalAdapter` records.
- Output is converted `QuillenLocalFactorSequenceCertificate` evidence and a `QuillenSuppliedEvidencePatchAssembly` whose `global_elementary_factors` multiply exactly to `A`.
- Verify original target matrix, ring, size, selected variable, Murthy certificate replay, local factor replay, factor order, denominator/local-witness metadata, products, and adapter replay metadata before assembly.
- Converted sequence certificates must verify and expose denominator provenance through `raw_denominators`, `product_denominator`, `verification.denominator_data`, and replay metadata.
- Localized denominator-cleared Murthy adapters remain staged: reject them before assembly with a clear `ArgumentError`.
- Delegate base-term handling to #218; missing or contradictory `A(0)` evidence must produce the existing staged #218 error.
- Keep the consumer expert/internal; do not export new names from `src/Suslin.jl`.
- Do not duplicate Murthy branch logic, accept raw factor vectors, broaden Murthy solving, implement the general `SL_3` route, ECP, recursive `SL_n`, Laurent support, or Steinberg optimization.
- Focused Murthy command is `julia --project=. -e 'include("test/expert/quillen_murthy_adapter_consumption.jl")'`.
- Existing route-adapter command is `julia --project=. -e 'include("test/expert/park_woodburn_quillen_route_adapter.jl")'`.
- Full package command is `julia --project=. -e 'using Pkg; Pkg.test()'`.

---

## File Structure

- Modify `src/algorithm/quillen_induction.jl`: add adapter-consumption records, replay/verify helpers, local sequence extraction, and supplied-evidence assembly wrapper.
- Create `test/expert/quillen_murthy_adapter_consumption.jl`: red/green positive assembly and negative controls.
- Modify `test/expert/park_woodburn_quillen_route_adapter.jl`: add a focused regression that an existing #211 ordinary adapter can feed the new consumer helper.
- Modify `test/runtests.jl`: register the new expert test after `expert/park_woodburn_quillen_route_adapter.jl`.

### Task 1: Add Red Adapter-Consumption Tests

**Files:**
- Create: `test/expert/quillen_murthy_adapter_consumption.jl`
- Modify: `test/expert/park_woodburn_quillen_route_adapter.jl`
- Modify: `test/runtests.jl`

**Interfaces:**
- Consumes future `Suslin.quillen_local_sequences_from_murthy_adapters`, `Suslin.assemble_quillen_patch_from_murthy_adapters`, `Suslin.replay_quillen_murthy_adapter_consumption`, and `Suslin.verify_quillen_murthy_adapter_consumption`.
- Produces failing coverage for adapter conversion, supplied-evidence assembly, provenance replay, and required negative controls.

- [ ] **Step 1: Write the new failing expert test**

Create `test/expert/quillen_murthy_adapter_consumption.jl` with:

```julia
using Test
using Suslin
using Oscar

const QMA_MURTHY_FIXTURES = joinpath(@__DIR__, "..", "fixtures", "sl3_murthy_gupta_cases.jl")
if !isdefined(Main, :SL3MurthyGuptaFixtureCatalog)
    include(QMA_MURTHY_FIXTURES)
end

function qma_product(factors, R, n::Int)
    product = identity_matrix(R, n)
    for factor in factors
        product *= factor
    end
    return product
end

function qma_rebuild(record; kwargs...)
    overrides = Dict{Symbol,Any}()
    for pair in kwargs
        overrides[pair.first] = pair.second
    end
    values = map(fieldnames(typeof(record))) do name
        get(overrides, name, getproperty(record, name))
    end
    return typeof(record)(values...)
end

function qma_fixture(id::AbstractString)
    catalog = SL3MurthyGuptaFixtureCatalog.catalog()
    return Dict(entry.id => entry for entry in catalog.cases)[id]
end

function qma_ordinary_adapter(; fixture_id = "mg-q0-unit-recursion")
    fixture = qma_fixture(fixture_id)
    certificate = Suslin.realize_sl3_local_certificate(
        fixture.entries.p,
        fixture.entries.q,
        fixture.entries.r,
        fixture.entries.s,
        fixture.variable,
    )
    adapter = Suslin._murthy_quillen_local_adapter(
        certificate,
        fixture.target,
        fixture.variable;
        witness_metadata = (;
            fixture_id = fixture.id,
            consumer_issue = 219,
        ),
    )
    return (; fixture, certificate, adapter)
end

function qma_localized_adapter(fixture_id)
    fixture = qma_fixture(fixture_id)
    context = Suslin.sl3_local_murthy_input_context(
        fixture.target,
        fixture.variable;
        witness = first(fixture.witnesses),
    )
    certificate = Suslin.realize_sl3_local_certificate(context)
    adapter = Suslin._murthy_quillen_local_adapter(
        certificate,
        fixture.target,
        fixture.variable;
        witness_metadata = (;
            fixture_id = fixture.id,
            consumer_issue = 219,
        ),
    )
    return (; fixture, certificate, adapter)
end

@testset "Murthy Quillen adapter consumption" begin
    ordinary = qma_ordinary_adapter()
    R = base_ring(ordinary.fixture.target)
    X = ordinary.fixture.variable

    sequences = Suslin.quillen_local_sequences_from_murthy_adapters(
        ordinary.fixture.target,
        X,
        [ordinary.adapter],
    )
    @test length(sequences) == 1
    @test all(Suslin.verify_quillen_local_factor_sequence_certificate, sequences)
    @test sequences[1] === ordinary.adapter.quillen_factor_sequence
    @test sequences[1].original_input == ordinary.fixture.target
    @test sequences[1].selected_variable == X
    @test sequences[1].raw_denominators == [one(R) for _ in ordinary.adapter.local_factor_replay.factors]
    @test sequences[1].product_denominator == one(R)
    @test sequences[1].local_product == ordinary.fixture.target
    @test sequences[1].local_correction == ordinary.fixture.target
    @test sequences[1].verification.denominator_data ==
          Suslin._quillen_denominator_data(sequences[1].normalized_local_contributions)
    @test all(
        provenance -> provenance.source == :murthy_local_sl3,
        [factor.provenance for factor in sequences[1].factors],
    )

    patch = Suslin.assemble_quillen_patch_from_murthy_adapters(
        ordinary.fixture.target,
        X,
        [ordinary.adapter];
        max_exponent = 1,
        base_term_policy = :already_handled,
        metadata = (; fixture_id = ordinary.fixture.id, consumer_issue = 219),
    )
    @test patch isa Suslin.QuillenSuppliedEvidencePatchAssembly
    @test Suslin.verify_quillen_patch(patch)
    @test patch.local_certificates == sequences
    @test patch.denominator_candidate.raw_denominators == [one(R)]
    @test patch.solver_result.coverage_sum == one(R)
    @test patch.product == ordinary.fixture.target
    @test qma_product(patch.global_elementary_factors, patch.ring, patch.size) ==
          ordinary.fixture.target
    @test patch.replay_metadata.metadata.source == :quillen_murthy_adapter_consumption

    consumption = Suslin.QuillenMurthyAdapterConsumption(
        ordinary.fixture.target,
        R,
        3,
        X,
        [ordinary.adapter],
        sequences,
        patch,
        (; source = :bad_metadata),
        nothing,
    )
    @test !Suslin.verify_quillen_murthy_adapter_consumption(consumption)

    result = Suslin.consume_murthy_quillen_adapters_for_patch(
        ordinary.fixture.target,
        X,
        [ordinary.adapter];
        max_exponent = 1,
        base_term_policy = :already_handled,
        metadata = (; fixture_id = ordinary.fixture.id, consumer_issue = 219),
    )
    @test result isa Suslin.QuillenMurthyAdapterConsumption
    @test Suslin.verify_quillen_murthy_adapter_consumption(result)
    @test result.local_sequence_certificates == sequences
    @test result.patch == patch
    @test result.replay_metadata.murthy_adapter_metadata[1] == ordinary.adapter.replay_metadata

    bad_factor_replay = qma_rebuild(
        ordinary.adapter.local_factor_replay;
        factors = reverse(ordinary.adapter.local_factor_replay.factors),
    )
    bad_factor_adapter = qma_rebuild(ordinary.adapter; local_factor_replay = bad_factor_replay)
    @test !Suslin._verify_murthy_quillen_local_adapter(bad_factor_adapter)
    @test_throws ArgumentError Suslin.quillen_local_sequences_from_murthy_adapters(
        ordinary.fixture.target,
        X,
        [bad_factor_adapter],
    )

    @test_throws ArgumentError Suslin.quillen_local_sequences_from_murthy_adapters(
        ordinary.fixture.target,
        one(R),
        [ordinary.adapter],
    )

    tampered_witness_factors = copy(ordinary.adapter.local_factor_replay.factors)
    tampered_record = tampered_witness_factors[1]
    tampered_witness_factors[1] = Suslin.SL3LocalElementaryFactor(
        tampered_record.R,
        tampered_record.n,
        tampered_record.row,
        tampered_record.col,
        tampered_record.numerator,
        tampered_record.denominator,
        tampered_record.selected_variable,
        (; tampered = true),
    )
    tampered_replay = qma_rebuild(ordinary.adapter.local_factor_replay; factors = tampered_witness_factors)
    tampered_adapter = qma_rebuild(ordinary.adapter; local_factor_replay = tampered_replay)
    @test !Suslin._verify_murthy_quillen_local_adapter(tampered_adapter)
    @test_throws ArgumentError Suslin.quillen_local_sequences_from_murthy_adapters(
        ordinary.fixture.target,
        X,
        [tampered_adapter],
    )

    wrong_target = ordinary.fixture.target * elementary_matrix(3, 1, 2, one(R), R)
    @test_throws ArgumentError Suslin.quillen_local_sequences_from_murthy_adapters(
        wrong_target,
        X,
        [ordinary.adapter],
    )

    localized = qma_localized_adapter("mg-local-q0-unit-at-u")
    @test localized.adapter.mode == :localized_replay_handoff
    @test_throws ArgumentError Suslin.quillen_local_sequences_from_murthy_adapters(
        localized.fixture.target,
        localized.fixture.variable,
        [localized.adapter],
    )

    @test_throws ArgumentError Suslin.assemble_quillen_patch_from_murthy_adapters(
        ordinary.fixture.target,
        X,
        [ordinary.adapter];
        max_exponent = 1,
    )
end
```

- [ ] **Step 2: Register the new expert test**

Add `"expert/quillen_murthy_adapter_consumption.jl"` to `test/runtests.jl` after the existing route adapter test.

- [ ] **Step 3: Add the existing route-adapter regression**

Append a focused assertion to `test/expert/park_woodburn_quillen_route_adapter.jl` inside the Murthy handoff test:

```julia
    consumed_sequence = only(Suslin.quillen_local_sequences_from_murthy_adapters(
        ordinary_fixture.target,
        ordinary_fixture.variable,
        [ordinary_adapter],
    ))
    @test consumed_sequence === ordinary_sequence
```

- [ ] **Step 4: Run the red tests**

Run:

```bash
julia --project=. -e 'include("test/expert/quillen_murthy_adapter_consumption.jl")'
julia --project=. -e 'include("test/expert/park_woodburn_quillen_route_adapter.jl")'
```

Expected: FAIL because the consumer entry points and records do not exist yet.

### Task 2: Implement Adapter Consumption

**Files:**
- Modify: `src/algorithm/quillen_induction.jl`

**Interfaces:**
- Produces `QuillenMurthyAdapterConsumption`, `QuillenMurthyAdapterConsumptionVerification`, `quillen_local_sequences_from_murthy_adapters`, `assemble_quillen_patch_from_murthy_adapters`, `consume_murthy_quillen_adapters_for_patch`, `replay_quillen_murthy_adapter_consumption`, and `verify_quillen_murthy_adapter_consumption`.
- Consumes existing `_verify_murthy_quillen_local_adapter`, `_murthy_quillen_local_factor_sequence_certificate`, `assemble_quillen_patch_from_local_evidence`, and #218 patch verification helpers.

- [ ] **Step 1: Add consumer records near `MurthyQuillenLocalAdapter`**

Add:

```julia
struct QuillenMurthyAdapterConsumptionVerification
    adapters_ok::Bool
    context_ok::Bool
    local_sequences_ok::Bool
    patch_ok::Bool
    patch_alignment_ok::Bool
    replay_metadata
    replay_metadata_ok::Bool
    overall_ok::Bool
end

struct QuillenMurthyAdapterConsumption
    original_input
    ring
    size::Int
    selected_variable
    murthy_adapters::Vector{MurthyQuillenLocalAdapter}
    local_sequence_certificates::Vector{QuillenLocalFactorSequenceCertificate}
    patch::QuillenSuppliedEvidencePatchAssembly
    replay_metadata
    verification::QuillenMurthyAdapterConsumptionVerification
end
```

- [ ] **Step 2: Add adapter validation and sequence extraction helpers**

Implement helpers that:

- collect a nonempty `Vector{MurthyQuillenLocalAdapter}`;
- verify `_verify_murthy_quillen_local_adapter`;
- compare ring, size, `original_input`, and selected variable against `A` and `X`;
- call `_murthy_quillen_local_factor_sequence_certificate(adapter)`;
- check `verify_quillen_local_factor_sequence_certificate(sequence)`;
- require `sequence.verification.denominator_data == _quillen_denominator_data(sequence.normalized_local_contributions)`;
- require every factor provenance has `source = :murthy_local_sl3`.

- [ ] **Step 3: Add the public-internal extraction function**

Implement:

```julia
function quillen_local_sequences_from_murthy_adapters(A, selected_variable, adapters)
    R = _require_quillen_denominator_cover_ring(base_ring(A))
    n = _require_square_matrix(A, "Murthy adapter consumption original input")
    n == 3 || throw(ArgumentError("Murthy adapter consumption requires a 3x3 input"))
    selected = _require_substitution_generator(R, selected_variable)
    collected = _quillen_murthy_adapter_vector(adapters)
    return _quillen_murthy_adapter_sequences(A, R, n, selected, collected)
end
```

- [ ] **Step 4: Add assembly and consumption wrappers**

Implement:

```julia
function assemble_quillen_patch_from_murthy_adapters(A, selected_variable, adapters; kwargs...)
    sequences = quillen_local_sequences_from_murthy_adapters(A, selected_variable, adapters)
    metadata = haskey(kwargs, :metadata) ? kwargs[:metadata] : (;)
    return assemble_quillen_patch_from_local_evidence(
        A,
        selected_variable,
        sequences;
        kwargs...,
    )
end
```

In actual Julia code, use explicit keyword arguments rather than splatting
`kwargs` from a named tuple so the method remains type-stable and compatible
with the existing #218 signature.

Implement `consume_murthy_quillen_adapters_for_patch` to return the full
`QuillenMurthyAdapterConsumption` certificate with replay metadata source
`:quillen_murthy_adapter_consumption`.

- [ ] **Step 5: Add replay and verify**

Replay recomputes local sequences from stored adapters, verifies the patch,
checks patch local certificates and original context, recomputes replay
metadata, and returns a `QuillenMurthyAdapterConsumptionVerification`.

- [ ] **Step 6: Run focused green tests**

Run:

```bash
julia --project=. -e 'include("test/expert/quillen_murthy_adapter_consumption.jl")'
julia --project=. -e 'include("test/expert/park_woodburn_quillen_route_adapter.jl")'
```

Expected: PASS.

### Task 3: Final Verification And Review

**Files:**
- Modify only files from Tasks 1 and 2.

**Interfaces:**
- Produces a branch ready for review when verification passes.

- [ ] **Step 1: Run required verification**

Run:

```bash
julia --project=. -e 'include("test/expert/quillen_murthy_adapter_consumption.jl")'
julia --project=. -e 'include("test/expert/park_woodburn_quillen_route_adapter.jl")'
julia --project=. -e 'using Pkg; Pkg.test()'
git diff --check
```

Expected: all commands exit 0.

- [ ] **Step 2: Review the diff**

Check:

```bash
git diff -- src/algorithm/quillen_induction.jl test/expert/quillen_murthy_adapter_consumption.jl test/expert/park_woodburn_quillen_route_adapter.jl test/runtests.jl docs/superpowers
```

Expected: no unrelated files or route-solver broadening.

- [ ] **Step 3: Commit and create PR**

Run:

```bash
git add docs/superpowers/specs/2026-06-30-issue-219-murthy-adapter-consumption-design.md \
        docs/superpowers/plans/2026-06-30-issue-219-murthy-adapter-consumption.md \
        src/algorithm/quillen_induction.jl \
        test/expert/quillen_murthy_adapter_consumption.jl \
        test/expert/park_woodburn_quillen_route_adapter.jl \
        test/runtests.jl
git commit -m "Implement #219: consume Murthy adapter evidence"
git push -u origin agent/issue-219-consume-murthy-quillen-adapter-outputs-in-the-ev-run-1
gh pr create --base main --head agent/issue-219-consume-murthy-quillen-adapter-outputs-in-the-ev-run-1 --title "Implement #219: consume Murthy adapter evidence" --body "## Summary

- consume verified Murthy Quillen-adapter records as #214 local sequence evidence
- delegate evidence-backed assembly to the #218 supplied-evidence patch path
- add adapter-provenance replay checks and negative controls

Closes #219"
```

Expected: PR URL is printed. If git writes or network are blocked by the
Agent Desk sandbox, report the exact failure in the final result.

## Plan Self-Review

- Every issue requirement maps to a validation, construction, or negative-test step.
- No Murthy solving or raw factor-vector path is introduced.
- The localized adapter path remains a staged error before assembly.
- The plan uses the recommended Subagent-Driven execution mode under the
  non-interactive Standing Answer Policy; if this sandbox blocks subagent git
  commits, execute the same tasks inline and record the deviation.
