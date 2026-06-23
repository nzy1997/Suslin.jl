# Issue 105 Constructive Quillen Induction Acceptance Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add the final #63 expert acceptance test proving the constructive Quillen local-to-global path works end to end through the existing #99 through #104 expert data and verifiers.

**Architecture:** Keep #105 as a focused expert acceptance harness. The new test builds cover certificates, local realization certificates, normalized local contributions, and a deterministic global patch from fixture catalog entries, then checks exact replay at every stage and registers the file in the expert test group.

**Tech Stack:** Julia, Oscar exact ordinary polynomial rings, existing Suslin Quillen fixture catalog and deterministic patch assembly helpers, Test stdlib.

## Global Constraints

- No `AGENTS.md`, `CLAUDE.md`, `GEMINI.md`, or `CODEX.md` file is present in this checkout.
- The worker branch is `agent/issue-105-add-the-final-constructive-quillen-induction-acc-run-1`.
- Dependencies #103 and #104 are merged.
- Reuse #99 fixture entries from `test/fixtures/quillen_patch_cases.jl`.
- Reuse #103 and #104 verifiers rather than duplicating unchecked construction logic.
- Keep the deterministic Quillen assembly surface expert/internal and do not export new names from `src/Suslin.jl`.
- Do not change `test/public/api_surface.jl`.
- Do not implement the final public `elementary_factorization` driver from #64.
- Do not solve arbitrary local realizability.
- Do not handle the ToricBuilder Laurent `GL_n` boundary from #38.
- Do not optimize or minimize the number of elementary factors.
- Focused command is `julia --project=. -e 'include("test/expert/quillen_induction_constructive.jl")'`.
- Full suite command is `julia --project=. test/runtests.jl all`.
- Full package command required by Agent Desk is `julia --project=. -e 'using Pkg; Pkg.test()'`.

---

## File Structure

- Create `test/expert/quillen_induction_constructive.jl`: final #63 expert acceptance test over two deterministic fixture entries, including replay assertions and negative controls.
- Modify `test/runtests.jl`: register `expert/quillen_induction_constructive.jl` immediately after `expert/quillen_patch_verification_hardening.jl`.
- Leave `src/Suslin.jl` and `test/public/api_surface.jl` unchanged.
- Modify `src/algorithm/quillen_induction.jl` only if the acceptance test exposes a real gap in the existing expert path.

---

### Task 1: Constructive Quillen Expert Acceptance Test

**Files:**
- Create: `test/expert/quillen_induction_constructive.jl`
- Modify: `test/runtests.jl`

**Interfaces:**
- Consumes: `test/fixtures/quillen_patch_cases.jl`, `Suslin.quillen_denominator_cover_certificate`, `Suslin.quillen_local_realization_certificate`, `Suslin.normalize_quillen_local_contributions`, `Suslin.assemble_deterministic_quillen_patch`, `Suslin.replay_deterministic_quillen_patch`, and `Suslin.verify_quillen_patch`.
- Produces: A focused expert acceptance command and expert runner registration for the full constructive Quillen path.

- [ ] **Step 1: Confirm the RED focused command**

Run:

```bash
julia --project=. -e 'include("test/expert/quillen_induction_constructive.jl")'
```

Expected before implementation: exit 1 with `SystemError: opening file ... test/expert/quillen_induction_constructive.jl: No such file or directory`.

- [ ] **Step 2: Add the acceptance test file**

Create `test/expert/quillen_induction_constructive.jl` with:

```julia
using Test
using Suslin
using Oscar

const QUILLEN_CONSTRUCTIVE_CATALOG_PATH =
    joinpath(@__DIR__, "..", "fixtures", "quillen_patch_cases.jl")
if !isdefined(Main, :QuillenPatchFixtureCatalog)
    include(QUILLEN_CONSTRUCTIVE_CATALOG_PATH)
end

function constructive_product(factors, R, n)
    product = identity_matrix(R, n)
    for factor in factors
        product *= factor
    end
    return product
end

function constructive_rebuild(record; kwargs...)
    overrides = Dict{Symbol,Any}()
    for pair in kwargs
        overrides[pair.first] = pair.second
    end
    values = map(fieldnames(typeof(record))) do name
        get(overrides, name, getproperty(record, name))
    end
    return typeof(record)(values...)
end

function constructive_local_certificate(local_factor)
    return Suslin.LocalCertificate(
        local_factor.certificate.indices,
        local_factor.certificate.denominators,
    )
end

function constructive_correction(local_factor)
    correction = local_factor.correction
    return Suslin.QuillenElementaryCorrection(
        correction.row,
        correction.col,
        correction.entry,
    )
end

function constructive_local_certificate_from_fixture(entry; local_index::Int)
    local_factor = entry.local_factors[local_index]
    return Suslin.quillen_local_realization_certificate(
        entry.target_matrix,
        entry.substitution_variable;
        local_certificate = constructive_local_certificate(local_factor),
        denominator = local_factor.denominator,
        coverage_multiplier = local_factor.coverage_multiplier,
        correction = constructive_correction(local_factor),
        factors = [local_factor.factor],
        patched_substitution_witness = entry.patched_substitution_witness,
        witness_metadata = (;
            fixture_id = entry.id,
            local_index = local_index,
            source_refs = entry.source_refs,
            consumer_issue_ids = entry.consumer_issue_ids,
        ),
    )
end

function constructive_cover(entry)
    denominators = [data.denominator for data in entry.denominator_data]
    multipliers = [data.coverage_multiplier for data in entry.denominator_data]
    return Suslin.quillen_denominator_cover_certificate(
        entry.ring.object,
        denominators,
        multipliers,
    )
end

function constructive_inputs(entry)
    cover = constructive_cover(entry)
    local_certificates = [
        constructive_local_certificate_from_fixture(entry; local_index = idx)
        for idx in eachindex(entry.local_factors)
    ]
    normalized = Suslin.normalize_quillen_local_contributions(
        local_certificates,
        cover;
        original_input = entry.target_matrix,
        selected_variable = entry.substitution_variable,
    )
    return cover, local_certificates, normalized
end

function constructive_patch(entry)
    cover, local_certificates, normalized = constructive_inputs(entry)
    patch = Suslin.assemble_deterministic_quillen_patch(
        entry.target_matrix,
        entry.substitution_variable,
        local_certificates,
        normalized,
        cover;
        target = entry.expected.global_correction,
    )
    return cover, local_certificates, normalized, patch
end

function constructive_rebuild_global_patch(patch; kwargs...)
    fields = merge((
        ring = patch.ring,
        size = patch.size,
        substitution_variable = patch.substitution_variable,
        original_input = patch.original_input,
        cover_certificate = patch.cover_certificate,
        denominator_data = patch.denominator_data,
        local_certificates = patch.local_certificates,
        normalized_local_contributions = patch.normalized_local_contributions,
        global_elementary_factors = patch.global_elementary_factors,
        patched_product = patch.patched_product,
        target = patch.target,
        replay_metadata = patch.replay_metadata,
        verification = patch.verification,
    ), kwargs)
    return Suslin.QuillenGlobalPatchAssembly(
        fields.ring,
        fields.size,
        fields.substitution_variable,
        fields.original_input,
        fields.cover_certificate,
        fields.denominator_data,
        fields.local_certificates,
        fields.normalized_local_contributions,
        fields.global_elementary_factors,
        fields.patched_product,
        fields.target,
        fields.replay_metadata,
        fields.verification,
    )
end

function assert_constructive_patch_replays(entry, cover, local_certificates, normalized, patch)
    R = entry.ring.object
    @test patch isa Suslin.QuillenGlobalPatchAssembly
    @test Suslin.verify_quillen_denominator_cover(cover)
    @test cover.coverage_sum == one(R)
    @test sum(
        data.coverage_multiplier * data.denominator
        for data in entry.denominator_data;
        init = zero(R),
    ) == one(R)
    @test [data.denominator for data in entry.denominator_data] == cover.denominators
    @test [data.coverage_multiplier for data in entry.denominator_data] ==
          cover.coverage_multipliers

    @test all(Suslin.verify_quillen_local_certificate, local_certificates)
    for (idx, certificate) in enumerate(local_certificates)
        local_replay = certificate.verification
        @test local_replay.overall_ok
        @test local_replay.denominator_ok
        @test local_replay.local_correction_ok
        @test local_replay.factors_ok
        @test local_replay.stored_product_ok
        @test local_replay.patched_substitution_ok
        @test certificate.local_product == entry.local_factors[idx].factor
        @test certificate.local_correction == entry.local_factors[idx].expected_correction
    end

    @test all(Suslin.verify_quillen_local_contribution_normalization, normalized)
    for item in normalized
        replay = Suslin.replay_quillen_local_contribution_normalization(item)
        @test replay.overall_ok
        @test replay.cover_certificate_ok
        @test replay.local_certificate_ok
        @test replay.cover_pair_ok
        @test replay.patched_substitution_ok
        @test replay.local_product_ok
        @test replay.local_correction_ok
        @test replay.local_contribution_ok
        @test replay.weighted_global_elementary_factors_ok
        @test replay.replay_metadata_ok
    end

    @test Suslin.verify_quillen_patch(patch)
    replay = Suslin.replay_deterministic_quillen_patch(patch)
    @test replay.overall_ok
    @test replay.cover_certificate_ok
    @test replay.local_certificates_ok
    @test replay.normalized_contributions_ok
    @test replay.local_count_ok
    @test replay.local_alignment_ok
    @test replay.cover_alignment_ok
    @test replay.normalized_input_ok
    @test replay.selected_variable_ok
    @test replay.denominator_data_ok
    @test replay.coverage_ok
    @test replay.coverage_sum == one(R)
    @test replay.global_elementary_factors_ok
    @test replay.product_ok
    @test replay.target_ok
    @test replay.replay_metadata_ok

    expected_factors = reduce(
        vcat,
        [item.weighted_global_elementary_factors for item in normalized];
        init = typeof(identity_matrix(R, entry.size))[],
    )
    @test patch.global_elementary_factors == expected_factors
    @test constructive_product(patch.global_elementary_factors, R, entry.size) ==
          entry.expected.global_correction
    @test patch.patched_product == entry.expected.global_correction
    @test patch.target == entry.expected.global_correction
    @test replay.product == entry.expected.global_correction
    @test replay.target == entry.expected.global_correction
end

@testset "constructive Quillen induction expert acceptance" begin
    entries = Main.QuillenPatchFixtureCatalog.cases_by_id()
    fixture_ids = [
        "quillen-patched-substitution-witness-qq",
        "quillen-constructive-acceptance-gf2",
    ]

    built = Dict{String,Any}()
    for fixture_id in fixture_ids
        entry = entries[fixture_id]
        cover, local_certificates, normalized, patch = constructive_patch(entry)
        assert_constructive_patch_replays(
            entry,
            cover,
            local_certificates,
            normalized,
            patch,
        )
        built[fixture_id] = (;
            entry = entry,
            cover = cover,
            local_certificates = local_certificates,
            normalized = normalized,
            patch = patch,
        )
    end

    witness_bundle = built["quillen-patched-substitution-witness-qq"]
    witness_certificate = witness_bundle.local_certificates[1]
    witness = witness_certificate.patched_substitution_witness
    tampered_witness = merge(witness, (; shift = witness.shift + one(witness_bundle.patch.ring)))
    tampered_certificate = constructive_rebuild(
        witness_certificate;
        patched_substitution_witness = tampered_witness,
    )
    tampered_witness_patch = constructive_rebuild_global_patch(
        witness_bundle.patch;
        local_certificates = [
            tampered_certificate,
            witness_bundle.local_certificates[2],
        ],
    )
    @test !Suslin.verify_quillen_local_certificate(tampered_certificate)
    @test !Suslin.verify_quillen_patch(tampered_witness_patch)

    gf2_bundle = built["quillen-constructive-acceptance-gf2"]
    bad_cover = Suslin.QuillenDenominatorCoverCertificate(
        gf2_bundle.cover.ring,
        gf2_bundle.cover.denominators,
        [
            gf2_bundle.cover.coverage_multipliers[1] + one(gf2_bundle.cover.ring),
            gf2_bundle.cover.coverage_multipliers[2],
        ],
        gf2_bundle.cover.coverage_sum + gf2_bundle.cover.denominators[1],
        gf2_bundle.cover.verification,
    )
    @test !Suslin.verify_quillen_denominator_cover(bad_cover)
    @test_throws ArgumentError Suslin.assemble_deterministic_quillen_patch(
        gf2_bundle.entry.target_matrix,
        gf2_bundle.entry.substitution_variable,
        gf2_bundle.local_certificates,
        gf2_bundle.normalized,
        bad_cover;
        target = gf2_bundle.entry.expected.global_correction,
    )
    @test !Suslin.verify_quillen_patch(
        constructive_rebuild_global_patch(gf2_bundle.patch; cover_certificate = bad_cover),
    )

    tampered_local_certificates = copy(gf2_bundle.local_certificates)
    tampered_factors = copy(tampered_local_certificates[1].factors)
    tampered_factors[1] =
        tampered_factors[1] *
        elementary_matrix(gf2_bundle.patch.size, 1, 3, one(gf2_bundle.patch.ring), gf2_bundle.patch.ring)
    tampered_local_certificates[1] = constructive_rebuild(
        tampered_local_certificates[1];
        factors = tampered_factors,
    )
    @test !Suslin.verify_quillen_local_certificate(tampered_local_certificates[1])
    @test !Suslin.verify_quillen_patch(
        constructive_rebuild_global_patch(
            gf2_bundle.patch;
            local_certificates = tampered_local_certificates,
        ),
    )
end
```

- [ ] **Step 3: Register the acceptance test**

In `test/runtests.jl`, add this line immediately after
`"expert/quillen_patch_verification_hardening.jl",`:

```julia
        "expert/quillen_induction_constructive.jl",
```

- [ ] **Step 4: Run focused GREEN verification**

Run:

```bash
julia --project=. -e 'include("test/expert/quillen_induction_constructive.jl")'
```

Expected: exit 0.

- [ ] **Step 5: Run full requested suite**

Run:

```bash
julia --project=. test/runtests.jl all
```

Expected: exit 0.

- [ ] **Step 6: Run required Agent Desk package test**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: exit 0.

- [ ] **Step 7: Commit**

Run:

```bash
git add docs/superpowers/specs/2026-06-23-issue-105-quillen-induction-constructive-design.md docs/superpowers/plans/2026-06-23-issue-105-quillen-induction-constructive.md test/expert/quillen_induction_constructive.jl test/runtests.jl
git commit -m "test: add constructive quillen induction acceptance"
```

## Self-Review

- Spec coverage: the plan covers the focused acceptance command, full-suite registration, #99 fixture reuse, #103/#104 replay checks, no public API changes, and negative controls.
- Placeholder scan: no incomplete markers remain.
- Type consistency: helper and fixture names match existing #103/#104 test style and qualified `Suslin.<name>` access.
