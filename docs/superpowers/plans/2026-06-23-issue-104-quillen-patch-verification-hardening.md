# Issue 104 Quillen Patch Verification Hardening Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add focused tamper-resistance coverage for deterministic Quillen patch replay and make the global patch verification summary explicitly record count and denominator-data replay checks.

**Architecture:** Keep #104 in the existing expert/internal Quillen assembly surface. The test constructs a valid #103 `QuillenGlobalPatchAssembly`, rebuilds one-field tampered copies through direct struct constructors, and expects `verify_quillen_patch` to reject each. The implementation adds explicit `local_count_ok` and `denominator_data_ok` fields to `QuillenGlobalPatchAssemblyVerification`, recomputes them during replay, compares them in stored verification metadata, and includes them in `overall_ok`.

**Tech Stack:** Julia, Oscar exact ordinary polynomial rings, existing Suslin Quillen fixture catalog and deterministic patch assembly helpers, Test stdlib.

## Global Constraints

- No `AGENTS.md`, `CLAUDE.md`, `GEMINI.md`, or `CODEX.md` file is present in this checkout.
- The worker branch is `agent/issue-104-harden-quillen-patch-verification-against-tamper-run-1`.
- Dependency #103 is merged; reuse `QuillenGlobalPatchAssembly`, `assemble_deterministic_quillen_patch`, `replay_deterministic_quillen_patch`, and `verify_quillen_patch(::QuillenGlobalPatchAssembly)`.
- Keep verifier-related #103 assembly names expert/internal and do not export them from `src/Suslin.jl`.
- Tests must use qualified `Suslin.<name>` access for unexported expert/internal names.
- Do not change `test/public/api_surface.jl`.
- Preserve `test/expert/quillen_patching_exact.jl`.
- Constructors continue to throw `ArgumentError` for invalid inputs; verifiers return `false` for malformed stored patches.
- Do not broaden the supported mathematical family.
- Do not optimize factor count.
- Do not change public `elementary_factorization` behavior.
- Focused command is `julia --project=. -e 'include("test/expert/quillen_patch_verification_hardening.jl")'`.
- Expert group command is `julia --project=. test/runtests.jl expert`.
- Full package command is `julia --project=. -e 'using Pkg; Pkg.test()'`.

---

## File Structure

- Create `test/expert/quillen_patch_verification_hardening.jl`: fixture-backed regression test for the #104 tamper surfaces and explicit replay metadata fields.
- Modify `test/runtests.jl`: register `expert/quillen_patch_verification_hardening.jl` immediately after `expert/quillen_global_patch_assembly.jl`.
- Modify `src/algorithm/quillen_induction.jl`: add `local_count_ok` and `denominator_data_ok` to `QuillenGlobalPatchAssemblyVerification`, replay construction, stored-verification comparison, and provisional constructor placeholder.
- Leave `src/Suslin.jl` and `test/public/api_surface.jl` unchanged.

---

### Task 1: RED Quillen Patch Verification Hardening Test

**Files:**
- Create: `test/expert/quillen_patch_verification_hardening.jl`
- Modify: `test/runtests.jl`

**Interfaces:**
- Consumes: `Suslin.assemble_deterministic_quillen_patch`, `Suslin.replay_deterministic_quillen_patch`, `Suslin.verify_quillen_patch`, direct constructors for expert/internal Quillen records, and `test/fixtures/quillen_patch_cases.jl`.
- Produces: A focused expert test that currently fails because replay summaries do not expose `local_count_ok` and `denominator_data_ok`, and that verifies every requested one-field tamper returns `false`.

- [ ] **Step 1: Add the new expert test file**

Create `test/expert/quillen_patch_verification_hardening.jl` with:

```julia
using Test
using Suslin
using Oscar

const QUILLEN_PATCH_HARDENING_CATALOG_PATH =
    joinpath(@__DIR__, "..", "fixtures", "quillen_patch_cases.jl")
if !isdefined(Main, :QuillenPatchFixtureCatalog)
    include(QUILLEN_PATCH_HARDENING_CATALOG_PATH)
end

function hardening_rebuild(record; kwargs...)
    overrides = Dict{Symbol,Any}()
    for pair in kwargs
        overrides[pair.first] = pair.second
    end
    values = map(fieldnames(typeof(record))) do name
        get(overrides, name, getproperty(record, name))
    end
    return typeof(record)(values...)
end

function hardening_extra_factor(patch)
    return elementary_matrix(patch.size, 1, patch.size, one(patch.ring), patch.ring)
end

function hardening_local_certificate(local_factor)
    return Suslin.LocalCertificate(
        local_factor.certificate.indices,
        local_factor.certificate.denominators,
    )
end

function hardening_correction(local_factor)
    correction = local_factor.correction
    return Suslin.QuillenElementaryCorrection(
        correction.row,
        correction.col,
        correction.entry,
    )
end

function hardening_local_certificate_from_fixture(entry; local_index::Int)
    local_factor = entry.local_factors[local_index]
    return Suslin.quillen_local_realization_certificate(
        entry.target_matrix,
        entry.substitution_variable;
        local_certificate = hardening_local_certificate(local_factor),
        denominator = local_factor.denominator,
        coverage_multiplier = local_factor.coverage_multiplier,
        correction = hardening_correction(local_factor),
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

function hardening_cover(entry)
    denominators = [data.denominator for data in entry.denominator_data]
    multipliers = [data.coverage_multiplier for data in entry.denominator_data]
    return Suslin.quillen_denominator_cover_certificate(
        entry.ring.object,
        denominators,
        multipliers,
    )
end

function hardening_inputs(entry)
    cover = hardening_cover(entry)
    local_certificates = [
        hardening_local_certificate_from_fixture(entry; local_index = idx)
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

function hardening_valid_patch()
    entries = Main.QuillenPatchFixtureCatalog.cases_by_id()
    entry = entries["quillen-patched-substitution-witness-qq"]
    cover, local_certificates, normalized = hardening_inputs(entry)
    patch = Suslin.assemble_deterministic_quillen_patch(
        entry.target_matrix,
        entry.substitution_variable,
        local_certificates,
        normalized,
        cover;
        target = entry.expected.global_correction,
    )
    return entry, patch
end

function hardening_tamper_cover_multiplier(patch)
    cover = patch.cover_certificate
    coverage_multipliers = copy(cover.coverage_multipliers)
    coverage_multipliers[1] += one(patch.ring)
    tampered_cover = hardening_rebuild(
        cover;
        coverage_multipliers = coverage_multipliers,
    )
    return hardening_rebuild(patch; cover_certificate = tampered_cover)
end

function hardening_tamper_local_certificate_factor(patch)
    local_certificates = copy(patch.local_certificates)
    factors = copy(local_certificates[1].factors)
    factors[1] = factors[1] * hardening_extra_factor(patch)
    local_certificates[1] = hardening_rebuild(local_certificates[1]; factors = factors)
    return hardening_rebuild(patch; local_certificates = local_certificates)
end

function hardening_tamper_patched_substitution_witness(patch)
    local_certificates = copy(patch.local_certificates)
    witness = local_certificates[1].patched_substitution_witness
    tampered_witness = merge(witness, (; shift = witness.shift + one(patch.ring)))
    local_certificates[1] = hardening_rebuild(
        local_certificates[1];
        patched_substitution_witness = tampered_witness,
    )
    return hardening_rebuild(patch; local_certificates = local_certificates)
end

function hardening_tamper_normalized_denominator(patch)
    normalized = copy(patch.normalized_local_contributions)
    normalized[1] = hardening_rebuild(
        normalized[1];
        denominator = normalized[1].denominator + one(patch.ring),
    )
    return hardening_rebuild(patch; normalized_local_contributions = normalized)
end

function hardening_tamper_global_factor(patch)
    factors = copy(patch.global_elementary_factors)
    factors[1] = factors[1] * hardening_extra_factor(patch)
    return hardening_rebuild(patch; global_elementary_factors = factors)
end

function hardening_tamper_stored_product(patch)
    return hardening_rebuild(
        patch;
        patched_product = patch.patched_product * hardening_extra_factor(patch),
    )
end

function hardening_tamper_verification_summary(patch)
    tampered_verification = hardening_rebuild(
        patch.verification;
        overall_ok = !patch.verification.overall_ok,
    )
    return hardening_rebuild(patch; verification = tampered_verification)
end

@testset "Quillen patch verifier rejects tampered replay data" begin
    _, patch = hardening_valid_patch()
    @test Suslin.verify_quillen_patch(patch)

    replay = Suslin.replay_deterministic_quillen_patch(patch)
    @test hasproperty(replay, :local_count_ok)
    @test hasproperty(replay, :denominator_data_ok)
    if hasproperty(replay, :local_count_ok)
        @test replay.local_count_ok
    end
    if hasproperty(replay, :denominator_data_ok)
        @test replay.denominator_data_ok
    end
    @test replay.cover_certificate_ok
    @test replay.local_certificates_ok
    @test replay.normalized_contributions_ok
    @test replay.global_elementary_factors_ok
    @test replay.product_ok
    @test replay.target_ok
    @test replay.replay_metadata_ok
    @test replay.overall_ok

    tampered_cases = [
        "cover multiplier" => hardening_tamper_cover_multiplier(patch),
        "local certificate factor" => hardening_tamper_local_certificate_factor(patch),
        "patched-substitution witness" => hardening_tamper_patched_substitution_witness(patch),
        "normalized contribution denominator" => hardening_tamper_normalized_denominator(patch),
        "global elementary factor" => hardening_tamper_global_factor(patch),
        "stored product" => hardening_tamper_stored_product(patch),
        "stored verification summary" => hardening_tamper_verification_summary(patch),
    ]

    for (label, tampered_patch) in tampered_cases
        @test !Suslin.verify_quillen_patch(tampered_patch)
    end

    upstream_corrupted = hardening_tamper_local_certificate_factor(patch)
    negative_control = hardening_rebuild(upstream_corrupted; patched_product = patch.target)
    @test negative_control.patched_product == negative_control.target
    @test !Suslin.verify_quillen_patch(negative_control)
end
```

- [ ] **Step 2: Register the expert test**

In `test/runtests.jl`, add this line immediately after
`"expert/quillen_global_patch_assembly.jl",`:

```julia
        "expert/quillen_patch_verification_hardening.jl",
```

- [ ] **Step 3: Run RED verification**

Run:

```bash
julia --project=. -e 'include("test/expert/quillen_patch_verification_hardening.jl")'
```

Expected before Task 2: exits nonzero with failed `hasproperty(replay, :local_count_ok)` and `hasproperty(replay, :denominator_data_ok)` assertions.

- [ ] **Step 4: Commit the RED test**

Run:

```bash
git add test/expert/quillen_patch_verification_hardening.jl test/runtests.jl
git commit -m "test: cover quillen patch verifier tamper rejection"
```

### Task 2: Explicit Global Patch Replay Metadata Checks

**Files:**
- Modify: `src/algorithm/quillen_induction.jl`

**Interfaces:**
- Consumes: existing global patch replay inputs and helpers in `src/algorithm/quillen_induction.jl`.
- Produces: `QuillenGlobalPatchAssemblyVerification.local_count_ok::Bool` and `QuillenGlobalPatchAssemblyVerification.denominator_data_ok::Bool`, both included in replay, stored verification comparison, and `overall_ok`.

- [ ] **Step 1: Add fields to the verification struct**

Modify `QuillenGlobalPatchAssemblyVerification` so the beginning of the struct is:

```julia
struct QuillenGlobalPatchAssemblyVerification
    cover_certificate_ok::Bool
    local_certificates_ok::Bool
    normalized_contributions_ok::Bool
    local_count_ok::Bool
    local_alignment_ok::Bool
    cover_alignment_ok::Bool
    normalized_input_ok::Bool
    selected_variable_ok::Bool
    denominator_data_ok::Bool
    coverage_sum
    coverage_ok::Bool
    global_elementary_factors::Vector
    global_elementary_factors_ok::Bool
    product
    product_ok::Bool
    target
    target_ok::Bool
    replay_metadata
    replay_metadata_ok::Bool
    overall_ok::Bool
end
```

- [ ] **Step 2: Compare the new stored verification fields**

In `_same_quillen_global_patch_verification`, add comparisons for the two new
fields:

```julia
           left.normalized_contributions_ok == right.normalized_contributions_ok &&
           left.local_count_ok == right.local_count_ok &&
           left.local_alignment_ok == right.local_alignment_ok &&
```

and:

```julia
           left.selected_variable_ok == right.selected_variable_ok &&
           left.denominator_data_ok == right.denominator_data_ok &&
           left.coverage_sum == right.coverage_sum &&
```

- [ ] **Step 3: Store count and denominator-data replay results**

In `replay_deterministic_quillen_patch`, rename the existing `count_ok`
binding to `local_count_ok`, reuse it for `local_alignment_ok`, and keep
`denominator_data_ok` as the explicit value returned in verification:

```julia
    local_count_ok = length(local_certificates) == length(normalized)

    local_alignment_ok = local_count_ok && all(eachindex(local_certificates)) do idx
        _same_quillen_local_certificate_data(
            normalized[idx].local_certificate,
            local_certificates[idx],
        )
    end
```

Keep the existing denominator-data calculation, but include
`denominator_data_ok` in `overall_ok`:

```julia
    overall_ok =
        cover_certificate_ok &&
        local_certificates_ok &&
        normalized_contributions_ok &&
        local_count_ok &&
        local_alignment_ok &&
        cover_alignment_ok &&
        normalized_input_ok &&
        selected_variable_ok &&
        denominator_data_ok &&
        coverage_ok &&
        global_elementary_factors_ok &&
        product_ok &&
        target_ok &&
        replay_metadata_ok
```

Update the `QuillenGlobalPatchAssemblyVerification(...)` return call so the
new arguments appear after `normalized_contributions_ok` and after
`selected_variable_ok`:

```julia
        normalized_contributions_ok,
        local_count_ok,
        local_alignment_ok,
```

and:

```julia
        selected_variable_ok,
        denominator_data_ok,
        coverage_sum,
```

- [ ] **Step 4: Update the provisional false verification placeholder**

In `assemble_deterministic_quillen_patch`, update the provisional
`QuillenGlobalPatchAssemblyVerification(...)` call to include the two new
`false` placeholders after `normalized_contributions_ok` and after
`selected_variable_ok`:

```julia
            false,
            false,
            false,
            false,
            false,
            false,
            false,
            false,
            false,
            zero(R),
```

The nine `false` values correspond to:
`cover_certificate_ok`, `local_certificates_ok`,
`normalized_contributions_ok`, `local_count_ok`, `local_alignment_ok`,
`cover_alignment_ok`, `normalized_input_ok`, `selected_variable_ok`, and
`denominator_data_ok`.

- [ ] **Step 5: Run GREEN verification**

Run:

```bash
julia --project=. -e 'include("test/expert/quillen_patch_verification_hardening.jl")'
```

Expected after Task 2: exits 0 and the new hardening test passes.

- [ ] **Step 6: Run adjacent regression coverage**

Run:

```bash
julia --project=. -e 'include("test/expert/quillen_global_patch_assembly.jl"); include("test/expert/quillen_patching_exact.jl")'
```

Expected: exits 0.

- [ ] **Step 7: Commit the implementation**

Run:

```bash
git add src/algorithm/quillen_induction.jl
git commit -m "feat: harden quillen patch replay metadata"
```

### Task 3: Verification Gate

**Files:**
- No source files.

**Interfaces:**
- Consumes: completed Task 1 and Task 2 commits.
- Produces: fresh verification evidence before publishing the PR.

- [ ] **Step 1: Run the focused #104 command**

Run:

```bash
julia --project=. -e 'include("test/expert/quillen_patch_verification_hardening.jl")'
```

Expected: exits 0.

- [ ] **Step 2: Run the registered expert group**

Run:

```bash
julia --project=. test/runtests.jl expert
```

Expected: exits 0.

- [ ] **Step 3: Run the Agent Desk required package command**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: exits 0.

- [ ] **Step 4: Confirm public API remains unchanged**

Run:

```bash
git diff -- src/Suslin.jl test/public/api_surface.jl
```

Expected: no output.

## Self-Review

- The plan covers every issue-required tamper case and the final-product negative control.
- The implementation scope is limited to expert/internal verification metadata.
- No public API file changes are planned.
- Each implementation task has a focused RED/GREEN command.
