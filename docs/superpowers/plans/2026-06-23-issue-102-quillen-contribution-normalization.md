# Issue 102 Quillen Contribution Normalization Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add replayable Quillen local contribution normalization records that bind #100 local certificates to #101 exact denominator covers.

**Architecture:** Keep the new normalization surface expert/internal in `src/algorithm/quillen_induction.jl`. A normalized record stores the verified local certificate, verified cover certificate, exact cover pair index, patched-substitution replay summary, local product, weighted global elementary factor, and replay metadata; verification recomputes each part from the source certificate and cover.

**Tech Stack:** Julia, Oscar exact polynomial rings, existing Suslin Quillen helpers, Test stdlib.

## Global Constraints

- No `AGENTS.md`, `CLAUDE.md`, or `GEMINI.md` file is present in this checkout.
- Dependencies #100 and #101 are closed and merged; reuse `QuillenLocalRealizationCertificate` and `QuillenDenominatorCoverCertificate`.
- Keep `QuillenLocalContribution` compatible with existing exact patching tests.
- Do not implement the final global patch constructor.
- Do not broaden to arbitrary local solver calls.
- Do not route through `elementary_factorization`.
- Keep new normalization names expert/internal and do not export them.
- Tests must use qualified `Suslin.<name>` access.
- Focused command is `julia --project=. -e 'include("test/expert/quillen_contribution_normalization.jl")'`.
- Expert group command is `julia --project=. test/runtests.jl expert`.
- Full package command is `julia --project=. -e 'using Pkg; Pkg.test()'`.

---

## File Structure

- Modify `src/algorithm/quillen_induction.jl`: add `QuillenLocalContributionNormalization`, `QuillenLocalContributionNormalizationVerification`, normalization constructors, replay, and verifier.
- Create `test/expert/quillen_contribution_normalization.jl`: normalize two fixture-backed local certificates against a verified cover and check tamper controls.
- Modify `test/runtests.jl`: register the focused test in the `expert` group after `expert/quillen_local_certificate.jl`.
- Leave `src/Suslin.jl` and `test/public/api_surface.jl` unchanged.

---

### Task 1: RED Expert Normalization Test

**Files:**
- Create: `test/expert/quillen_contribution_normalization.jl`
- Modify: `test/runtests.jl`

**Interfaces:**
- Consumes: #99 fixture catalog entries, #100 `Suslin.quillen_local_realization_certificate`, and #101 `Suslin.quillen_denominator_cover_certificate`.
- Produces: RED coverage for `Suslin.normalize_quillen_local_contributions`, `Suslin.replay_quillen_local_contribution_normalization`, and `Suslin.verify_quillen_local_contribution_normalization`.

- [ ] **Step 1: Write the failing test**

Create `test/expert/quillen_contribution_normalization.jl` with:

```julia
using Test
using Suslin
using Oscar

const QUILLEN_PATCH_CATALOG_PATH = joinpath(@__DIR__, "..", "fixtures", "quillen_patch_cases.jl")
if !isdefined(Main, :QuillenPatchFixtureCatalog)
    include(QUILLEN_PATCH_CATALOG_PATH)
end

function contribution_normalization_test_product(factors, R, n)
    product = identity_matrix(R, n)
    for factor in factors
        product *= factor
    end
    return product
end

function contribution_normalization_local_certificate(local_factor)
    return Suslin.LocalCertificate(
        local_factor.certificate.indices,
        local_factor.certificate.denominators,
    )
end

function contribution_normalization_correction(local_factor)
    correction = local_factor.correction
    return Suslin.QuillenElementaryCorrection(
        correction.row,
        correction.col,
        correction.entry,
    )
end

function contribution_normalization_local_certificate_from_fixture(entry; local_index::Int)
    local_factor = entry.local_factors[local_index]
    return Suslin.quillen_local_realization_certificate(
        entry.target_matrix,
        entry.substitution_variable;
        local_certificate = contribution_normalization_local_certificate(local_factor),
        denominator = local_factor.denominator,
        coverage_multiplier = local_factor.coverage_multiplier,
        correction = contribution_normalization_correction(local_factor),
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

function contribution_normalization_cover(entry)
    denominators = [data.denominator for data in entry.denominator_data]
    multipliers = [data.coverage_multiplier for data in entry.denominator_data]
    return Suslin.quillen_denominator_cover_certificate(entry.ring.object, denominators, multipliers)
end

function contribution_normalization_expected_factor(entry, local_index::Int)
    local_factor = entry.local_factors[local_index]
    correction = local_factor.correction
    weighted_entry =
        local_factor.coverage_multiplier * local_factor.denominator * correction.entry
    return elementary_matrix(
        entry.size,
        correction.row,
        correction.col,
        weighted_entry,
        entry.ring.object,
    )
end

function rebuild_normalized_contribution(normalized; kwargs...)
    fields = merge((
        local_certificate = normalized.local_certificate,
        cover_certificate = normalized.cover_certificate,
        original_input = normalized.original_input,
        selected_variable = normalized.selected_variable,
        denominator = normalized.denominator,
        coverage_multiplier = normalized.coverage_multiplier,
        cover_index = normalized.cover_index,
        patched_substitution_witness = normalized.patched_substitution_witness,
        patched_substitution = normalized.patched_substitution,
        local_product = normalized.local_product,
        local_correction = normalized.local_correction,
        local_contribution = normalized.local_contribution,
        weighted_global_elementary_factors =
            normalized.weighted_global_elementary_factors,
        replay_metadata = normalized.replay_metadata,
        verification = normalized.verification,
    ), kwargs)
    return Suslin.QuillenLocalContributionNormalization(
        fields.local_certificate,
        fields.cover_certificate,
        fields.original_input,
        fields.selected_variable,
        fields.denominator,
        fields.coverage_multiplier,
        fields.cover_index,
        fields.patched_substitution_witness,
        fields.patched_substitution,
        fields.local_product,
        fields.local_correction,
        fields.local_contribution,
        fields.weighted_global_elementary_factors,
        fields.replay_metadata,
        fields.verification,
    )
end

@testset "Quillen local contribution normalization" begin
    entries = Main.QuillenPatchFixtureCatalog.cases_by_id()
    entry = entries["quillen-patched-substitution-witness-qq"]
    cover = contribution_normalization_cover(entry)
    local_certificates = [
        contribution_normalization_local_certificate_from_fixture(entry; local_index = idx)
        for idx in eachindex(entry.local_factors)
    ]

    normalized = Suslin.normalize_quillen_local_contributions(
        local_certificates,
        cover;
        original_input = entry.target_matrix,
        selected_variable = entry.substitution_variable,
    )

    @test length(normalized) == 2
    @test all(item -> item isa Suslin.QuillenLocalContributionNormalization, normalized)
    @test all(Suslin.verify_quillen_local_contribution_normalization, normalized)

    for (idx, item) in enumerate(normalized)
        replay = Suslin.replay_quillen_local_contribution_normalization(item)
        @test replay.overall_ok
        @test replay.local_certificate_ok
        @test replay.cover_certificate_ok
        @test replay.cover_pair_ok
        @test replay.local_product_ok
        @test replay.local_correction_ok
        @test replay.patched_substitution_ok
        @test replay.weighted_global_elementary_factors_ok
        @test replay.replay_metadata_ok

        local_certificate = local_certificates[idx]
        @test item.denominator == cover.denominators[item.cover_index]
        @test item.coverage_multiplier == cover.coverage_multipliers[item.cover_index]
        @test item.local_product == local_certificate.local_product
        @test item.local_correction == local_certificate.local_correction
        @test contribution_normalization_test_product(
            local_certificate.factors,
            local_certificate.ring,
            local_certificate.size,
        ) == item.local_product

        witness = item.patched_substitution_witness
        @test witness !== nothing
        @test item.patched_substitution.actual_matrix == witness.expected_matrix
        @test Suslin.patched_substitution(
            witness.matrix,
            witness.variable,
            witness.denominator,
            witness.exponent,
            witness.shift,
        ) == witness.expected_matrix

        expected_factor = contribution_normalization_expected_factor(entry, idx)
        @test item.weighted_global_elementary_factors == [expected_factor]
        @test replay.weighted_global_elementary_factors == [expected_factor]
        @test item.replay_metadata.cover_index == item.cover_index
        @test item.replay_metadata.local_witness_metadata.local_index == idx
    end

    base = normalized[1]
    R = base.local_certificate.ring
    replacement_variable = first(filter(gen -> gen != base.selected_variable, collect(gens(R))))
    @test_throws ArgumentError Suslin.normalize_quillen_local_contribution(
        base.local_certificate,
        base.cover_certificate;
        selected_variable = replacement_variable,
    )
    @test !Suslin.verify_quillen_local_contribution_normalization(
        rebuild_normalized_contribution(base; selected_variable = replacement_variable),
    )

    bad_exponent_witness = merge(
        base.patched_substitution_witness,
        (; exponent = base.patched_substitution_witness.exponent + 1),
    )
    @test_throws ArgumentError Suslin.normalize_quillen_local_contribution(
        base.local_certificate,
        base.cover_certificate;
        patched_substitution_witness = bad_exponent_witness,
    )
    @test !Suslin.verify_quillen_local_contribution_normalization(
        rebuild_normalized_contribution(base; patched_substitution_witness = bad_exponent_witness),
    )

    bad_shift_witness = merge(
        base.patched_substitution_witness,
        (; shift = base.patched_substitution_witness.shift + one(R)),
    )
    @test_throws ArgumentError Suslin.normalize_quillen_local_contribution(
        base.local_certificate,
        base.cover_certificate;
        patched_substitution_witness = bad_shift_witness,
    )
    @test !Suslin.verify_quillen_local_contribution_normalization(
        rebuild_normalized_contribution(base; patched_substitution_witness = bad_shift_witness),
    )

    @test !Suslin.verify_quillen_local_contribution_normalization(
        rebuild_normalized_contribution(
            base;
            coverage_multiplier = base.coverage_multiplier + one(R),
        ),
    )

    tampered_factor = base.weighted_global_elementary_factors[1] *
        elementary_matrix(base.local_certificate.size, 2, 3, one(R), R)
    @test !Suslin.verify_quillen_local_contribution_normalization(
        rebuild_normalized_contribution(
            base;
            weighted_global_elementary_factors = [tampered_factor],
        ),
    )
end
```

- [ ] **Step 2: Register the expert test**

In `test/runtests.jl`, add:

```julia
        "expert/quillen_contribution_normalization.jl",
```

immediately after `expert/quillen_local_certificate.jl`.

- [ ] **Step 3: Run RED verification**

Run:

```bash
julia --project=. -e 'include("test/expert/quillen_contribution_normalization.jl")'
```

Expected: FAIL with `UndefVarError` for the new normalization helper.

- [ ] **Step 4: Commit RED test**

Commit the test and runner registration:

```bash
git add test/expert/quillen_contribution_normalization.jl test/runtests.jl
git commit -m "test: cover quillen contribution normalization replay"
```

### Task 2: Normalization Implementation

**Files:**
- Modify: `src/algorithm/quillen_induction.jl`

**Interfaces:**
- Consumes: `_quillen_local_certificate_replay_summary`, `_quillen_local_patched_witness_summary`, `verify_quillen_local_certificate`, `verify_quillen_denominator_cover`, `_normalize_quillen_contribution`, `_quillen_factors`, `_same_quillen_factors`.
- Produces:
  - `QuillenLocalContributionNormalization`
  - `QuillenLocalContributionNormalizationVerification`
  - `normalize_quillen_local_contribution`
  - `normalize_quillen_local_contributions`
  - `replay_quillen_local_contribution_normalization`
  - `verify_quillen_local_contribution_normalization`

- [ ] **Step 1: Add normalization structs and cover pairing helpers**

In `src/algorithm/quillen_induction.jl`, after
`QuillenLocalRealizationCertificate`, add:

```julia
struct QuillenLocalContributionNormalizationVerification
    local_certificate_ok::Bool
    cover_certificate_ok::Bool
    original_input_ok::Bool
    selected_variable_ok::Bool
    cover_ring_ok::Bool
    cover_index_ok::Bool
    cover_pair
    cover_pair_ok::Bool
    patched_substitution
    patched_substitution_ok::Bool
    local_product
    local_product_ok::Bool
    local_correction
    local_correction_ok::Bool
    weighted_global_elementary_factors::Vector
    weighted_global_elementary_factors_ok::Bool
    replay_metadata
    replay_metadata_ok::Bool
    overall_ok::Bool
end

struct QuillenLocalContributionNormalization
    local_certificate::QuillenLocalRealizationCertificate
    cover_certificate::QuillenDenominatorCoverCertificate
    original_input
    selected_variable
    denominator
    coverage_multiplier
    cover_index::Int
    patched_substitution_witness
    patched_substitution
    local_product
    local_correction
    local_contribution::QuillenLocalContribution
    weighted_global_elementary_factors::Vector
    replay_metadata
    verification::QuillenLocalContributionNormalizationVerification
end

function _quillen_cover_pair_index(cover::QuillenDenominatorCoverCertificate, denominator, coverage_multiplier)
    matches = Int[]
    for idx in eachindex(cover.denominators)
        if cover.denominators[idx] == denominator &&
           cover.coverage_multipliers[idx] == coverage_multiplier
            push!(matches, idx)
        end
    end
    isempty(matches) &&
        throw(ArgumentError("local contribution denominator and coverage multiplier must match an exact cover pair"))
    length(matches) == 1 ||
        throw(ArgumentError("local contribution denominator and coverage multiplier match multiple cover pairs"))
    return only(matches)
end
```

- [ ] **Step 2: Add replay metadata and verification equality helpers**

Add:

```julia
function _quillen_normalization_replay_metadata(certificate, cover, cover_index::Int)
    return (;
        local_witness_metadata = certificate.witness_metadata,
        cover_index = cover_index,
        cover_denominator_count = length(cover.denominators),
        cover_coverage_sum = cover.coverage_sum,
    )
end

function _same_quillen_normalization_verification(
    left::QuillenLocalContributionNormalizationVerification,
    right::QuillenLocalContributionNormalizationVerification,
)::Bool
    return left.local_certificate_ok == right.local_certificate_ok &&
           left.cover_certificate_ok == right.cover_certificate_ok &&
           left.original_input_ok == right.original_input_ok &&
           left.selected_variable_ok == right.selected_variable_ok &&
           left.cover_ring_ok == right.cover_ring_ok &&
           left.cover_index_ok == right.cover_index_ok &&
           left.cover_pair == right.cover_pair &&
           left.cover_pair_ok == right.cover_pair_ok &&
           left.patched_substitution == right.patched_substitution &&
           left.patched_substitution_ok == right.patched_substitution_ok &&
           left.local_product == right.local_product &&
           left.local_product_ok == right.local_product_ok &&
           left.local_correction == right.local_correction &&
           left.local_correction_ok == right.local_correction_ok &&
           _same_quillen_factors(
               left.weighted_global_elementary_factors,
               right.weighted_global_elementary_factors,
           ) &&
           left.weighted_global_elementary_factors_ok ==
               right.weighted_global_elementary_factors_ok &&
           left.replay_metadata == right.replay_metadata &&
           left.replay_metadata_ok == right.replay_metadata_ok &&
           left.overall_ok == right.overall_ok
end
```

- [ ] **Step 3: Add replay function**

Add:

```julia
function replay_quillen_local_contribution_normalization(
    normalized::QuillenLocalContributionNormalization,
)
    certificate = normalized.local_certificate
    cover = normalized.cover_certificate
    R = _require_supported_quillen_ring(certificate.ring)
    n = certificate.size
    selected_variable = _require_substitution_generator(R, normalized.selected_variable)
    local_replay = _quillen_local_certificate_replay_summary(certificate)

    local_certificate_ok = verify_quillen_local_certificate(certificate)
    cover_certificate_ok = verify_quillen_denominator_cover(cover)
    original_input_ok =
        normalized.original_input == certificate.original_input &&
        local_replay.original_input == normalized.original_input
    selected_variable_ok =
        selected_variable == certificate.selected_variable &&
        local_replay.selected_variable == selected_variable
    cover_ring_ok = cover.ring == R
    cover_index_ok = 1 <= normalized.cover_index <= length(cover.denominators)
    cover_pair = cover_index_ok ?
        QuillenDenominatorData(
            cover.denominators[normalized.cover_index],
            cover.coverage_multipliers[normalized.cover_index],
        ) :
        nothing
    cover_pair_ok =
        cover_ring_ok &&
        cover_index_ok &&
        cover_pair.denominator == normalized.denominator &&
        cover_pair.coverage_multiplier == normalized.coverage_multiplier &&
        certificate.denominator == normalized.denominator &&
        certificate.coverage_multiplier == normalized.coverage_multiplier

    patched_substitution = _quillen_local_patched_witness_summary(
        normalized.patched_substitution_witness,
        R,
        n,
        selected_variable,
    )
    patched_substitution_ok =
        normalized.patched_substitution_witness == certificate.patched_substitution_witness &&
        patched_substitution.ok &&
        normalized.patched_substitution == patched_substitution

    replayed_contribution = _normalize_quillen_contribution(
        QuillenLocalContribution(
            certificate.local_certificate,
            normalized.denominator,
            normalized.coverage_multiplier,
            certificate.correction,
        ),
        R,
        n,
    )
    weighted_global_elementary_factors = _quillen_factors(R, n, [replayed_contribution])
    weighted_global_elementary_factors_ok = _same_quillen_factors(
        normalized.weighted_global_elementary_factors,
        weighted_global_elementary_factors,
    )
    local_product_ok = normalized.local_product == local_replay.local_product
    local_correction_ok = normalized.local_correction == local_replay.local_correction
    replay_metadata = _quillen_normalization_replay_metadata(certificate, cover, normalized.cover_index)
    replay_metadata_ok = normalized.replay_metadata == replay_metadata
    overall_ok =
        local_certificate_ok &&
        cover_certificate_ok &&
        original_input_ok &&
        selected_variable_ok &&
        cover_ring_ok &&
        cover_pair_ok &&
        patched_substitution_ok &&
        local_product_ok &&
        local_correction_ok &&
        weighted_global_elementary_factors_ok &&
        replay_metadata_ok

    return QuillenLocalContributionNormalizationVerification(
        local_certificate_ok,
        cover_certificate_ok,
        original_input_ok,
        selected_variable_ok,
        cover_ring_ok,
        cover_index_ok,
        cover_pair,
        cover_pair_ok,
        patched_substitution,
        patched_substitution_ok,
        local_replay.local_product,
        local_product_ok,
        local_replay.local_correction,
        local_correction_ok,
        weighted_global_elementary_factors,
        weighted_global_elementary_factors_ok,
        replay_metadata,
        replay_metadata_ok,
        overall_ok,
    )
end
```

- [ ] **Step 4: Add constructors and verifier**

Add:

```julia
function normalize_quillen_local_contribution(
    certificate::QuillenLocalRealizationCertificate,
    cover::QuillenDenominatorCoverCertificate;
    original_input = certificate.original_input,
    selected_variable = certificate.selected_variable,
    patched_substitution_witness = certificate.patched_substitution_witness,
)
    verify_quillen_local_certificate(certificate) ||
        throw(ArgumentError("Quillen local realization certificate does not replay"))
    verify_quillen_denominator_cover(cover) ||
        throw(ArgumentError("Quillen denominator cover certificate does not replay"))
    R = _require_supported_quillen_ring(certificate.ring)
    n = certificate.size
    cover.ring == R ||
        throw(ArgumentError("Quillen local contribution normalization requires cover and local certificate rings to match"))
    selected = _require_substitution_generator(R, selected_variable)
    selected == certificate.selected_variable ||
        throw(ArgumentError("Quillen local contribution normalization requires the selected variable recorded by the local certificate"))
    original_input == certificate.original_input ||
        throw(ArgumentError("Quillen local contribution normalization requires the original input recorded by the local certificate"))
    patched_substitution_witness == certificate.patched_substitution_witness ||
        throw(ArgumentError("Quillen local contribution normalization requires the patched-substitution witness recorded by the local certificate"))

    cover_index = _quillen_cover_pair_index(cover, certificate.denominator, certificate.coverage_multiplier)
    local_contribution = _normalize_quillen_contribution(
        QuillenLocalContribution(
            certificate.local_certificate,
            certificate.denominator,
            certificate.coverage_multiplier,
            certificate.correction,
        ),
        R,
        n,
    )
    weighted_global_elementary_factors = _quillen_factors(R, n, [local_contribution])
    patched_substitution = _quillen_local_patched_witness_summary(
        patched_substitution_witness,
        R,
        n,
        selected,
    )
    patched_substitution.ok ||
        throw(ArgumentError("Quillen local contribution normalization patched substitution does not replay"))
    replay_metadata = _quillen_normalization_replay_metadata(certificate, cover, cover_index)
    provisional = QuillenLocalContributionNormalization(
        certificate,
        cover,
        original_input,
        selected,
        certificate.denominator,
        certificate.coverage_multiplier,
        cover_index,
        patched_substitution_witness,
        patched_substitution,
        certificate.local_product,
        certificate.local_correction,
        local_contribution,
        weighted_global_elementary_factors,
        replay_metadata,
        QuillenLocalContributionNormalizationVerification(
            false,
            false,
            false,
            false,
            false,
            false,
            nothing,
            false,
            (; present = false, ok = false),
            false,
            certificate.local_product,
            false,
            certificate.local_correction,
            false,
            Any[],
            false,
            replay_metadata,
            false,
            false,
        ),
    )
    verification = replay_quillen_local_contribution_normalization(provisional)
    verification.overall_ok ||
        throw(ArgumentError("Quillen local contribution normalization data does not replay"))
    return QuillenLocalContributionNormalization(
        provisional.local_certificate,
        provisional.cover_certificate,
        provisional.original_input,
        provisional.selected_variable,
        provisional.denominator,
        provisional.coverage_multiplier,
        provisional.cover_index,
        provisional.patched_substitution_witness,
        provisional.patched_substitution,
        provisional.local_product,
        provisional.local_correction,
        provisional.local_contribution,
        provisional.weighted_global_elementary_factors,
        provisional.replay_metadata,
        verification,
    )
end

function normalize_quillen_local_contributions(
    certificates,
    cover::QuillenDenominatorCoverCertificate;
    original_input = nothing,
    selected_variable = nothing,
)
    return [
        normalize_quillen_local_contribution(
            certificate,
            cover;
            original_input = original_input === nothing ? certificate.original_input : original_input,
            selected_variable = selected_variable === nothing ? certificate.selected_variable : selected_variable,
        )
        for certificate in collect(certificates)
    ]
end

function verify_quillen_local_contribution_normalization(normalized)::Bool
    try
        replay = replay_quillen_local_contribution_normalization(normalized)
        return replay.overall_ok &&
               _same_quillen_normalization_verification(normalized.verification, replay)
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end
```

- [ ] **Step 5: Run GREEN focused verification**

Run:

```bash
julia --project=. -e 'include("test/expert/quillen_contribution_normalization.jl")'
```

Expected: PASS.

- [ ] **Step 6: Commit implementation**

Commit the implementation:

```bash
git add src/algorithm/quillen_induction.jl
git commit -m "feat: normalize quillen local contributions"
```

### Task 3: Verification and Review

**Files:**
- No planned source edits unless verification or review identifies a defect.

**Interfaces:**
- Consumes: Task 1 and Task 2.
- Produces: verified branch ready for PR.

- [ ] **Step 1: Run focused command**

Run:

```bash
julia --project=. -e 'include("test/expert/quillen_contribution_normalization.jl")'
```

Expected: PASS.

- [ ] **Step 2: Run expert group**

Run:

```bash
julia --project=. test/runtests.jl expert
```

Expected: PASS.

- [ ] **Step 3: Run full package command**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: PASS.

- [ ] **Step 4: Review changed files**

Run:

```bash
git status --short
git diff --stat origin/main..HEAD
git diff --name-status origin/main..HEAD
```

Expected: only the design/plan docs, `src/algorithm/quillen_induction.jl`,
`test/expert/quillen_contribution_normalization.jl`, and `test/runtests.jl`
are changed.

---

## Self-Review

- Every requirement from the design has a task.
- The plan keeps the API expert/internal and leaves public exports unchanged.
- The RED test covers two local certificate entries, patched substitution,
  local products, weighted factors, cover pairings, and the required negative
  controls.
- No broad local solver, final patch constructor, or `elementary_factorization`
  route is introduced.
