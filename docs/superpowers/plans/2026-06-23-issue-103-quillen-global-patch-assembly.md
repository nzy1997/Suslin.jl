# Issue 103 Quillen Global Patch Assembly Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add deterministic Quillen global patch assembly from verified local certificates, verified denominator covers, and #102 normalized local contributions.

**Architecture:** Keep the new assembly surface expert/internal in `src/algorithm/quillen_induction.jl`. The assembly record stores the cover, local certificates, normalized contributions, global elementary factors, exact product, target, replay metadata, and a replayed verification object; `verify_quillen_patch` gets a method for the new record while the public API surface remains unchanged.

**Tech Stack:** Julia, Oscar exact ordinary polynomial rings, existing Suslin Quillen certificate and normalization helpers, Test stdlib.

## Global Constraints

- No `AGENTS.md`, `CLAUDE.md`, or `GEMINI.md` file is present in this checkout.
- The worker branch is `agent/issue-103-assemble-deterministic-quillen-patches-into-glob-run-1`.
- Dependency #102 is merged; reuse `QuillenLocalContributionNormalization` and `verify_quillen_local_contribution_normalization`.
- Dependency #101 is merged; reuse `QuillenDenominatorCoverCertificate` and `verify_quillen_denominator_cover`.
- Dependency #100 is merged; reuse `QuillenLocalRealizationCertificate` and `verify_quillen_local_certificate`.
- Keep new assembly names expert/internal and do not export them from `src/Suslin.jl`.
- Tests must use qualified `Suslin.<name>` access for new names.
- Do not change `test/public/api_surface.jl`.
- Do not implement the final public Park-Woodburn `elementary_factorization` driver.
- Do not solve arbitrary local certificates inside this issue.
- Do not handle Laurent `GL_n` determinant correction.
- Focused command is `julia --project=. -e 'include("test/expert/quillen_global_patch_assembly.jl")'`.
- Expert group command is `julia --project=. test/runtests.jl expert`.
- Full package command is `julia --project=. -e 'using Pkg; Pkg.test()'`.

---

## File Structure

- Modify `src/algorithm/quillen_induction.jl`: add `QuillenGlobalPatchAssembly`, `QuillenGlobalPatchAssemblyVerification`, exact comparison helpers, replay, constructor, and `verify_quillen_patch` method.
- Create `test/expert/quillen_global_patch_assembly.jl`: assemble two deterministic ordinary-polynomial patches from fixture catalog entries and cover negative controls.
- Modify `test/runtests.jl`: register `expert/quillen_global_patch_assembly.jl` immediately after `expert/quillen_contribution_normalization.jl`.
- Leave `src/Suslin.jl` and `test/public/api_surface.jl` unchanged.

---

### Task 1: RED Expert Global Assembly Test

**Files:**
- Create: `test/expert/quillen_global_patch_assembly.jl`
- Modify: `test/runtests.jl`

**Interfaces:**
- Consumes: fixture catalog entries, #100 `Suslin.quillen_local_realization_certificate`, #101 `Suslin.quillen_denominator_cover_certificate`, and #102 `Suslin.normalize_quillen_local_contributions`.
- Produces: RED coverage for `Suslin.assemble_deterministic_quillen_patch`, `Suslin.replay_deterministic_quillen_patch`, `Suslin.QuillenGlobalPatchAssembly`, and `verify_quillen_patch(::Suslin.QuillenGlobalPatchAssembly)`.

- [ ] **Step 1: Write the failing test**

Create `test/expert/quillen_global_patch_assembly.jl` with:

```julia
using Test
using Suslin
using Oscar

const QUILLEN_GLOBAL_PATCH_CATALOG_PATH =
    joinpath(@__DIR__, "..", "fixtures", "quillen_patch_cases.jl")
if !isdefined(Main, :QuillenPatchFixtureCatalog)
    include(QUILLEN_GLOBAL_PATCH_CATALOG_PATH)
end

function global_patch_test_product(factors, R, n)
    product = identity_matrix(R, n)
    for factor in factors
        product *= factor
    end
    return product
end

function global_patch_local_certificate(local_factor)
    return Suslin.LocalCertificate(
        local_factor.certificate.indices,
        local_factor.certificate.denominators,
    )
end

function global_patch_correction(local_factor)
    correction = local_factor.correction
    return Suslin.QuillenElementaryCorrection(
        correction.row,
        correction.col,
        correction.entry,
    )
end

function global_patch_local_certificate_from_fixture(entry; local_index::Int)
    local_factor = entry.local_factors[local_index]
    return Suslin.quillen_local_realization_certificate(
        entry.target_matrix,
        entry.substitution_variable;
        local_certificate = global_patch_local_certificate(local_factor),
        denominator = local_factor.denominator,
        coverage_multiplier = local_factor.coverage_multiplier,
        correction = global_patch_correction(local_factor),
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

function global_patch_cover(entry)
    denominators = [data.denominator for data in entry.denominator_data]
    multipliers = [data.coverage_multiplier for data in entry.denominator_data]
    return Suslin.quillen_denominator_cover_certificate(
        entry.ring.object,
        denominators,
        multipliers,
    )
end

function global_patch_inputs(entry)
    cover = global_patch_cover(entry)
    local_certificates = [
        global_patch_local_certificate_from_fixture(entry; local_index = idx)
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

function rebuild_global_patch(patch; kwargs...)
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

@testset "deterministic Quillen global patch assembly" begin
    entries = Main.QuillenPatchFixtureCatalog.cases_by_id()
    fixture_ids = [
        "quillen-patched-substitution-witness-qq",
        "quillen-nontrivial-multipliers-qq",
    ]

    for fixture_id in fixture_ids
        entry = entries[fixture_id]
        cover, local_certificates, normalized = global_patch_inputs(entry)
        patch = Suslin.assemble_deterministic_quillen_patch(
            entry.target_matrix,
            entry.substitution_variable,
            local_certificates,
            normalized,
            cover;
            target = entry.expected.global_correction,
        )

        @test patch isa Suslin.QuillenGlobalPatchAssembly
        @test verify_quillen_patch(patch)
        @test patch.ring == entry.ring.object
        @test patch.size == entry.size
        @test patch.substitution_variable == entry.substitution_variable
        @test patch.original_input == entry.target_matrix
        @test patch.cover_certificate == cover
        @test patch.local_certificates == local_certificates
        @test patch.normalized_local_contributions == normalized
        @test length(patch.global_elementary_factors) == length(entry.local_factors)
        @test patch.global_elementary_factors ==
            reduce(vcat, [item.weighted_global_elementary_factors for item in normalized])
        @test global_patch_test_product(
            patch.global_elementary_factors,
            patch.ring,
            patch.size,
        ) == entry.expected.global_correction
        @test patch.patched_product == entry.expected.global_correction
        @test patch.target == entry.expected.global_correction

        replay = Suslin.replay_deterministic_quillen_patch(patch)
        @test replay.overall_ok
        @test replay.cover_certificate_ok
        @test replay.local_certificates_ok
        @test replay.normalized_contributions_ok
        @test replay.local_alignment_ok
        @test replay.cover_alignment_ok
        @test replay.normalized_input_ok
        @test replay.selected_variable_ok
        @test replay.global_elementary_factors_ok
        @test replay.product_ok
        @test replay.replay_metadata_ok
        @test replay.coverage_sum == one(entry.ring.object)
        @test replay.product == entry.expected.global_correction
        @test replay.target == entry.expected.global_correction
        @test patch.replay_metadata.fixture_count == length(local_certificates)
    end

    entry = entries["quillen-patched-substitution-witness-qq"]
    cover, local_certificates, normalized = global_patch_inputs(entry)

    wrong_local_order = copy(local_certificates)
    wrong_local_order[1], wrong_local_order[2] = wrong_local_order[2], wrong_local_order[1]
    @test_throws ArgumentError Suslin.assemble_deterministic_quillen_patch(
        entry.target_matrix,
        entry.substitution_variable,
        wrong_local_order,
        normalized,
        cover;
        target = entry.expected.global_correction,
    )

    patch = Suslin.assemble_deterministic_quillen_patch(
        entry.target_matrix,
        entry.substitution_variable,
        local_certificates,
        normalized,
        cover;
        target = entry.expected.global_correction,
    )
    tampered_local_patch = rebuild_global_patch(patch; local_certificates = wrong_local_order)
    @test !verify_quillen_patch(tampered_local_patch)

    bad_cover = Suslin.QuillenDenominatorCoverCertificate(
        cover.ring,
        cover.denominators,
        [cover.coverage_multipliers[1] + one(cover.ring), cover.coverage_multipliers[2]],
        cover.coverage_sum + cover.denominators[1],
        cover.verification,
    )
    @test !Suslin.verify_quillen_denominator_cover(bad_cover)
    @test_throws ArgumentError Suslin.assemble_deterministic_quillen_patch(
        entry.target_matrix,
        entry.substitution_variable,
        local_certificates,
        normalized,
        bad_cover;
        target = entry.expected.global_correction,
    )
    tampered_cover_patch = rebuild_global_patch(patch; cover_certificate = bad_cover)
    @test !verify_quillen_patch(tampered_cover_patch)
end
```

- [ ] **Step 2: Register the expert test**

In `test/runtests.jl`, add:

```julia
        "expert/quillen_global_patch_assembly.jl",
```

immediately after `expert/quillen_contribution_normalization.jl`.

- [ ] **Step 3: Run RED verification**

Run:

```bash
julia --project=. -e 'include("test/expert/quillen_global_patch_assembly.jl")'
```

Expected: FAIL with `UndefVarError` for `assemble_deterministic_quillen_patch` or `QuillenGlobalPatchAssembly`.

- [ ] **Step 4: Commit RED test**

Commit the test and runner registration:

```bash
git add test/expert/quillen_global_patch_assembly.jl test/runtests.jl
git commit -m "test: cover deterministic quillen global patch assembly"
```

### Task 2: Global Assembly Implementation

**Files:**
- Modify: `src/algorithm/quillen_induction.jl`

**Interfaces:**
- Consumes: `verify_quillen_local_certificate`, `verify_quillen_denominator_cover`, `verify_quillen_local_contribution_normalization`, `replay_quillen_local_contribution_normalization`, `_quillen_denominator_data`, `_quillen_product`, `_same_quillen_denominator_data`, `_same_quillen_factors`, `_same_quillen_denominator_cover_verification`, and `_same_quillen_normalization_verification`.
- Produces: `QuillenGlobalPatchAssembly`, `QuillenGlobalPatchAssemblyVerification`, `assemble_deterministic_quillen_patch`, `replay_deterministic_quillen_patch`, and `verify_quillen_patch(::QuillenGlobalPatchAssembly)`.

- [ ] **Step 1: Add assembly structs and equality helpers**

Append this code to `src/algorithm/quillen_induction.jl` after
`verify_quillen_local_contribution_normalization`:

```julia
struct QuillenGlobalPatchAssemblyVerification
    cover_certificate_ok::Bool
    local_certificates_ok::Bool
    normalized_contributions_ok::Bool
    local_alignment_ok::Bool
    cover_alignment_ok::Bool
    normalized_input_ok::Bool
    selected_variable_ok::Bool
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

struct QuillenGlobalPatchAssembly
    ring
    size::Int
    substitution_variable
    original_input
    cover_certificate::QuillenDenominatorCoverCertificate
    denominator_data::Vector{QuillenDenominatorData}
    local_certificates::Vector{QuillenLocalRealizationCertificate}
    normalized_local_contributions::Vector{QuillenLocalContributionNormalization}
    global_elementary_factors::Vector
    patched_product
    target
    replay_metadata
    verification::QuillenGlobalPatchAssemblyVerification
end

function _same_quillen_local_certificate_data(
    left::QuillenLocalRealizationCertificate,
    right::QuillenLocalRealizationCertificate,
)::Bool
    return left.original_input == right.original_input &&
           left.ring == right.ring &&
           left.size == right.size &&
           left.selected_variable == right.selected_variable &&
           left.local_certificate.indices == right.local_certificate.indices &&
           left.local_certificate.denominators == right.local_certificate.denominators &&
           left.denominator == right.denominator &&
           left.coverage_multiplier == right.coverage_multiplier &&
           left.correction == right.correction &&
           _same_quillen_factors(left.factors, right.factors) &&
           left.local_product == right.local_product &&
           left.local_correction == right.local_correction &&
           left.patched_substitution_witness == right.patched_substitution_witness &&
           left.witness_metadata == right.witness_metadata &&
           left.verification == right.verification
end

function _same_quillen_cover_certificate_data(
    left::QuillenDenominatorCoverCertificate,
    right::QuillenDenominatorCoverCertificate,
)::Bool
    return left.ring == right.ring &&
           left.denominators == right.denominators &&
           left.coverage_multipliers == right.coverage_multipliers &&
           left.coverage_sum == right.coverage_sum &&
           _same_quillen_denominator_cover_verification(left.verification, right.verification)
end

function _quillen_global_patch_replay_metadata(
    cover::QuillenDenominatorCoverCertificate,
    local_certificates,
    normalized_contributions,
)
    return (;
        fixture_count = length(local_certificates),
        normalized_count = length(normalized_contributions),
        cover_denominator_count = length(cover.denominators),
        cover_coverage_sum = cover.coverage_sum,
        local_witness_metadata = [certificate.witness_metadata for certificate in local_certificates],
        normalized_cover_indices = [normalized.cover_index for normalized in normalized_contributions],
    )
end

function _same_quillen_global_patch_verification(
    left::QuillenGlobalPatchAssemblyVerification,
    right::QuillenGlobalPatchAssemblyVerification,
)::Bool
    return left.cover_certificate_ok == right.cover_certificate_ok &&
           left.local_certificates_ok == right.local_certificates_ok &&
           left.normalized_contributions_ok == right.normalized_contributions_ok &&
           left.local_alignment_ok == right.local_alignment_ok &&
           left.cover_alignment_ok == right.cover_alignment_ok &&
           left.normalized_input_ok == right.normalized_input_ok &&
           left.selected_variable_ok == right.selected_variable_ok &&
           left.coverage_sum == right.coverage_sum &&
           left.coverage_ok == right.coverage_ok &&
           _same_quillen_factors(
               left.global_elementary_factors,
               right.global_elementary_factors,
           ) &&
           left.global_elementary_factors_ok == right.global_elementary_factors_ok &&
           left.product == right.product &&
           left.product_ok == right.product_ok &&
           left.target == right.target &&
           left.target_ok == right.target_ok &&
           left.replay_metadata == right.replay_metadata &&
           left.replay_metadata_ok == right.replay_metadata_ok &&
           left.overall_ok == right.overall_ok
end
```

- [ ] **Step 2: Add replay and verifier**

Append:

```julia
function replay_deterministic_quillen_patch(patch::QuillenGlobalPatchAssembly)
    R = _require_supported_quillen_ring(patch.ring)
    n = patch.size
    n >= 2 || throw(ArgumentError("patch size must be at least 2"))
    selected = _require_substitution_generator(R, patch.substitution_variable)
    target = _quillen_global_target_matrix(
        patch.target;
        ring = R,
        size = n,
        label = "target",
    )

    cover = patch.cover_certificate
    local_certificates = patch.local_certificates
    normalized = patch.normalized_local_contributions
    cover_certificate_ok = verify_quillen_denominator_cover(cover)
    local_certificates_ok = all(verify_quillen_local_certificate, local_certificates)
    normalized_contributions_ok =
        all(verify_quillen_local_contribution_normalization, normalized)
    count_ok = length(local_certificates) == length(normalized)

    local_alignment_ok = count_ok && all(eachindex(local_certificates)) do idx
        _same_quillen_local_certificate_data(
            normalized[idx].local_certificate,
            local_certificates[idx],
        )
    end
    cover_alignment_ok =
        _same_quillen_cover_certificate_data(cover, patch.cover_certificate) &&
        all(normalized) do item
            _same_quillen_cover_certificate_data(item.cover_certificate, cover)
        end
    normalized_input_ok = all(normalized) do item
        item.original_input == patch.original_input &&
            item.local_certificate.original_input == patch.original_input
    end
    selected_variable_ok =
        all(normalized) do item
            item.selected_variable == selected &&
                item.local_certificate.selected_variable == selected
        end

    expected_denominator_data = _quillen_denominator_data([
        item.local_contribution for item in normalized
    ])
    denominator_data_ok =
        _same_quillen_denominator_data(patch.denominator_data, expected_denominator_data)
    coverage_sum = _quillen_coverage_sum(R, patch.denominator_data)
    coverage_ok =
        denominator_data_ok && cover_certificate_ok && coverage_sum == one(R) &&
        coverage_sum == cover.coverage_sum
    global_elementary_factors =
        reduce(vcat, [item.weighted_global_elementary_factors for item in normalized]; init = Any[])
    global_elementary_factors_ok = _same_quillen_factors(
        patch.global_elementary_factors,
        global_elementary_factors,
    )
    product = _quillen_product(R, n, patch.global_elementary_factors)
    product_ok = global_elementary_factors_ok && patch.patched_product == product
    target_ok = product == target
    replay_metadata =
        _quillen_global_patch_replay_metadata(cover, local_certificates, normalized)
    replay_metadata_ok = patch.replay_metadata == replay_metadata
    overall_ok =
        cover_certificate_ok &&
        local_certificates_ok &&
        normalized_contributions_ok &&
        local_alignment_ok &&
        cover_alignment_ok &&
        normalized_input_ok &&
        selected_variable_ok &&
        coverage_ok &&
        global_elementary_factors_ok &&
        product_ok &&
        target_ok &&
        replay_metadata_ok

    return QuillenGlobalPatchAssemblyVerification(
        cover_certificate_ok,
        local_certificates_ok,
        normalized_contributions_ok,
        local_alignment_ok,
        cover_alignment_ok,
        normalized_input_ok,
        selected_variable_ok,
        coverage_sum,
        coverage_ok,
        global_elementary_factors,
        global_elementary_factors_ok,
        product,
        product_ok,
        target,
        target_ok,
        replay_metadata,
        replay_metadata_ok,
        overall_ok,
    )
end

function verify_quillen_patch(patch::QuillenGlobalPatchAssembly)::Bool
    try
        replay = replay_deterministic_quillen_patch(patch)
        return replay.overall_ok &&
               _same_quillen_global_patch_verification(patch.verification, replay)
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end
```

- [ ] **Step 3: Add target coercion and constructor**

Append this helper before `replay_deterministic_quillen_patch` or above the
constructor:

```julia
function _quillen_global_target_matrix(target; ring, size::Int, label::AbstractString)
    if target isa QuillenElementaryCorrection
        correction = _normalize_quillen_contribution(
            QuillenLocalContribution(
                LocalCertificate([target.row, target.col], [one(ring), one(ring)]),
                one(ring),
                one(ring),
                target,
            ),
            ring,
            size,
        ).correction
        return elementary_matrix(size, correction.row, correction.col, correction.entry, ring)
    end
    _quillen_local_require_factor_matrix(target, ring, size, label)
    return target
end

function assemble_deterministic_quillen_patch(
    original_input,
    selected_variable,
    local_certificates,
    normalized_local_contributions,
    cover::QuillenDenominatorCoverCertificate;
    target = original_input,
    ring = nothing,
    size = nothing,
)
    certificates = collect(local_certificates)
    normalized = collect(normalized_local_contributions)
    isempty(certificates) &&
        throw(ArgumentError("deterministic Quillen patch assembly requires at least one local certificate"))
    length(certificates) == length(normalized) ||
        throw(ArgumentError("deterministic Quillen patch assembly requires matching local certificate and normalized contribution counts"))
    verify_quillen_denominator_cover(cover) ||
        throw(ArgumentError("Quillen denominator cover certificate does not replay"))
    all(verify_quillen_local_certificate, certificates) ||
        throw(ArgumentError("Quillen local realization certificate does not replay"))
    all(verify_quillen_local_contribution_normalization, normalized) ||
        throw(ArgumentError("Quillen local contribution normalization does not replay"))

    R, n = _quillen_local_input_ring_size(original_input; ring, size)
    cover.ring == R ||
        throw(ArgumentError("deterministic Quillen patch assembly requires cover and target rings to match"))
    selected = _require_substitution_generator(R, selected_variable)
    target_matrix = _quillen_global_target_matrix(target; ring = R, size = n, label = "target")

    all(certificate -> certificate.ring == R && certificate.size == n, certificates) ||
        throw(ArgumentError("deterministic Quillen patch assembly requires local certificate rings and sizes to match the target"))
    all(certificate -> certificate.selected_variable == selected, certificates) ||
        throw(ArgumentError("deterministic Quillen patch assembly requires local certificate variables to match"))
    all(certificate -> certificate.original_input == original_input, certificates) ||
        throw(ArgumentError("deterministic Quillen patch assembly requires local certificate original inputs to match"))

    for idx in eachindex(certificates)
        _same_quillen_local_certificate_data(normalized[idx].local_certificate, certificates[idx]) ||
            throw(ArgumentError("deterministic Quillen patch assembly requires normalized contribution local certificates to match input order"))
        _same_quillen_cover_certificate_data(normalized[idx].cover_certificate, cover) ||
            throw(ArgumentError("deterministic Quillen patch assembly requires normalized contribution cover certificates to match"))
        normalized[idx].original_input == original_input ||
            throw(ArgumentError("deterministic Quillen patch assembly requires normalized contribution original inputs to match"))
        normalized[idx].selected_variable == selected ||
            throw(ArgumentError("deterministic Quillen patch assembly requires normalized contribution selected variables to match"))
    end

    denominator_data = _quillen_denominator_data([
        item.local_contribution for item in normalized
    ])
    coverage_sum = _quillen_coverage_sum(R, denominator_data)
    coverage_sum == one(R) ||
        throw(ArgumentError("deterministic Quillen patch assembly denominator coverage must sum to one"))

    factor_type = typeof(identity_matrix(R, n))
    global_elementary_factors = factor_type[]
    for item in normalized
        append!(global_elementary_factors, item.weighted_global_elementary_factors)
    end
    patched_product = _quillen_product(R, n, global_elementary_factors)
    patched_product == target_matrix ||
        throw(ArgumentError("deterministic Quillen patch assembly product does not equal the target"))

    replay_metadata =
        _quillen_global_patch_replay_metadata(cover, certificates, normalized)
    provisional = QuillenGlobalPatchAssembly(
        R,
        n,
        selected,
        original_input,
        cover,
        denominator_data,
        certificates,
        normalized,
        global_elementary_factors,
        patched_product,
        target_matrix,
        replay_metadata,
        QuillenGlobalPatchAssemblyVerification(
            false,
            false,
            false,
            false,
            false,
            false,
            false,
            zero(R),
            false,
            Any[],
            false,
            identity_matrix(R, n),
            false,
            target_matrix,
            false,
            replay_metadata,
            false,
            false,
        ),
    )
    verification = replay_deterministic_quillen_patch(provisional)
    verification.overall_ok ||
        throw(ArgumentError("deterministic Quillen patch assembly data does not replay"))
    return QuillenGlobalPatchAssembly(
        provisional.ring,
        provisional.size,
        provisional.substitution_variable,
        provisional.original_input,
        provisional.cover_certificate,
        provisional.denominator_data,
        provisional.local_certificates,
        provisional.normalized_local_contributions,
        provisional.global_elementary_factors,
        provisional.patched_product,
        provisional.target,
        provisional.replay_metadata,
        verification,
    )
end
```

- [ ] **Step 4: Run GREEN verification**

Run:

```bash
julia --project=. -e 'include("test/expert/quillen_global_patch_assembly.jl")'
```

Expected: PASS.

- [ ] **Step 5: Run expert group verification**

Run:

```bash
julia --project=. test/runtests.jl expert
```

Expected: PASS.

- [ ] **Step 6: Commit implementation**

Commit source changes:

```bash
git add src/algorithm/quillen_induction.jl
git commit -m "feat: assemble deterministic quillen global patches"
```
