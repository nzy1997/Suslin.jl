# Issue 218 Quillen Supplied Evidence Patch Assembly Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an internal supplied-evidence constructor that assembles a verified global Quillen patch from checked local factor sequence evidence.

**Architecture:** Add a new supplied-evidence assembly certificate in `src/algorithm/quillen_induction.jl` without changing public exports. The constructor verifies #214 local sequences, extracts a #215 candidate, solves a #216 cover, builds or validates a #217 substitution chain, expands each local sequence factor under its solved cover term, records explicit base-term handling, and verifies exact multiplication to the input matrix.

**Tech Stack:** Julia, Oscar exact ordinary polynomial rings and matrices, existing Suslin Quillen sequence/candidate/solver/chain helpers, `Test`.

## Global Constraints

- Input is an ordinary-polynomial matrix `A`, selected variable `X`, nonempty verified `QuillenLocalFactorSequenceCertificate` evidence for the same input, and either supplied base-term factor evidence or explicit `base_term_policy = :trivial` / `:already_handled`.
- Output is a `QuillenSuppliedEvidencePatchAssembly` whose `global_elementary_factors` multiply exactly to `A`.
- The constructor must build the denominator-cover candidate, solver result, cover certificate, substitution chain, sequence expansions, base-term boundary record, product, target, and replay metadata from verified inputs.
- Every generated global elementary factor must come from a replayed local sequence factor and the solved cover term `g_i * r_i^l`; do not collapse a sequence into one opaque correction.
- The verifier must reject incomplete local evidence, corrupted local sequence factors, unproven cover multipliers, non-replaying substitution chains, missing base-term evidence, and tampered assembly records.
- Keep the constructor expert/internal; do not export new names from `src/Suslin.jl`.
- Do not implement Murthy local solving, the general `SL_3` driver, ECP, route-boundary wiring, public APIs, or a base-term factorization algorithm.
- Focused supplied-evidence command is `julia --project=. -e 'include("test/expert/quillen_supplied_evidence_patch_assembly.jl")'`.
- Constructive regression command is `julia --project=. -e 'include("test/expert/quillen_induction_constructive.jl")'`.
- Full package command is `julia --project=. -e 'using Pkg; Pkg.test()'`.

---

## File Structure

- Modify `src/algorithm/quillen_induction.jl`: add sequence expansion records, supplied assembly records, replay/verify helpers, base-term handling, and `assemble_quillen_patch_from_local_evidence`.
- Create `test/expert/quillen_supplied_evidence_patch_assembly.jl`: positive supplied-evidence assembly plus negative controls.
- Modify `test/expert/quillen_induction_constructive.jl`: add one supplied-evidence regression in the existing constructive test.
- Modify `test/runtests.jl`: register the new expert test after `expert/quillen_patch_substitution_chain.jl`.

### Task 1: Add Red Supplied-Evidence Tests

**Files:**
- Create: `test/expert/quillen_supplied_evidence_patch_assembly.jl`
- Modify: `test/runtests.jl`

**Interfaces:**
- Consumes future `Suslin.QuillenSuppliedEvidencePatchAssembly`, `Suslin.assemble_quillen_patch_from_local_evidence`, `Suslin.replay_quillen_supplied_evidence_patch`, and `Suslin.verify_quillen_patch(::QuillenSuppliedEvidencePatchAssembly)`.
- Produces failing coverage for supplied local factor sequence assembly and required negative controls.

- [ ] **Step 1: Write the failing expert test**

Create `test/expert/quillen_supplied_evidence_patch_assembly.jl`:

```julia
using Test
using Suslin
using Oscar

const QUILLEN_SUPPLIED_EVIDENCE_CATALOG_PATH =
    joinpath(@__DIR__, "..", "fixtures", "quillen_patch_cases.jl")
if !isdefined(Main, :QuillenPatchFixtureCatalog)
    include(QUILLEN_SUPPLIED_EVIDENCE_CATALOG_PATH)
end

function qse_product(factors, R, n::Int)
    product = identity_matrix(R, n)
    for factor in factors
        product *= factor
    end
    return product
end

function qse_local_certificate(local_factor)
    return Suslin.LocalCertificate(
        local_factor.certificate.indices,
        local_factor.certificate.denominators,
    )
end

function qse_correction(local_factor)
    correction = local_factor.correction
    return Suslin.QuillenElementaryCorrection(
        correction.row,
        correction.col,
        correction.entry,
    )
end

function qse_sequence_certificate(entry, index::Int)
    local_factor = entry.local_factors[index]
    realization = Suslin.quillen_local_realization_certificate(
        entry.target_matrix,
        entry.substitution_variable;
        local_certificate = qse_local_certificate(local_factor),
        denominator = local_factor.denominator,
        coverage_multiplier = local_factor.coverage_multiplier,
        correction = qse_correction(local_factor),
        factors = [local_factor.factor],
        patched_substitution_witness = entry.patched_substitution_witness,
        ring = entry.ring.object,
        size = entry.size,
    )
    return Suslin.quillen_local_factor_sequence_certificate(
        realization;
        factor_provenance = (;
            factor_index = 1,
            sequence_index = 1,
            local_index = 1,
            fixture_id = entry.id,
            source = :supplied_evidence_patch_test,
        ),
        metadata = (;
            source_refs = entry.source_refs,
            consumer_issue_ids = entry.consumer_issue_ids,
        ),
    )
end

function qse_sequence_certificates(entry)
    return [
        qse_sequence_certificate(entry, index)
        for index in eachindex(entry.local_factors)
    ]
end

function qse_rebuild(record; kwargs...)
    overrides = Dict{Symbol,Any}()
    for pair in kwargs
        overrides[pair.first] = pair.second
    end
    values = map(fieldnames(typeof(record))) do name
        get(overrides, name, getproperty(record, name))
    end
    return typeof(record)(values...)
end

function qse_rebuild_factor(factor; kwargs...)
    fields = merge((
        row = factor.row,
        col = factor.col,
        numerator = factor.numerator,
        denominator = factor.denominator,
        coverage_multiplier = factor.coverage_multiplier,
        local_certificate = factor.local_certificate,
        provenance = factor.provenance,
        metadata = factor.metadata,
    ), NamedTuple(kwargs))
    return Suslin.QuillenLocalElementaryFactor(
        fields.row,
        fields.col,
        fields.numerator,
        fields.denominator,
        fields.coverage_multiplier,
        fields.local_certificate,
        fields.provenance,
        fields.metadata,
    )
end

function qse_rebuild_sequence_certificate(cert; kwargs...)
    fields = merge((
        original_input = cert.original_input,
        ring = cert.ring,
        size = cert.size,
        selected_variable = cert.selected_variable,
        factors = cert.factors,
        raw_denominators = cert.raw_denominators,
        product_denominator = cert.product_denominator,
        local_product = cert.local_product,
        local_correction = cert.local_correction,
        normalized_local_contributions = cert.normalized_local_contributions,
        normalized_global_elementary_factors = cert.normalized_global_elementary_factors,
        patched_substitution_witness = cert.patched_substitution_witness,
        chain_witness = cert.chain_witness,
        witness_metadata = cert.witness_metadata,
        replay_metadata = cert.replay_metadata,
        verification = cert.verification,
    ), NamedTuple(kwargs))
    return Suslin.QuillenLocalFactorSequenceCertificate(
        fields.original_input,
        fields.ring,
        fields.size,
        fields.selected_variable,
        fields.factors,
        fields.raw_denominators,
        fields.product_denominator,
        fields.local_product,
        fields.local_correction,
        fields.normalized_local_contributions,
        fields.normalized_global_elementary_factors,
        fields.patched_substitution_witness,
        fields.chain_witness,
        fields.witness_metadata,
        fields.replay_metadata,
        fields.verification,
    )
end

function qse_rebuild_chain(chain; kwargs...)
    fields = merge((
        original_matrix = chain.original_matrix,
        ring = chain.ring,
        size = chain.size,
        selected_variable = chain.selected_variable,
        sign_convention = chain.sign_convention,
        solver_result = chain.solver_result,
        cumulative_coefficients = chain.cumulative_coefficients,
        intermediate_matrices = chain.intermediate_matrices,
        steps = chain.steps,
        bracket_matrices = chain.bracket_matrices,
        base_term = chain.base_term,
        metadata = chain.metadata,
        replay_metadata = chain.replay_metadata,
        verification = chain.verification,
    ), NamedTuple(kwargs))
    return Suslin.QuillenPatchSubstitutionChain(
        fields.original_matrix,
        fields.ring,
        fields.size,
        fields.selected_variable,
        fields.sign_convention,
        fields.solver_result,
        fields.cumulative_coefficients,
        fields.intermediate_matrices,
        fields.steps,
        fields.bracket_matrices,
        fields.base_term,
        fields.metadata,
        fields.replay_metadata,
        fields.verification,
    )
end

function qse_expected_sequence_factors(patch)
    factor_type = typeof(identity_matrix(patch.ring, patch.size))
    factors = factor_type[]
    for expansion in patch.sequence_expansions
        append!(factors, expansion.global_elementary_factors)
    end
    return factors
end

@testset "Quillen supplied local evidence patch assembly" begin
    entries = Main.QuillenPatchFixtureCatalog.cases_by_id()
    entry = entries["quillen-two-open-cover-qq"]
    certificates = qse_sequence_certificates(entry)

    patch = Suslin.assemble_quillen_patch_from_local_evidence(
        entry.target_matrix,
        entry.substitution_variable,
        certificates;
        max_exponent = 2,
        base_term_policy = :already_handled,
        metadata = (; fixture_id = entry.id, consumer_issue_id = "#218"),
    )

    @test patch isa Suslin.QuillenSuppliedEvidencePatchAssembly
    @test all(Suslin.verify_quillen_local_factor_sequence_certificate, certificates)
    @test Suslin.verify_quillen_denominator_cover_candidate(patch.denominator_candidate)
    @test Suslin.verify_quillen_denominator_cover_solver_result(patch.solver_result)
    @test Suslin.verify_quillen_denominator_cover(patch.cover_certificate)
    @test Suslin.verify_quillen_patch_substitution_chain(patch.substitution_chain)
    @test Suslin.verify_quillen_patch(patch)

    replay = Suslin.replay_quillen_supplied_evidence_patch(patch)
    @test replay.overall_ok
    @test replay.local_certificates_ok
    @test replay.denominator_candidate_ok
    @test replay.denominator_candidate_matches
    @test replay.solver_result_ok
    @test replay.solver_source_candidate_ok
    @test replay.cover_certificate_ok
    @test replay.substitution_chain_ok
    @test replay.substitution_chain_matches
    @test replay.base_term_ok
    @test replay.sequence_expansions_ok
    @test replay.global_elementary_factors_ok
    @test replay.product_ok
    @test replay.target_ok
    @test replay.replay_metadata_ok

    @test patch.local_certificates == certificates
    @test patch.denominator_candidate.raw_denominators ==
          [certificate.product_denominator for certificate in certificates]
    @test patch.solver_result.coverage_sum == one(patch.ring)
    @test patch.cover_certificate.denominators == patch.solver_result.powered_denominators
    @test patch.substitution_chain.original_matrix == entry.target_matrix
    @test patch.substitution_chain.verification.telescope_ok
    @test patch.base_term_policy == :already_handled
    @test isempty(patch.base_term_factors)
    @test patch.base_term == patch.substitution_chain.base_term
    @test patch.replay_metadata.metadata == (; fixture_id = entry.id, consumer_issue_id = "#218")

    expected_sequence_factors = qse_expected_sequence_factors(patch)
    @test patch.sequence_elementary_factors == expected_sequence_factors
    @test patch.global_elementary_factors == expected_sequence_factors
    @test qse_product(patch.global_elementary_factors, patch.ring, patch.size) ==
          entry.target_matrix
    @test patch.product == entry.target_matrix
    @test patch.target == entry.target_matrix

    for (local_index, expansion) in enumerate(patch.sequence_expansions)
        certificate = certificates[local_index]
        @test Suslin.verify_quillen_sequence_contribution_expansion(expansion)
        @test expansion.local_certificate == certificate
        @test expansion.local_index == local_index
        @test expansion.coverage_multiplier == patch.solver_result.coverage_multipliers[local_index]
        @test expansion.powered_denominator == patch.solver_result.powered_denominators[local_index]
        @test expansion.cover_term == patch.solver_result.coverage_terms[local_index]
        @test length(expansion.global_elementary_factors) == length(certificate.factors)
        for (factor_index, factor) in enumerate(certificate.factors)
            expected = elementary_matrix(
                patch.size,
                factor.row,
                factor.col,
                expansion.cover_term * factor.numerator,
                patch.ring,
            )
            @test expansion.global_elementary_factors[factor_index] == expected
        end
    end

    @test_throws ArgumentError Suslin.assemble_quillen_patch_from_local_evidence(
        entry.target_matrix,
        entry.substitution_variable,
        certificates[1:1];
        max_exponent = 2,
        base_term_policy = :already_handled,
    )

    bad_certificates = copy(certificates)
    bad_factors = copy(bad_certificates[1].factors)
    bad_factors[1] = qse_rebuild_factor(
        bad_factors[1];
        numerator = bad_factors[1].numerator + one(entry.ring.object),
    )
    bad_certificates[1] = qse_rebuild_sequence_certificate(
        bad_certificates[1];
        factors = bad_factors,
    )
    @test !Suslin.verify_quillen_local_factor_sequence_certificate(bad_certificates[1])
    @test_throws ArgumentError Suslin.assemble_quillen_patch_from_local_evidence(
        entry.target_matrix,
        entry.substitution_variable,
        bad_certificates;
        max_exponent = 2,
        base_term_policy = :already_handled,
    )

    @test_throws ArgumentError Suslin.assemble_quillen_patch_from_local_evidence(
        entry.target_matrix,
        entry.substitution_variable,
        certificates;
        max_exponent = 1,
        coverage_multipliers = [
            one(entry.ring.object) + entry.ring.generators[2],
            one(entry.ring.object),
        ],
        base_term_policy = :already_handled,
    )

    tampered_chain = qse_rebuild_chain(
        patch.substitution_chain;
        sign_convention = :park_woodburn_plus,
    )
    @test !Suslin.verify_quillen_patch_substitution_chain(tampered_chain)
    @test_throws ArgumentError Suslin.assemble_quillen_patch_from_local_evidence(
        entry.target_matrix,
        entry.substitution_variable,
        certificates;
        max_exponent = 2,
        base_term_policy = :already_handled,
        substitution_chain = tampered_chain,
    )

    @test_throws ArgumentError Suslin.assemble_quillen_patch_from_local_evidence(
        entry.target_matrix,
        entry.substitution_variable,
        certificates;
        max_exponent = 2,
    )

    tampered_patch = qse_rebuild(
        patch;
        product = patch.product * elementary_matrix(patch.size, 1, 2, one(patch.ring), patch.ring),
    )
    @test !Suslin.verify_quillen_patch(tampered_patch)
end
```

- [ ] **Step 2: Register the red expert test**

In `test/runtests.jl`, add:

```julia
        "expert/quillen_supplied_evidence_patch_assembly.jl",
```

immediately after:

```julia
        "expert/quillen_patch_substitution_chain.jl",
```

- [ ] **Step 3: Run the focused red test**

Run:

```bash
julia --project=. -e 'include("test/expert/quillen_supplied_evidence_patch_assembly.jl")'
```

Expected: FAIL with `UndefVarError` for `assemble_quillen_patch_from_local_evidence` or `QuillenSuppliedEvidencePatchAssembly`.

- [ ] **Step 4: Commit the red tests**

```bash
git add test/expert/quillen_supplied_evidence_patch_assembly.jl test/runtests.jl
git commit -m "test: cover supplied evidence quillen assembly"
```

### Task 2: Add Sequence Expansion Replay

**Files:**
- Modify: `src/algorithm/quillen_induction.jl`

**Interfaces:**
- Consumes existing `QuillenLocalFactorSequenceCertificate` and `QuillenDenominatorCoverSolverResult`.
- Produces `QuillenSequenceContributionExpansion`, `replay_quillen_sequence_contribution_expansion`, and `verify_quillen_sequence_contribution_expansion`.

- [ ] **Step 1: Add sequence expansion records**

Insert after `QuillenLocalContributionNormalization`:

```julia
struct QuillenSequenceContributionExpansionVerification
    local_certificate_ok::Bool
    solver_result_ok::Bool
    local_index_ok::Bool
    solver_context_ok::Bool
    powered_denominator
    coverage_multiplier
    cover_term
    cover_term_ok::Bool
    factor_provenance::Vector
    global_elementary_factors::Vector
    global_elementary_factors_ok::Bool
    replay_metadata
    replay_metadata_ok::Bool
    overall_ok::Bool
end

struct QuillenSequenceContributionExpansion
    local_certificate::QuillenLocalFactorSequenceCertificate
    solver_result::QuillenDenominatorCoverSolverResult
    local_index::Int
    powered_denominator
    coverage_multiplier
    cover_term
    factor_provenance::Vector
    global_elementary_factors::Vector
    replay_metadata
    verification::QuillenSequenceContributionExpansionVerification
end
```

- [ ] **Step 2: Add expansion helpers**

Add helpers before `QuillenGlobalPatchAssemblyVerification`:

```julia
function _quillen_sequence_expansion_metadata(
    certificate::QuillenLocalFactorSequenceCertificate,
    solver_result::QuillenDenominatorCoverSolverResult,
    local_index::Int,
    factor_provenance,
)
    return (;
        source = :quillen_supplied_local_sequence_expansion,
        local_index = local_index,
        factor_count = length(certificate.factors),
        raw_denominator = solver_result.raw_denominators[local_index],
        powered_denominator = solver_result.powered_denominators[local_index],
        coverage_multiplier = solver_result.coverage_multipliers[local_index],
        cover_term = solver_result.coverage_terms[local_index],
        factor_provenance = factor_provenance,
        local_replay_metadata = certificate.replay_metadata,
    )
end

function _quillen_sequence_expansion_factors(
    certificate::QuillenLocalFactorSequenceCertificate,
    cover_term,
)
    R = _require_supported_quillen_ring(certificate.ring)
    n = certificate.size
    factor_type = typeof(identity_matrix(R, n))
    factors = factor_type[]
    for factor in certificate.factors
        push!(
            factors,
            elementary_matrix(
                n,
                factor.row,
                factor.col,
                _coerce_into_ring(R, cover_term * factor.numerator, "sequence expansion entry"),
                R,
            ),
        )
    end
    return factors
end
```

- [ ] **Step 3: Add replay, constructor, and verifier**

Add:

```julia
function replay_quillen_sequence_contribution_expansion(
    expansion::QuillenSequenceContributionExpansion,
)
    certificate = expansion.local_certificate
    solver_result = expansion.solver_result
    R = _require_supported_quillen_ring(certificate.ring)
    n = certificate.size
    local_certificate_ok = verify_quillen_local_factor_sequence_certificate(certificate)
    solver_result_ok = verify_quillen_denominator_cover_solver_result(solver_result)
    local_index_ok = 1 <= expansion.local_index <= length(solver_result.raw_denominators)
    solver_context_ok =
        solver_result_ok &&
        solver_result.ring == R &&
        local_index_ok &&
        solver_result.raw_denominators[expansion.local_index] == certificate.product_denominator
    powered_denominator = local_index_ok ? solver_result.powered_denominators[expansion.local_index] : zero(R)
    coverage_multiplier = local_index_ok ? solver_result.coverage_multipliers[expansion.local_index] : zero(R)
    cover_term = local_index_ok ? solver_result.coverage_terms[expansion.local_index] : zero(R)
    cover_term_ok =
        expansion.powered_denominator == powered_denominator &&
        expansion.coverage_multiplier == coverage_multiplier &&
        expansion.cover_term == cover_term &&
        cover_term == coverage_multiplier * powered_denominator
    factor_provenance = _quillen_local_sequence_factor_provenance(certificate.factors)
    global_elementary_factors = _quillen_sequence_expansion_factors(certificate, cover_term)
    global_elementary_factors_ok =
        _same_quillen_factors(expansion.global_elementary_factors, global_elementary_factors)
    replay_metadata = _quillen_sequence_expansion_metadata(
        certificate,
        solver_result,
        expansion.local_index,
        factor_provenance,
    )
    replay_metadata_ok = expansion.replay_metadata == replay_metadata
    overall_ok =
        local_certificate_ok &&
        solver_result_ok &&
        local_index_ok &&
        solver_context_ok &&
        cover_term_ok &&
        global_elementary_factors_ok &&
        replay_metadata_ok
    return QuillenSequenceContributionExpansionVerification(
        local_certificate_ok,
        solver_result_ok,
        local_index_ok,
        solver_context_ok,
        powered_denominator,
        coverage_multiplier,
        cover_term,
        cover_term_ok,
        factor_provenance,
        global_elementary_factors,
        global_elementary_factors_ok,
        replay_metadata,
        replay_metadata_ok,
        overall_ok,
    )
end

function _same_quillen_sequence_expansion_verification(
    left::QuillenSequenceContributionExpansionVerification,
    right::QuillenSequenceContributionExpansionVerification,
)::Bool
    return left.local_certificate_ok == right.local_certificate_ok &&
           left.solver_result_ok == right.solver_result_ok &&
           left.local_index_ok == right.local_index_ok &&
           left.solver_context_ok == right.solver_context_ok &&
           left.powered_denominator == right.powered_denominator &&
           left.coverage_multiplier == right.coverage_multiplier &&
           left.cover_term == right.cover_term &&
           left.cover_term_ok == right.cover_term_ok &&
           left.factor_provenance == right.factor_provenance &&
           _same_quillen_factors(left.global_elementary_factors, right.global_elementary_factors) &&
           left.global_elementary_factors_ok == right.global_elementary_factors_ok &&
           left.replay_metadata == right.replay_metadata &&
           left.replay_metadata_ok == right.replay_metadata_ok &&
           left.overall_ok == right.overall_ok
end

function quillen_sequence_contribution_expansion(
    certificate::QuillenLocalFactorSequenceCertificate,
    solver_result::QuillenDenominatorCoverSolverResult,
    local_index::Int,
)
    verify_quillen_local_factor_sequence_certificate(certificate) ||
        throw(ArgumentError("Quillen supplied evidence assembly requires verified local sequence certificates"))
    verify_quillen_denominator_cover_solver_result(solver_result) ||
        throw(ArgumentError("Quillen supplied evidence assembly requires a verified denominator-cover solver result"))
    1 <= local_index <= length(solver_result.raw_denominators) ||
        throw(ArgumentError("Quillen supplied evidence assembly local index is outside the solver result"))
    solver_result.ring == certificate.ring &&
        solver_result.raw_denominators[local_index] == certificate.product_denominator ||
        throw(ArgumentError("Quillen supplied evidence assembly solver denominator must match local sequence provenance"))
    factor_provenance = _quillen_local_sequence_factor_provenance(certificate.factors)
    powered_denominator = solver_result.powered_denominators[local_index]
    coverage_multiplier = solver_result.coverage_multipliers[local_index]
    cover_term = solver_result.coverage_terms[local_index]
    global_elementary_factors = _quillen_sequence_expansion_factors(certificate, cover_term)
    replay_metadata = _quillen_sequence_expansion_metadata(
        certificate,
        solver_result,
        local_index,
        factor_provenance,
    )
    provisional = QuillenSequenceContributionExpansion(
        certificate,
        solver_result,
        Int(local_index),
        powered_denominator,
        coverage_multiplier,
        cover_term,
        factor_provenance,
        global_elementary_factors,
        replay_metadata,
        QuillenSequenceContributionExpansionVerification(
            false,
            false,
            false,
            false,
            powered_denominator,
            coverage_multiplier,
            cover_term,
            false,
            factor_provenance,
            typeof(identity_matrix(certificate.ring, certificate.size))[],
            false,
            replay_metadata,
            false,
            false,
        ),
    )
    verification = replay_quillen_sequence_contribution_expansion(provisional)
    verification.overall_ok ||
        throw(ArgumentError("Quillen supplied evidence sequence expansion does not replay"))
    return QuillenSequenceContributionExpansion(
        provisional.local_certificate,
        provisional.solver_result,
        provisional.local_index,
        provisional.powered_denominator,
        provisional.coverage_multiplier,
        provisional.cover_term,
        provisional.factor_provenance,
        provisional.global_elementary_factors,
        provisional.replay_metadata,
        verification,
    )
end

function verify_quillen_sequence_contribution_expansion(expansion)::Bool
    try
        replay = replay_quillen_sequence_contribution_expansion(expansion)
        return replay.overall_ok &&
               _same_quillen_sequence_expansion_verification(expansion.verification, replay)
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end
```

- [ ] **Step 4: Run red-to-green focused test**

Run:

```bash
julia --project=. -e 'include("test/expert/quillen_supplied_evidence_patch_assembly.jl")'
```

Expected: still FAIL, now on missing supplied assembly type or constructor.

- [ ] **Step 5: Commit sequence expansion**

```bash
git add src/algorithm/quillen_induction.jl
git commit -m "feat: replay quillen sequence expansions"
```

### Task 3: Add Supplied Evidence Assembly

**Files:**
- Modify: `src/algorithm/quillen_induction.jl`

**Interfaces:**
- Consumes Task 2 sequence expansion helpers plus existing #215/#216/#217 records.
- Produces `QuillenSuppliedEvidencePatchAssembly`, `assemble_quillen_patch_from_local_evidence`, `replay_quillen_supplied_evidence_patch`, and `verify_quillen_patch(::QuillenSuppliedEvidencePatchAssembly)`.

- [ ] **Step 1: Add supplied assembly records**

Insert after sequence expansion helpers:

```julia
struct QuillenSuppliedEvidencePatchAssemblyVerification
    local_certificates_ok::Bool
    denominator_candidate_ok::Bool
    denominator_candidate_matches::Bool
    solver_result_ok::Bool
    solver_source_candidate_ok::Bool
    cover_certificate_ok::Bool
    substitution_chain_ok::Bool
    substitution_chain_matches::Bool
    base_term_ok::Bool
    sequence_expansions_ok::Bool
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

struct QuillenSuppliedEvidencePatchAssembly
    ring
    size::Int
    substitution_variable
    original_input
    local_certificates::Vector{QuillenLocalFactorSequenceCertificate}
    denominator_candidate::QuillenDenominatorCoverCandidate
    solver_result::QuillenDenominatorCoverSolverResult
    cover_certificate::QuillenDenominatorCoverCertificate
    substitution_chain::QuillenPatchSubstitutionChain
    base_term_policy::Symbol
    base_term
    base_term_factors::Vector
    base_term_product
    sequence_expansions::Vector{QuillenSequenceContributionExpansion}
    sequence_elementary_factors::Vector
    global_elementary_factors::Vector
    product
    target
    replay_metadata
    verification::QuillenSuppliedEvidencePatchAssemblyVerification
end
```

- [ ] **Step 2: Add metadata and comparison helpers**

Add:

```julia
function _quillen_supplied_patch_metadata(
    candidate::QuillenDenominatorCoverCandidate,
    solver_result::QuillenDenominatorCoverSolverResult,
    substitution_chain::QuillenPatchSubstitutionChain,
    base_term_policy::Symbol,
    sequence_expansions,
    metadata,
)
    return (;
        source = :quillen_supplied_local_evidence_patch_assembly,
        local_count = length(candidate.local_certificates),
        raw_denominators = candidate.raw_denominators,
        exponent = solver_result.exponent,
        powered_denominators = solver_result.powered_denominators,
        coverage_multipliers = solver_result.coverage_multipliers,
        coverage_sum = solver_result.coverage_sum,
        substitution_chain_replay_metadata = substitution_chain.replay_metadata,
        base_term_policy = base_term_policy,
        sequence_expansion_metadata = [expansion.replay_metadata for expansion in sequence_expansions],
        metadata = metadata,
    )
end

function _same_quillen_sequence_expansions(left, right)::Bool
    length(left) == length(right) || return false
    for idx in eachindex(left)
        left[idx].local_index == right[idx].local_index || return false
        left[idx].cover_term == right[idx].cover_term || return false
        _same_quillen_factors(left[idx].global_elementary_factors, right[idx].global_elementary_factors) || return false
        left[idx].replay_metadata == right[idx].replay_metadata || return false
        _same_quillen_sequence_expansion_verification(left[idx].verification, right[idx].verification) || return false
    end
    return true
end

function _same_quillen_supplied_patch_verification(
    left::QuillenSuppliedEvidencePatchAssemblyVerification,
    right::QuillenSuppliedEvidencePatchAssemblyVerification,
)::Bool
    return left.local_certificates_ok == right.local_certificates_ok &&
           left.denominator_candidate_ok == right.denominator_candidate_ok &&
           left.denominator_candidate_matches == right.denominator_candidate_matches &&
           left.solver_result_ok == right.solver_result_ok &&
           left.solver_source_candidate_ok == right.solver_source_candidate_ok &&
           left.cover_certificate_ok == right.cover_certificate_ok &&
           left.substitution_chain_ok == right.substitution_chain_ok &&
           left.substitution_chain_matches == right.substitution_chain_matches &&
           left.base_term_ok == right.base_term_ok &&
           left.sequence_expansions_ok == right.sequence_expansions_ok &&
           _same_quillen_factors(left.global_elementary_factors, right.global_elementary_factors) &&
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

- [ ] **Step 3: Add base-term helpers**

Add:

```julia
function _quillen_supplied_base_term_policy(base_term_policy, base_term_factors)
    if base_term_policy === nothing
        base_term_factors === nothing &&
            throw(ArgumentError("Quillen supplied evidence patch assembly requires supplied A(0) factors or base_term_policy = :trivial or :already_handled"))
        return :supplied
    end
    policy = Symbol(base_term_policy)
    policy in (:supplied, :trivial, :already_handled) ||
        throw(ArgumentError("Quillen supplied evidence patch assembly has unsupported base-term policy"))
    return policy
end

function _quillen_supplied_base_term_factors(R, n::Int, base_term_factors)
    factor_type = typeof(identity_matrix(R, n))
    base_term_factors === nothing && return factor_type[]
    return [
        _quillen_local_require_factor_matrix(factor, R, n, "base-term factor")
        for factor in collect(base_term_factors)
    ]
end

function _quillen_supplied_base_term_ok(policy::Symbol, base_term, factors, product, R, n::Int)
    policy == :supplied && return product == base_term
    policy == :trivial && return isempty(factors) && base_term == identity_matrix(R, n)
    policy == :already_handled && return isempty(factors)
    return false
end
```

- [ ] **Step 4: Add replay and verifier**

Add:

```julia
function replay_quillen_supplied_evidence_patch(
    patch::QuillenSuppliedEvidencePatchAssembly,
)
    R = _require_quillen_denominator_cover_ring(patch.ring)
    n = patch.size
    _quillen_local_require_factor_matrix(patch.original_input, R, n, "supplied evidence original input")
    selected = _require_substitution_generator(R, patch.substitution_variable)
    local_certificates = patch.local_certificates
    local_certificates_ok = !isempty(local_certificates) &&
        all(verify_quillen_local_factor_sequence_certificate, local_certificates)
    expected_candidate = extract_quillen_denominator_cover_candidate(local_certificates)
    denominator_candidate_ok = verify_quillen_denominator_cover_candidate(patch.denominator_candidate)
    denominator_candidate_matches =
        denominator_candidate_ok &&
        patch.denominator_candidate.original_input == expected_candidate.original_input &&
        patch.denominator_candidate.raw_denominators == expected_candidate.raw_denominators &&
        _same_quillen_local_denominator_supports(
            patch.denominator_candidate.local_supports,
            expected_candidate.local_supports,
        )
    solver_result_ok = verify_quillen_denominator_cover_solver_result(patch.solver_result)
    solver_source_candidate_ok =
        solver_result_ok &&
        patch.solver_result.source_candidate isa QuillenDenominatorCoverCandidate &&
        patch.solver_result.source_candidate.raw_denominators == patch.denominator_candidate.raw_denominators &&
        patch.solver_result.source_candidate.original_input == patch.denominator_candidate.original_input
    cover_certificate_ok =
        verify_quillen_denominator_cover(patch.cover_certificate) &&
        _same_quillen_cover_certificate_data(patch.cover_certificate, patch.solver_result.cover_certificate)
    substitution_chain_ok = verify_quillen_patch_substitution_chain(patch.substitution_chain)
    substitution_chain_matches =
        substitution_chain_ok &&
        patch.substitution_chain.original_matrix == patch.original_input &&
        patch.substitution_chain.selected_variable == selected &&
        patch.substitution_chain.solver_result.raw_denominators == patch.solver_result.raw_denominators &&
        patch.substitution_chain.solver_result.exponent == patch.solver_result.exponent &&
        patch.substitution_chain.solver_result.coverage_multipliers == patch.solver_result.coverage_multipliers
    base_term_factors = _quillen_supplied_base_term_factors(R, n, patch.base_term_factors)
    base_term_product = _quillen_product(R, n, base_term_factors)
    base_term_ok =
        patch.base_term == patch.substitution_chain.base_term &&
        patch.base_term_product == base_term_product &&
        _quillen_supplied_base_term_ok(
            patch.base_term_policy,
            patch.base_term,
            base_term_factors,
            base_term_product,
            R,
            n,
        )
    sequence_expansions = [
        quillen_sequence_contribution_expansion(certificate, patch.solver_result, index)
        for (index, certificate) in enumerate(local_certificates)
    ]
    sequence_expansions_ok =
        all(verify_quillen_sequence_contribution_expansion, patch.sequence_expansions) &&
        _same_quillen_sequence_expansions(patch.sequence_expansions, sequence_expansions)
    factor_type = typeof(identity_matrix(R, n))
    sequence_elementary_factors = factor_type[]
    for expansion in sequence_expansions
        append!(sequence_elementary_factors, expansion.global_elementary_factors)
    end
    global_elementary_factors = copy(base_term_factors)
    append!(global_elementary_factors, sequence_elementary_factors)
    global_elementary_factors_ok =
        _same_quillen_factors(patch.sequence_elementary_factors, sequence_elementary_factors) &&
        _same_quillen_factors(patch.global_elementary_factors, global_elementary_factors)
    product = _quillen_product(R, n, global_elementary_factors)
    target = _quillen_local_require_factor_matrix(patch.target, R, n, "supplied evidence target")
    product_ok = global_elementary_factors_ok && patch.product == product && product == target
    target_ok = target == patch.original_input
    replay_metadata = _quillen_supplied_patch_metadata(
        patch.denominator_candidate,
        patch.solver_result,
        patch.substitution_chain,
        patch.base_term_policy,
        sequence_expansions,
        patch.replay_metadata.metadata,
    )
    replay_metadata_ok = patch.replay_metadata == replay_metadata
    overall_ok =
        local_certificates_ok &&
        denominator_candidate_ok &&
        denominator_candidate_matches &&
        solver_result_ok &&
        solver_source_candidate_ok &&
        cover_certificate_ok &&
        substitution_chain_ok &&
        substitution_chain_matches &&
        base_term_ok &&
        sequence_expansions_ok &&
        global_elementary_factors_ok &&
        product_ok &&
        target_ok &&
        replay_metadata_ok
    return QuillenSuppliedEvidencePatchAssemblyVerification(
        local_certificates_ok,
        denominator_candidate_ok,
        denominator_candidate_matches,
        solver_result_ok,
        solver_source_candidate_ok,
        cover_certificate_ok,
        substitution_chain_ok,
        substitution_chain_matches,
        base_term_ok,
        sequence_expansions_ok,
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

function verify_quillen_patch(patch::QuillenSuppliedEvidencePatchAssembly)::Bool
    try
        replay = replay_quillen_supplied_evidence_patch(patch)
        return replay.overall_ok &&
               _same_quillen_supplied_patch_verification(patch.verification, replay)
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end
```

- [ ] **Step 5: Add constructor**

Add:

```julia
function assemble_quillen_patch_from_local_evidence(
    A,
    selected_variable,
    local_certificates;
    max_exponent::Integer = 4,
    exponent = nothing,
    coverage_multipliers = nothing,
    supplied_multipliers = nothing,
    substitution_chain = nothing,
    base_term_policy = nothing,
    base_term_factors = nothing,
    metadata = (;),
)
    R = _require_quillen_denominator_cover_ring(base_ring(A))
    n = _require_square_matrix(A, "supplied evidence original input")
    selected = _require_substitution_generator(R, selected_variable)
    certificates = collect(local_certificates)
    isempty(certificates) &&
        throw(ArgumentError("Quillen supplied evidence patch assembly requires at least one local sequence certificate"))
    all(verify_quillen_local_factor_sequence_certificate, certificates) ||
        throw(ArgumentError("Quillen supplied evidence patch assembly requires verified local sequence certificates"))
    all(certificate -> certificate.original_input == A, certificates) ||
        throw(ArgumentError("Quillen supplied evidence patch assembly requires local evidence for the input matrix"))
    candidate = extract_quillen_denominator_cover_candidate(certificates)
    candidate.ring == R && candidate.size == n && candidate.selected_variable == selected ||
        throw(ArgumentError("Quillen supplied evidence patch assembly candidate context does not match the input"))
    solver_result = solve_quillen_denominator_cover(
        candidate;
        max_exponent,
        exponent,
        coverage_multipliers,
        supplied_multipliers,
    )
    chain = substitution_chain === nothing ?
        quillen_patch_substitution_chain(
            A,
            selected,
            solver_result;
            metadata = merge((; source = :quillen_supplied_evidence_patch_assembly), metadata),
        ) :
        substitution_chain
    verify_quillen_patch_substitution_chain(chain) ||
        throw(ArgumentError("Quillen supplied evidence patch assembly substitution chain does not replay"))
    chain.original_matrix == A && chain.selected_variable == selected ||
        throw(ArgumentError("Quillen supplied evidence patch assembly substitution chain context does not match"))
    chain.solver_result.raw_denominators == solver_result.raw_denominators &&
        chain.solver_result.exponent == solver_result.exponent &&
        chain.solver_result.coverage_multipliers == solver_result.coverage_multipliers ||
        throw(ArgumentError("Quillen supplied evidence patch assembly substitution chain solver data does not match"))
    policy = _quillen_supplied_base_term_policy(base_term_policy, base_term_factors)
    normalized_base_factors = _quillen_supplied_base_term_factors(R, n, base_term_factors)
    base_product = _quillen_product(R, n, normalized_base_factors)
    _quillen_supplied_base_term_ok(policy, chain.base_term, normalized_base_factors, base_product, R, n) ||
        throw(ArgumentError("Quillen supplied evidence patch assembly base-term evidence is missing or does not replay"))
    sequence_expansions = [
        quillen_sequence_contribution_expansion(certificate, solver_result, index)
        for (index, certificate) in enumerate(certificates)
    ]
    factor_type = typeof(identity_matrix(R, n))
    sequence_elementary_factors = factor_type[]
    for expansion in sequence_expansions
        append!(sequence_elementary_factors, expansion.global_elementary_factors)
    end
    global_elementary_factors = copy(normalized_base_factors)
    append!(global_elementary_factors, sequence_elementary_factors)
    product = _quillen_product(R, n, global_elementary_factors)
    product == A ||
        throw(ArgumentError("Quillen supplied evidence patch assembly factors do not multiply to the input matrix"))
    replay_metadata = _quillen_supplied_patch_metadata(
        candidate,
        solver_result,
        chain,
        policy,
        sequence_expansions,
        metadata,
    )
    provisional = QuillenSuppliedEvidencePatchAssembly(
        R,
        n,
        selected,
        A,
        certificates,
        candidate,
        solver_result,
        solver_result.cover_certificate,
        chain,
        policy,
        chain.base_term,
        normalized_base_factors,
        base_product,
        sequence_expansions,
        sequence_elementary_factors,
        global_elementary_factors,
        product,
        A,
        replay_metadata,
        QuillenSuppliedEvidencePatchAssemblyVerification(
            false,
            false,
            false,
            false,
            false,
            false,
            false,
            false,
            false,
            false,
            factor_type[],
            false,
            identity_matrix(R, n),
            false,
            A,
            false,
            replay_metadata,
            false,
            false,
        ),
    )
    verification = replay_quillen_supplied_evidence_patch(provisional)
    verification.overall_ok ||
        throw(ArgumentError("Quillen supplied evidence patch assembly data does not replay"))
    return QuillenSuppliedEvidencePatchAssembly(
        provisional.ring,
        provisional.size,
        provisional.substitution_variable,
        provisional.original_input,
        provisional.local_certificates,
        provisional.denominator_candidate,
        provisional.solver_result,
        provisional.cover_certificate,
        provisional.substitution_chain,
        provisional.base_term_policy,
        provisional.base_term,
        provisional.base_term_factors,
        provisional.base_term_product,
        provisional.sequence_expansions,
        provisional.sequence_elementary_factors,
        provisional.global_elementary_factors,
        provisional.product,
        provisional.target,
        provisional.replay_metadata,
        verification,
    )
end
```

- [ ] **Step 6: Run focused test**

Run:

```bash
julia --project=. -e 'include("test/expert/quillen_supplied_evidence_patch_assembly.jl")'
```

Expected: PASS.

- [ ] **Step 7: Commit assembly implementation**

```bash
git add src/algorithm/quillen_induction.jl
git commit -m "feat: assemble quillen patches from supplied evidence"
```

### Task 4: Add Constructive Regression And Verify

**Files:**
- Modify: `test/expert/quillen_induction_constructive.jl`

**Interfaces:**
- Consumes `assemble_quillen_patch_from_local_evidence` from Task 3.
- Produces a regression that the existing constructive acceptance fixture can use supplied evidence assembly.

- [ ] **Step 1: Add sequence-certificate helpers to the constructive test**

In `test/expert/quillen_induction_constructive.jl`, add helpers near `constructive_local_certificate_from_fixture`:

```julia
function constructive_sequence_certificate_from_fixture(entry; local_index::Int)
    local_factor = entry.local_factors[local_index]
    realization = constructive_local_certificate_from_fixture(entry; local_index)
    return Suslin.quillen_local_factor_sequence_certificate(
        realization;
        factor_provenance = (;
            factor_index = 1,
            sequence_index = 1,
            local_index = 1,
            fixture_id = entry.id,
            source = :constructive_supplied_evidence_regression,
        ),
        metadata = (;
            source_refs = entry.source_refs,
            consumer_issue_ids = entry.consumer_issue_ids,
        ),
    )
end

function constructive_sequence_certificates(entry)
    return [
        constructive_sequence_certificate_from_fixture(entry; local_index = idx)
        for idx in eachindex(entry.local_factors)
    ]
end
```

- [ ] **Step 2: Add the regression inside the existing testset**

After the positive loop that fills `built`, add:

```julia
    constructive_entry = entries["quillen-constructive-acceptance-gf2"]
    supplied_sequence_certificates = constructive_sequence_certificates(constructive_entry)
    supplied_patch = Suslin.assemble_quillen_patch_from_local_evidence(
        constructive_entry.target_matrix,
        constructive_entry.substitution_variable,
        supplied_sequence_certificates;
        max_exponent = 2,
        base_term_policy = :already_handled,
        metadata = (; fixture_id = constructive_entry.id, consumer_issue_id = "#218"),
    )
    @test supplied_patch isa Suslin.QuillenSuppliedEvidencePatchAssembly
    @test Suslin.verify_quillen_patch(supplied_patch)
    @test supplied_patch.product == constructive_entry.target_matrix
    @test supplied_patch.global_elementary_factors ==
          reduce(
              vcat,
              [expansion.global_elementary_factors for expansion in supplied_patch.sequence_expansions];
              init = typeof(identity_matrix(supplied_patch.ring, supplied_patch.size))[],
          )
```

- [ ] **Step 3: Run required focused commands**

Run:

```bash
julia --project=. -e 'include("test/expert/quillen_supplied_evidence_patch_assembly.jl")'
julia --project=. -e 'include("test/expert/quillen_induction_constructive.jl")'
```

Expected: both PASS.

- [ ] **Step 4: Commit regression**

```bash
git add test/expert/quillen_induction_constructive.jl
git commit -m "test: add constructive supplied evidence assembly regression"
```

### Task 5: Final Verification And Review

**Files:**
- No planned edits unless verification reveals defects.

**Interfaces:**
- Consumes all implementation tasks.
- Produces final verification evidence for PR creation.

- [ ] **Step 1: Run diff hygiene**

Run:

```bash
git diff --check
```

Expected: exit 0.

- [ ] **Step 2: Run required focused tests**

Run:

```bash
julia --project=. -e 'include("test/expert/quillen_supplied_evidence_patch_assembly.jl")'
julia --project=. -e 'include("test/expert/quillen_induction_constructive.jl")'
```

Expected: both exit 0.

- [ ] **Step 3: Run package entry point**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: exit 0.

- [ ] **Step 4: Request final code review**

Use `superpowers:requesting-code-review` with a diff package from the merge base to `HEAD`. Fix Critical and Important findings, then rerun the covering verification commands.

- [ ] **Step 5: Finish branch**

Use `superpowers:verification-before-completion` and `superpowers:finishing-a-development-branch`. Choose "Push and create a Pull Request" under the Standing Answer Policy.

## Automatic Execution Choice

Plan complete. Under the Standing Answer Policy, choose the recommended execution option:

1. Subagent-Driven (recommended) - use `superpowers:subagent-driven-development`.
2. Inline Execution - use `superpowers:executing-plans`.

Automatic choice: Subagent-Driven, because it is marked recommended by the Superpowers plan handoff and this run is non-interactive.

## Plan Self-Review

- Every issue requirement maps to a task: local verification, denominator extraction, cover solving, substitution-chain replay, explicit sequence expansion, base-term boundary handling, exact product verification, negative controls, focused commands, and `Pkg.test()`.
- No `TBD` or `TODO` placeholders remain.
- Type and function names are consistent across tasks.
- The plan keeps new names internal and does not add public exports.
