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
