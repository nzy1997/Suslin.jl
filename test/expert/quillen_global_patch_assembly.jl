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
