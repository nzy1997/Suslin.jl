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
