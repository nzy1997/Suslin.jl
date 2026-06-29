using Test
using Suslin
using Oscar

const QUILLEN_MAINLINE_CATALOG_PATH = joinpath(@__DIR__, "..", "fixtures", "quillen_mainline_cases.jl")
if !isdefined(Main, :QuillenMainlineFixtureCatalog)
    include(QUILLEN_MAINLINE_CATALOG_PATH)
end

function qde_local_certificate(local_factor)
    return Suslin.LocalCertificate(
        local_factor.certificate.indices,
        local_factor.certificate.denominators,
    )
end

function qde_local_correction(local_factor)
    correction = local_factor.correction
    return Suslin.QuillenElementaryCorrection(
        correction.row,
        correction.col,
        correction.entry,
    )
end

function qde_factor_record(entry, index::Int)
    local_factor = entry.patch_case.local_factors[index]
    evidence = entry.local_evidence.records[index]
    return (;
        row = local_factor.correction.row,
        col = local_factor.correction.col,
        numerator = local_factor.correction.entry,
        denominator = local_factor.denominator,
        coverage_multiplier = local_factor.coverage_multiplier,
        local_certificate = qde_local_certificate(local_factor),
        provenance = (;
            factor_index = 1,
            sequence_index = 1,
            local_index = index,
            fixture_id = entry.id,
            source = evidence.source,
            sequence_status = evidence.sequence_status,
        ),
        metadata = (;
            source_refs = entry.source_refs,
            consumer_issue_ids = entry.consumer_issue_ids,
        ),
    )
end

function qde_sequence_certificate(entry, index::Int; selected_variable = entry.patch_case.substitution_variable)
    local_factor = entry.patch_case.local_factors[index]
    record = qde_factor_record(entry, index)
    # #214 sequence certificates are one-factor artifacts, so their positional provenance
    # is normalized to 1 even when the surrounding mainline fixture index is larger.
    factor_provenance = merge(record.provenance, (; factor_index = 1, sequence_index = 1, local_index = 1))
    realization = Suslin.quillen_local_realization_certificate(
        entry.patch_case.target_matrix,
        selected_variable;
        local_certificate = qde_local_certificate(local_factor),
        denominator = local_factor.denominator,
        coverage_multiplier = local_factor.coverage_multiplier,
        correction = qde_local_correction(local_factor),
        factors = [local_factor.factor],
        patched_substitution_witness = entry.patch_case.patched_substitution_witness,
        ring = entry.patch_case.ring.object,
        size = entry.patch_case.size,
    )
    return Suslin.quillen_local_factor_sequence_certificate(
        realization;
        factor_provenance = factor_provenance,
        metadata = record.metadata,
    )
end

function qde_sequence_certificates(entry)
    return [
        qde_sequence_certificate(entry, index)
        for index in eachindex(entry.patch_case.local_factors)
    ]
end

function qde_rebuild_factor(factor; kwargs...)
    overrides = NamedTuple(kwargs)
    fields = merge((
        row = factor.row,
        col = factor.col,
        numerator = factor.numerator,
        denominator = factor.denominator,
        coverage_multiplier = factor.coverage_multiplier,
        local_certificate = factor.local_certificate,
        provenance = factor.provenance,
        metadata = factor.metadata,
    ), overrides)
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

function qde_rebuild_sequence_certificate(cert; kwargs...)
    overrides = NamedTuple(kwargs)
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
    ), overrides)
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

function qde_rebuild_support(support; kwargs...)
    overrides = NamedTuple(kwargs)
    fields = merge((
        local_index = support.local_index,
        support_denominator = support.support_denominator,
        support_kind = support.support_kind,
        factor_denominators = support.factor_denominators,
        factor_entries = support.factor_entries,
        factor_provenance = support.factor_provenance,
        replayed_denominator = support.replayed_denominator,
        replay_equality = support.replay_equality,
        replay_ok = support.replay_ok,
    ), overrides)
    return Suslin.QuillenLocalDenominatorSupport(
        fields.local_index,
        fields.support_denominator,
        fields.support_kind,
        fields.factor_denominators,
        fields.factor_entries,
        fields.factor_provenance,
        fields.replayed_denominator,
        fields.replay_equality,
        fields.replay_ok,
    )
end

function qde_rebuild_candidate(candidate; kwargs...)
    overrides = NamedTuple(kwargs)
    fields = merge((
        original_input = candidate.original_input,
        ring = candidate.ring,
        size = candidate.size,
        selected_variable = candidate.selected_variable,
        local_certificates = candidate.local_certificates,
        raw_denominators = candidate.raw_denominators,
        local_supports = candidate.local_supports,
        replay_metadata = candidate.replay_metadata,
        verification = candidate.verification,
    ), overrides)
    return Suslin.QuillenDenominatorCoverCandidate(
        fields.original_input,
        fields.ring,
        fields.size,
        fields.selected_variable,
        fields.local_certificates,
        fields.raw_denominators,
        fields.local_supports,
        fields.replay_metadata,
        fields.verification,
    )
end

function qde_trivial_sequence_certificate(R, selected_variable, n::Int)
    correction = Suslin.QuillenElementaryCorrection(1, 2, zero(R))
    local_certificate = Suslin.LocalCertificate([1, 2], [one(R), one(R)])
    factor = elementary_matrix(n, 1, 2, zero(R), R)
    realization = Suslin.quillen_local_realization_certificate(
        correction,
        selected_variable;
        local_certificate = local_certificate,
        denominator = one(R),
        coverage_multiplier = one(R),
        correction = correction,
        factors = [factor],
        ring = R,
        size = n,
    )
    return Suslin.quillen_local_factor_sequence_certificate(
        realization;
        factor_provenance = (;
            factor_index = 1,
            sequence_index = 1,
            local_index = 1,
            source = :trivial_verified_sequence,
        ),
        metadata = (; source = :trivial_verified_sequence),
    )
end

@testset "Quillen denominator extraction" begin
    entries = Main.QuillenMainlineFixtureCatalog.cases_by_id()
    entry = entries["quillen-two-open-cover-qq"]
    nontrivial_entry = entries["quillen-nontrivial-multipliers-qq"]
    constructive_entry = entries["quillen-constructive-acceptance-gf2"]
    certificates = qde_sequence_certificates(entry)
    nontrivial_certificates = qde_sequence_certificates(nontrivial_entry)
    constructive_certificates = qde_sequence_certificates(constructive_entry)

    @test all(Suslin.verify_quillen_local_factor_sequence_certificate, certificates)
    @test all(Suslin.verify_quillen_local_factor_sequence_certificate, nontrivial_certificates)
    @test all(Suslin.verify_quillen_local_factor_sequence_certificate, constructive_certificates)

    candidate = Suslin.extract_quillen_denominator_cover_candidate(certificates)

    @test candidate isa Suslin.QuillenDenominatorCoverCandidate
    @test Suslin.verify_quillen_denominator_cover_candidate(candidate)
    @test candidate.raw_denominators == [data.denominator for data in entry.patch_case.denominator_data]
    @test candidate.selected_variable == entry.patch_case.substitution_variable
    @test candidate.original_input == entry.patch_case.target_matrix
    @test length(candidate.local_supports) == 2
    @test all(support -> support.support_kind == :product, candidate.local_supports)
    @test all(support -> support.replay_ok, candidate.local_supports)
    @test [only(support.factor_denominators) for support in candidate.local_supports] == candidate.raw_denominators
    @test candidate.verification.local_certificates_ok
    @test candidate.verification.local_supports_ok
    @test candidate.verification.raw_denominators_ok
    @test candidate.verification.overall_ok

    @test_throws ArgumentError Suslin.extract_quillen_denominator_cover_candidate(
        Suslin.QuillenLocalFactorSequenceCertificate[],
    )

    dropped_candidate = qde_rebuild_candidate(
        candidate;
        local_certificates = candidate.local_certificates[1:1],
        raw_denominators = candidate.raw_denominators[1:1],
        local_supports = candidate.local_supports[1:1],
    )
    @test !Suslin.verify_quillen_denominator_cover_candidate(dropped_candidate)

    edited_raw_denominator_candidate = qde_rebuild_candidate(
        candidate;
        raw_denominators = [
            candidate.raw_denominators[1] + one(candidate.ring),
            candidate.raw_denominators[2],
        ],
    )
    @test !Suslin.verify_quillen_denominator_cover_candidate(edited_raw_denominator_candidate)

    edited_support_denominator_candidate = qde_rebuild_candidate(
        candidate;
        local_supports = begin
            modified = collect(candidate.local_supports)
            modified[1] = qde_rebuild_support(
                modified[1];
                support_denominator = modified[1].support_denominator + one(candidate.ring),
            )
            modified
        end,
    )
    @test !Suslin.verify_quillen_denominator_cover_candidate(edited_support_denominator_candidate)

    edited_factor_denominator_candidate = qde_rebuild_candidate(
        candidate;
        local_supports = begin
            modified = collect(candidate.local_supports)
            factor_denominators = copy(modified[1].factor_denominators)
            factor_denominators[1] += one(candidate.ring)
            modified[1] = qde_rebuild_support(
                modified[1];
                factor_denominators = factor_denominators,
            )
            modified
        end,
    )
    @test !Suslin.verify_quillen_denominator_cover_candidate(edited_factor_denominator_candidate)

    mixed_variable_certificates = collect(certificates)
    mixed_variable_certificates[2] = qde_sequence_certificate(
        entry,
        2;
        selected_variable = first(filter(
            gen -> gen != certificates[1].selected_variable,
            collect(gens(certificates[1].ring)),
        )),
    )
    @test_throws ArgumentError Suslin.extract_quillen_denominator_cover_candidate(mixed_variable_certificates)

    unverified_local_sequence_certificates = collect(certificates)
    unverified_local_sequence_certificates[2] = qde_rebuild_sequence_certificate(
        certificates[2];
        factors = begin
            modified = collect(certificates[2].factors)
            modified[1] = qde_rebuild_factor(
                modified[1];
                numerator = modified[1].numerator + one(certificates[2].ring),
            )
            modified
        end,
    )
    @test !Suslin.verify_quillen_local_factor_sequence_certificate(unverified_local_sequence_certificates[2])
    @test_throws ArgumentError Suslin.extract_quillen_denominator_cover_candidate(unverified_local_sequence_certificates)

    mixed_original_input_certificates = [
        certificates[1],
        nontrivial_certificates[1],
    ]
    @test certificates[1].original_input != nontrivial_certificates[1].original_input
    @test certificates[1].ring == nontrivial_certificates[1].ring
    @test certificates[1].size == nontrivial_certificates[1].size
    @test certificates[1].selected_variable == nontrivial_certificates[1].selected_variable
    @test_throws ArgumentError Suslin.extract_quillen_denominator_cover_candidate(mixed_original_input_certificates)

    mixed_ring_certificates = [
        certificates[1],
        constructive_certificates[1],
    ]
    @test all(Suslin.verify_quillen_local_factor_sequence_certificate, mixed_ring_certificates)
    @test mixed_ring_certificates[1].ring != mixed_ring_certificates[2].ring
    @test_throws ArgumentError Suslin.extract_quillen_denominator_cover_candidate(mixed_ring_certificates)

    larger_size_certificate = qde_trivial_sequence_certificate(
        certificates[1].ring,
        certificates[1].selected_variable,
        certificates[1].size + 1,
    )
    mixed_size_certificates = [
        certificates[1],
        larger_size_certificate,
    ]
    @test all(Suslin.verify_quillen_local_factor_sequence_certificate, mixed_size_certificates)
    @test mixed_size_certificates[1].ring == mixed_size_certificates[2].ring
    @test mixed_size_certificates[1].selected_variable == mixed_size_certificates[2].selected_variable
    @test mixed_size_certificates[1].size != mixed_size_certificates[2].size
    @test_throws ArgumentError Suslin.extract_quillen_denominator_cover_candidate(mixed_size_certificates)
end
