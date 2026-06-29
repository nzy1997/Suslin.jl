using Test
using Suslin
using Oscar

include(joinpath(@__DIR__, "..", "fixtures", "quillen_mainline_cases.jl"))

mutable struct QLFSFlakyProvenance
    isempty_calls::Int
end

function Base.isempty(provenance::QLFSFlakyProvenance)
    provenance.isempty_calls += 1
    return provenance.isempty_calls > 1
end

function qlfs_product(factors, R, n)
    product = identity_matrix(R, n)
    for factor in factors
        product *= factor
    end
    return product
end

function qlfs_local_certificate(local_factor)
    return Suslin.LocalCertificate(
        local_factor.certificate.indices,
        local_factor.certificate.denominators,
    )
end

function qlfs_factor_record(entry, index::Int)
    local_factor = entry.patch_case.local_factors[index]
    evidence = entry.local_evidence.records[index]
    return (;
        row = local_factor.correction.row,
        col = local_factor.correction.col,
        numerator = local_factor.correction.entry,
        denominator = local_factor.denominator,
        coverage_multiplier = local_factor.coverage_multiplier,
        local_certificate = qlfs_local_certificate(local_factor),
        provenance = (;
            factor_index = index,
            sequence_index = index,
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

function qlfs_local_certificate_with_denominator(factor, denominator)
    indices = factor.local_certificate.indices
    denominators = [
        (index == factor.row || index == factor.col) ? denominator : factor.local_certificate.denominators[position]
        for (position, index) in enumerate(indices)
    ]
    return Suslin.LocalCertificate(indices, denominators)
end

function qlfs_rebuild_factor(factor; kwargs...)
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

function qlfs_replay_factor(factor, R, n)
    weighted_entry = factor.coverage_multiplier * factor.denominator * factor.numerator
    return elementary_matrix(n, factor.row, factor.col, weighted_entry, R)
end

function qlfs_sequence_fields(factors, R, n)
    raw_denominators = [factor.denominator for factor in factors]
    normalized_local_contributions = [
        Suslin.QuillenLocalContribution(
            factor.local_certificate,
            factor.denominator,
            factor.coverage_multiplier,
            Suslin.QuillenElementaryCorrection(factor.row, factor.col, factor.numerator),
        )
        for factor in factors
    ]
    normalized_global_elementary_factors =
        [qlfs_replay_factor(factor, R, n) for factor in factors]
    local_product = qlfs_product(normalized_global_elementary_factors, R, n)
    return (;
        raw_denominators = raw_denominators,
        product_denominator = prod(raw_denominators; init = one(R)),
        normalized_local_contributions = normalized_local_contributions,
        normalized_global_elementary_factors = normalized_global_elementary_factors,
        local_product = local_product,
    )
end

function qlfs_denominator_data(contributions)
    return [
        Suslin.QuillenDenominatorData(
            contribution.denominator,
            contribution.coverage_multiplier,
        )
        for contribution in contributions
    ]
end

function qlfs_replay_metadata(fields)
    return (;
        factor_count = length(fields.factors),
        raw_denominators = [factor.denominator for factor in fields.factors],
        denominator_data = qlfs_denominator_data(fields.normalized_local_contributions),
        factor_provenance = [factor.provenance for factor in fields.factors],
        factor_metadata = [factor.metadata for factor in fields.factors],
        factor_local_certificates = [factor.local_certificate for factor in fields.factors],
        has_patched_substitution_witness = fields.patched_substitution_witness !== nothing,
        has_chain_witness = fields.chain_witness !== nothing,
        chain_witness = fields.chain_witness,
        witness_metadata = fields.witness_metadata,
    )
end

function qlfs_rebuild_sequence_certificate(cert; recompute_derived::Bool = false, kwargs...)
    overrides = NamedTuple(kwargs)
    derived = recompute_derived ? qlfs_sequence_fields(cert.factors, cert.ring, cert.size) : (;)
    fields = merge((
        original_input = cert.original_input,
        ring = cert.ring,
        size = cert.size,
        selected_variable = cert.selected_variable,
        factors = cert.factors,
        raw_denominators = cert.raw_denominators,
        product_denominator = cert.product_denominator,
        normalized_local_contributions = cert.normalized_local_contributions,
        normalized_global_elementary_factors = cert.normalized_global_elementary_factors,
        local_product = cert.local_product,
        local_correction = cert.local_correction,
        patched_substitution_witness = cert.patched_substitution_witness,
        chain_witness = cert.chain_witness,
        witness_metadata = cert.witness_metadata,
        replay_metadata = cert.replay_metadata,
        verification = cert.verification,
    ), derived, overrides)
    if recompute_derived && !(:factors in keys(overrides))
        recomputed = qlfs_sequence_fields(fields.factors, fields.ring, fields.size)
        fields = merge(fields, recomputed)
    elseif :factors in keys(overrides) &&
           !(:raw_denominators in keys(overrides)) &&
           !(:product_denominator in keys(overrides)) &&
           !(:normalized_local_contributions in keys(overrides)) &&
           !(:normalized_global_elementary_factors in keys(overrides)) &&
           !(:local_product in keys(overrides))
        recomputed = qlfs_sequence_fields(fields.factors, fields.ring, fields.size)
        fields = merge(fields, recomputed)
    end
    if recompute_derived && !(:replay_metadata in keys(overrides))
        fields = merge(fields, (; replay_metadata = qlfs_replay_metadata(fields)))
    end
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

@testset "Quillen local factor sequence certificates" begin
    entries = QuillenMainlineFixtureCatalog.cases_by_id()
    entry = entries["quillen-patched-substitution-witness-qq"]
    factors = [qlfs_factor_record(entry, index) for index in eachindex(entry.patch_case.local_factors)]

    cert = Suslin.quillen_local_factor_sequence_certificate(
        entry.patch_case.target_matrix,
        entry.patch_case.substitution_variable;
        factors = factors,
        ring = entry.patch_case.ring.object,
        size = entry.patch_case.size,
        local_evidence = entry.local_evidence,
        provenance = (;
            fixture_id = entry.id,
            source_refs = entry.source_refs,
            consumer_issue_ids = entry.consumer_issue_ids,
        ),
    )

    @test cert isa Suslin.QuillenLocalFactorSequenceCertificate
    @test Suslin.verify_quillen_local_factor_sequence_certificate(cert)
    replay = Suslin.replay_quillen_local_factor_sequence(cert)
    @test length(cert.factors) == 2
    @test cert.raw_denominators == [factor.denominator for factor in cert.factors]
    @test cert.product_denominator == prod(cert.raw_denominators; init = one(cert.ring))
    @test cert.normalized_global_elementary_factors == collect(entry.local_evidence.factors)
    @test cert.local_product == qlfs_product(cert.normalized_global_elementary_factors, cert.ring, cert.size)
    @test cert.local_correction == entry.local_evidence.expected_product
    @test replay.raw_denominators == cert.raw_denominators
    @test replay.product_denominator == cert.product_denominator
    @test replay.normalized_global_elementary_factors == cert.normalized_global_elementary_factors
    @test replay.local_product == cert.local_product
    @test cert.verification.factor_provenance_ok
    @test cert.verification.product_denominator_ok
    @test cert.verification.overall_ok

    bad_numerator = qlfs_rebuild_sequence_certificate(cert; factors = begin
        modified = collect(cert.factors)
        modified[1] = qlfs_rebuild_factor(
            modified[1];
            numerator = modified[1].numerator + one(cert.ring),
        )
        modified
    end)
    @test !Suslin.verify_quillen_local_factor_sequence_certificate(bad_numerator)

    bad_denominator = qlfs_rebuild_sequence_certificate(cert; factors = begin
        modified = collect(cert.factors)
        new_denominator = modified[1].denominator + one(cert.ring)
        modified[1] = qlfs_rebuild_factor(
            modified[1];
            denominator = new_denominator,
            local_certificate = qlfs_local_certificate_with_denominator(modified[1], new_denominator),
        )
        modified
    end)
    bad_denominator_replay = Suslin.replay_quillen_local_factor_sequence(bad_denominator)
    @test bad_denominator.raw_denominators != cert.raw_denominators
    @test bad_denominator.normalized_global_elementary_factors != cert.normalized_global_elementary_factors
    @test bad_denominator_replay.raw_denominators == bad_denominator.raw_denominators
    @test bad_denominator_replay.product_denominator == bad_denominator.product_denominator
    @test bad_denominator_replay.normalized_global_elementary_factors ==
          bad_denominator.normalized_global_elementary_factors
    @test bad_denominator_replay.local_product == bad_denominator.local_product
    @test !bad_denominator_replay.overall_ok
    @test !Suslin.verify_quillen_local_factor_sequence_certificate(bad_denominator)

    bad_variable = qlfs_rebuild_sequence_certificate(cert; selected_variable = first(filter(
        gen -> gen != cert.selected_variable,
        collect(gens(cert.ring)),
    )))
    @test !Suslin.verify_quillen_local_factor_sequence_certificate(bad_variable)

    bad_provenance = qlfs_rebuild_sequence_certificate(cert; factors = begin
        modified = collect(cert.factors)
        modified[1] = qlfs_rebuild_factor(
            modified[1];
            provenance = merge(modified[1].provenance, (; factor_index = 99)),
        )
        modified
    end, recompute_derived = true)
    bad_provenance = qlfs_rebuild_sequence_certificate(
        bad_provenance;
        verification = Suslin.replay_quillen_local_factor_sequence(bad_provenance),
    )
    bad_provenance_replay = Suslin.replay_quillen_local_factor_sequence(bad_provenance)
    @test !bad_provenance_replay.factor_provenance_ok
    @test bad_provenance_replay.product_denominator_ok
    @test bad_provenance_replay.replay_metadata_ok
    @test !bad_provenance_replay.overall_ok
    @test !Suslin.verify_quillen_local_factor_sequence_certificate(bad_provenance)

    reversed_order = qlfs_rebuild_sequence_certificate(
        cert;
        factors = reverse(collect(cert.factors)),
        recompute_derived = true,
    )
    reversed_order = qlfs_rebuild_sequence_certificate(
        reversed_order;
        verification = Suslin.replay_quillen_local_factor_sequence(reversed_order),
    )
    reversed_order_replay = Suslin.replay_quillen_local_factor_sequence(reversed_order)
    @test reversed_order.factors == reverse(collect(cert.factors))
    @test reversed_order.normalized_global_elementary_factors ==
          reverse(cert.normalized_global_elementary_factors)
    @test reversed_order_replay.raw_denominators == reversed_order.raw_denominators
    @test reversed_order_replay.normalized_global_elementary_factors ==
          reversed_order.normalized_global_elementary_factors
    @test reversed_order_replay.local_product == reversed_order.local_product
    @test !reversed_order_replay.factor_provenance_ok
    @test reversed_order_replay.product_denominator_ok
    @test reversed_order_replay.replay_metadata_ok
    @test !reversed_order_replay.overall_ok
    @test !Suslin.verify_quillen_local_factor_sequence_certificate(reversed_order)

    dict_provenance_factor = begin
        source = cert.factors[1]
        (;
            row = source.row,
            col = source.col,
            numerator = source.numerator,
            denominator = source.denominator,
            coverage_multiplier = source.coverage_multiplier,
            provenance = Dict(
                :factor_index => 1,
                :sequence_index => 1,
                :local_index => 1,
                :fixture_id => entry.id,
            ),
        )
    end
    dict_default_cert = Suslin.quillen_local_factor_sequence_certificate(
        entry.patch_case.target_matrix,
        entry.patch_case.substitution_variable;
        factors = [dict_provenance_factor],
        ring = entry.patch_case.ring.object,
        size = entry.patch_case.size,
    )
    @test Suslin.verify_quillen_local_factor_sequence_certificate(dict_default_cert)
    @test dict_default_cert.local_correction == dict_default_cert.local_product
    @test dict_default_cert.factors[1].metadata == (; factor_index = 1)
    @test dict_default_cert.factors[1].local_certificate.indices ==
          [dict_default_cert.factors[1].row, dict_default_cert.factors[1].col]
    @test dict_default_cert.factors[1].local_certificate.denominators ==
          [dict_default_cert.factors[1].denominator, dict_default_cert.factors[1].denominator]
    @test dict_default_cert.verification.factor_provenance_ok

    empty_vector_provenance_factor = merge(
        dict_provenance_factor,
        (; provenance = Symbol[]),
    )
    @test_throws ArgumentError Suslin.quillen_local_factor_sequence_certificate(
        entry.patch_case.target_matrix,
        entry.patch_case.substitution_variable;
        factors = [empty_vector_provenance_factor],
        ring = entry.patch_case.ring.object,
        size = entry.patch_case.size,
    )

    legacy_source = qlfs_rebuild_sequence_certificate(
        cert;
        patched_substitution_witness = nothing,
        chain_witness = nothing,
        witness_metadata = (;),
        recompute_derived = true,
    )
    legacy_source = qlfs_rebuild_sequence_certificate(
        legacy_source;
        verification = Suslin.replay_quillen_local_factor_sequence(legacy_source),
    )
    legacy_cert = Suslin.QuillenLocalFactorSequenceCertificate(
        legacy_source.ring,
        legacy_source.size,
        legacy_source.selected_variable,
        legacy_source.factors,
        legacy_source.raw_denominators,
        legacy_source.product_denominator,
        legacy_source.normalized_global_elementary_factors,
        legacy_source.local_product,
        legacy_source.local_correction,
        legacy_source.verification,
    )
    @test Suslin.verify_quillen_local_factor_sequence_certificate(legacy_cert)
    @test legacy_cert.original_input == legacy_source.original_input
    @test legacy_cert.replay_metadata == legacy_source.replay_metadata

    incomplete_certificate_factor = qlfs_rebuild_factor(
        legacy_source.factors[1];
        local_certificate = Suslin.LocalCertificate(
            [legacy_source.factors[1].row],
            [legacy_source.factors[1].denominator],
        ),
    )
    incomplete_legacy_cert = Suslin.QuillenLocalFactorSequenceCertificate(
        legacy_source.ring,
        legacy_source.size,
        legacy_source.selected_variable,
        [incomplete_certificate_factor],
        [incomplete_certificate_factor.denominator],
        incomplete_certificate_factor.denominator,
        legacy_source.normalized_global_elementary_factors,
        legacy_source.local_product,
        legacy_source.local_correction,
        legacy_source.verification,
    )
    @test isempty(incomplete_legacy_cert.normalized_local_contributions)

    flaky_provenance_factor = qlfs_rebuild_factor(
        cert.factors[1];
        provenance = QLFSFlakyProvenance(0),
    )
    flaky_provenance = qlfs_rebuild_sequence_certificate(
        cert;
        factors = [flaky_provenance_factor],
        recompute_derived = true,
    )
    flaky_provenance_replay = Suslin.replay_quillen_local_factor_sequence(flaky_provenance)
    @test !flaky_provenance_replay.factor_provenance_ok
    @test !flaky_provenance_replay.overall_ok

    malformed_size = qlfs_rebuild_sequence_certificate(cert; size = 1)
    @test !Suslin.verify_quillen_local_factor_sequence_certificate(malformed_size)
end
