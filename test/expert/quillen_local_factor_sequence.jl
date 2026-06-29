using Test
using Suslin
using Oscar

include(joinpath(@__DIR__, "..", "fixtures", "quillen_mainline_cases.jl"))

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
    normalized_global_elementary_factors =
        [qlfs_replay_factor(factor, R, n) for factor in factors]
    local_product = qlfs_product(normalized_global_elementary_factors, R, n)
    return (;
        raw_denominators = raw_denominators,
        product_denominator = prod(raw_denominators; init = one(R)),
        normalized_global_elementary_factors = normalized_global_elementary_factors,
        local_product = local_product,
    )
end

function qlfs_rebuild_sequence_certificate(cert; recompute_derived::Bool = false, kwargs...)
    overrides = NamedTuple(kwargs)
    derived = recompute_derived ? qlfs_sequence_fields(cert.factors, cert.ring, cert.size) : (;)
    fields = merge((
        ring = cert.ring,
        size = cert.size,
        selected_variable = cert.selected_variable,
        factors = cert.factors,
        raw_denominators = cert.raw_denominators,
        product_denominator = cert.product_denominator,
        normalized_global_elementary_factors = cert.normalized_global_elementary_factors,
        local_product = cert.local_product,
        local_correction = cert.local_correction,
        verification = cert.verification,
    ), derived, overrides)
    if recompute_derived && !(:factors in keys(overrides))
        recomputed = qlfs_sequence_fields(fields.factors, fields.ring, fields.size)
        fields = merge(fields, recomputed)
    elseif :factors in keys(overrides) &&
           !(:raw_denominators in keys(overrides)) &&
           !(:product_denominator in keys(overrides)) &&
           !(:normalized_global_elementary_factors in keys(overrides)) &&
           !(:local_product in keys(overrides))
        recomputed = qlfs_sequence_fields(fields.factors, fields.ring, fields.size)
        fields = merge(fields, recomputed)
    end
    return Suslin.QuillenLocalFactorSequenceCertificate(
        fields.ring,
        fields.size,
        fields.selected_variable,
        fields.factors,
        fields.raw_denominators,
        fields.product_denominator,
        fields.normalized_global_elementary_factors,
        fields.local_product,
        fields.local_correction,
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
            provenance = merge(modified[1].provenance, (; source = :tampered_source)),
        )
        modified
    end)
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
end
