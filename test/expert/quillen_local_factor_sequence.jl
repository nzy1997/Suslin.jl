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

function qlfs_factor_record(entry, index::Int)
    local_factor = entry.patch_case.local_factors[index]
    evidence = entry.local_evidence.records[index]
    return (;
        row = local_factor.correction.row,
        col = local_factor.correction.col,
        numerator = local_factor.correction.entry,
        denominator = local_factor.denominator,
        coverage_multiplier = local_factor.coverage_multiplier,
        local_certificate = Suslin.LocalCertificate(
            local_factor.certificate.indices,
            local_factor.certificate.denominators,
        ),
        provenance = (;
            factor_index = index,
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

function qlfs_rebuild_sequence_certificate(cert; kwargs...)
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
    ), kwargs)
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
    @test length(cert.factors) == 2
    @test cert.raw_denominators == [factor.denominator for factor in cert.factors]
    @test cert.product_denominator == prod(cert.raw_denominators; init = one(cert.ring))
    @test cert.normalized_global_elementary_factors == collect(entry.local_evidence.factors)
    @test cert.local_product == qlfs_product(cert.normalized_global_elementary_factors, cert.ring, cert.size)
    @test cert.local_correction == entry.local_evidence.expected_product
    @test cert.verification.factor_provenance_ok
    @test cert.verification.product_denominator_ok
    @test cert.verification.overall_ok

    bad_numerator = qlfs_rebuild_sequence_certificate(cert; factors = begin
        modified = collect(cert.factors)
        modified[1] = merge(modified[1], (; numerator = modified[1].numerator + one(cert.ring)))
        modified
    end)
    @test !Suslin.verify_quillen_local_factor_sequence_certificate(bad_numerator)

    bad_denominator = qlfs_rebuild_sequence_certificate(cert; factors = begin
        modified = collect(cert.factors)
        modified[1] = merge(modified[1], (; denominator = modified[1].denominator + one(cert.ring)))
        modified
    end)
    @test !Suslin.verify_quillen_local_factor_sequence_certificate(bad_denominator)

    bad_variable = qlfs_rebuild_sequence_certificate(cert; selected_variable = first(filter(
        gen -> gen != cert.selected_variable,
        collect(gens(cert.ring)),
    )))
    @test !Suslin.verify_quillen_local_factor_sequence_certificate(bad_variable)

    bad_provenance = qlfs_rebuild_sequence_certificate(cert; factors = begin
        modified = collect(cert.factors)
        modified[1] = merge(modified[1], (; provenance = merge(modified[1].provenance, (; source = :tampered_source))))
        modified
    end)
    @test !Suslin.verify_quillen_local_factor_sequence_certificate(bad_provenance)

    reversed_order = qlfs_rebuild_sequence_certificate(cert; factors = reverse(collect(cert.factors)))
    @test !Suslin.verify_quillen_local_factor_sequence_certificate(reversed_order)
end
