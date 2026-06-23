using Test
using Suslin
using Oscar

include(joinpath(@__DIR__, "..", "fixtures", "quillen_patch_cases.jl"))

function quillen_local_test_product(factors, R, n)
    product = identity_matrix(R, n)
    for factor in factors
        product *= factor
    end
    return product
end

function local_certificate_from_fixture(local_factor)
    return Suslin.LocalCertificate(
        local_factor.certificate.indices,
        local_factor.certificate.denominators,
    )
end

function local_correction_from_fixture(local_factor)
    correction = local_factor.correction
    return Suslin.QuillenElementaryCorrection(
        correction.row,
        correction.col,
        correction.entry,
    )
end

function realization_certificate_from_fixture(entry; local_index::Int = 1)
    local_factor = entry.local_factors[local_index]
    return Suslin.quillen_local_realization_certificate(
        entry.target_matrix,
        entry.substitution_variable;
        local_certificate = local_certificate_from_fixture(local_factor),
        denominator = local_factor.denominator,
        coverage_multiplier = local_factor.coverage_multiplier,
        correction = local_correction_from_fixture(local_factor),
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

function correction_realization_certificate_from_fixture(entry; local_index::Int = 1)
    local_factor = entry.local_factors[local_index]
    correction = local_correction_from_fixture(local_factor)
    return Suslin.quillen_local_realization_certificate(
        correction,
        entry.substitution_variable;
        local_certificate = local_certificate_from_fixture(local_factor),
        denominator = local_factor.denominator,
        coverage_multiplier = local_factor.coverage_multiplier,
        correction = correction,
        local_correction = local_factor.factor,
        ring = entry.ring.object,
        size = entry.size,
        witness_metadata = (;
            fixture_id = entry.id,
            local_index = local_index,
            input_kind = :elementary_correction,
        ),
    )
end

function rebuild_local_certificate(cert; kwargs...)
    fields = merge((
        original_input = cert.original_input,
        ring = cert.ring,
        size = cert.size,
        selected_variable = cert.selected_variable,
        local_certificate = cert.local_certificate,
        denominator = cert.denominator,
        coverage_multiplier = cert.coverage_multiplier,
        correction = cert.correction,
        factors = cert.factors,
        local_product = cert.local_product,
        local_correction = cert.local_correction,
        patched_substitution_witness = cert.patched_substitution_witness,
        witness_metadata = cert.witness_metadata,
        verification = cert.verification,
    ), kwargs)
    return Suslin.QuillenLocalRealizationCertificate(
        fields.original_input,
        fields.ring,
        fields.size,
        fields.selected_variable,
        fields.local_certificate,
        fields.denominator,
        fields.coverage_multiplier,
        fields.correction,
        fields.factors,
        fields.local_product,
        fields.local_correction,
        fields.patched_substitution_witness,
        fields.witness_metadata,
        fields.verification,
    )
end

@testset "Quillen local realization certificates" begin
    entries = QuillenPatchFixtureCatalog.cases_by_id()
    fixture_ids = [
        "quillen-two-open-cover-qq",
        "quillen-supplied-local-certificate-gf2",
        "quillen-patched-substitution-witness-qq",
    ]

    certs = [realization_certificate_from_fixture(entries[id]) for id in fixture_ids]
    @test all(cert -> cert isa Suslin.QuillenLocalRealizationCertificate, certs)

    for cert in certs
        @test Suslin.verify_quillen_local_certificate(cert)
        @test quillen_local_test_product(cert.factors, cert.ring, cert.size) == cert.local_correction
        @test cert.local_product == cert.local_correction
        @test cert.verification.local_product == cert.local_correction
        @test cert.verification.local_correction_ok
        @test cert.verification.factors_ok
        @test cert.verification.overall_ok
        @test Suslin._coerce_into_ring(cert.ring, cert.denominator, "test denominator") ==
              cert.denominator
        for denominator in cert.local_certificate.denominators
            @test Suslin._coerce_into_ring(
                cert.ring,
                denominator,
                "test certificate denominator",
            ) == denominator
        end
    end

    patched_cert = certs[3]
    @test patched_cert.patched_substitution_witness !== nothing
    @test patched_cert.verification.patched_substitution_ok
    witness = patched_cert.patched_substitution_witness
    @test Suslin.patched_substitution(
        witness.matrix,
        witness.variable,
        witness.denominator,
        witness.exponent,
        witness.shift,
    ) == witness.expected_matrix

    correction_cert = correction_realization_certificate_from_fixture(entries[fixture_ids[1]])
    @test Suslin.verify_quillen_local_certificate(correction_cert)
    @test correction_cert.original_input == correction_cert.correction
    @test correction_cert.factors == [correction_cert.local_correction]
    @test correction_cert.local_product == correction_cert.local_correction
    @test correction_cert.verification.original_input == correction_cert.correction
    @test correction_cert.verification.witness_metadata.input_kind == :elementary_correction

    local_factor = entries[fixture_ids[1]].local_factors[1]
    correction_input = local_correction_from_fixture(local_factor)
    @test_throws ArgumentError Suslin.quillen_local_realization_certificate(
        correction_input,
        entries[fixture_ids[1]].substitution_variable;
        local_certificate = local_certificate_from_fixture(local_factor),
        denominator = local_factor.denominator,
        coverage_multiplier = local_factor.coverage_multiplier,
        correction = correction_input,
        local_correction = local_factor.factor,
        size = entries[fixture_ids[1]].size,
    )
    @test_throws ArgumentError Suslin.quillen_local_realization_certificate(
        correction_input,
        entries[fixture_ids[1]].substitution_variable;
        local_certificate = local_certificate_from_fixture(local_factor),
        denominator = local_factor.denominator,
        coverage_multiplier = local_factor.coverage_multiplier,
        correction = correction_input,
        local_correction = local_factor.factor,
        ring = entries[fixture_ids[1]].ring.object,
    )
    @test_throws ArgumentError Suslin.quillen_local_realization_certificate(
        correction_input,
        entries[fixture_ids[1]].substitution_variable;
        local_certificate = local_certificate_from_fixture(local_factor),
        denominator = local_factor.denominator,
        coverage_multiplier = local_factor.coverage_multiplier,
        correction = correction_input,
        local_correction = local_factor.factor,
        ring = entries[fixture_ids[1]].ring.object,
        size = 1,
    )
    @test_throws ArgumentError Suslin.quillen_local_realization_certificate(
        entries[fixture_ids[1]].target_matrix,
        entries[fixture_ids[1]].substitution_variable;
        local_certificate = local_certificate_from_fixture(local_factor),
        denominator = local_factor.denominator,
        coverage_multiplier = local_factor.coverage_multiplier,
        correction = correction_input,
    )

    base_cert = certs[1]
    R = base_cert.ring
    bad_factor = base_cert.factors[1] *
        elementary_matrix(base_cert.size, 1, 3, one(R), R)
    @test !Suslin.verify_quillen_local_certificate(rebuild_local_certificate(base_cert; factors = [bad_factor]))

    @test !Suslin.verify_quillen_local_certificate(rebuild_local_certificate(
        base_cert;
        original_input = identity_matrix(R, base_cert.size),
    ))

    @test !Suslin.verify_quillen_local_certificate(rebuild_local_certificate(
        base_cert;
        denominator = base_cert.denominator + one(R),
    ))

    generators = collect(gens(R))
    replacement_variable = first(filter(gen -> gen != base_cert.selected_variable, generators))
    @test !Suslin.verify_quillen_local_certificate(rebuild_local_certificate(
        base_cert;
        selected_variable = replacement_variable,
    ))

    patched_cert = certs[3]
    bad_witness = merge(
        patched_cert.patched_substitution_witness,
        (;
            exponent = patched_cert.patched_substitution_witness.exponent + 1,
        ),
    )
    @test !Suslin.verify_quillen_local_certificate(rebuild_local_certificate(
        patched_cert;
        patched_substitution_witness = bad_witness,
    ))

    @test !Suslin.verify_quillen_local_certificate(rebuild_local_certificate(
        base_cert;
        witness_metadata = merge(base_cert.witness_metadata, (; local_index = 99)),
    ))
end
