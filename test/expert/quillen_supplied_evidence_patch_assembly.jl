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
