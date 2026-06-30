using Test
using Suslin
using Oscar

const QMA_MURTHY_FIXTURES =
    joinpath(@__DIR__, "..", "fixtures", "sl3_murthy_gupta_cases.jl")

if !isdefined(Main, :SL3MurthyGuptaFixtureCatalog)
    include(QMA_MURTHY_FIXTURES)
end

function qma_product(factors, R, n::Int)
    product = identity_matrix(R, n)
    for factor in factors
        product *= factor
    end
    return product
end

function qma_rebuild(record; kwargs...)
    overrides = Dict{Symbol,Any}()
    for pair in kwargs
        overrides[pair.first] = pair.second
    end
    values = map(fieldnames(typeof(record))) do name
        get(overrides, name, getproperty(record, name))
    end
    return typeof(record)(values...)
end

function qma_fixture(id::AbstractString)
    catalog = SL3MurthyGuptaFixtureCatalog.catalog()
    return Dict(entry.id => entry for entry in catalog.cases)[id]
end

function qma_ordinary_adapter()
    fixture = qma_fixture("mg-q0-unit-recursion")
    certificate = Suslin.realize_sl3_local_certificate(
        fixture.entries.p,
        fixture.entries.q,
        fixture.entries.r,
        fixture.entries.s,
        fixture.variable,
    )
    adapter = Suslin._murthy_quillen_local_adapter(
        certificate,
        fixture.target,
        fixture.variable;
        witness_metadata = (;
            fixture_id = fixture.id,
            consumer_issue = 219,
        ),
    )
    return (; fixture, certificate, adapter)
end

function qma_localized_adapter(fixture_id::AbstractString)
    fixture = qma_fixture(fixture_id)
    context = Suslin.sl3_local_murthy_input_context(
        fixture.target,
        fixture.variable;
        witness = first(fixture.witnesses),
    )
    certificate = Suslin.realize_sl3_local_certificate(context)
    adapter = Suslin._murthy_quillen_local_adapter(
        certificate,
        fixture.target,
        fixture.variable;
        witness_metadata = (;
            fixture_id = fixture.id,
            consumer_issue = 219,
        ),
    )
    return (; fixture, certificate, adapter)
end

function qma_refactor_sequence_factors(sequence; leading_zero::Bool = false)
    R = sequence.ring
    refactored = Suslin.QuillenLocalElementaryFactor[]
    offset = leading_zero ? 1 : 0
    if leading_zero
        push!(
            refactored,
            Suslin.QuillenLocalElementaryFactor(
                1,
                3,
                zero(R),
                one(R),
                one(R),
                Suslin.LocalCertificate([1, 3], [one(R), one(R)]),
                (;
                    source = :murthy_local_sl3,
                    factor_index = 1,
                    murthy_denominator = one(R),
                    murthy_local_unit_witness = nothing,
                ),
                (; source = :forged_murthy_adapter_sequence),
            ),
        )
    end
    for (index, factor) in enumerate(sequence.factors)
        push!(
            refactored,
            Suslin.QuillenLocalElementaryFactor(
                factor.row,
                factor.col,
                factor.numerator,
                factor.denominator,
                factor.coverage_multiplier,
                factor.local_certificate,
                merge(factor.provenance, (; factor_index = index + offset)),
                factor.metadata,
            ),
        )
    end
    return refactored
end

@testset "Murthy Quillen adapter consumption" begin
    ordinary = qma_ordinary_adapter()
    R = base_ring(ordinary.fixture.target)
    X = ordinary.fixture.variable

    sequences = Suslin.quillen_local_sequences_from_murthy_adapters(
        ordinary.fixture.target,
        X,
        [ordinary.adapter],
    )
    @test length(sequences) == 1
    @test all(Suslin.verify_quillen_local_factor_sequence_certificate, sequences)
    @test Suslin._same_quillen_local_elementary_factors(
        sequences[1].factors,
        ordinary.adapter.quillen_factor_sequence.factors,
    )
    @test sequences[1].original_input == ordinary.fixture.target
    @test sequences[1].selected_variable == X
    @test sequences[1].raw_denominators ==
          [one(R) for _ in ordinary.adapter.local_factor_replay.factors]
    @test sequences[1].product_denominator == one(R)
    @test sequences[1].local_product == ordinary.fixture.target
    @test sequences[1].local_correction == ordinary.fixture.target
    @test sequences[1].verification.denominator_data ==
          Suslin._quillen_denominator_data(sequences[1].normalized_local_contributions)
    @test all(
        provenance -> provenance.source == :murthy_local_sl3,
        [factor.provenance for factor in sequences[1].factors],
    )

    patch = Suslin.assemble_quillen_patch_from_murthy_adapters(
        ordinary.fixture.target,
        X,
        [ordinary.adapter];
        max_exponent = 1,
        base_term_policy = :already_handled,
        metadata = (; fixture_id = ordinary.fixture.id, consumer_issue = 219),
    )
    @test patch isa Suslin.QuillenSuppliedEvidencePatchAssembly
    @test Suslin.verify_quillen_patch(patch)
    @test Suslin._same_quillen_local_factor_sequence_certificates(
        patch.local_certificates,
        sequences,
    )
    @test patch.denominator_candidate.raw_denominators == [one(R)]
    @test patch.solver_result.coverage_sum == one(R)
    @test patch.product == ordinary.fixture.target
    @test qma_product(patch.global_elementary_factors, patch.ring, patch.size) ==
          ordinary.fixture.target
    @test patch.replay_metadata.metadata.source == :quillen_murthy_adapter_consumption
    @test patch.replay_metadata.metadata.metadata ==
          (; fixture_id = ordinary.fixture.id, consumer_issue = 219)

    result = Suslin.consume_murthy_quillen_adapters_for_patch(
        ordinary.fixture.target,
        X,
        [ordinary.adapter];
        max_exponent = 1,
        base_term_policy = :already_handled,
        metadata = (; fixture_id = ordinary.fixture.id, consumer_issue = 219),
    )
    @test result isa Suslin.QuillenMurthyAdapterConsumption
    @test Suslin.verify_quillen_murthy_adapter_consumption(result)
    @test Suslin._same_quillen_local_factor_sequence_certificates(
        result.local_sequence_certificates,
        sequences,
    )
    @test Suslin._same_quillen_local_factor_sequence_certificates(
        result.patch.local_certificates,
        sequences,
    )
    @test result.patch.product == ordinary.fixture.target
    @test result.replay_metadata.murthy_adapter_metadata[1] ==
          ordinary.adapter.replay_metadata
    @test result.replay_metadata.local_sequence_metadata[1].factor_count ==
          sequences[1].replay_metadata.factor_count
    @test result.replay_metadata.local_sequence_metadata[1].raw_denominators ==
          sequences[1].replay_metadata.raw_denominators
    @test length(result.replay_metadata.local_sequence_metadata[1].denominator_data) ==
          length(sequences[1].replay_metadata.denominator_data)
    @test all(eachindex(sequences[1].replay_metadata.denominator_data)) do idx
        left = result.replay_metadata.local_sequence_metadata[1].denominator_data[idx]
        right = sequences[1].replay_metadata.denominator_data[idx]
        left.denominator == right.denominator &&
            left.coverage_multiplier == right.coverage_multiplier
    end
    @test result.replay_metadata.local_sequence_metadata[1].factor_provenance ==
          sequences[1].replay_metadata.factor_provenance
    @test result.replay_metadata.local_sequence_metadata[1].witness_metadata ==
          sequences[1].replay_metadata.witness_metadata
    @test result.replay_metadata.patch_metadata == result.patch.replay_metadata

    tampered_consumption = qma_rebuild(
        result;
        replay_metadata = (; source = :bad_metadata),
    )
    @test !Suslin.verify_quillen_murthy_adapter_consumption(tampered_consumption)
    bad_patch_metadata = merge(
        result.patch.replay_metadata,
        (; metadata = (; source = :not_murthy_adapter_consumption)),
    )
    bad_patch = qma_rebuild(result.patch; replay_metadata = bad_patch_metadata)
    bad_patch = qma_rebuild(
        bad_patch;
        verification = Suslin.replay_quillen_supplied_evidence_patch(bad_patch),
    )
    @test Suslin.verify_quillen_patch(bad_patch)
    matching_bad_replay_metadata = merge(
        result.replay_metadata,
        (; patch_metadata = bad_patch.replay_metadata),
    )
    matching_bad_patch_consumption = qma_rebuild(
        result;
        patch = bad_patch,
        replay_metadata = matching_bad_replay_metadata,
    )
    matching_bad_patch_consumption = qma_rebuild(
        matching_bad_patch_consumption;
        verification = Suslin.replay_quillen_murthy_adapter_consumption(
            matching_bad_patch_consumption,
        ),
    )
    @test !Suslin.verify_quillen_murthy_adapter_consumption(
        matching_bad_patch_consumption,
    )

    forged_sequence = Suslin.quillen_local_factor_sequence_certificate(
        ordinary.fixture.target,
        X;
        factors = qma_refactor_sequence_factors(
            ordinary.adapter.quillen_factor_sequence;
            leading_zero = true,
        ),
        local_correction = ordinary.fixture.target,
        witness_metadata = ordinary.adapter.witness_metadata,
        local_evidence = (;
            source = :forged_cached_quillen_sequence,
            expected_product = ordinary.fixture.target,
        ),
    )
    @test Suslin.verify_quillen_local_factor_sequence_certificate(forged_sequence)
    @test length(forged_sequence.factors) ==
          length(ordinary.adapter.quillen_factor_sequence.factors) + 1
    forged_sequence_adapter = qma_rebuild(
        ordinary.adapter;
        quillen_factor_sequence = forged_sequence,
    )
    forged_sequence_adapter = qma_rebuild(
        forged_sequence_adapter;
        verification = Suslin._murthy_quillen_local_adapter_summary(
            forged_sequence_adapter,
        ),
    )
    @test Suslin._verify_murthy_quillen_local_adapter(forged_sequence_adapter)
    @test_throws ArgumentError Suslin.quillen_local_sequences_from_murthy_adapters(
        ordinary.fixture.target,
        X,
        [forged_sequence_adapter],
    )
    forged_consumption = qma_rebuild(
        result;
        murthy_adapters = [forged_sequence_adapter],
    )
    forged_consumption_replay =
        Suslin.replay_quillen_murthy_adapter_consumption(forged_consumption)
    @test !forged_consumption_replay.local_sequences_ok
    @test !forged_consumption_replay.overall_ok
    @test !Suslin.verify_quillen_murthy_adapter_consumption(
        qma_rebuild(result; size = 2),
    )

    bad_factor_replay = qma_rebuild(
        ordinary.adapter.local_factor_replay;
        factors = reverse(ordinary.adapter.local_factor_replay.factors),
    )
    bad_factor_adapter = qma_rebuild(
        ordinary.adapter;
        local_factor_replay = bad_factor_replay,
    )
    @test !Suslin._verify_murthy_quillen_local_adapter(bad_factor_adapter)
    @test_throws ArgumentError Suslin.quillen_local_sequences_from_murthy_adapters(
        ordinary.fixture.target,
        X,
        [bad_factor_adapter],
    )

    @test_throws ArgumentError Suslin.quillen_local_sequences_from_murthy_adapters(
        ordinary.fixture.target,
        one(R),
        [ordinary.adapter],
    )

    tampered_witness_factors = copy(ordinary.adapter.local_factor_replay.factors)
    tampered_record = tampered_witness_factors[1]
    tampered_witness_factors[1] = Suslin.SL3LocalElementaryFactor(
        tampered_record.R,
        tampered_record.n,
        tampered_record.row,
        tampered_record.col,
        tampered_record.numerator,
        tampered_record.denominator,
        tampered_record.selected_variable,
        (; tampered = true),
    )
    tampered_replay = qma_rebuild(
        ordinary.adapter.local_factor_replay;
        factors = tampered_witness_factors,
    )
    tampered_adapter = qma_rebuild(
        ordinary.adapter;
        local_factor_replay = tampered_replay,
    )
    @test !Suslin._verify_murthy_quillen_local_adapter(tampered_adapter)
    @test_throws ArgumentError Suslin.quillen_local_sequences_from_murthy_adapters(
        ordinary.fixture.target,
        X,
        [tampered_adapter],
    )

    wrong_target =
        ordinary.fixture.target * elementary_matrix(3, 1, 2, one(R), R)
    @test_throws ArgumentError Suslin.quillen_local_sequences_from_murthy_adapters(
        wrong_target,
        X,
        [ordinary.adapter],
    )

    localized = qma_localized_adapter("mg-local-q0-unit-at-u")
    @test localized.adapter.mode == :localized_replay_handoff
    @test_throws ArgumentError Suslin.quillen_local_sequences_from_murthy_adapters(
        localized.fixture.target,
        localized.fixture.variable,
        [localized.adapter],
    )

    @test_throws ArgumentError Suslin.assemble_quillen_patch_from_murthy_adapters(
        ordinary.fixture.target,
        X,
        [ordinary.adapter];
        max_exponent = 1,
    )
end
