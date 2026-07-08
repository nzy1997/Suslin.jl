using Test
using Suslin
using Oscar

if !(@isdefined case008_d14_laurent_link_witness_search_report)
    include(joinpath(@__DIR__, "..", "expert", "case008_d14_laurent_link_witness_search.jl"))
end

function _internal_link_without_field(value::NamedTuple, field::Symbol)
    kept = tuple((name for name in keys(value) if name != field)...)
    return NamedTuple{kept}(tuple((getproperty(value, name) for name in kept)...))
end

@testset "internal Laurent link-witness helpers" begin
    runtests = read(joinpath(@__DIR__, "..", "runtests.jl"), String)
    @test occursin("\"internal/laurent_link_witness_helpers.jl\"", runtests)

    R, (u, v) = Suslin.suslin_laurent_polynomial_ring(GF(2), ["u", "v"])
    column = [one(R), u + v]
    context = (;
        case_id = "synthetic",
        dimension = 2,
        ring_generators = ("u", "v"),
        status = :link_witness_context,
    )
    witness = (;
        family = :two_entry_laurent_combination,
        pivot_index = 2,
        partner_index = 1,
        coefficient = 1,
        exponent = (1, 0),
        ring_generators = ("u", "v"),
    )

    @test Suslin._laurent_link_witness_status(witness, length(column), R) == :ok
    candidate = Suslin._laurent_link_witness_candidate_from_replay(
        column,
        R,
        witness;
        case_id = context.case_id,
    )
    @test Suslin._verify_laurent_link_witness_candidate(column, R, candidate)
    @test candidate.source_endpoint.leading_exponent == (1, 0)
    @test candidate.target_endpoint.leading_exponent == (0, 1)
    @test candidate.measure_relation == :strict_decrease

    cert = Suslin._laurent_link_witness_certificate_from_replay(
        context,
        witness,
        column,
        R,
    )
    @test cert.case_id == "synthetic"
    @test cert.next_boundary == :laurent_endpoint_reduction
    @test cert.status == :link_witness_certificate
    @test Suslin._validate_laurent_link_witness_certificate(cert, column, R) == :ok

    missing_field = _internal_link_without_field(witness, :family)
    @test Suslin._laurent_link_witness_status(missing_field, length(column), R) ==
          :malformed_witness
    missing_field_cert = merge(cert, (; witness = missing_field))
    @test Suslin._validate_laurent_link_witness_certificate(
        missing_field_cert,
        column,
        R,
    ) == :identity_replay_failed

    equal_indices = merge(witness, (; partner_index = 2))
    @test Suslin._laurent_link_witness_status(equal_indices, length(column), R) ==
          :pivot_partner_index_equality
    @test !Suslin._verify_laurent_link_witness_candidate(
        column,
        R,
        merge(candidate, (; witness = equal_indices)),
    )

    wrong_generators = merge(witness, (; ring_generators = ("v", "u")))
    @test Suslin._laurent_link_witness_status(wrong_generators, length(column), R) ==
          :wrong_ring_generators

    stale_source = merge(
        cert,
        (;
            source_endpoint = merge(
                cert.source_endpoint,
                (; term_count = cert.source_endpoint.term_count + 1),
            ),
        ),
    )
    @test Suslin._validate_laurent_link_witness_certificate(
        stale_source,
        column,
        R,
    ) == :stale_source_endpoint

    stale_target = merge(
        cert,
        (;
            target_endpoint = merge(
                cert.target_endpoint,
                (; term_count = cert.target_endpoint.term_count + 1),
            ),
        ),
    )
    @test Suslin._validate_laurent_link_witness_certificate(
        stale_target,
        column,
        R,
    ) == :stale_target_endpoint

    coefficient_tampering = merge(
        cert,
        (; witness = merge(witness, (; coefficient = 0))),
    )
    @test Suslin._validate_laurent_link_witness_certificate(
        coefficient_tampering,
        column,
        R,
    ) == :identity_replay_failed

    nonidentity = merge(
        cert,
        (; witness = merge(witness, (; exponent = (0, 0)))),
    )
    @test Suslin._validate_laurent_link_witness_certificate(
        nonidentity,
        column,
        R,
    ) == :identity_replay_failed
end

@testset "internal d14 Laurent link-witness candidate" begin
    fixture = ToricBuilderCase008D14ColumnBoundary.boundary_fixture()
    source = _case008_d14_link_witness_source_data(fixture)
    report = case008_d14_laurent_link_witness_search_report(fixture)
    @test validate_case008_d14_laurent_link_witness_search_report(report, fixture) == :ok
    @test report.status == :candidate_found
    @test report.next_boundary == :laurent_link_witness_certificate

    candidate = first(report.candidates)
    witness = candidate.witness
    @test candidate.source_endpoint.case_id == "case_008"
    @test candidate.target_endpoint.case_id == "case_008"
    @test candidate.replay_status == :ok
    @test candidate.identity_status == :verified
    @test candidate.measure_relation == :strict_decrease
    @test witness.pivot_index == 10
    @test witness.partner_index == 1
    @test witness.coefficient == 1
    @test witness.exponent == (1, -1)

    @test Suslin._verify_laurent_link_witness_candidate(
        source.replay.after_column,
        fixture.ring,
        candidate,
    )
    cert = Suslin._laurent_link_witness_certificate_from_replay(
        source.context,
        witness,
        source.replay.after_column,
        fixture.ring,
    )
    @test cert.case_id == "case_008"
    @test cert.next_boundary == :laurent_endpoint_reduction
    @test Suslin._validate_laurent_link_witness_certificate(
        cert,
        source.replay.after_column,
        fixture.ring,
    ) == :ok
end
