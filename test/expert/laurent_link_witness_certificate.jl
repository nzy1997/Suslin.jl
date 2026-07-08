using Test
using Suslin
using Oscar

if !(@isdefined case008_d14_laurent_link_witness_search_report)
    include(joinpath(@__DIR__, "case008_d14_laurent_link_witness_search.jl"))
end

@testset "Laurent link-witness certificate shell" begin
    runtests = read(joinpath(@__DIR__, "..", "runtests.jl"), String)
    @test occursin(
        "\"expert/laurent_link_witness_certificate.jl\"",
        runtests,
    )

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

    cert = laurent_link_witness_certificate(context, witness, column, R)
    @test cert.case_id == "synthetic"
    @test cert.dimension == 2
    @test cert.ring_generators == ("u", "v")
    @test cert.context_status == :link_witness_context
    @test cert.witness == witness
    @test cert.source_endpoint.leading_exponent == (1, 0)
    @test cert.target_endpoint.leading_exponent == (0, 1)
    @test cert.replay_status == :ok
    @test cert.identity_status == :verified
    @test cert.next_boundary == :laurent_endpoint_reduction
    @test cert.status == :link_witness_certificate
    @test validate_laurent_link_witness_certificate(cert, column, R) == :ok
end

@testset "case_008 d=14 Laurent link-witness certificate summary" begin
    summary = case008_d14_laurent_link_witness_certificate_summary()
    @test summary.case_id == "case_008"
    @test summary.dimension == 14
    @test summary.search_status in (:candidate_found, :exhausted)
    @test summary.d14_status in (:link_witness_certificate, :no_link_witness_candidate)

    if summary.search_status == :candidate_found
        @test summary.d14_status == :link_witness_certificate
        @test summary.certificate.case_id == "case_008"
        @test validate_laurent_link_witness_certificate(
            summary.certificate,
            summary.source_column,
            summary.ring,
        ) == :ok
    else
        @test summary.d14_status == :no_link_witness_candidate
        @test !hasproperty(summary, :certificate)
    end
end
