using Test
using Suslin
using Oscar

if !(@isdefined case008_d14_laurent_link_witness_context)
    include(joinpath(@__DIR__, "case008_d14_laurent_link_witness_context.jl"))
end

@testset "case_008 d=14 Laurent link-witness search report" begin
    runtests = read(joinpath(@__DIR__, "..", "runtests.jl"), String)
    @test occursin(
        "\"expert/case008_d14_laurent_link_witness_search.jl\"",
        runtests,
    )

    report = case008_d14_laurent_link_witness_search_report()
    @test report.case_id == "case_008"
    @test report.dimension == 14
    @test report.context_status == :link_witness_context
    @test report.witness_families == (:two_entry_laurent_combination,)
    @test report.witness_semantics == :pivot_plus_shifted_partner_endpoint
    @test report.pivot_index == 10
    @test report.partner_indices ==
          (1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 14)
    @test report.exponent_radius == 1
    @test report.exponent_vectors ==
          ((-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 0), (0, 1), (1, -1), (1, 0), (1, 1))
    @test report.coefficient_family == (1,)
    @test report.checked_candidate_count == 117
    @test report.status in (:candidate_found, :exhausted)
    @test report.candidate_count == length(report.candidates)
    @test report.replay_verified_count == report.candidate_count
    @test validate_case008_d14_laurent_link_witness_search_report(report) == :ok

    if report.status == :candidate_found
        @test report.candidate_count > 0
        @test report.next_boundary == :laurent_link_witness_certificate
        @test all(candidate -> candidate.replay_status == :ok, report.candidates)
        @test all(candidate -> candidate.identity_status == :verified, report.candidates)
    else
        @test report.candidate_count == 0
        @test report.replay_verified_count == 0
        @test isempty(report.candidates)
        @test report.next_boundary == :laurent_link_witness_search_expansion
    end
end
