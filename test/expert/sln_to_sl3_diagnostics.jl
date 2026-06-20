using Test
using Suslin
using Oscar

include("../fixtures/toricbuilder_issue38_cases.jl")
include("../fixtures/laurent_large_acceptance_cases.jl")

function _issue41_captured_error(f)
    try
        f()
        return nothing
    catch err
        return err
    end
end

function _issue41_failure_diagnostics(diagnostic)
    return filter(block -> block.status == :failure, diagnostic.block_diagnostics)
end

function _issue41_assert_issue38_diagnostic(core, label::AbstractString)
    diagnostic = diagnose_sln_to_sl3_reduction(core)
    @test diagnostic isa SLNToSL3ReductionDiagnostic
    @test diagnostic.status == :failure
    @test diagnostic.failure_code == :local_shape_failure
    @test diagnostic.determinant_status == :determinant_one
    @test diagnostic.determinant_classification == :one

    failures = _issue41_failure_diagnostics(diagnostic)
    @test length(failures) >= 1
    failure = first(failures)
    @test failure isa SL3LocalReductionDiagnostic
    @test failure.block_location == [1, 2, 3]
    @test failure.failure_code == :local_shape_failure
    @test failure.local_shape_reason == :not_embedded_2x2_with_trailing_identity
    @test failure.solver_status == :not_attempted

    @test diagnostic.partition_search.searched
    @test diagnostic.partition_search.status == :no_success
    @test diagnostic.partition_search.attempted_count == 10
    @test isempty(diagnostic.partition_search.successful_partitions)

    err = _issue41_captured_error(() -> reduce_sln_to_sl3(core))
    @test err isa ArgumentError
    message = sprint(showerror, err)
    @test occursin("staged SL_n to local SL_3 reduction failure", message)
    @test occursin("failed to solve local SL_3 obligation on block [1, 2, 3]", message)

    return diagnostic
end

@testset "Issue 38 Laurent core diagnostics" begin
    entry = only(ToricBuilderIssue38Cases.catalog().cases)
    row_core = entry.normalizations.row.core
    column_core = entry.normalizations.column.core

    row_diagnostic = _issue41_assert_issue38_diagnostic(row_core, "row")
    column_diagnostic = _issue41_assert_issue38_diagnostic(column_core, "column")

    @test row_diagnostic.message !== nothing
    @test column_diagnostic.message !== nothing
end

@testset "Issue 38 partition search remains disabled when requested off" begin
    entry = only(ToricBuilderIssue38Cases.catalog().cases)
    issue38_core = entry.normalizations.row.core

    diagnostic = diagnose_sln_to_sl3_reduction(issue38_core; search_partitions = false)
    @test diagnostic.status == :failure
    @test diagnostic.partition_search.searched == false
    @test diagnostic.partition_search.status == :disabled
    @test diagnostic.partition_search.attempted_count == 0
end

@testset "Reassembly failure does not trigger partition search" begin
    R, (X,) = Oscar.polynomial_ring(QQ, ["X"])
    A = identity_matrix(R, 6)
    A[1, 4] = X

    diagnostic = diagnose_sln_to_sl3_reduction(A)
    @test diagnostic.status == :failure
    @test diagnostic.failure_code == :reassembly_failure
    @test diagnostic.partition_search.searched == false
    @test diagnostic.partition_search.status == :not_applicable
    @test diagnostic.partition_search.attempted_count == 0
end

@testset "Supported Laurent block-local diagnostics" begin
    catalog = LaurentLargeAcceptanceCases.acceptance_catalog()
    case = only(filter(entry -> entry.id == "laurent-block-local-40x40", catalog.cases))

    diagnostic = diagnose_sln_to_sl3_reduction(case.matrix; block_locations = case.block_locations)
    @test diagnostic isa SLNToSL3ReductionDiagnostic
    @test diagnostic.status == :success
    @test diagnostic.failure_code === nothing
    @test diagnostic.determinant_status == :determinant_one
    @test isempty(_issue41_failure_diagnostics(diagnostic))
    @test !diagnostic.partition_search.searched
    @test diagnostic.partition_search.status == :not_applicable
end
