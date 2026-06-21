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

struct Issue41DetErrorMatrix
    ring
    error
end

Oscar.base_ring(A::Issue41DetErrorMatrix) = A.ring
Oscar.det(A::Issue41DetErrorMatrix) = throw(A.error)

function _issue41_assert_issue38_diagnostic(core)
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

    row_diagnostic = _issue41_assert_issue38_diagnostic(row_core)
    column_diagnostic = _issue41_assert_issue38_diagnostic(column_core)

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

@testset "Polynomial determinant precondition failure uses stable diagnostic enums" begin
    R, (X,) = Oscar.polynomial_ring(QQ, ["X"])
    A = diagonal_matrix(R, [X, one(R), one(R)])

    diagnostic = diagnose_sln_to_sl3_reduction(A)
    @test diagnostic isa SLNToSL3ReductionDiagnostic
    @test diagnostic.status == :failure
    @test diagnostic.failure_code == :determinant_failure
    @test diagnostic.determinant_status == :determinant_not_one
    @test diagnostic.block_diagnostics == SL3LocalReductionDiagnostic[]
    @test diagnostic.partition_search.status == :not_applicable
end

@testset "Laurent unsupported determinant normalization returns diagnostic instead of throwing" begin
    L, (t,) = suslin_laurent_polynomial_ring(QQ, ["t"])
    A = diagonal_matrix(L, [one(L) + t, one(L), one(L)])

    diagnostic = diagnose_sln_to_sl3_reduction(A)
    @test diagnostic isa SLNToSL3ReductionDiagnostic
    @test diagnostic.status == :failure
    @test diagnostic.failure_code == :determinant_failure
    @test diagnostic.determinant_status == :determinant_not_one
    @test diagnostic.determinant_classification == :non_unit
    @test isempty(diagnostic.block_diagnostics)
    @test diagnostic.partition_search.searched == false
    @test diagnostic.partition_search.status == :not_applicable
    @test diagnostic.partition_search.attempted_count == 0
    @test diagnostic.message !== nothing
    @test occursin("unsupported Laurent GL_n determinant", diagnostic.message)
end

@testset "Laurent determinant correction diagnostic uses stable failure code" begin
    L, (t,) = suslin_laurent_polynomial_ring(QQ, ["t"])
    A = diagonal_matrix(L, [t, one(L), one(L)])

    diagnostic = diagnose_sln_to_sl3_reduction(A)
    @test diagnostic.status == :failure
    @test diagnostic.failure_code == :determinant_failure
    @test diagnostic.determinant_status == :determinant_requires_correction
    @test diagnostic.determinant_classification == :laurent_monomial_unit
    @test occursin("Laurent determinant correction", diagnostic.message)
end

@testset "Diagnostic determinant fallback helpers are stable" begin
    R, (X,) = Oscar.polynomial_ring(QQ, ["X"])

    status, classification = Suslin._sln_determinant_status(
        Issue41DetErrorMatrix(R, ArgumentError("issue 41 determinant sentinel")),
        :polynomial,
        nothing,
    )
    @test status == :determinant_check_failed
    @test classification === nothing

    @test Suslin._sl3_local_determinant_status(
        Issue41DetErrorMatrix(R, ArgumentError("issue 41 local determinant sentinel")),
        R,
    ) == :determinant_check_failed

    @test Suslin._diagnostic_message_for_determinant_failure(:polynomial, :determinant_check_failed, nothing) ==
          "determinant check failed for the staged SL_n to local SL_3 reduction path"
    @test Suslin._determinant_status_from_laurent_classification(:one) == :determinant_one
    @test Suslin._determinant_status_from_laurent_classification(:laurent_monomial_unit) == :determinant_requires_correction

    @test_throws InterruptException Suslin._sln_determinant_status(
        Issue41DetErrorMatrix(R, InterruptException()),
        :polynomial,
        nothing,
    )
    @test_throws InterruptException Suslin._sl3_local_determinant_status(
        Issue41DetErrorMatrix(R, InterruptException()),
        R,
    )
end

@testset "Local determinant failure reports local solver diagnostic" begin
    R, (X,) = Oscar.polynomial_ring(QQ, ["X"])
    A = identity_matrix(R, 6)
    A[1:3, 1:3] = matrix(R, [
        X zero(R) zero(R);
        zero(R) one(R) zero(R);
        zero(R) zero(R) one(R)
    ])

    diagnostic = Suslin._diagnose_sl3_local_obligation(A, R, [1, 2, 3], X, :polynomial)
    @test diagnostic.status == :failure
    @test diagnostic.failure_code == :local_solver_failure
    @test diagnostic.determinant_status == :determinant_not_one
    @test diagnostic.local_shape_reason == :embedded_2x2_with_trailing_identity
    @test diagnostic.solver_status == :not_attempted
end

@testset "Embedded q-unit local solver success uses stable diagnostic enums" begin
    R, (X,) = Oscar.polynomial_ring(QQ, ["X"])
    A = identity_matrix(R, 6)
    A[1:3, 1:3] = matrix(R, [
        X -one(R) zero(R);
        one(R) zero(R) zero(R);
        zero(R) zero(R) one(R)
    ])

    diagnostic = diagnose_sln_to_sl3_reduction(A; search_partitions = false)
    @test diagnostic isa SLNToSL3ReductionDiagnostic
    @test diagnostic.status == :success
    @test diagnostic.failure_code === nothing
    @test diagnostic.determinant_status == :determinant_one
    @test diagnostic.partition_search.status == :not_applicable

    failures = _issue41_failure_diagnostics(diagnostic)
    @test isempty(failures)
    @test length(diagnostic.block_diagnostics) == 1
    success = only(diagnostic.block_diagnostics)
    @test success.failure_code === nothing
    @test success.local_shape_reason == :embedded_2x2_with_trailing_identity
    @test success.solver_status == :success
    @test success.message === nothing
end

@testset "Embedded q0-nonunit local solver failure uses stable diagnostic enums" begin
    R, (X,) = Oscar.polynomial_ring(QQ, ["X"])
    A = identity_matrix(R, 6)
    A[1:3, 1:3] = matrix(R, [
        X^2 + one(R) X zero(R);
        X^2 + X + one(R) X + one(R) zero(R);
        zero(R) zero(R) one(R)
    ])

    diagnostic = diagnose_sln_to_sl3_reduction(A; search_partitions = false)
    @test diagnostic isa SLNToSL3ReductionDiagnostic
    @test diagnostic.status == :failure
    @test diagnostic.failure_code == :local_solver_failure
    @test diagnostic.determinant_status == :determinant_one
    @test diagnostic.partition_search.status == :disabled

    failures = _issue41_failure_diagnostics(diagnostic)
    @test length(failures) == 1
    failure = only(failures)
    @test failure.failure_code == :local_solver_failure
    @test failure.local_shape_reason == :embedded_2x2_with_trailing_identity
    @test failure.solver_status == :failure
    @test failure.message !== nothing
    @test occursin("Murthy q(0)-nonunit", failure.message)
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
