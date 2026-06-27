using Dates
using Test

const TORICBUILDER_CACHE_STATUS_REPORT_SCRIPT =
    joinpath(@__DIR__, "..", "..", "scripts", "report_toricbuilder_cache_q_blocks.jl")
const TORICBUILDER_CACHE_STATUS_REPORT_PATH =
    joinpath(@__DIR__, "..", "..", "docs", "audits", "2026-06-24-toricbuilder-cache-q-block-status.md")

struct FakeBoundedEntry
    id::String
    dimensions::NamedTuple{(:matrix,), Tuple{Tuple{Int, Int}}}
    sparse_entry_count::Int
    expected_test_level::Symbol
    mode::Symbol
end

FakeBoundedEntry(id::AbstractString; matrix = (5, 5), sparse_entry_count = 25, mode::Symbol) =
    FakeBoundedEntry(
        String(id),
        (; matrix),
        sparse_entry_count,
        :default_contract,
        mode,
    )

const FORCE_WAIT_FOR_EXIT_FAILURE = Ref(false)

function _matching_stage_names(row, status::Symbol)
    return [
        stage for stage in ToricBuilderCacheQBlockStatusReport.STAGE_NAMES if
        getproperty(row.stage_timings, stage).status == status
    ]
end

function _stage_has_numeric_elapsed(row, stage::Symbol)
    return getproperty(row.stage_timings, stage).elapsed_seconds isa Number
end

function _stage_has_stable_error(row, stage::Symbol)
    details = getproperty(row.stage_timings, stage).error_details
    return details isa AbstractString && details != "none" && details != "not_run"
end

function _bounded_route_row_is_structured(row; timeout_seconds)
    if row.route_status == :gl_certificate_pass
        return row.verified == true && row.error_details == "none"
    elseif row.route_status == :certified_algorithm_boundary
        stages = _matching_stage_names(row, :certified_algorithm_boundary)
        length(stages) == 1 || return false
        stage = only(stages)
        return _stage_has_numeric_elapsed(row, stage) &&
            _stage_has_stable_error(row, stage) &&
            occursin(string(stage), row.evidence)
    elseif row.route_status == :timed_out
        stages = _matching_stage_names(row, :timed_out)
        length(stages) == 1 || return false
        stage = only(stages)
        budget_text = "timed out after $(ToricBuilderCacheQBlockStatusReport._runtime_text(timeout_seconds)) seconds"
        timing = getproperty(row.stage_timings, stage)
        return _stage_has_numeric_elapsed(row, stage) &&
            timing.error_details isa AbstractString &&
            occursin(budget_text, timing.error_details) &&
            row.error_details isa AbstractString &&
            occursin(budget_text, row.error_details) &&
            occursin(string(stage), row.error_details)
    end
    return false
end

function _bounded_route_row_has_issue147_peel_evidence(row; timeout_seconds)
    _bounded_route_row_is_structured(row; timeout_seconds) || return false
    row.route_status == :timed_out || return true
    stages = _matching_stage_names(row, :timed_out)
    length(stages) == 1 || return false
    stage = only(stages)
    stage != :certificate_construction && return true
    details = row.error_details
    details isa AbstractString || return false
    return occursin("peel progress:", details) &&
        occursin("current d=", details) &&
        occursin("completed steps=", details) &&
        (occursin("last-column nnz=", details) || occursin("max entry terms=", details))
end

@testset "ToricBuilder cache Q-block status report" begin
    @test isfile(TORICBUILDER_CACHE_STATUS_REPORT_SCRIPT)

    include(TORICBUILDER_CACHE_STATUS_REPORT_SCRIPT)
    @test !isdefined(ToricBuilderCacheQBlockStatusReport, :_peel_column_stats)
    @test !isdefined(ToricBuilderCacheQBlockStatusReport, :_last_column_or_nothing)
    @test !isdefined(ToricBuilderCacheQBlockStatusReport, :_worker_laurent_column_peel_recursive)
    @test !isdefined(ToricBuilderCacheQBlockStatusReport, :_worker_factor_laurent_sl_column_peel)
    @test !isdefined(ToricBuilderCacheQBlockStatusReport, :_worker_laurent_gl_factorization_certificate)
    @eval ToricBuilderCacheQBlockStatusReport begin
        const _TEST_FORCE_WAIT_FOR_EXIT_FAILURE = Main.FORCE_WAIT_FOR_EXIT_FAILURE

        function _worker_command(entry::Main.FakeBoundedEntry, progress_path, result_path = nothing)
            project_path = dirname(Base.active_project())
            script = if entry.mode == :invalid_stdout
                """
                write(stderr, "fake worker stderr")
                write(stdout, "not a serialized row")
                """
            elseif entry.mode == :slow_timeout
                """
                open($(repr(progress_path)), "w") do io
                    println(io, "current_stage=normalization")
                    println(io, "stage_started_at=", time())
                end
                sleep(10)
                """
            elseif entry.mode == :peel_progress_timeout
                """
                open($(repr(progress_path)), "w") do io
                    println(io, "current_stage=certificate_construction")
                    println(io, "stage_started_at=", time())
                    println(io, "stage.determinant_classification.status=pass")
                    println(io, "stage.determinant_classification.elapsed_seconds=0.001")
                    println(io, "stage.determinant_classification.error_details=none")
                    println(io, "stage.normalization.status=pass")
                    println(io, "stage.normalization.elapsed_seconds=0.002")
                    println(io, "stage.normalization.error_details=none")
                    println(io, "peel.current_dimension=21")
                    println(io, "peel.completed_steps=9")
                    println(io, "peel.last_completed_dimension=22")
                    println(io, "peel.last_completed_elapsed_seconds=0.875")
                    println(io, "peel.last_completed_left_factors=7")
                    println(io, "peel.last_completed_right_factors=11")
                    println(io, "peel.last_column_nnz=20")
                    println(io, "peel.max_entry_terms=26")
                end
                sleep(10)
                """
            elseif entry.mode == :certificate_timeout_without_peel
                """
                open($(repr(progress_path)), "w") do io
                    println(io, "current_stage=certificate_construction")
                    println(io, "stage_started_at=", time())
                    println(io, "stage.determinant_classification.status=pass")
                    println(io, "stage.determinant_classification.elapsed_seconds=0.012")
                    println(io, "stage.determinant_classification.error_details=none")
                    println(io, "stage.normalization.status=pass")
                    println(io, "stage.normalization.elapsed_seconds=0.034")
                    println(io, "stage.normalization.error_details=none")
                end
                sleep(10)
                """
            elseif entry.mode == :partial_peel_progress_timeout
                """
                open($(repr(progress_path)), "w") do io
                    println(io, "current_stage=certificate_construction")
                    println(io, "stage_started_at=", time())
                    println(io, "stage.determinant_classification.status=pass")
                    println(io, "stage.determinant_classification.elapsed_seconds=0.001")
                    println(io, "stage.determinant_classification.error_details=none")
                    println(io, "stage.normalization.status=pass")
                    println(io, "stage.normalization.elapsed_seconds=0.002")
                    println(io, "stage.normalization.error_details=none")
                    println(io, "peel.current_dimension=21")
                    println(io, "peel.completed_steps=9")
                    println(io, "peel.last_column_nnz=20")
                end
                sleep(10)
                """
            elseif entry.mode == :serialized_success
                result_writer = result_path === nothing ?
                    """
                    serialize(stdout, row)
                    """ :
                    """
                    open($(repr(result_path)), "w") do io
                        serialize(io, row)
                    end
                    """
                """
                using Serialization
                row = (
                    case_id = $(repr(entry.id)),
                    matrix_size = (5, 5),
                    sparse_entry_count = 25,
                    expected_test_level = :default_contract,
                    route_status = :gl_certificate_pass,
                    public_elementary_status = :not_run,
                    determinant_class = :laurent_monomial_unit,
                    determinant = "1",
                    normalization_status = :pass,
                    gl_certificate_status = :pass,
                    verified = true,
                    factor_count = 3,
                    decomposed_base_matrix_count = 3,
                    runtime_seconds = 0.01,
                    error_details = "none",
                    evidence = "fake bounded worker success",
                    stage_timings = (
                        determinant_classification = (status = :pass, elapsed_seconds = 0.001, error_details = "none"),
                        normalization = (status = :pass, elapsed_seconds = 0.002, error_details = "none"),
                        certificate_construction = (status = :pass, elapsed_seconds = 0.003, error_details = "none"),
                        verification = (status = :pass, elapsed_seconds = 0.004, error_details = "none"),
                    ),
                )
                $result_writer
                """
            elseif entry.mode == :noisy_serialized_success
                result_writer = result_path === nothing ?
                    """
                    serialize(stdout, row)
                    """ :
                    """
                    open($(repr(result_path)), "w") do io
                        serialize(io, row)
                    end
                    """
                """
                using Serialization
                write(stdout, "diagnostic stdout before serialized worker row")
                row = (
                    case_id = $(repr(entry.id)),
                    matrix_size = (5, 5),
                    sparse_entry_count = 25,
                    expected_test_level = :default_contract,
                    route_status = :gl_certificate_pass,
                    public_elementary_status = :not_run,
                    determinant_class = :laurent_monomial_unit,
                    determinant = "1",
                    normalization_status = :pass,
                    gl_certificate_status = :pass,
                    verified = true,
                    factor_count = 3,
                    decomposed_base_matrix_count = 3,
                    runtime_seconds = 0.01,
                    error_details = "none",
                    evidence = "fake noisy bounded worker success",
                    stage_timings = (
                        determinant_classification = (status = :pass, elapsed_seconds = 0.001, error_details = "none"),
                        normalization = (status = :pass, elapsed_seconds = 0.002, error_details = "none"),
                        certificate_construction = (status = :pass, elapsed_seconds = 0.003, error_details = "none"),
                        verification = (status = :pass, elapsed_seconds = 0.004, error_details = "none"),
                    ),
                )
                $result_writer
                """
            else
                error("unsupported fake bounded worker mode $(entry.mode)")
            end
            return `$(Base.julia_cmd()) --project=$(project_path) -e $(script)`
        end

        function _wait_for_exit_after_kill(proc::Base.Process; grace_seconds = 1.0, poll_seconds = 0.05)
            exited = invoke(
                _wait_for_exit_after_kill,
                Tuple{Any},
                proc;
                grace_seconds = grace_seconds,
                poll_seconds = poll_seconds,
            )
            return _TEST_FORCE_WAIT_FOR_EXIT_FAILURE[] ? false : exited
        end
    end
    report = ToricBuilderCacheQBlockStatusReport.build_report()

    @test length(report.rows) == 12
    @test report.source_fixture == "test/fixtures/toricbuilder_cache_q_blocks.jl"

    by_id = Dict(row.case_id => row for row in report.rows)
    @test Set(keys(by_id)) == Set(["case_$(lpad(string(idx), 3, "0"))" for idx in 1:12])

    @test by_id["case_001"].route_status == :gl_certificate_pass
    @test by_id["case_001"].public_elementary_status == :staged_boundary
    @test by_id["case_001"].determinant_class == :laurent_monomial_unit
    @test by_id["case_001"].verified == true
    @test by_id["case_001"].factor_count > 0
    @test by_id["case_001"].runtime_seconds > 0

    @test by_id["case_002"].route_status == :gl_certificate_pass
    @test by_id["case_002"].public_elementary_status == :staged_boundary
    @test by_id["case_002"].verified == true
    @test by_id["case_002"].decomposed_base_matrix_count == 186
    @test by_id["case_002"].runtime_seconds > 0

    @test by_id["case_003"].route_status == :gl_certificate_pass
    @test by_id["case_003"].public_elementary_status == :staged_boundary
    @test by_id["case_003"].verified == true
    @test by_id["case_003"].factor_count > 0
    @test by_id["case_003"].runtime_seconds > 0

    @test by_id["case_004"].route_status == :gl_certificate_pass
    @test by_id["case_004"].decomposed_base_matrix_count == 189
    @test by_id["case_005"].route_status == :gl_certificate_pass
    @test by_id["case_005"].decomposed_base_matrix_count == 168
    @test by_id["case_006"].route_status == :gl_certificate_pass
    @test by_id["case_006"].decomposed_base_matrix_count == 212

    @test by_id["case_010"].route_status == :gl_certificate_pass
    @test by_id["case_010"].public_elementary_status == :staged_boundary
    @test by_id["case_010"].determinant_class == :laurent_monomial_unit
    @test by_id["case_010"].verified == true
    @test by_id["case_010"].decomposed_base_matrix_count > 0
    @test by_id["case_010"].runtime_seconds > 0
    @test by_id["case_010"].error_details == "none"

    @test by_id["case_011"].matrix_size == (288, 288)
    @test by_id["case_011"].route_status == :not_exercised_in_default_report
    @test by_id["case_011"].public_elementary_status == :not_run
    @test by_id["case_011"].decomposed_base_matrix_count == :not_run
    @test by_id["case_011"].runtime_seconds == :not_run

    @test by_id["case_009"].matrix_size == (62, 62)
    @test by_id["case_009"].route_status == :not_exercised_in_default_report
    @test by_id["case_009"].public_elementary_status == :not_run
    @test by_id["case_009"].determinant_class == :not_run
    @test by_id["case_009"].verified == false
    @test by_id["case_009"].decomposed_base_matrix_count == :not_run
    @test by_id["case_009"].runtime_seconds == :not_run
    @test by_id["case_009"].error_details == "not_run"

    markdown = ToricBuilderCacheQBlockStatusReport.render_markdown(report)
    @test occursin("# ToricBuilder Cache Q-Block Status Report", markdown)
    @test occursin("Matrix size", markdown)
    @test occursin("Decomposed base matrices", markdown)
    @test occursin("Runtime seconds", markdown)
    @test occursin(r"\| case_001 \| 6x6 \| 30 \| default_contract \| gl_certificate_pass \| staged_boundary \| laurent_monomial_unit \| 50 \| [0-9]+\.[0-9]{3} \|", markdown)
    @test occursin(r"\| case_002 \| 14x14 \| 79 \| default_contract \| gl_certificate_pass \| staged_boundary \| laurent_monomial_unit \| 186 \| [0-9]+\.[0-9]{3} \|", markdown)
    @test occursin(r"\| case_003 \| 6x6 \| 27 \| default_contract \| gl_certificate_pass \| staged_boundary \| laurent_monomial_unit \| 49 \| [0-9]+\.[0-9]{3} \|", markdown)
    @test occursin(r"\| case_004 \| 18x18 \| 73 \| default_contract \| gl_certificate_pass \| staged_boundary \| laurent_monomial_unit \| 189 \| [0-9]+\.[0-9]{3} \|", markdown)
    @test occursin(r"\| case_005 \| 14x14 \| 90 \| default_contract \| gl_certificate_pass \| staged_boundary \| laurent_monomial_unit \| 168 \| [0-9]+\.[0-9]{3} \|", markdown)
    @test occursin(r"\| case_006 \| 18x18 \| 99 \| default_contract \| gl_certificate_pass \| staged_boundary \| laurent_monomial_unit \| 212 \| [0-9]+\.[0-9]{3} \|", markdown)
    @test occursin(r"\| case_010 \| 6x6 \| 34 \| default_contract \| gl_certificate_pass \| staged_boundary \| laurent_monomial_unit \| [1-9][0-9]* \| [0-9]+\.[0-9]{3} \|", markdown)
    @test occursin("| case_011 | 288x288 | 14713 | optional_slow | not_exercised_in_default_report | not_run | not_run | not_run | not_run |", markdown)
    @test occursin("| case_009 | 62x62 | 739 | default_contract | not_exercised_in_default_report | not_run | not_run | not_run | not_run |", markdown)
    @test hasproperty(by_id["case_001"], :stage_timings)
    @test by_id["case_001"].stage_timings.determinant_classification.status == :pass
    @test by_id["case_001"].stage_timings.normalization.status == :pass
    @test by_id["case_001"].stage_timings.certificate_construction.status == :pass
    @test by_id["case_001"].stage_timings.verification.status == :pass
    @test by_id["case_001"].stage_timings.determinant_classification.elapsed_seconds >= 0

    @test hasproperty(by_id["case_007"], :stage_timings)
    @test by_id["case_007"].stage_timings.determinant_classification.status == :not_run
    @test by_id["case_007"].stage_timings.normalization.status == :not_run
    @test by_id["case_007"].stage_timings.certificate_construction.status == :not_run
    @test by_id["case_007"].stage_timings.verification.status == :not_run

    @test hasproperty(by_id["case_009"], :stage_timings)
    @test by_id["case_009"].stage_timings.determinant_classification.status == :not_run
    @test by_id["case_009"].stage_timings.normalization.status == :not_run
    @test by_id["case_009"].stage_timings.certificate_construction.status == :not_run
    @test by_id["case_009"].stage_timings.verification.status == :not_run

    @test occursin("## Stage Timing Details", markdown)
    @test occursin("Determinant classification", markdown)
    @test occursin("Certificate construction", markdown)
    @test !occursin("## Route Error Details", markdown)
    @test !occursin("unsupported exact unimodular column reduction", markdown)
    @test occursin("julia --project=. scripts/report_toricbuilder_cache_q_blocks.jl", markdown)

    parsed_timeout = ToricBuilderCacheQBlockStatusReport._parse_args([
        "--exercise=case_007",
        "--timeout-seconds=1.5",
        "--output=/tmp/qblock-timeout.md",
    ])
    @test parsed_timeout.exercised == ["case_007"]
    @test parsed_timeout.timeout_seconds == 1.5

    @test_throws ArgumentError ToricBuilderCacheQBlockStatusReport._parse_args([
        "--timeout-seconds=0",
    ])
    @test_throws ArgumentError ToricBuilderCacheQBlockStatusReport._parse_args([
        "--timeout-seconds=not-a-number",
    ])

    @test_throws ArgumentError ToricBuilderCacheQBlockStatusReport.build_report(;
        exercised_case_ids = ("case_999",),
    )

    unknown_output = tempname()
    @test_throws ArgumentError ToricBuilderCacheQBlockStatusReport.main([
        "--exercise=case_999",
        "--output=$(unknown_output)",
    ])

    timeout_report = ToricBuilderCacheQBlockStatusReport.build_report(;
        exercised_case_ids = ("case_007",),
        timeout_seconds = 1.0,
    )
    timeout_by_id = Dict(row.case_id => row for row in timeout_report.rows)
    @test timeout_by_id["case_007"].route_status == :timed_out
    @test timeout_by_id["case_007"].runtime_seconds >= 1.0
    @test timeout_by_id["case_007"].runtime_seconds < 20.0
    timeout_stage_statuses = [
        getproperty(timeout_by_id["case_007"].stage_timings, stage).status for
        stage in ToricBuilderCacheQBlockStatusReport.STAGE_NAMES
    ]
    @test :timed_out in timeout_stage_statuses
    @test !all(==(:not_run), timeout_stage_statuses)
    @test timeout_by_id["case_008"].route_status == :not_exercised_in_default_report

    case008_timeout_report = ToricBuilderCacheQBlockStatusReport.build_report(;
        exercised_case_ids = ("case_008",),
        timeout_seconds = 1.0,
    )
    case008_timeout_by_id = Dict(row.case_id => row for row in case008_timeout_report.rows)
    case008_timeout_row = case008_timeout_by_id["case_008"]
    @test case008_timeout_row.route_status != :not_exercised_in_default_report
    @test case008_timeout_row.route_status != :route_error
    @test _bounded_route_row_is_structured(case008_timeout_row; timeout_seconds = 1.0)

    raw_route_error_row = merge(case008_timeout_row, (;
        route_status = :route_error,
        error_details = "timeout",
        stage_timings = ToricBuilderCacheQBlockStatusReport._not_run_stage_timings(),
    ))
    @test !_bounded_route_row_is_structured(raw_route_error_row; timeout_seconds = 1.0)

    plain_timeout_row = merge(case008_timeout_row, (;
        route_status = :timed_out,
        error_details = "timeout",
        stage_timings = ToricBuilderCacheQBlockStatusReport._not_run_stage_timings(),
    ))
    @test !_bounded_route_row_is_structured(plain_timeout_row; timeout_seconds = 1.0)

    timeout_output = tempname()
    timeout_path = ToricBuilderCacheQBlockStatusReport.main([
        "--exercise=case_007",
        "--timeout-seconds=1",
        "--output=$(timeout_output)",
    ])
    @test timeout_path == timeout_output
    timeout_markdown = read(timeout_output, String)
    @test occursin("| case_007 | 42x42 | 546 | default_contract | timed_out |", timeout_markdown)
    @test occursin("## Stage Timing Details", timeout_markdown)
    rm(timeout_output; force = true)

    @testset "bounded worker helpers" begin
        worker = run(`$(Base.julia_cmd()) --startup-file=no -e "exit()"`; wait = false)
        wait(worker)
        start_time = time()
        exited = ToricBuilderCacheQBlockStatusReport._wait_for_exit_after_kill(
            worker;
            grace_seconds = 0.5,
            poll_seconds = 0.01,
        )
        @test exited
        @test !process_running(worker)
        @test time() - start_time < 1.0

        entry = (;
            id = "case_route_error",
            dimensions = (; matrix = (5, 5)),
            sparse_entry_count = 25,
            expected_test_level = :default_contract,
        )
        runtime_seconds = 1.234
        stderr_text = "bounded worker stderr"
        route_error_row = ToricBuilderCacheQBlockStatusReport._worker_route_error_row(
            entry,
            runtime_seconds,
            stderr_text,
        )
        @test route_error_row.route_status == :route_error
        @test route_error_row.error_details == stderr_text
        @test route_error_row.stage_timings.determinant_classification.status == :route_error
        @test route_error_row.stage_timings.determinant_classification.elapsed_seconds == runtime_seconds
        @test route_error_row.stage_timings.determinant_classification.error_details == stderr_text

        timeout_entry = (;
            id = "case_timeout",
            dimensions = (; matrix = (4, 4)),
            sparse_entry_count = 16,
            expected_test_level = :default_contract,
        )
        timed_out_row = ToricBuilderCacheQBlockStatusReport._timed_out_row(
            timeout_entry,
            1.0,
            1.25,
            (;
                current_stage = :normalization,
                stage_started_at = time() - 0.05,
                timings = Dict{Symbol, Any}(),
            ),
        )
        @test timed_out_row.route_status == :timed_out
        @test timed_out_row.error_details == "timed out after 1.000 seconds while running normalization"
        @test timed_out_row.stage_timings.normalization.status == :timed_out
        @test timed_out_row.stage_timings.normalization.elapsed_seconds < 1.0

        fallback_timed_out_row = ToricBuilderCacheQBlockStatusReport._timed_out_row(
            timeout_entry,
            1.0,
            1.25,
            (; current_stage = :verification, timings = Dict{Symbol, Any}()),
        )
        @test fallback_timed_out_row.stage_timings.verification.elapsed_seconds == 1.25

        progress_path = tempname()
        worker_row = ToricBuilderCacheQBlockStatusReport._worker_exercised_row("case_010", progress_path)
        @test worker_row.case_id == "case_010"
        @test worker_row.route_status == :gl_certificate_pass
        @test worker_row.public_elementary_status == :not_run
        @test worker_row.error_details == "none"
        @test worker_row.stage_timings.determinant_classification.status == :pass
        @test worker_row.stage_timings.verification.status == :pass
        for path in (progress_path, string(progress_path, ".tmp"))
            isfile(path) && rm(path; force = true)
        end

        invalid_stdout_row = ToricBuilderCacheQBlockStatusReport._bounded_exercised_row(
            FakeBoundedEntry("case_invalid_stdout"; mode = :invalid_stdout),
            2.0,
        )
        @test invalid_stdout_row.route_status == :route_error
        @test occursin("fake worker stderr", invalid_stdout_row.error_details)
        @test occursin("deserialize", invalid_stdout_row.error_details)
        @test invalid_stdout_row.stage_timings.determinant_classification.status == :route_error

        FORCE_WAIT_FOR_EXIT_FAILURE[] = true
        timeout_row = try
            ToricBuilderCacheQBlockStatusReport._bounded_exercised_row(
                FakeBoundedEntry("case_kill_grace"; mode = :slow_timeout),
                0.2,
            )
        finally
            FORCE_WAIT_FOR_EXIT_FAILURE[] = false
        end
        @test timeout_row.route_status == :timed_out
        @test occursin("did not exit after kill grace", timeout_row.error_details)
        @test occursin("did not exit after kill grace", timeout_row.evidence)
        kill_grace_stage_statuses = [
            getproperty(timeout_row.stage_timings, stage).status for
            stage in ToricBuilderCacheQBlockStatusReport.STAGE_NAMES
        ]
        @test :timed_out in kill_grace_stage_statuses

        peel_timeout_row = ToricBuilderCacheQBlockStatusReport._bounded_exercised_row(
            FakeBoundedEntry("case_peel_progress"; mode = :peel_progress_timeout),
            0.2,
        )
        @test peel_timeout_row.route_status == :timed_out
        @test occursin("peel progress:", peel_timeout_row.error_details)
        @test occursin("current d=21", peel_timeout_row.error_details)
        @test occursin("completed steps=9", peel_timeout_row.error_details)
        @test occursin("last completed d=22", peel_timeout_row.error_details)
        @test occursin("last-column nnz=20", peel_timeout_row.error_details)
        @test occursin("max entry terms=26", peel_timeout_row.error_details)
        @test occursin("peel progress:", peel_timeout_row.evidence)
        peel_timeout_markdown = ToricBuilderCacheQBlockStatusReport.render_markdown((;
            title = "Synthetic Timeout Report",
            generated_on = Dates.today(),
            source_fixture = "synthetic",
            exercised_case_ids = ("case_peel_progress",),
            timeout_seconds = 0.2,
            rows = [peel_timeout_row],
        ))
        @test occursin("| case_peel_progress | 5x5 | 25 | default_contract | timed_out |", peel_timeout_markdown)
        @test occursin("peel progress:", peel_timeout_markdown)
        @test occursin("current d=21", peel_timeout_markdown)
        @test occursin("completed steps=9", peel_timeout_markdown)
        @test occursin("last completed d=22", peel_timeout_markdown)
        @test occursin("last-column nnz=20", peel_timeout_markdown)
        @test occursin("max entry terms=26", peel_timeout_markdown)
        @test _bounded_route_row_has_issue147_peel_evidence(peel_timeout_row; timeout_seconds = 0.2)

        missing_peel_timeout_row = ToricBuilderCacheQBlockStatusReport._bounded_exercised_row(
            FakeBoundedEntry("case_missing_peel_progress"; matrix = (30, 30), sparse_entry_count = 477, mode = :certificate_timeout_without_peel),
            0.2,
        )
        @test missing_peel_timeout_row.route_status == :timed_out
        @test occursin("certificate_construction", missing_peel_timeout_row.error_details)
        @test _bounded_route_row_is_structured(missing_peel_timeout_row; timeout_seconds = 0.2)
        @test !_bounded_route_row_has_issue147_peel_evidence(missing_peel_timeout_row; timeout_seconds = 0.2)

        partial_peel_timeout_row = ToricBuilderCacheQBlockStatusReport._bounded_exercised_row(
            FakeBoundedEntry("case_partial_peel_progress"; mode = :partial_peel_progress_timeout),
            0.2,
        )
        @test partial_peel_timeout_row.route_status == :timed_out
        @test occursin("peel progress:", partial_peel_timeout_row.error_details)
        @test occursin("current d=21", partial_peel_timeout_row.error_details)
        @test occursin("completed steps=9", partial_peel_timeout_row.error_details)
        @test occursin("last-column nnz=20", partial_peel_timeout_row.error_details)
        @test !occursin("nothing", partial_peel_timeout_row.error_details)
        @test occursin("peel progress:", partial_peel_timeout_row.evidence)
        partial_peel_timeout_markdown = ToricBuilderCacheQBlockStatusReport.render_markdown((;
            title = "Synthetic Partial Timeout Report",
            generated_on = Dates.today(),
            source_fixture = "synthetic",
            exercised_case_ids = ("case_partial_peel_progress",),
            timeout_seconds = 0.2,
            rows = [partial_peel_timeout_row],
        ))
        @test occursin(
            "| case_partial_peel_progress | 5x5 | 25 | default_contract | timed_out |",
            partial_peel_timeout_markdown,
        )
        @test occursin("peel progress:", partial_peel_timeout_markdown)
        @test occursin("current d=21", partial_peel_timeout_markdown)
        @test occursin("completed steps=9", partial_peel_timeout_markdown)
        @test occursin("last-column nnz=20", partial_peel_timeout_markdown)
        @test !occursin("nothing", partial_peel_timeout_markdown)

        bounded_worker_row = ToricBuilderCacheQBlockStatusReport._bounded_exercised_row(
            FakeBoundedEntry("case_success"; mode = :serialized_success),
            5.0,
        )
        @test bounded_worker_row.case_id == "case_success"
        @test bounded_worker_row.route_status == :gl_certificate_pass
        @test bounded_worker_row.error_details == "none"

        noisy_worker_row = ToricBuilderCacheQBlockStatusReport._bounded_exercised_row(
            FakeBoundedEntry("case_noisy_success"; mode = :noisy_serialized_success),
            5.0,
        )
        @test noisy_worker_row.case_id == "case_noisy_success"
        @test noisy_worker_row.route_status == :gl_certificate_pass
        @test noisy_worker_row.error_details == "none"
    end

    @testset "_record_stage! failure paths" begin
        entry = (;
            id = "case_failure",
            dimensions = (; matrix = (3, 3)),
            sparse_entry_count = 9,
            expected_test_level = :default_contract,
        )

        start_ns = time_ns()
        boundary_timings = Dict{Symbol, Any}()
        boundary_result = ToricBuilderCacheQBlockStatusReport._record_stage!(
            boundary_timings,
            :determinant_classification,
            () -> throw(ArgumentError("unsupported Laurent GL_n determinant test boundary")),
        )
        @test boundary_result.status == :certified_algorithm_boundary
        boundary_row = ToricBuilderCacheQBlockStatusReport._stage_failure_row(
            entry,
            boundary_timings,
            boundary_result,
            start_ns,
        )
        @test boundary_row.route_status == :certified_algorithm_boundary
        @test boundary_row.stage_timings.determinant_classification.status ==
              :certified_algorithm_boundary
        @test occursin("determinant_classification", boundary_row.evidence)

        exact_boundary_timings = Dict{Symbol, Any}()
        exact_boundary_result = ToricBuilderCacheQBlockStatusReport._record_stage!(
            exact_boundary_timings,
            :certificate_construction,
            () -> throw(ArgumentError(
                "unsupported exact unimodular column reduction for Laurent-normalized column of length 21: no supported unit, witness-unit, monicity-normalized, or 3-entry block reduction stage applies",
            )),
        )
        @test exact_boundary_result.status == :certified_algorithm_boundary

        start_ns = time_ns()
        route_timings = Dict{Symbol, Any}()
        route_result = ToricBuilderCacheQBlockStatusReport._record_stage!(
            route_timings,
            :determinant_classification,
            () -> throw(ErrorException("unexpected qblock test failure")),
        )
        @test route_result.status == :route_error
        route_row = ToricBuilderCacheQBlockStatusReport._stage_failure_row(
            entry,
            route_timings,
            route_result,
            start_ns,
        )
        @test route_row.route_status == :route_error
        @test route_row.stage_timings.determinant_classification.status == :route_error
    end
    @test !isfile(unknown_output)
    @test isfile(TORICBUILDER_CACHE_STATUS_REPORT_PATH)
end
