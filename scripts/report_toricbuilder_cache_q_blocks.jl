module ToricBuilderCacheQBlockStatusReport

using Dates
using Serialization
using Oscar
using Printf
using Suslin

include(joinpath(@__DIR__, "..", "test", "fixtures", "toricbuilder_cache_q_blocks.jl"))

const DEFAULT_EXERCISED_CASE_IDS =
    ("case_001", "case_002", "case_003", "case_004", "case_005", "case_006", "case_010")
const DEFAULT_REPORT_PATH =
    joinpath(@__DIR__, "..", "docs", "audits", "2026-06-24-toricbuilder-cache-q-block-status.md")
const SOURCE_FIXTURE = "test/fixtures/toricbuilder_cache_q_blocks.jl"

function _symbol_text(value)
    text = string(value)
    return startswith(text, ":") ? text[2:end] : text
end

function _size_text(size_tuple)
    return string(size_tuple[1], "x", size_tuple[2])
end

function _runtime_text(value)
    value isa Number && return @sprintf("%.3f", value)
    return _symbol_text(value)
end

function _elapsed_seconds(start_ns::UInt64)
    return round((time_ns() - start_ns) / 1.0e9; digits = 3)
end

const STAGE_NAMES = (
    :determinant_classification,
    :normalization,
    :certificate_construction,
    :verification,
)

function _stage_timing(status::Symbol; elapsed_seconds = :not_run, error_details = "none")
    return (; status, elapsed_seconds, error_details)
end

function _not_run_stage_timings()
    return (;
        determinant_classification = _stage_timing(:not_run; error_details = "not_run"),
        normalization = _stage_timing(:not_run; error_details = "not_run"),
        certificate_construction = _stage_timing(:not_run; error_details = "not_run"),
        verification = _stage_timing(:not_run; error_details = "not_run"),
    )
end

function _stage_timings_from_dict(timings::Dict{Symbol, Any})
    return (;
        determinant_classification = get(
            timings,
            :determinant_classification,
            _stage_timing(:not_run; error_details = "not_run"),
        ),
        normalization = get(timings, :normalization, _stage_timing(:not_run; error_details = "not_run")),
        certificate_construction = get(
            timings,
            :certificate_construction,
            _stage_timing(:not_run; error_details = "not_run"),
        ),
        verification = get(timings, :verification, _stage_timing(:not_run; error_details = "not_run")),
    )
end

function _write_worker_progress(progress_path, current_stage::Symbol, stage_started_at::Float64, timings)
    progress_path === nothing && return nothing
    tmp = string(progress_path, ".tmp")
    open(tmp, "w") do io
        println(io, "current_stage=", current_stage)
        println(io, "stage_started_at=", stage_started_at)
        for stage in STAGE_NAMES
            timing = get(timings, stage, nothing)
            timing === nothing && continue
            println(io, "stage.", stage, ".status=", timing.status)
            println(io, "stage.", stage, ".elapsed_seconds=", timing.elapsed_seconds)
            println(io, "stage.", stage, ".error_details=", replace(timing.error_details, '\n' => "\\n"))
        end
    end
    mv(tmp, progress_path; force = true)
    return nothing
end

function _read_worker_progress(progress_path)
    progress_path === nothing && return (; current_stage = :determinant_classification, timings = Dict{Symbol, Any}())
    isfile(progress_path) || return (; current_stage = :determinant_classification, timings = Dict{Symbol, Any}())
    data = Dict{String, String}()
    for line in eachline(progress_path)
        key, value = split(line, "="; limit = 2)
        data[key] = value
    end
    timings = Dict{Symbol, Any}()
    for stage in STAGE_NAMES
        prefix = string("stage.", stage, ".")
        haskey(data, string(prefix, "status")) || continue
        elapsed_raw = data[string(prefix, "elapsed_seconds")]
        elapsed = elapsed_raw == "not_run" ? :not_run : parse(Float64, elapsed_raw)
        timings[stage] = _stage_timing(
            Symbol(data[string(prefix, "status")]);
            elapsed_seconds = elapsed,
            error_details = replace(get(data, string(prefix, "error_details"), "none"), "\\n" => "\n"),
        )
    end
    return (;
        current_stage = Symbol(get(data, "current_stage", "determinant_classification")),
        stage_started_at = parse(Float64, get(data, "stage_started_at", string(time()))),
        timings,
    )
end

function _record_stage!(timings::Dict{Symbol, Any}, stage::Symbol, f)
    start_ns = time_ns()
    try
        value = f()
        timings[stage] = _stage_timing(:pass; elapsed_seconds = _elapsed_seconds(start_ns))
        return (; status = :pass, value)
    catch err
        err isa InterruptException && rethrow()
        error_details = sprint(showerror, err)
        status = _staged_argument_error(err) ? :certified_algorithm_boundary : :route_error
        timings[stage] = _stage_timing(status; elapsed_seconds = _elapsed_seconds(start_ns), error_details)
        return (; status, error = err, error_details)
    end
end

function _record_worker_stage!(f, timings::Dict{Symbol, Any}, stage::Symbol, progress_path)
    stage_started_at = time()
    _write_worker_progress(progress_path, stage, stage_started_at, timings)
    result = _record_stage!(timings, stage, f)
    _write_worker_progress(progress_path, stage, stage_started_at, timings)
    return result
end

function _stage_timing_text(timing)
    status = _symbol_text(timing.status)
    timing.elapsed_seconds isa Number && return string(status, " (", _runtime_text(timing.elapsed_seconds), "s)")
    return status
end

function _print_stage_timing_details(io, rows)
    println(io, "## Stage Timing Details")
    println(io)
    println(io, "| Case | Determinant classification | Normalization | Certificate construction | Verification |")
    println(io, "| --- | --- | --- | --- | --- |")
    for row in rows
        timings = row.stage_timings
        println(
            io,
            "| $(row.case_id) | $(_stage_timing_text(timings.determinant_classification)) | " *
            "$(_stage_timing_text(timings.normalization)) | " *
            "$(_stage_timing_text(timings.certificate_construction)) | " *
            "$(_stage_timing_text(timings.verification)) |",
        )
    end
    println(io)
end

function _parse_timeout_seconds(raw::AbstractString)
    value = try
        parse(Float64, raw)
    catch err
        err isa ArgumentError || rethrow()
        throw(ArgumentError("--timeout-seconds must be a positive number, got $(repr(raw))"))
    end
    isfinite(value) && value > 0 ||
        throw(ArgumentError("--timeout-seconds must be positive, got $(raw)"))
    return value
end

function _staged_argument_error(err)
    err isa ArgumentError || return false
    message = sprint(showerror, err)
    return occursin("staged", message) ||
        occursin("unsupported Laurent GL_n determinant", message) ||
        occursin("normalization boundary", message)
end

function _public_elementary_status(A)
    try
        factors = elementary_factorization(A)
        return verify_factorization(A, factors) ? :elementary_factorization_pass : :factorization_unverified
    catch err
        err isa InterruptException && rethrow()
        return _staged_argument_error(err) ? :staged_boundary : :route_error
    end
end

function _pending_row(entry)
    return (;
        case_id = entry.id,
        matrix_size = entry.dimensions.matrix,
        sparse_entry_count = entry.sparse_entry_count,
        expected_test_level = entry.expected_test_level,
        route_status = :not_exercised_in_default_report,
        public_elementary_status = :not_run,
        determinant_class = :not_run,
        determinant = "not_run",
        normalization_status = :not_run,
        gl_certificate_status = :not_run,
        verified = false,
        factor_count = 0,
        decomposed_base_matrix_count = :not_run,
        runtime_seconds = :not_run,
        error_details = "not_run",
        evidence = "Fixture recorded; full Suslin route audit not exercised in the default report run.",
        stage_timings = _not_run_stage_timings(),
    )
end

function _stage_failure_row(
    entry,
    timings,
    result,
    start_ns::UInt64;
    public_status = :not_run,
    profile = nothing,
)
    determinant_class = profile === nothing ? result.status : profile.classification
    determinant = profile === nothing ? string(result.status) : string(profile.determinant)
    return (;
        case_id = entry.id,
        matrix_size = entry.dimensions.matrix,
        sparse_entry_count = entry.sparse_entry_count,
        expected_test_level = entry.expected_test_level,
        route_status = result.status,
        public_elementary_status = public_status,
        determinant_class,
        determinant,
        normalization_status = get(timings, :normalization, _stage_timing(:not_run)).status,
        gl_certificate_status = get(timings, :certificate_construction, _stage_timing(:not_run)).status,
        verified = false,
        factor_count = 0,
        decomposed_base_matrix_count = 0,
        runtime_seconds = _elapsed_seconds(start_ns),
        error_details = result.error_details,
        evidence = "Bounded certificate route stopped at $(result.status); see Route Error Details.",
        stage_timings = _stage_timings_from_dict(timings),
    )
end

function _exercised_row(entry)
    start_ns = time_ns()
    A = ToricBuilderCacheQBlocks.materialize_matrix(entry)
    public_status = _public_elementary_status(A)
    stage_timings = Dict{Symbol, Any}()

    profile_result = _record_stage!(stage_timings, :determinant_classification, () -> classify_laurent_determinant(A))
    profile_result.status == :pass || return _stage_failure_row(entry, stage_timings, profile_result, start_ns)
    profile = profile_result.value

    try
        normalization_result =
            _record_stage!(stage_timings, :normalization, () -> normalize_laurent_gl_matrix(A))
        normalization_result.status == :pass ||
            return _stage_failure_row(entry, stage_timings, normalization_result, start_ns; profile = profile)
        normalization = normalization_result.value

        certificate_result = _record_stage!(
            stage_timings,
            :certificate_construction,
            () -> laurent_gl_factorization_certificate(A),
        )
        certificate_result.status == :pass ||
            return _stage_failure_row(entry, stage_timings, certificate_result, start_ns; profile = profile)
        certificate = certificate_result.value

        verification_result = _record_stage!(
            stage_timings,
            :verification,
            () -> verify_laurent_gl_factorization_certificate(certificate),
        )
        verification_result.status == :pass ||
            return _stage_failure_row(entry, stage_timings, verification_result, start_ns; profile = profile)
        verified = verification_result.value
        return (;
            case_id = entry.id,
            matrix_size = entry.dimensions.matrix,
            sparse_entry_count = entry.sparse_entry_count,
            expected_test_level = entry.expected_test_level,
            route_status = verified ? :gl_certificate_pass : :gl_certificate_fail,
            public_elementary_status = public_status,
            determinant_class = profile.classification,
            determinant = string(profile.determinant),
            normalization_status = :normalization_pass,
            gl_certificate_status = verified ? :gl_certificate_pass : :gl_certificate_fail,
            verified,
            factor_count = length(certificate.core_factors),
            decomposed_base_matrix_count = length(certificate.core_factors),
            runtime_seconds = _elapsed_seconds(start_ns),
            error_details = "none",
            evidence = "normalize_laurent_gl_matrix and laurent_gl_factorization_certificate exercised; normalized determinant is $(det(normalization.normalized_matrix)).",
            stage_timings = _stage_timings_from_dict(stage_timings),
        )
    catch err
        err isa InterruptException && rethrow()
        error_details = sprint(showerror, err)
        return (;
            case_id = entry.id,
            matrix_size = entry.dimensions.matrix,
            sparse_entry_count = entry.sparse_entry_count,
            expected_test_level = entry.expected_test_level,
            route_status = _staged_argument_error(err) ? :unsupported_staged : :route_error,
            public_elementary_status = public_status,
            determinant_class = profile.classification,
            determinant = string(profile.determinant),
            normalization_status = :normalization_or_certificate_fail,
            gl_certificate_status = :gl_certificate_fail,
            verified = false,
            factor_count = 0,
            decomposed_base_matrix_count = 0,
            runtime_seconds = _elapsed_seconds(start_ns),
            error_details,
            evidence = "Route probe failed; see Route Error Details.",
            stage_timings = _stage_timings_from_dict(stage_timings),
        )
    end
end

function _worker_exercised_row(case_id::AbstractString, progress_path)
    catalog = ToricBuilderCacheQBlocks.catalog()
    _validate_exercised_case_ids!(Set([String(case_id)]), catalog)
    entry = only(filter(entry -> entry.id == case_id, catalog.cases))

    start_ns = time_ns()
    timings = Dict{Symbol, Any}()
    A = ToricBuilderCacheQBlocks.materialize_matrix(entry)

    profile_result = _record_worker_stage!(timings, :determinant_classification, progress_path) do
        classify_laurent_determinant(A)
    end
    if profile_result.status != :pass
        return _stage_failure_row(entry, timings, profile_result, start_ns)
    end
    profile = profile_result.value

    normalization_result = _record_worker_stage!(timings, :normalization, progress_path) do
        normalize_laurent_gl_matrix(A)
    end
    normalization_result.status == :pass ||
        return _stage_failure_row(entry, timings, normalization_result, start_ns; profile)
    normalization = normalization_result.value

    certificate_result = _record_worker_stage!(timings, :certificate_construction, progress_path) do
        laurent_gl_factorization_certificate(A)
    end
    certificate_result.status == :pass ||
        return _stage_failure_row(entry, timings, certificate_result, start_ns; profile)
    certificate = certificate_result.value

    verification_result = _record_worker_stage!(timings, :verification, progress_path) do
        verify_laurent_gl_factorization_certificate(certificate)
    end
    verification_result.status == :pass ||
        return _stage_failure_row(entry, timings, verification_result, start_ns; profile)
    verified = verification_result.value

    return (;
        case_id = entry.id,
        matrix_size = entry.dimensions.matrix,
        sparse_entry_count = entry.sparse_entry_count,
        expected_test_level = entry.expected_test_level,
        route_status = verified ? :gl_certificate_pass : :gl_certificate_fail,
        public_elementary_status = :not_run,
        determinant_class = profile.classification,
        determinant = string(profile.determinant),
        normalization_status = :normalization_pass,
        gl_certificate_status = verified ? :gl_certificate_pass : :gl_certificate_fail,
        verified,
        factor_count = length(certificate.core_factors),
        decomposed_base_matrix_count = length(certificate.core_factors),
        runtime_seconds = _elapsed_seconds(start_ns),
        error_details = "none",
        evidence = "Bounded certificate route exercised; normalized determinant is $(det(normalization.normalized_matrix)).",
        stage_timings = _stage_timings_from_dict(timings),
    )
end

function _worker_main(case_id::AbstractString, progress_path)
    row = _worker_exercised_row(case_id, progress_path)
    serialize(stdout, row)
    return nothing
end

function _catalog_case_ids(catalog)
    return Set(entry.id for entry in catalog.cases)
end

function _validate_exercised_case_ids!(exercised::Set{String}, catalog)
    known = _catalog_case_ids(catalog)
    unknown = sort!(collect(setdiff(exercised, known)))
    if !isempty(unknown)
        known_text = join(sort!(collect(known)), ", ")
        throw(ArgumentError(
            "unknown --exercise case ID(s): $(join(unknown, ", ")); known case IDs: $(known_text)",
        ))
    end
    return exercised
end

function _timed_out_row(entry, timeout_seconds, runtime_seconds, progress)
    timings = copy(progress.timings)
    current_stage = progress.current_stage
    stage_elapsed =
        hasproperty(progress, :stage_started_at) ? round(time() - progress.stage_started_at; digits = 3) :
        runtime_seconds
    timings[current_stage] = _stage_timing(
        :timed_out;
        elapsed_seconds = max(stage_elapsed, timeout_seconds),
        error_details = "timed out after $(_runtime_text(timeout_seconds)) seconds",
    )
    return (;
        case_id = entry.id,
        matrix_size = entry.dimensions.matrix,
        sparse_entry_count = entry.sparse_entry_count,
        expected_test_level = entry.expected_test_level,
        route_status = :timed_out,
        public_elementary_status = :not_run,
        determinant_class = :timed_out,
        determinant = "timed_out",
        normalization_status = current_stage == :normalization ? :timed_out : get(
            timings,
            :normalization,
            _stage_timing(:not_run),
        ).status,
        gl_certificate_status = current_stage in (:certificate_construction, :verification) ? :timed_out : :not_run,
        verified = false,
        factor_count = 0,
        decomposed_base_matrix_count = 0,
        runtime_seconds,
        error_details = "timed out after $(_runtime_text(timeout_seconds)) seconds while running $(current_stage)",
        evidence = "Bounded worker exceeded the configured process-level budget.",
        stage_timings = _stage_timings_from_dict(timings),
    )
end

function _worker_command(entry, progress_path)
    project_path = dirname(Base.active_project())
    return `$(Base.julia_cmd()) --project=$(project_path) $(abspath(@__FILE__)) --bounded-worker=$(entry.id) --worker-progress=$(progress_path)`
end

function _worker_route_error_row(entry, runtime_seconds, stderr_text)
    error_details = isempty(stderr_text) ? "bounded worker exited nonzero" : stderr_text
    timings = Dict{Symbol, Any}(
        :determinant_classification =>
            _stage_timing(:route_error; elapsed_seconds = runtime_seconds, error_details),
    )
    return (;
        case_id = entry.id,
        matrix_size = entry.dimensions.matrix,
        sparse_entry_count = entry.sparse_entry_count,
        expected_test_level = entry.expected_test_level,
        route_status = :route_error,
        public_elementary_status = :not_run,
        determinant_class = :route_error,
        determinant = "route_error",
        normalization_status = :not_run,
        gl_certificate_status = :not_run,
        verified = false,
        factor_count = 0,
        decomposed_base_matrix_count = 0,
        runtime_seconds,
        error_details,
        evidence = "Bounded worker exited before producing a report row.",
        stage_timings = _stage_timings_from_dict(timings),
    )
end

function _wait_for_exit_after_kill(proc; grace_seconds = 1.0, poll_seconds = 0.05)
    deadline = time() + grace_seconds
    while process_running(proc)
        remaining = deadline - time()
        remaining <= 0 && break
        sleep(min(poll_seconds, remaining))
    end
    return !process_running(proc)
end

function _bounded_exercised_row(entry, timeout_seconds::Float64)
    stdout_path = tempname()
    stderr_path = tempname()
    progress_path = tempname()
    start_time = time()
    proc = run(
        pipeline(_worker_command(entry, progress_path); stdout = stdout_path, stderr = stderr_path),
        wait = false,
    )

    try
        while process_running(proc)
            runtime_seconds = round(time() - start_time; digits = 3)
            if runtime_seconds >= timeout_seconds
                try
                    kill(proc)
                catch
                end
                _wait_for_exit_after_kill(proc)
                progress = _read_worker_progress(progress_path)
                return _timed_out_row(entry, timeout_seconds, runtime_seconds, progress)
            end
            sleep(0.05)
        end

        wait(proc)
        runtime_seconds = round(time() - start_time; digits = 3)
        if success(proc)
            return open(deserialize, stdout_path)
        end
        stderr_text = isfile(stderr_path) ? read(stderr_path, String) : ""
        return _worker_route_error_row(entry, runtime_seconds, stderr_text)
    finally
        for path in (stdout_path, stderr_path, progress_path, string(progress_path, ".tmp"))
            isfile(path) && rm(path; force = true)
        end
    end
end

function build_report(;
    exercised_case_ids = DEFAULT_EXERCISED_CASE_IDS,
    generated_on = Dates.today(),
    timeout_seconds = nothing,
)
    exercised = Set(string.(exercised_case_ids))
    catalog = ToricBuilderCacheQBlocks.catalog()
    _validate_exercised_case_ids!(exercised, catalog)
    rows = [
        entry.id in exercised ?
        (timeout_seconds === nothing ? _exercised_row(entry) : _bounded_exercised_row(entry, timeout_seconds)) :
        _pending_row(entry)
        for entry in catalog.cases
    ]
    return (;
        title = "ToricBuilder Cache Q-Block Status Report",
        generated_on,
        source_fixture = SOURCE_FIXTURE,
        exercised_case_ids = Tuple(sort!(collect(exercised))),
        timeout_seconds,
        rows,
    )
end

function render_markdown(report)
    io = IOBuffer()
    println(io, "# $(report.title)")
    println(io)
    println(io, "Date: $(report.generated_on)")
    println(io, "Source fixture: `$(report.source_fixture)`")
    println(io)
    println(io, "This report records the current Suslin status for the checked-in ToricBuilder cache Q-block fixtures.")
    println(io, "The default run fully exercises `$(join(report.exercised_case_ids, "`, `"))`; every other case is saved fixture data marked `not_exercised_in_default_report`.")
    println(io)
    println(io, "## Summary")
    println(io)
    exercised_count = count(row -> row.public_elementary_status != :not_run, report.rows)
    pass_count = count(row -> row.route_status == :gl_certificate_pass, report.rows)
    error_count = count(row -> row.route_status == :route_error, report.rows)
    not_exercised_count = count(row -> row.route_status == :not_exercised_in_default_report, report.rows)
    println(io, "- Total checked-in Q-block cases: $(length(report.rows))")
    println(io, "- Fully exercised cases: $(exercised_count)")
    println(io, "- GL certificate passes: $(pass_count)")
    println(io, "- Route errors: $(error_count)")
    println(io, "- Not exercised in default report: $(not_exercised_count)")
    println(io)
    println(io, "## Case Table")
    println(io)
    println(io, "| Case | Matrix size | Sparse nnz | Test level | Route status | Public elementary status | Determinant class | Decomposed base matrices | Runtime seconds |")
    println(io, "| --- | ---: | ---: | --- | --- | --- | --- | ---: | ---: |")
    for row in report.rows
        println(
            io,
            "| $(row.case_id) | $(_size_text(row.matrix_size)) | $(row.sparse_entry_count) | " *
            "$(_symbol_text(row.expected_test_level)) | $(_symbol_text(row.route_status)) | " *
            "$(_symbol_text(row.public_elementary_status)) | $(_symbol_text(row.determinant_class)) | " *
            "$(_symbol_text(row.decomposed_base_matrix_count)) | $(_runtime_text(row.runtime_seconds)) |",
        )
    end
    println(io)
    println(io, "## Exercised Evidence")
    println(io)
    for row in report.rows
        row.public_elementary_status != :not_run || continue
        println(io, "- `$(row.case_id)`: determinant `$(row.determinant)`, route `$(row.route_status)`, public `$(row.public_elementary_status)`, decomposed base matrices `$(row.decomposed_base_matrix_count)`, verified `$(row.verified)`. $(row.evidence)")
    end
    println(io)
    _print_stage_timing_details(io, report.rows)
    error_rows = filter(row -> row.error_details != "none" && row.error_details != "not_run", report.rows)
    if !isempty(error_rows)
        println(io, "## Route Error Details")
        println(io)
        for row in error_rows
            println(io, "- `$(row.case_id)`: $(row.error_details)")
        end
        println(io)
    end
    println(io, "## Not Exercised Boundary")
    println(io)
    println(io, "The non-exercised rows are real stored matrices, not placeholders. They are intentionally not routed through the full factorization stack in the default report because large Laurent GL certification can be much slower than fixture/schema validation.")
    println(io, "Use `--exercise=case_001,case_003,...` to expand the set of cases that the report probes.")
    println(io)
    println(io, "## Reproduction")
    println(io)
    println(io, "```text")
    println(io, "julia --project=. scripts/report_toricbuilder_cache_q_blocks.jl")
    println(io, "julia --project=. test/internal/toricbuilder_cache_status_report.jl")
    println(io, "```")
    return String(take!(io))
end

function write_report(path::AbstractString = DEFAULT_REPORT_PATH; report = build_report())
    mkpath(dirname(path))
    open(path, "w") do io
        write(io, render_markdown(report))
    end
    return path
end

function _parse_args(args)
    output = DEFAULT_REPORT_PATH
    exercised = collect(DEFAULT_EXERCISED_CASE_IDS)
    timeout_seconds = nothing

    for arg in args
        if startswith(arg, "--output=")
            output = arg[length("--output=")+1:end]
        elseif startswith(arg, "--exercise=")
            raw = split(arg[length("--exercise=")+1:end], ",")
            exercised = [strip(case_id) for case_id in raw if !isempty(strip(case_id))]
        elseif startswith(arg, "--timeout-seconds=")
            timeout_seconds = _parse_timeout_seconds(arg[length("--timeout-seconds=")+1:end])
        else
            throw(ArgumentError("unsupported argument: $(arg)"))
        end
    end

    return (; output, exercised, timeout_seconds)
end

function _worker_arg_value(args, prefix::AbstractString)
    matches = [arg[length(prefix)+1:end] for arg in args if startswith(arg, prefix)]
    return isempty(matches) ? nothing : only(matches)
end

function main(args = ARGS)
    worker_case = _worker_arg_value(args, "--bounded-worker=")
    if worker_case !== nothing
        progress_path = _worker_arg_value(args, "--worker-progress=")
        _worker_main(worker_case, progress_path)
        return nothing
    end

    options = _parse_args(args)
    report = build_report(;
        exercised_case_ids = options.exercised,
        timeout_seconds = options.timeout_seconds,
    )
    path = write_report(options.output; report)
    println("wrote ToricBuilder cache Q-block status report to $(path)")
    return path
end

end

if abspath(PROGRAM_FILE) == @__FILE__
    ToricBuilderCacheQBlockStatusReport.main()
end
