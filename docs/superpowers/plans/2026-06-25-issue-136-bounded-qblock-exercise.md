# Issue 136 Bounded Q-Block Exercise Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add explicit bounded exercise mode to the ToricBuilder Q-block report script, with case-ID validation, process-level timeout enforcement, stable bounded statuses, and stage timing metadata.

**Architecture:** Preserve the existing default report path and default exercised case set. Add validation and timeout parsing to the report script, then evaluate timed exercised rows in a single-case Julia subprocess that the parent can kill after the configured budget while preserving partial stage progress for the final report row.

**Tech Stack:** Julia, Oscar, Suslin, Test stdlib, Serialization stdlib, local Julia subprocesses.

## Global Constraints

- Repository has no `AGENTS.md`, `CLAUDE.md`, `GEMINI.md`, or `CONVENTIONS.md`; follow `README.md` test commands and existing Suslin style.
- Extend `scripts/report_toricbuilder_cache_q_blocks.jl`; do not change the default exercised set.
- Keep optional slow cases out of routine unbounded tests unless the command is run explicitly.
- Validate requested `--exercise` IDs explicitly and throw `ArgumentError` for unknown IDs such as `case_999`.
- `--timeout-seconds` must be a positive numeric process-level budget.
- A bounded slow row must use stable route statuses such as `:gl_certificate_pass`, `:certified_algorithm_boundary`, or `:timed_out`.
- Stage timing metadata must cover determinant classification, normalization, certificate construction, and verification.
- Do not make any slow case pass in this issue.
- Required issue command: `julia --project=. scripts/report_toricbuilder_cache_q_blocks.jl --exercise=case_007 --timeout-seconds=1 --output=/tmp/qblock-case007-timeout.md`.
- Required internal report command: `julia --project=. -e 'include("test/internal/toricbuilder_cache_status_report.jl")'`.
- Required package verification command: `julia --project=. -e 'using Pkg; Pkg.test()'`.

---

## File Structure

- Modify `scripts/report_toricbuilder_cache_q_blocks.jl`: add stage timing row fields, exercise validation, timeout argument parsing, worker subprocess mode, parent timeout handling, and markdown timing details.
- Modify `test/internal/toricbuilder_cache_status_report.jl`: add regression tests for validation, timeout parsing, default non-exercise behavior, stage timing metadata, markdown timing output, and the bounded `case_007` command path.
- Modify `docs/audits/2026-06-24-toricbuilder-cache-q-block-status.md`: refresh the generated default report after the script output format changes.
- Keep `docs/superpowers/specs/2026-06-25-issue-136-bounded-qblock-exercise-design.md` and this plan in the PR.

---

### Task 1: Validate Exercise Inputs And Parse Timeout

**Files:**
- Modify: `test/internal/toricbuilder_cache_status_report.jl`
- Modify: `scripts/report_toricbuilder_cache_q_blocks.jl`

**Interfaces:**
- Consumes: `ToricBuilderCacheQBlockStatusReport._parse_args(args)`, `ToricBuilderCacheQBlockStatusReport.build_report(; exercised_case_ids)`, `ToricBuilderCacheQBlockStatusReport.main(args)`.
- Produces: `_parse_timeout_seconds(raw::AbstractString)`, `timeout_seconds` in parsed options, and catalog-backed exercise validation inside `build_report`.

- [ ] **Step 1: Write failing validation and parsing tests**

Append these assertions near the end of the existing `@testset` in
`test/internal/toricbuilder_cache_status_report.jl`, before the final `end`:

```julia
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
    @test !isfile(unknown_output)
```

- [ ] **Step 2: Run the status-report test to verify RED**

Run:

```bash
julia --project=. -e 'include("test/internal/toricbuilder_cache_status_report.jl")'
```

Expected: FAIL because `_parse_args` does not accept `--timeout-seconds` yet and
`build_report` still accepts unknown exercise IDs.

- [ ] **Step 3: Implement timeout parsing and exercise validation**

In `scripts/report_toricbuilder_cache_q_blocks.jl`, add this helper after
`_elapsed_seconds`:

```julia
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
```

Add catalog validation helpers before `build_report`:

```julia
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
```

Change `build_report` to accept `timeout_seconds = nothing`, load the catalog
before rows are built, and validate the exercised set:

```julia
function build_report(;
    exercised_case_ids = DEFAULT_EXERCISED_CASE_IDS,
    generated_on = Dates.today(),
    timeout_seconds = nothing,
)
    exercised = Set(string.(exercised_case_ids))
    catalog = ToricBuilderCacheQBlocks.catalog()
    _validate_exercised_case_ids!(exercised, catalog)
    rows = [
        entry.id in exercised ? _exercised_row(entry) : _pending_row(entry)
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
```

Change `_parse_args` to initialize and parse `timeout_seconds`:

```julia
    timeout_seconds = nothing
```

and inside the argument loop:

```julia
        elseif startswith(arg, "--timeout-seconds=")
            timeout_seconds = _parse_timeout_seconds(arg[length("--timeout-seconds=")+1:end])
```

Return `(; output, exercised, timeout_seconds)`, and pass the option from
`main`:

```julia
report = build_report(; exercised_case_ids = options.exercised, timeout_seconds = options.timeout_seconds)
```

- [ ] **Step 4: Run the status-report test to verify GREEN**

Run:

```bash
julia --project=. -e 'include("test/internal/toricbuilder_cache_status_report.jl")'
```

Expected: PASS for the new argument parsing and unknown-ID validation tests.

- [ ] **Step 5: Commit Task 1**

Run:

```bash
git add scripts/report_toricbuilder_cache_q_blocks.jl test/internal/toricbuilder_cache_status_report.jl
git commit -m "test: validate qblock exercise arguments"
```

---

### Task 2: Add Stage Timing Metadata To Report Rows

**Files:**
- Modify: `test/internal/toricbuilder_cache_status_report.jl`
- Modify: `scripts/report_toricbuilder_cache_q_blocks.jl`

**Interfaces:**
- Consumes: row named tuples returned by `_pending_row(entry)` and `_exercised_row(entry)`.
- Produces: row field `stage_timings`, `_stage_timing(status; elapsed_seconds, error_details)`, `_not_run_stage_timings()`, `_stage_timing_text(timing)`, and `Stage Timing Details` markdown section.

- [ ] **Step 1: Write failing stage timing and markdown tests**

Append these assertions near the existing markdown assertions in
`test/internal/toricbuilder_cache_status_report.jl`:

```julia
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

    @test occursin("## Stage Timing Details", markdown)
    @test occursin("Determinant classification", markdown)
    @test occursin("Certificate construction", markdown)
```

- [ ] **Step 2: Run the status-report test to verify RED**

Run:

```bash
julia --project=. -e 'include("test/internal/toricbuilder_cache_status_report.jl")'
```

Expected: FAIL because rows do not yet expose `stage_timings` and markdown has
no stage timing section.

- [ ] **Step 3: Implement stage timing helpers**

In `scripts/report_toricbuilder_cache_q_blocks.jl`, add `Serialization` to the
imports now because the next task's worker will use it:

```julia
using Serialization
```

Add these helpers after `_elapsed_seconds`:

```julia
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
```

- [ ] **Step 4: Attach stage timings to pending and exercised rows**

Add `stage_timings = _not_run_stage_timings(),` to `_pending_row`.

In `_exercised_row`, replace direct calls to `classify_laurent_determinant`,
`normalize_laurent_gl_matrix`, `laurent_gl_factorization_certificate`, and
`verify_laurent_gl_factorization_certificate` with `_record_stage!` calls using
a local `stage_timings = Dict{Symbol, Any}()`. The successful path should set:

```julia
profile_result = _record_stage!(stage_timings, :determinant_classification) do
    classify_laurent_determinant(A)
end
profile_result.status == :pass || return _stage_failure_row(entry, stage_timings, profile_result, start_ns)
profile = profile_result.value
```

Use the same pattern for normalization, certificate construction, and
verification. Add `stage_timings = _stage_timings_from_dict(stage_timings),` to
both success and failure rows.

Add `_stage_failure_row(entry, timings, result, start_ns; public_status = :not_run, profile = nothing)`
before `_exercised_row`:

```julia
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
```

- [ ] **Step 5: Add markdown timing details**

Add helper functions before `render_markdown`:

```julia
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
```

Call `_print_stage_timing_details(io, report.rows)` after the exercised
evidence section and before route error details.

- [ ] **Step 6: Run the status-report test to verify GREEN**

Run:

```bash
julia --project=. -e 'include("test/internal/toricbuilder_cache_status_report.jl")'
```

Expected: PASS, with default slow rows still not exercised.

- [ ] **Step 7: Commit Task 2**

Run:

```bash
git add scripts/report_toricbuilder_cache_q_blocks.jl test/internal/toricbuilder_cache_status_report.jl
git commit -m "feat: record qblock report stage timings"
```

---

### Task 3: Enforce Bounded Worker Timeout

**Files:**
- Modify: `test/internal/toricbuilder_cache_status_report.jl`
- Modify: `scripts/report_toricbuilder_cache_q_blocks.jl`

**Interfaces:**
- Consumes: `timeout_seconds` parsed by Task 1 and stage timing helpers from Task 2.
- Produces: worker mode arguments `--bounded-worker=<case_id>` and `--worker-progress=<path>`, `_bounded_exercised_row(entry, timeout_seconds)`, `_worker_exercised_row(case_id, progress_path)`, and timeout rows with `route_status == :timed_out`.

- [ ] **Step 1: Write failing bounded timeout tests**

Append these assertions near the end of the existing testset:

```julia
    timeout_report = ToricBuilderCacheQBlockStatusReport.build_report(;
        exercised_case_ids = ("case_007",),
        timeout_seconds = 1.0,
    )
    timeout_by_id = Dict(row.case_id => row for row in timeout_report.rows)
    @test timeout_by_id["case_007"].route_status == :timed_out
    @test timeout_by_id["case_007"].runtime_seconds >= 1.0
    @test timeout_by_id["case_007"].runtime_seconds < 20.0
    @test timeout_by_id["case_007"].stage_timings.determinant_classification.status in
        (:timed_out, :pass, :not_run)
    @test timeout_by_id["case_008"].route_status == :not_exercised_in_default_report

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
```

- [ ] **Step 2: Run the status-report test to verify RED**

Run:

```bash
julia --project=. -e 'include("test/internal/toricbuilder_cache_status_report.jl")'
```

Expected: FAIL because `timeout_seconds` is parsed but not yet used to run a
bounded subprocess.

- [ ] **Step 3: Add worker progress file helpers**

In `scripts/report_toricbuilder_cache_q_blocks.jl`, add:

```julia
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
    progress_path === nothing || isfile(progress_path) || return (; current_stage = :determinant_classification, timings = Dict{Symbol, Any}())
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
```

- [ ] **Step 4: Add worker-side measured execution**

Add `_record_worker_stage!`:

```julia
function _record_worker_stage!(timings::Dict{Symbol, Any}, stage::Symbol, progress_path, f)
    stage_started_at = time()
    _write_worker_progress(progress_path, stage, stage_started_at, timings)
    result = _record_stage!(timings, stage, f)
    _write_worker_progress(progress_path, stage, stage_started_at, timings)
    return result
end
```

Add `_worker_exercised_row(case_id, progress_path)`:

```julia
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
```

Add `_worker_main(case_id, progress_path)`:

```julia
function _worker_main(case_id::AbstractString, progress_path)
    row = _worker_exercised_row(case_id, progress_path)
    serialize(stdout, row)
    return nothing
end
```

- [ ] **Step 5: Add parent-side bounded subprocess execution**

Add `_timed_out_row(entry, timeout_seconds, runtime_seconds, progress)`:

```julia
function _timed_out_row(entry, timeout_seconds, runtime_seconds, progress)
    timings = copy(progress.timings)
    current_stage = progress.current_stage
    stage_elapsed = hasproperty(progress, :stage_started_at) ? round(time() - progress.stage_started_at; digits = 3) : runtime_seconds
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
        normalization_status = current_stage == :normalization ? :timed_out : get(timings, :normalization, _stage_timing(:not_run)).status,
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
```

Add parent process helpers:

```julia
function _worker_command(entry, progress_path)
    project_path = dirname(Base.active_project())
    return `$(Base.julia_cmd()) --project=$(project_path) $(abspath(@__FILE__)) --bounded-worker=$(entry.id) --worker-progress=$(progress_path)`
end

function _worker_route_error_row(entry, runtime_seconds, stderr_text)
    timings = Dict{Symbol, Any}(
        :determinant_classification =>
            _stage_timing(:route_error; elapsed_seconds = runtime_seconds, error_details = stderr_text),
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
        error_details = isempty(stderr_text) ? "bounded worker exited nonzero" : stderr_text,
        evidence = "Bounded worker exited before producing a report row.",
        stage_timings = _stage_timings_from_dict(timings),
    )
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
                kill(proc)
                try
                    wait(proc)
                catch
                end
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
```

Change `build_report` row selection to:

```julia
rows = [
    entry.id in exercised ?
    (timeout_seconds === nothing ? _exercised_row(entry) : _bounded_exercised_row(entry, timeout_seconds)) :
    _pending_row(entry)
    for entry in catalog.cases
]
```

- [ ] **Step 6: Route worker arguments before normal CLI parsing**

Add helper parsers:

```julia
function _worker_arg_value(args, prefix::AbstractString)
    matches = [arg[length(prefix)+1:end] for arg in args if startswith(arg, prefix)]
    return isempty(matches) ? nothing : only(matches)
end
```

At the top of `main(args = ARGS)`, add:

```julia
    worker_case = _worker_arg_value(args, "--bounded-worker=")
    if worker_case !== nothing
        progress_path = _worker_arg_value(args, "--worker-progress=")
        _worker_main(worker_case, progress_path)
        return nothing
    end
```

- [ ] **Step 7: Run the status-report test to verify GREEN**

Run:

```bash
julia --project=. -e 'include("test/internal/toricbuilder_cache_status_report.jl")'
```

Expected: PASS, and the bounded `case_007` assertions complete in less than 20
seconds.

- [ ] **Step 8: Commit Task 3**

Run:

```bash
git add scripts/report_toricbuilder_cache_q_blocks.jl test/internal/toricbuilder_cache_status_report.jl
git commit -m "feat: bound slow qblock exercise rows"
```

---

### Task 4: Refresh Report And Run Verification

**Files:**
- Modify: `docs/audits/2026-06-24-toricbuilder-cache-q-block-status.md`

**Interfaces:**
- Consumes: final report script behavior from Tasks 1-3.
- Produces: refreshed default audit report with the `Stage Timing Details` section.

- [ ] **Step 1: Refresh the default report markdown**

Run:

```bash
julia --project=. scripts/report_toricbuilder_cache_q_blocks.jl
```

Expected: command exits zero and refreshes
`docs/audits/2026-06-24-toricbuilder-cache-q-block-status.md` with default
exercise rows only.

- [ ] **Step 2: Run the required bounded slow command**

Run:

```bash
julia --project=. scripts/report_toricbuilder_cache_q_blocks.jl --exercise=case_007 --timeout-seconds=1 --output=/tmp/qblock-case007-timeout.md
```

Expected: command exits zero, writes `/tmp/qblock-case007-timeout.md`, and the
row for `case_007` contains `timed_out` plus `Stage Timing Details`.

- [ ] **Step 3: Run the required unknown-ID negative control**

Run:

```bash
rm -f /tmp/qblock-case999.md
julia --project=. scripts/report_toricbuilder_cache_q_blocks.jl --exercise=case_999 --timeout-seconds=1 --output=/tmp/qblock-case999.md
```

Expected: command exits nonzero with `ArgumentError: unknown --exercise case
ID(s): case_999`, and `/tmp/qblock-case999.md` does not exist.

- [ ] **Step 4: Run the required internal report test**

Run:

```bash
julia --project=. -e 'include("test/internal/toricbuilder_cache_status_report.jl")'
```

Expected: PASS.

- [ ] **Step 5: Run the required package test**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: PASS for the default package test entry point.

- [ ] **Step 6: Commit Task 4**

Run:

```bash
git add docs/audits/2026-06-24-toricbuilder-cache-q-block-status.md
git commit -m "docs: refresh qblock status timing report"
```
