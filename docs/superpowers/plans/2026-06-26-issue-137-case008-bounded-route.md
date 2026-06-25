# Issue 137 Case 008 Bounded Route Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make explicit bounded `case_008` report runs produce a stable structured route classification instead of a raw worker serialization route error.

**Architecture:** Preserve the #136 subprocess timeout model and default exercised set. Route the worker's serialized row through a dedicated result file while keeping stdout/stderr as diagnostics, then add internal validator coverage that accepts only the three issue-approved bounded route statuses for `case_008`.

**Tech Stack:** Julia, Oscar, Suslin, Test stdlib, Serialization stdlib, local Julia subprocesses.

## Global Constraints

- Repository has no `AGENTS.md`, `CLAUDE.md`, `GEMINI.md`, or `CONVENTIONS.md`; follow `README.md` test commands and existing Suslin style.
- Do not add `case_008` to `DEFAULT_EXERCISED_CASE_IDS`.
- Do not require `case_008` to pass the Laurent GL certificate route.
- Do not add a new determinant or normalization algorithm.
- Preserve #136's process-level timeout behavior and stage timing keys: `:determinant_classification`, `:normalization`, `:certificate_construction`, and `:verification`.
- A bounded `case_008` row must have one of `:gl_certificate_pass`, `:certified_algorithm_boundary`, or `:timed_out`; it must not remain `:not_exercised_in_default_report`, `:route_error`, or a raw unstructured timeout.
- `:timed_out` rows must include the timed-out stage, elapsed time, and configured timeout budget text.
- `:certified_algorithm_boundary` rows must include a stable failure code, stage name, and stage timing.
- The negative control is a raw route error or plain timeout row for `case_008`; the status-report validator must reject it.
- Required issue command: `julia --project=. scripts/report_toricbuilder_cache_q_blocks.jl --exercise=case_008 --timeout-seconds=120 --output=/tmp/qblock-case008.md`.
- Required internal report command: `julia --project=. -e 'include("test/internal/toricbuilder_cache_status_report.jl")'`.
- Required package verification command: `julia --project=. -e 'using Pkg; Pkg.test()'`.

---

## File Structure

- Modify `scripts/report_toricbuilder_cache_q_blocks.jl`: add the `--worker-result` worker argument, serialize worker rows to a dedicated result file, and deserialize that file in the bounded parent path.
- Modify `test/internal/toricbuilder_cache_status_report.jl`: add fake-worker coverage for stdout noise plus result-file serialization, and add the `case_008` structured route validator and negative controls.
- Keep `docs/superpowers/specs/2026-06-26-issue-137-case008-bounded-route-design.md` and this plan in the PR.

---

### Task 1: Serialize Bounded Worker Rows Through A Result File

**Files:**
- Modify: `test/internal/toricbuilder_cache_status_report.jl`
- Modify: `scripts/report_toricbuilder_cache_q_blocks.jl`

**Interfaces:**
- Consumes: `_bounded_exercised_row(entry, timeout_seconds)`, `_worker_command(entry, progress_path)`, `_worker_main(case_id, progress_path)`, and `_deserialize_worker_row(stdout_path)`.
- Produces: `_worker_command(entry, progress_path, result_path)`, `_worker_main(case_id, progress_path, result_path)`, `--worker-result=<path>`, and result-file deserialization that ignores unrelated worker stdout bytes.

- [ ] **Step 1: Write the failing fake-worker regression test**

In `test/internal/toricbuilder_cache_status_report.jl`, change the test-only
worker override signature:

```julia
        function _worker_command(entry::Main.FakeBoundedEntry, progress_path, result_path = nothing)
```

Then add this branch after the existing `:serialized_success` branch:

```julia
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
```

Add this assertion in the existing `"bounded worker helpers"` testset after the
`:serialized_success` assertion:

```julia
        noisy_worker_row = ToricBuilderCacheQBlockStatusReport._bounded_exercised_row(
            FakeBoundedEntry("case_noisy_success"; mode = :noisy_serialized_success),
            5.0,
        )
        @test noisy_worker_row.case_id == "case_noisy_success"
        @test noisy_worker_row.route_status == :gl_certificate_pass
        @test noisy_worker_row.error_details == "none"
```

- [ ] **Step 2: Run the focused test to verify RED**

Run:

```bash
julia --project=. -e 'include("test/internal/toricbuilder_cache_status_report.jl")'
```

Expected: FAIL on the new noisy fake worker because the parent still deserializes
captured stdout and sees non-serialized diagnostic bytes before the row.

- [ ] **Step 3: Implement the worker result-file protocol**

In `scripts/report_toricbuilder_cache_q_blocks.jl`, change `_worker_command` to:

```julia
function _worker_command(entry, progress_path, result_path)
    project_path = dirname(Base.active_project())
    return `$(Base.julia_cmd()) --project=$(project_path) $(abspath(@__FILE__)) --bounded-worker=$(entry.id) --worker-progress=$(progress_path) --worker-result=$(result_path)`
end
```

Change `_worker_main` to:

```julia
function _worker_main(case_id::AbstractString, progress_path, result_path)
    row = _worker_exercised_row(case_id, progress_path)
    if result_path === nothing
        serialize(stdout, row)
    else
        open(result_path, "w") do io
            serialize(io, row)
        end
    end
    return nothing
end
```

Change `_bounded_exercised_row` so it creates a result path, passes it to the
worker, deserializes that result file, and cleans it up:

```julia
    result_path = tempname()
```

Use the result path in the process command:

```julia
        pipeline(_worker_command(entry, progress_path, result_path); stdout = stdout_path, stderr = stderr_path),
```

Use the result path on success:

```julia
                return _deserialize_worker_row(result_path)
```

Add `result_path` to the cleanup tuple:

```julia
        for path in (stdout_path, stderr_path, progress_path, string(progress_path, ".tmp"), result_path)
```

Change the worker-mode branch in `main` to pass the optional result path:

```julia
        result_path = _worker_arg_value(args, "--worker-result=")
        _worker_main(worker_case, progress_path, result_path)
```

- [ ] **Step 4: Run the focused test to verify GREEN**

Run:

```bash
julia --project=. -e 'include("test/internal/toricbuilder_cache_status_report.jl")'
```

Expected: PASS, including the noisy fake worker result-file regression.

- [ ] **Step 5: Commit Task 1**

Run:

```bash
git add scripts/report_toricbuilder_cache_q_blocks.jl test/internal/toricbuilder_cache_status_report.jl
git commit -m "fix: isolate qblock worker result serialization"
```

---

### Task 2: Validate Case 008 Structured Bounded Status

**Files:**
- Modify: `test/internal/toricbuilder_cache_status_report.jl`

**Interfaces:**
- Consumes: `ToricBuilderCacheQBlockStatusReport.build_report(; exercised_case_ids, timeout_seconds)`, `ToricBuilderCacheQBlockStatusReport.STAGE_NAMES`, and `ToricBuilderCacheQBlockStatusReport._runtime_text(timeout_seconds)`.
- Produces: `_bounded_route_row_is_structured(row; timeout_seconds)` test validator and `case_008` bounded-route acceptance/negative-control coverage.

- [ ] **Step 1: Write failing case_008 validator tests**

In `test/internal/toricbuilder_cache_status_report.jl`, add these assertions
near the existing bounded `case_007` timeout assertions:

```julia
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
```

- [ ] **Step 2: Run the focused test to verify RED**

Run:

```bash
julia --project=. -e 'include("test/internal/toricbuilder_cache_status_report.jl")'
```

Expected: FAIL with `UndefVarError: _bounded_route_row_is_structured not defined`.

- [ ] **Step 3: Implement the test validator**

Add these helper functions near the top of
`test/internal/toricbuilder_cache_status_report.jl`, after
`const FORCE_WAIT_FOR_EXIT_FAILURE = Ref(false)`:

```julia
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
```

- [ ] **Step 4: Run the focused test to verify GREEN**

Run:

```bash
julia --project=. -e 'include("test/internal/toricbuilder_cache_status_report.jl")'
```

Expected: PASS, and the new `case_008` row is a structured bounded status.

- [ ] **Step 5: Commit Task 2**

Run:

```bash
git add test/internal/toricbuilder_cache_status_report.jl
git commit -m "test: validate case008 bounded route status"
```

---

## Final Verification

After both tasks complete and reviews are clean, run:

```bash
julia --project=. scripts/report_toricbuilder_cache_q_blocks.jl --exercise=case_008 --timeout-seconds=120 --output=/tmp/qblock-case008.md
julia --project=. -e 'include("test/internal/toricbuilder_cache_status_report.jl")'
julia --project=. -e 'using Pkg; Pkg.test()'
```

Inspect `/tmp/qblock-case008.md` and confirm the `case_008` row is not
`not_exercised_in_default_report`, not `route_error`, and not an unstructured
timeout.

## Plan Self-Review

- Spec coverage: Task 1 covers the observed worker payload corruption. Task 2
  covers `case_008` structured status acceptance and negative controls. Final
  verification covers the exact issue commands.
- Placeholder scan: no placeholders, TBDs, or deferred implementation steps.
- Type consistency: task interfaces use the same worker-result and validator
  names across tests, implementation, and final verification.
