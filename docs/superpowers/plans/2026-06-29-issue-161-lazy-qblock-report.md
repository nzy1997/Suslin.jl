# Issue 161 Lazy Q-Block Report Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add explicit lazy determinant certificate reporting to the bounded ToricBuilder cache Q-block status report without changing the conservative default case scope.

**Architecture:** Keep eager report behavior as the default. Thread explicit determinant strategy options through CLI parsing, row construction, bounded worker execution, worker progress, and markdown rendering. Tests prove that lazy rows always carry determinant strategy, correction side, determinant source, and stage timing, and that pass claims require exact verification plus deferred determinant evidence.

**Tech Stack:** Julia, Oscar, Suslin Laurent GL certificate APIs, Serialization worker rows, Test stdlib.

## Global Constraints

- Keep `determinant_strategy = :eager` as the report default.
- Accepted report strategy values are exactly `:eager` and `:lazy`.
- Accepted lazy correction sides are exactly `:row` and `:column`.
- `--correction-side` without `--determinant-strategy=lazy` must throw `ArgumentError`.
- Do not add cases to `DEFAULT_EXERCISED_CASE_IDS`.
- Do not update `docs/audits/2026-06-24-toricbuilder-cache-q-block-status.md`.
- Do not commit `Manifest.toml`.
- Every report row must expose `determinant_strategy`, `correction_side`, and `determinant_source`.
- Non-exercised rows must report `determinant_strategy = :not_run`, `correction_side = :not_run`, and `determinant_source = :not_run`.
- Eager exercised rows must report `determinant_strategy = :eager`, `correction_side = :not_run`, and `determinant_source = :not_run`.
- Lazy exercised rows must report `determinant_strategy = :lazy` and `correction_side in (:row, :column)`.
- Lazy rows must report `determinant_source = :deferred_submatrix` only after deferred determinant work is reached; otherwise report `:not_reached`.
- A lazy row with `route_status = :gl_certificate_pass` must have `verified = true` and `determinant_source = :deferred_submatrix`.
- A timeout must not be reported as a pass.
- Preserve stage timing details in generated markdown.
- Required focused verification command is `julia --project=. -e 'include("test/internal/toricbuilder_cache_status_report.jl")'`.
- Required lazy case_009 command is `julia --project=. scripts/report_toricbuilder_cache_q_blocks.jl --exercise=case_009 --timeout-seconds=180 --determinant-strategy=lazy --output=/tmp/qblock-case009-lazy.md`.
- Required package verification command is `julia --project=. -e 'using Pkg; Pkg.test()'`.

---

## File Structure

- Modify `test/internal/toricbuilder_cache_status_report.jl`: add synthetic lazy worker modes, row-truthfulness helpers, lazy option parse tests, lazy markdown tests, and one real lazy `case_010` worker fixture.
- Modify `scripts/report_toricbuilder_cache_q_blocks.jl`: add report option parsing, row metadata fields, lazy worker route, determinant-source progress persistence, and markdown metadata output.
- Create or update `docs/superpowers/plans/2026-06-29-issue-161-lazy-qblock-report-task-*.md` only if subagent reports need durable handoff notes.

---

### Task 1: Add Failing Lazy Report Tests

**Files:**
- Modify: `test/internal/toricbuilder_cache_status_report.jl`

**Interfaces:**
- Consumes existing `ToricBuilderCacheQBlockStatusReport._parse_args(args)`.
- Consumes planned `build_report(; determinant_strategy, correction_side)`.
- Consumes planned `_worker_exercised_row(case_id, progress_path; determinant_strategy, correction_side)`.
- Consumes planned `_bounded_exercised_row(entry, timeout_seconds; determinant_strategy, correction_side)`.
- Produces RED tests that fail before implementation because report rows and parser results do not expose lazy determinant metadata.

- [ ] **Step 1: Extend fake worker modes**

Update `FakeBoundedEntry` mode cases inside the test module override of
`_worker_command` so it accepts the planned keyword arguments:

```julia
function _worker_command(
    entry::Main.FakeBoundedEntry,
    progress_path,
    result_path = nothing;
    determinant_strategy = :eager,
    correction_side = :not_run,
)
```

Add metadata lines to existing synthetic serialized rows:

```julia
determinant_strategy = determinant_strategy,
correction_side = determinant_strategy == :lazy ? correction_side : :not_run,
determinant_source = determinant_strategy == :lazy ? :deferred_submatrix : :not_run,
```

Add a lazy progress timeout mode that has not reached deferred determinant work:

```julia
elseif entry.mode == :lazy_timeout_not_reached
    """
    open($(repr(progress_path)), "w") do io
        println(io, "current_stage=certificate_construction")
        println(io, "stage_started_at=", time())
        println(io, "determinant_strategy=lazy")
        println(io, "correction_side=$(correction_side)")
        println(io, "determinant_source=not_reached")
        println(io, "peel.current_dimension=62")
        println(io, "peel.completed_steps=0")
        println(io, "peel.last_column_nnz=19")
    end
    sleep(10)
    """
```

Add a lazy progress timeout mode that has reached deferred determinant work:

```julia
elseif entry.mode == :lazy_timeout_deferred
    """
    open($(repr(progress_path)), "w") do io
        println(io, "current_stage=certificate_construction")
        println(io, "stage_started_at=", time())
        println(io, "determinant_strategy=lazy")
        println(io, "correction_side=$(correction_side)")
        println(io, "determinant_source=deferred_submatrix")
        println(io, "peel.current_dimension=61")
        println(io, "peel.completed_steps=1")
        println(io, "peel.last_completed_dimension=62")
        println(io, "peel.last_column_nnz=23")
    end
    sleep(10)
    """
```

- [ ] **Step 2: Add row truthfulness helpers**

Add these helper functions near `_bounded_route_row_is_structured`:

```julia
function _row_lazy_metadata_is_structured(row)
    hasproperty(row, :determinant_strategy) || return false
    hasproperty(row, :correction_side) || return false
    hasproperty(row, :determinant_source) || return false
    row.determinant_strategy == :lazy || return true
    row.correction_side in (:row, :column) || return false
    row.determinant_source in (:deferred_submatrix, :not_reached, :not_run) || return false
    if row.route_status == :gl_certificate_pass
        row.verified == true || return false
        row.determinant_source == :deferred_submatrix || return false
    end
    row.route_status == :timed_out && row.verified == true && return false
    return true
end
```

Then update `_bounded_route_row_is_structured` so its first line is:

```julia
_row_lazy_metadata_is_structured(row) || return false
```

- [ ] **Step 3: Assert default rows expose non-lazy metadata**

After the default `by_id` dictionary is created, add:

```julia
@test by_id["case_001"].determinant_strategy == :eager
@test by_id["case_001"].correction_side == :not_run
@test by_id["case_001"].determinant_source == :not_run
@test by_id["case_009"].determinant_strategy == :not_run
@test by_id["case_009"].correction_side == :not_run
@test by_id["case_009"].determinant_source == :not_run
```

- [ ] **Step 4: Add option parser tests**

Add parser checks near the timeout parser tests:

```julia
parsed_lazy = ToricBuilderCacheQBlockStatusReport._parse_args([
    "--exercise=case_009",
    "--timeout-seconds=1.5",
    "--determinant-strategy=lazy",
    "--correction-side=column",
    "--output=/tmp/qblock-lazy.md",
])
@test parsed_lazy.determinant_strategy == :lazy
@test parsed_lazy.correction_side == :column

@test_throws ArgumentError ToricBuilderCacheQBlockStatusReport._parse_args([
    "--determinant-strategy=unsupported",
])
@test_throws ArgumentError ToricBuilderCacheQBlockStatusReport._parse_args([
    "--determinant-strategy=eager",
    "--correction-side=row",
])
@test_throws ArgumentError ToricBuilderCacheQBlockStatusReport._parse_args([
    "--determinant-strategy=lazy",
    "--correction-side=diagonal",
])
```

- [ ] **Step 5: Add synthetic lazy worker checks**

Inside the `"bounded worker helpers"` testset, add:

```julia
lazy_success_row = ToricBuilderCacheQBlockStatusReport._bounded_exercised_row(
    FakeBoundedEntry("case_lazy_success"; mode = :serialized_success),
    5.0;
    determinant_strategy = :lazy,
    correction_side = :column,
)
@test lazy_success_row.route_status == :gl_certificate_pass
@test lazy_success_row.verified
@test lazy_success_row.determinant_strategy == :lazy
@test lazy_success_row.correction_side == :column
@test lazy_success_row.determinant_source == :deferred_submatrix
@test _bounded_route_row_is_structured(lazy_success_row; timeout_seconds = 5.0)

lazy_not_reached_timeout = ToricBuilderCacheQBlockStatusReport._bounded_exercised_row(
    FakeBoundedEntry("case_lazy_not_reached"; mode = :lazy_timeout_not_reached),
    synthetic_timeout_seconds;
    determinant_strategy = :lazy,
    correction_side = :row,
)
@test lazy_not_reached_timeout.route_status == :timed_out
@test !lazy_not_reached_timeout.verified
@test lazy_not_reached_timeout.determinant_strategy == :lazy
@test lazy_not_reached_timeout.correction_side == :row
@test lazy_not_reached_timeout.determinant_source == :not_reached
@test _bounded_route_row_is_structured(lazy_not_reached_timeout; timeout_seconds = synthetic_timeout_seconds)

lazy_deferred_timeout = ToricBuilderCacheQBlockStatusReport._bounded_exercised_row(
    FakeBoundedEntry("case_lazy_deferred"; mode = :lazy_timeout_deferred),
    synthetic_timeout_seconds;
    determinant_strategy = :lazy,
    correction_side = :column,
)
@test lazy_deferred_timeout.route_status == :timed_out
@test lazy_deferred_timeout.determinant_source == :deferred_submatrix
@test _bounded_route_row_is_structured(lazy_deferred_timeout; timeout_seconds = synthetic_timeout_seconds)
```

- [ ] **Step 6: Add negative controls**

Add:

```julia
false_pass_row = merge(lazy_success_row, (; verified = false))
@test !_bounded_route_row_is_structured(false_pass_row; timeout_seconds = 5.0)

missing_lazy_metadata_row = merge(lazy_success_row, (; determinant_source = :not_run))
@test !_bounded_route_row_is_structured(missing_lazy_metadata_row; timeout_seconds = 5.0)

not_reached_pass_row = merge(lazy_success_row, (; determinant_source = :not_reached))
@test !_bounded_route_row_is_structured(not_reached_pass_row; timeout_seconds = 5.0)

timeout_claimed_pass_row = merge(lazy_not_reached_timeout, (; route_status = :gl_certificate_pass))
@test !_bounded_route_row_is_structured(timeout_claimed_pass_row; timeout_seconds = synthetic_timeout_seconds)
```

- [ ] **Step 7: Add real lazy fixture and markdown checks**

Add:

```julia
lazy_progress_path = tempname()
lazy_worker_row = ToricBuilderCacheQBlockStatusReport._worker_exercised_row(
    "case_010",
    lazy_progress_path;
    determinant_strategy = :lazy,
    correction_side = :row,
)
@test lazy_worker_row.case_id == "case_010"
@test lazy_worker_row.route_status == :gl_certificate_pass
@test lazy_worker_row.verified
@test lazy_worker_row.determinant_strategy == :lazy
@test lazy_worker_row.correction_side == :row
@test lazy_worker_row.determinant_source == :deferred_submatrix
@test lazy_worker_row.stage_timings.certificate_construction.status == :pass
@test lazy_worker_row.stage_timings.verification.status == :pass
for path in (lazy_progress_path, string(lazy_progress_path, ".tmp"))
    isfile(path) && rm(path; force = true)
end

lazy_report = ToricBuilderCacheQBlockStatusReport.build_report(;
    exercised_case_ids = ("case_010",),
    determinant_strategy = :lazy,
    correction_side = :column,
)
lazy_markdown = ToricBuilderCacheQBlockStatusReport.render_markdown(lazy_report)
@test occursin("determinant_strategy", lazy_markdown)
@test occursin("correction_side", lazy_markdown)
@test occursin("determinant_source", lazy_markdown)
@test occursin("| case_010 | lazy | column | deferred_submatrix |", lazy_markdown)
```

- [ ] **Step 8: Run focused test to verify RED**

Run:

```bash
julia --project=. -e 'include("test/internal/toricbuilder_cache_status_report.jl")'
```

Expected: FAIL because `_parse_args` has no determinant options and rows do not have the new metadata fields.

- [ ] **Step 9: Commit the RED test changes**

Commit only the test file:

```bash
git add test/internal/toricbuilder_cache_status_report.jl
git commit -m "test: cover lazy qblock report metadata"
```

---

### Task 2: Implement Lazy Report Metadata

**Files:**
- Modify: `scripts/report_toricbuilder_cache_q_blocks.jl`

**Interfaces:**
- Consumes test expectations from Task 1.
- Produces `_parse_args(args).determinant_strategy::Symbol`.
- Produces `_parse_args(args).correction_side::Union{Symbol, Nothing}`.
- Produces row fields `determinant_strategy`, `correction_side`, and `determinant_source`.
- Produces worker helper methods accepting `determinant_strategy` and `correction_side` keywords.

- [ ] **Step 1: Add option parsing helpers**

Add below `_parse_timeout_seconds`:

```julia
function _parse_choice_symbol(raw::AbstractString, allowed, option_name::AbstractString)
    value = Symbol(raw)
    value in allowed && return value
    allowed_text = join(string.(":", allowed), " or ")
    throw(ArgumentError("$(option_name) must be $(allowed_text), got $(repr(raw))"))
end

function _normalize_report_route_options(determinant_strategy::Symbol, correction_side)
    determinant_strategy in (:eager, :lazy) ||
        throw(ArgumentError("--determinant-strategy must be :eager or :lazy"))
    if determinant_strategy == :eager
        correction_side === nothing ||
            throw(ArgumentError("--correction-side is supported only with --determinant-strategy=lazy"))
        return (; determinant_strategy, correction_side = :not_run)
    end
    resolved_side = correction_side === nothing ? :row : correction_side
    resolved_side in (:row, :column) ||
        throw(ArgumentError("--correction-side must be :row or :column"))
    return (; determinant_strategy, correction_side = resolved_side)
end
```

- [ ] **Step 2: Thread metadata into progress files**

Update `_write_worker_progress` to accept keywords:

```julia
determinant_strategy = :eager,
correction_side = :not_run,
determinant_source = :not_run,
```

Write these lines after `stage_started_at`:

```julia
println(io, "determinant_strategy=", determinant_strategy)
println(io, "correction_side=", correction_side)
println(io, "determinant_source=", determinant_source)
```

Update `_read_worker_progress` to return those symbols with defaults:

```julia
determinant_strategy = Symbol(get(data, "determinant_strategy", "eager")),
correction_side = Symbol(get(data, "correction_side", "not_run")),
determinant_source = Symbol(get(data, "determinant_source", "not_run")),
```

Update `_write_peel_worker_progress` to pass through the same metadata keywords.

- [ ] **Step 3: Add row metadata helper and update row constructors**

Add:

```julia
function _row_route_metadata(; determinant_strategy = :eager, correction_side = :not_run, determinant_source = :not_run)
    return (; determinant_strategy, correction_side, determinant_source)
end
```

Merge this helper into every row `NamedTuple` constructor:

```julia
merge((; existing_fields...), _row_route_metadata(...))
```

Use these values:

- `_pending_row`: all `:not_run`.
- Eager pass/failure rows: `determinant_strategy = :eager`, `correction_side = :not_run`, `determinant_source = :not_run`.
- Lazy pass rows: `determinant_strategy = :lazy`, resolved `correction_side`, and certificate `determinant_source`.
- Lazy failure rows before deferred determinant work: `determinant_source = :not_reached`.
- Lazy timeout rows: use metadata from `_read_worker_progress`.
- Worker route errors: use the requested metadata and `determinant_source = :not_reached` for lazy.

- [ ] **Step 4: Add lazy worker route**

Add helper:

```julia
function _lazy_gl_certificate_with_progress(A; correction_side = :row, progress_path = nothing, stage_started_at = time(), timings = Dict{Symbol, Any}())
    determinant_source = Ref(:not_reached)
    progress_callback = progress -> _write_peel_worker_progress(
        progress_path,
        stage_started_at,
        timings,
        progress;
        determinant_strategy = :lazy,
        correction_side,
        determinant_source = determinant_source[],
    )
    deferred_certificate = Suslin._laurent_determinant_deferred_peel_certificate(
        A;
        progress_callback,
    )
    metadata = Suslin._normalize_laurent_determinant_deferred_submatrix(deferred_certificate)
    determinant_source[] = metadata.determinant_source
    _write_worker_progress(
        progress_path,
        :certificate_construction,
        stage_started_at,
        timings;
        determinant_strategy = :lazy,
        correction_side,
        determinant_source = determinant_source[],
    )
    return Suslin._laurent_gl_lazy_deferred_correction_certificate(
        metadata;
        correction_side,
        progress_callback,
    )
end
```

Update `_worker_exercised_row` to accept keywords. For eager, keep the existing staged determinant classification, normalization, certificate construction, and verification path. For lazy, skip eager determinant classification and normalization, run only certificate construction and verification, and return rows with lazy metadata.

- [ ] **Step 5: Thread options through bounded workers and report build**

Update signatures:

```julia
function _worker_command(entry, progress_path, result_path; determinant_strategy = :eager, correction_side = :not_run)
function _worker_route_error_row(entry, runtime_seconds, stderr_text; determinant_strategy = :eager, correction_side = :not_run)
function _bounded_exercised_row(entry, timeout_seconds::Float64; determinant_strategy = :eager, correction_side = :not_run)
function build_report(; exercised_case_ids = DEFAULT_EXERCISED_CASE_IDS, generated_on = Dates.today(), timeout_seconds = nothing, determinant_strategy = :eager, correction_side = nothing)
```

Pass normalized options into row construction. Add worker command flags:

```julia
--determinant-strategy=$(_symbol_text(determinant_strategy))
```

and for lazy:

```julia
--correction-side=$(_symbol_text(correction_side))
```

Update worker `main` to parse those flags before calling `_worker_main`.

- [ ] **Step 6: Render markdown metadata**

Add function:

```julia
function _print_determinant_route_metadata(io, rows)
    println(io, "## Determinant Route Metadata")
    println(io)
    println(io, "| Case | determinant_strategy | correction_side | determinant_source |")
    println(io, "| --- | --- | --- | --- |")
    for row in rows
        println(io, "| $(row.case_id) | $(_symbol_text(row.determinant_strategy)) | $(_symbol_text(row.correction_side)) | $(_symbol_text(row.determinant_source)) |")
    end
    println(io)
end
```

Call it after the case table and before exercised evidence. Add metadata to the exercised evidence line:

```julia
determinant_strategy `$(_symbol_text(row.determinant_strategy))`, correction_side `$(_symbol_text(row.correction_side))`, determinant_source `$(_symbol_text(row.determinant_source))`
```

- [ ] **Step 7: Run focused test to verify GREEN**

Run:

```bash
julia --project=. -e 'include("test/internal/toricbuilder_cache_status_report.jl")'
```

Expected: PASS.

- [ ] **Step 8: Commit implementation**

Commit only the report script:

```bash
git add scripts/report_toricbuilder_cache_q_blocks.jl
git commit -m "feat: report lazy qblock determinant metadata"
```

---

### Task 3: Verify Required Lazy Case and Package Tests

**Files:**
- Modify only if verification exposes a defect in `scripts/report_toricbuilder_cache_q_blocks.jl` or `test/internal/toricbuilder_cache_status_report.jl`.

**Interfaces:**
- Consumes implemented CLI options and lazy metadata.
- Produces verified `/tmp/qblock-case009-lazy.md`.

- [ ] **Step 1: Run focused internal verification**

Run:

```bash
julia --project=. -e 'include("test/internal/toricbuilder_cache_status_report.jl")'
```

Expected: PASS.

- [ ] **Step 2: Run required lazy case_009 command**

Run:

```bash
julia --project=. scripts/report_toricbuilder_cache_q_blocks.jl --exercise=case_009 --timeout-seconds=180 --determinant-strategy=lazy --output=/tmp/qblock-case009-lazy.md
```

Expected: exit 0. `/tmp/qblock-case009-lazy.md` contains:

```text
| case_009 |
determinant_strategy
correction_side
determinant_source
## Stage Timing Details
```

The case_009 row may pass, time out, or hit a certified boundary. It must not claim `gl_certificate_pass` unless `verified` is true and `determinant_source` is `deferred_submatrix`.

- [ ] **Step 3: Inspect the generated case_009 report**

Run:

```bash
rg -n "case_009|determinant_strategy|correction_side|determinant_source|Stage Timing Details|gl_certificate_pass|timed_out|certified_algorithm_boundary" /tmp/qblock-case009-lazy.md
```

Expected: output shows a structured case_009 row, the determinant metadata table, and stage timing details.

- [ ] **Step 4: Run package verification**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: PASS.

- [ ] **Step 5: Commit verification-only fixes if needed**

If any verification uncovers a code/test defect, fix it with a failing focused test first when possible, re-run the relevant command, then commit:

```bash
git add scripts/report_toricbuilder_cache_q_blocks.jl test/internal/toricbuilder_cache_status_report.jl
git commit -m "fix: keep lazy qblock report evidence truthful"
```

If no fixes are needed, do not create a verification-only commit.

---

## Self-Review

- Spec coverage: this plan covers explicit CLI options, row metadata, worker timeout metadata, markdown output, synthetic negative controls, one real lazy fixture, case_009 bounded verification, and full package verification.
- Placeholder scan: no placeholder tasks or unresolved option names remain.
- Type consistency: option fields are `Symbol` values throughout; `correction_side` is `nothing` only at parse/build boundaries and normalized to a row symbol before worker execution.
