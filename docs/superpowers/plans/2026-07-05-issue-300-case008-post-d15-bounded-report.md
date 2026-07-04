# Issue 300 Case008 Post-D15 Bounded Report Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Harden the explicit bounded `case_008` status-report acceptance check so it rejects the old d15 unsupported Laurent-column boundary after issue #299.

**Architecture:** Keep `scripts/report_toricbuilder_cache_q_blocks.jl` as the only report runner. Add a focused post-d15 predicate in `test/internal/toricbuilder_cache_status_report.jl` that builds on the existing structured bounded-row validator and checks route/stage detail text for old d15 boundary evidence, smaller-dimension boundary evidence, and completed-d15 timeout evidence. Leave the default exercised case set unchanged.

**Tech Stack:** Julia, `Test`, existing ToricBuilder bounded Q-block report helpers.

## Global Constraints

- Do not add `case_008` to `DEFAULT_EXERCISED_CASE_IDS`.
- Do not add a new report runner.
- Accepted explicit `case_008` outcomes are a Laurent `GL_n` certificate pass, a structured later algorithm boundary below dimension 15, or a structured timeout whose progress metadata proves the old d15 reducer boundary was passed.
- Do not accept route details containing `unsupported exact unimodular column reduction for Laurent-normalized column of length 15`.
- Do not accept `current d=15` with `unsupported_laurent_column_family`.
- Do not accept a timeout at or above d15 unless progress metadata proves the d15 reduction itself is complete.
- Do not require `case_007`, `case_009`, `case_011`, or `case_012` to pass.
- Do not add Steinberg factor-count optimization.

---

## File Structure

- Modify: `test/internal/toricbuilder_cache_status_report.jl`
  - Add local validator helper functions beside the existing bounded-row helper predicates.
  - Add synthetic negative and positive controls in the existing `case_008` bounded-report assertions.
- No production source file changes are expected unless the real bounded report exposes a generator-side defect during verification.

### Task 1: Post-D15 Case008 Bounded Validator

**Files:**
- Modify: `test/internal/toricbuilder_cache_status_report.jl`

**Interfaces:**
- Consumes: `_bounded_route_row_is_structured(row; timeout_seconds)`, `_matching_stage_names(row, status)`, `ToricBuilderCacheQBlockStatusReport.STAGE_NAMES`, and row fields `route_status`, `error_details`, `evidence`, and `stage_timings`.
- Produces: `_case008_bounded_route_row_passed_old_d15_boundary(row; timeout_seconds)::Bool`.

- [ ] **Step 1: Write failing post-d15 validator tests**

In `test/internal/toricbuilder_cache_status_report.jl`, after the existing `plain_timeout_row` negative control near the one-second `case_008` bounded report assertions, add:

```julia
    old_d15_boundary_details =
        "unsupported exact unimodular column reduction for Laurent-normalized column of length 15: current d=15 unsupported_laurent_column_family"
    boundary_stage_timings = merge(case008_timeout_row.stage_timings, (;
        certificate_construction = (;
            status = :certified_algorithm_boundary,
            elapsed_seconds = 0.125,
            error_details = old_d15_boundary_details,
        ),
    ))
    old_d15_boundary_row = merge(case008_timeout_row, (;
        route_status = :certified_algorithm_boundary,
        error_details = old_d15_boundary_details,
        evidence = "Bounded certificate route stopped at certified_algorithm_boundary during certificate_construction; see Route Error Details.",
        stage_timings = boundary_stage_timings,
    ))
    @test _bounded_route_row_is_structured(old_d15_boundary_row; timeout_seconds = 1.0)
    @test !_case008_bounded_route_row_passed_old_d15_boundary(
        old_d15_boundary_row;
        timeout_seconds = 1.0,
    )

    current_d15_boundary_details =
        "peel progress: current d=15, completed steps=15; unsupported_laurent_column_family"
    current_d15_boundary_row = merge(old_d15_boundary_row, (;
        error_details = current_d15_boundary_details,
        stage_timings = merge(case008_timeout_row.stage_timings, (;
            certificate_construction = (;
                status = :certified_algorithm_boundary,
                elapsed_seconds = 0.125,
                error_details = current_d15_boundary_details,
            ),
        )),
    ))
    @test _bounded_route_row_is_structured(current_d15_boundary_row; timeout_seconds = 1.0)
    @test !_case008_bounded_route_row_passed_old_d15_boundary(
        current_d15_boundary_row;
        timeout_seconds = 1.0,
    )

    current_d14_boundary_details =
        "peel progress: current d=14, completed steps=16; unsupported_laurent_column_family"
    current_d14_boundary_row = merge(old_d15_boundary_row, (;
        error_details = current_d14_boundary_details,
        stage_timings = merge(case008_timeout_row.stage_timings, (;
            certificate_construction = (;
                status = :certified_algorithm_boundary,
                elapsed_seconds = 0.125,
                error_details = current_d14_boundary_details,
            ),
        )),
    ))
    @test _bounded_route_row_is_structured(current_d14_boundary_row; timeout_seconds = 1.0)
    @test _case008_bounded_route_row_passed_old_d15_boundary(
        current_d14_boundary_row;
        timeout_seconds = 1.0,
    )

    timeout_at_d15_details =
        "timed out after 1.000 seconds; peel progress: current d=15, completed steps=15, last-column nnz=20 while running certificate_construction"
    timeout_at_d15_row = merge(case008_timeout_row, (;
        route_status = :timed_out,
        error_details = timeout_at_d15_details,
        stage_timings = merge(case008_timeout_row.stage_timings, (;
            certificate_construction = (;
                status = :timed_out,
                elapsed_seconds = 1.0,
                error_details = "timed out after 1.000 seconds; peel progress: current d=15, completed steps=15, last-column nnz=20",
            ),
        )),
    ))
    @test _bounded_route_row_is_structured(timeout_at_d15_row; timeout_seconds = 1.0)
    @test !_case008_bounded_route_row_passed_old_d15_boundary(
        timeout_at_d15_row;
        timeout_seconds = 1.0,
    )

    timeout_below_d15_details =
        "timed out after 1.000 seconds; peel progress: current d=14, completed steps=16, last-column nnz=20 while running certificate_construction"
    timeout_below_d15_row = merge(timeout_at_d15_row, (;
        error_details = timeout_below_d15_details,
        stage_timings = merge(case008_timeout_row.stage_timings, (;
            certificate_construction = (;
                status = :timed_out,
                elapsed_seconds = 1.0,
                error_details = "timed out after 1.000 seconds; peel progress: current d=14, completed steps=16, last-column nnz=20",
            ),
        )),
    ))
    @test _bounded_route_row_is_structured(timeout_below_d15_row; timeout_seconds = 1.0)
    @test _case008_bounded_route_row_passed_old_d15_boundary(
        timeout_below_d15_row;
        timeout_seconds = 1.0,
    )

    timeout_completed_d15_details =
        "timed out after 1.000 seconds; peel progress: current d=15, completed steps=16, last completed d=15, last-column nnz=20 while running certificate_construction"
    timeout_completed_d15_row = merge(timeout_at_d15_row, (;
        error_details = timeout_completed_d15_details,
        stage_timings = merge(case008_timeout_row.stage_timings, (;
            certificate_construction = (;
                status = :timed_out,
                elapsed_seconds = 1.0,
                error_details = "timed out after 1.000 seconds; peel progress: current d=15, completed steps=16, last completed d=15, last-column nnz=20",
            ),
        )),
    ))
    @test _bounded_route_row_is_structured(timeout_completed_d15_row; timeout_seconds = 1.0)
    @test _case008_bounded_route_row_passed_old_d15_boundary(
        timeout_completed_d15_row;
        timeout_seconds = 1.0,
    )
```

- [ ] **Step 2: Run test to verify it fails for the expected reason**

Run:

```bash
julia --project=. -e 'include("test/internal/toricbuilder_cache_status_report.jl")'
```

Expected: FAIL with `UndefVarError: _case008_bounded_route_row_passed_old_d15_boundary not defined`. If it fails earlier for an unrelated syntax or fixture error, fix the test setup before implementing the helper.

- [ ] **Step 3: Implement the post-d15 helper**

Near `_bounded_route_row_has_issue147_peel_evidence`, add:

```julia
const CASE008_OLD_D15_LAURENT_BOUNDARY_TEXT =
    "unsupported exact unimodular column reduction for Laurent-normalized column of length 15"
const CASE008_D15_BOUNDARY_DIMENSION = 15

function _row_route_detail_text(row)
    details = String[]
    row.error_details isa AbstractString && push!(details, row.error_details)
    hasproperty(row, :evidence) && row.evidence isa AbstractString && push!(details, row.evidence)
    if hasproperty(row, :stage_timings)
        for stage in ToricBuilderCacheQBlockStatusReport.STAGE_NAMES
            timing = getproperty(row.stage_timings, stage)
            timing.error_details isa AbstractString && push!(details, timing.error_details)
        end
    end
    return join(details, "\n")
end

function _dimension_mentions(details::AbstractString, pattern::Regex)
    return Int[parse(Int, only(match.captures)) for match in eachmatch(pattern, details)]
end

function _case008_details_have_old_d15_boundary(details::AbstractString)
    occursin(CASE008_OLD_D15_LAURENT_BOUNDARY_TEXT, details) && return true
    return occursin("current d=$(CASE008_D15_BOUNDARY_DIMENSION)", details) &&
        occursin("unsupported_laurent_column_family", details)
end

function _case008_details_prove_post_d15_progress(details::AbstractString)
    current_dimensions = _dimension_mentions(details, r"current d=([0-9]+)")
    any(<(CASE008_D15_BOUNDARY_DIMENSION), current_dimensions) && return true

    column_lengths = _dimension_mentions(
        details,
        r"Laurent-normalized column of length ([0-9]+)",
    )
    any(<(CASE008_D15_BOUNDARY_DIMENSION), column_lengths) && return true

    completed_dimensions = _dimension_mentions(details, r"last completed d=([0-9]+)")
    return any(==(CASE008_D15_BOUNDARY_DIMENSION), completed_dimensions)
end

function _case008_bounded_route_row_passed_old_d15_boundary(row; timeout_seconds)
    _bounded_route_row_is_structured(row; timeout_seconds) || return false
    details = _row_route_detail_text(row)
    _case008_details_have_old_d15_boundary(details) && return false
    row.route_status == :gl_certificate_pass && return true
    row.route_status == :certified_algorithm_boundary &&
        return _case008_details_prove_post_d15_progress(details)
    if row.route_status == :timed_out
        stages = _matching_stage_names(row, :timed_out)
        length(stages) == 1 || return false
        stage = only(stages)
        stage == :verification && return true
        stage == :certificate_construction || return false
        return _case008_details_prove_post_d15_progress(details)
    end
    return false
end
```

- [ ] **Step 4: Run focused internal test to verify it passes**

Run:

```bash
julia --project=. -e 'include("test/internal/toricbuilder_cache_status_report.jl")'
```

Expected: PASS, including the new post-d15 negative controls.

- [ ] **Step 5: Commit task changes**

Run:

```bash
git add test/internal/toricbuilder_cache_status_report.jl
git commit -m "test: gate case008 bounded report past d15"
```

### Task 2: Post-D15 Report Evidence Verification

**Files:**
- Modify only if verification exposes a defect in `scripts/report_toricbuilder_cache_q_blocks.jl` or `test/internal/toricbuilder_cache_status_report.jl`.

**Interfaces:**
- Consumes: `scripts/report_toricbuilder_cache_q_blocks.jl --exercise=case_008 --timeout-seconds=180`.
- Produces: `/tmp/qblock-case008-after-d15.md` evidence and full-suite verification.

- [ ] **Step 1: Run the issue report command**

Run:

```bash
julia --project=. scripts/report_toricbuilder_cache_q_blocks.jl --exercise=case_008 --timeout-seconds=180 --output=/tmp/qblock-case008-after-d15.md
```

Expected: exit 0.

- [ ] **Step 2: Inspect the generated report for forbidden old-d15 evidence**

Run:

```bash
rg -n "unsupported exact unimodular column reduction for Laurent-normalized column of length 15|current d=15.*unsupported_laurent_column_family|unsupported_laurent_column_family.*current d=15" /tmp/qblock-case008-after-d15.md
```

Expected: exit 1 with no matches.

Run:

```bash
rg -n "case_008|Route Error Details|peel progress|current d=|last completed d=|certified_algorithm_boundary|timed_out|gl_certificate_pass" /tmp/qblock-case008-after-d15.md
```

Expected: output includes the `case_008` row and either a pass, a smaller later boundary, or timeout progress proving post-d15 movement.

- [ ] **Step 3: Run focused and full verification**

Run:

```bash
julia --project=. -e 'include("test/internal/toricbuilder_cache_status_report.jl")'
julia --project=. test/runtests.jl all
julia --project=. -e 'using Pkg; Pkg.test()'
git diff --check
```

Expected: all commands exit 0.

- [ ] **Step 4: Commit any verification-driven fixes**

If Task 2 required edits, run:

```bash
git add scripts/report_toricbuilder_cache_q_blocks.jl test/internal/toricbuilder_cache_status_report.jl
git commit -m "fix: refresh case008 post-d15 report evidence"
```

If no edits were required, do not create an empty commit.

## Self-Review

- Spec coverage: Task 1 covers the validator and required negative controls; Task 2 covers the explicit report command, report inspection, focused internal test, full suite, package test, and diff hygiene.
- Red-flag scan: this plan contains no incomplete markers or unspecified implementation steps.
- Type consistency: the produced helper name is `_case008_bounded_route_row_passed_old_d15_boundary(row; timeout_seconds)::Bool`, and tests call the same helper.
