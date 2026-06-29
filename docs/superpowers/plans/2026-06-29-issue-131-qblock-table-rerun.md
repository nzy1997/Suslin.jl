# Issue 131 Q-Block Table Rerun Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Refresh the ToricBuilder Q-block status report table with source monomials, one combined status column, determinant values, and elementary-factor complexity metrics.

**Architecture:** Keep the existing `scripts/report_toricbuilder_cache_q_blocks.jl` generator as the single source of report rows. Add static ToricBuilder `a,b` monomial metadata for the checked-in cache cases, compute factor metrics at the same point certificates are produced, and render the updated markdown table. Preserve existing route metadata and stage timing sections.

**Tech Stack:** Julia, Oscar, Suslin Laurent GL certificates, existing Test stdlib report tests.

## Global Constraints

- Do not widen the default exercised case set beyond the current issue131 defaults.
- Remove the main-table `Determinant class` column.
- Combine `Route status` and `Public elementary status` into one `Status` column.
- Add ToricBuilder definition monomials as two columns before the status/metric evidence.
- Add main-table columns for determinant, `max_elementary_factor_monomial_degree(factors)`, and `total_elementary_factor_offdiagonal_monomials(factors)`.
- Use `not_run` for metrics when a case is not exercised or fails before factors are available.
- Regenerate `docs/audits/2026-06-24-toricbuilder-cache-q-block-status.md` from the current repository state.

---

### Task 1: Lock The Markdown Shape With Tests

**Files:**
- Modify: `test/internal/toricbuilder_cache_status_report.jl`

**Interfaces:**
- Consumes: `ToricBuilderCacheQBlockStatusReport.build_report()`
- Consumes: `ToricBuilderCacheQBlockStatusReport.render_markdown(report)`
- Produces: failing assertions for source monomial fields, combined status rendering, determinant rendering, and factor metric rendering.

- [ ] **Step 1: Add assertions for `case_001` and `case_002` monomials and factor metrics.**
- [ ] **Step 2: Replace old table-regex assertions with the new table layout.**
- [ ] **Step 3: Run `julia --project=. -e 'include("test/internal/toricbuilder_cache_status_report.jl")'` and confirm the failures identify missing fields/table columns.**

### Task 2: Implement Row Metadata And Metrics

**Files:**
- Modify: `scripts/report_toricbuilder_cache_q_blocks.jl`

**Interfaces:**
- Produces: row fields `toricbuilder_a`, `toricbuilder_b`, `factor_max_monomial_degree`, and `factor_total_offdiagonal_monomials`.
- Produces: updated main markdown table with the requested columns.

- [ ] **Step 1: Add a static mapping for ToricBuilder source monomials for `case_001` through `case_012`.**
- [ ] **Step 2: Add helpers for source monomial lookup, combined status text, and factor metric extraction.**
- [ ] **Step 3: Populate the new fields in pending, failure, timeout, route-error, eager success, lazy success, and worker success rows.**
- [ ] **Step 4: Render the requested main table columns and keep detailed route/stage metadata sections intact.**

### Task 3: Verify And Rerun The Report

**Files:**
- Modify: `docs/audits/2026-06-24-toricbuilder-cache-q-block-status.md`

**Interfaces:**
- Consumes: `scripts/report_toricbuilder_cache_q_blocks.jl`
- Produces: regenerated issue131 markdown audit.

- [ ] **Step 1: Run the focused status-report test.**
- [ ] **Step 2: Run `julia --project=. scripts/report_toricbuilder_cache_q_blocks.jl`.**
- [ ] **Step 3: Re-run the focused status-report test after regeneration.**
- [ ] **Step 4: Inspect `git diff -- scripts/report_toricbuilder_cache_q_blocks.jl test/internal/toricbuilder_cache_status_report.jl docs/audits/2026-06-24-toricbuilder-cache-q-block-status.md`.**
