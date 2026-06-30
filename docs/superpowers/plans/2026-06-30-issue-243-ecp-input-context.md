# Issue 243 ECP Input Context Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an internal replay-verified `ECPInputContext` for ordinary polynomial ECP columns.

**Architecture:** Mirror the existing ECP certificate pattern in `src/algorithm/column_reduction.jl`: construct a provisional context, replay all stored fields into verification metadata, store the verification, and expose `verify_ecp_input_context` for tamper checks. Keep Laurent columns on their existing path by rejecting them in this ordinary-polynomial context constructor.

**Tech Stack:** Julia, Oscar exact polynomial rings, existing Suslin ECP diagnostics and fixture catalogs, Julia `Test`.

## Global Constraints

- Do not produce elementary factors, extract link witnesses, call #184 `SL_3` routes, or change public reducer behavior.
- Keep the new API internal: do not add exports to `src/Suslin.jl`.
- Build on `_validated_unimodular_column`, `_column_reduction_ring_profile`, `_unimodular_witness`, and `diagnose_unimodular_column_reduction`.
- Keep Laurent columns on the existing separate path by rejecting Laurent rings in the context constructor.
- The verifier must recompute ring metadata, one-based indexing, length, unimodularity, selected-variable membership, and staged diagnostic metadata.
- Focused verification command is `julia --project=. -e 'include("test/expert/ecp_input_context.jl")'`.
- Package verification command is `julia --project=. -e 'using Pkg; Pkg.test()'`.

---

## File Structure

- Modify `src/algorithm/column_reduction.jl`: add `ECPInputContext`, constructor helpers, replay summary, and `verify_ecp_input_context`.
- Create `test/expert/ecp_input_context.jl`: focused positive and negative context tests.
- Modify `test/runtests.jl`: register the expert test file.
- Add workflow docs in `docs/superpowers/specs/2026-06-30-issue-243-ecp-input-context-design.md` and this plan.

### Task 1: Add Red ECP Input Context Tests

**Files:**
- Create: `test/expert/ecp_input_context.jl`
- Modify: `test/runtests.jl`

**Interfaces:**
- Consumes future `Suslin.ecp_input_context`, `Suslin.ECPInputContext`, and `Suslin.verify_ecp_input_context`.
- Produces a focused test contract for constructor validation and replay verification.

- [ ] **Step 1: Write the failing tests**

Create `test/expert/ecp_input_context.jl` with helper functions that load `test/fixtures/ecp_column_cases.jl` and `test/fixtures/ecp_mainline_cases.jl`, build fixture columns, and rebuild a tampered `ECPInputContext` by replacing one struct field.

The positive tests must construct contexts for:

```julia
mainline_cases["ecp-mainline-gf2-hard-slice"]
column_cases["ecp-unsupported-unimodular-gf2"]
mainline_cases["ecp-mainline-length4-coupled-qq"]
```

For each context, assert the coerced column, ring profile, variables, normalized variable order, length, witness identity, selected-variable membership, support classification, staged failure reason, staged diagnostic equality, and `verify_ecp_input_context(ctx)`.

The negative tests must assert:

```julia
@test_throws ArgumentError Suslin.ecp_input_context(non_unimodular, R)
@test_throws ArgumentError Suslin.ecp_input_context(length_two_unimodular, R)
@test_throws ArgumentError Suslin.ecp_input_context(column, R; selected_variable = x + y)
@test_throws ArgumentError Suslin.ecp_input_context(column, R; unimodularity_witness = bad_witness)
@test !Suslin.verify_ecp_input_context(tampered_witness_context)
```

- [ ] **Step 2: Register the expert test**

Add `"expert/ecp_input_context.jl"` to the expert file list in `test/runtests.jl` near the other ECP expert tests.

- [ ] **Step 3: Run RED verification**

Run:

```bash
julia --project=. -e 'include("test/expert/ecp_input_context.jl")'
```

Expected: FAIL because `Suslin.ecp_input_context` and `Suslin.verify_ecp_input_context` do not exist yet.

### Task 2: Implement ECP Input Context And Replay Verifier

**Files:**
- Modify: `src/algorithm/column_reduction.jl`
- Test: `test/expert/ecp_input_context.jl`

**Interfaces:**
- Produces `ECPInputContext`, `ecp_input_context(v, R; ...)`, outer `ECPInputContext(v, R; ...)`, `_ecp_input_context_replay_summary(context)`, and `verify_ecp_input_context(context)::Bool`.

- [ ] **Step 1: Add the internal context struct**

Add `struct ECPInputContext` beside the existing ECP certificate structs with fields for column, ring, ring profile, variables, variable order, column length, unimodularity witness, selected-variable index/value, support classification, staged failure reason, staged diagnostic, and verification metadata.

- [ ] **Step 2: Add constructor helpers**

Implement helpers to reject Laurent rings, validate supplied witness hints, compute the canonical witness with `_unimodular_witness`, normalize variable order, resolve selected variables, and derive staged metadata from `diagnose_unimodular_column_reduction`.

- [ ] **Step 3: Add replay verifier**

Implement `_ecp_input_context_replay_summary(context)` and `verify_ecp_input_context(context)::Bool` so the verifier recomputes every stored field and returns `false` for malformed or tampered objects.

- [ ] **Step 4: Run GREEN focused verification**

Run:

```bash
julia --project=. -e 'include("test/expert/ecp_input_context.jl")'
```

Expected: PASS.

### Task 3: Run Package Verification And Commit

**Files:**
- Modify: `src/algorithm/column_reduction.jl`
- Modify: `test/runtests.jl`
- Create: `test/expert/ecp_input_context.jl`
- Create: `docs/superpowers/specs/2026-06-30-issue-243-ecp-input-context-design.md`
- Create: `docs/superpowers/plans/2026-06-30-issue-243-ecp-input-context.md`

**Interfaces:**
- Consumes the finished implementation from Task 2.
- Produces a verified branch ready for PR.

- [ ] **Step 1: Run focused verification**

Run:

```bash
julia --project=. -e 'include("test/expert/ecp_input_context.jl")'
```

Expected: PASS.

- [ ] **Step 2: Run package verification**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: PASS.

- [ ] **Step 3: Review git status and commit**

Run:

```bash
git status --short
git diff --stat
git add docs/superpowers/specs/2026-06-30-issue-243-ecp-input-context-design.md docs/superpowers/plans/2026-06-30-issue-243-ecp-input-context.md src/algorithm/column_reduction.jl test/expert/ecp_input_context.jl test/runtests.jl
git commit -m "Add checked ECP input context"
```

Expected: commit succeeds with only the intended files.

## Plan Self-Review

- Every issue requirement maps to a task and test assertion.
- The plan preserves TDD by adding the failing focused test before implementation.
- The implementation is internal and does not alter public reducer behavior.
- Laurent behavior is explicitly rejected in the new ordinary-polynomial context.
- No placeholder markers remain.
