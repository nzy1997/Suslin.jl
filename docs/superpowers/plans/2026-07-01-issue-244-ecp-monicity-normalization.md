# Issue 244 ECP Monicity Normalization Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a replay-verified context-level ECP monicity normalization record that makes a selected entry first and monic.

**Architecture:** Extend `src/algorithm/column_reduction.jl` beside the existing #243 context and monicity-search code. The new record reuses the deterministic search and substitution helpers, adds concrete elementary coordinate moves, and verifies exact replay for both the standalone record and the existing `:monicity_normalization` certificate stage.

**Tech Stack:** Julia, Oscar polynomial rings and matrices, existing Suslin ECP helpers, Julia `Test`.

## Global Constraints

- Do not extract Park-Woodburn link witnesses, realize `SL_3` link obligations, or dispatch through `reduce_unimodular_column`.
- Keep the new API internal: do not add exports to `src/Suslin.jl`.
- Input is a #243 `ECPInputContext`; verify it before constructing a normalization record.
- Reuse `_deterministic_ecp_monicity_search`, `_is_monic_in_variable`, and existing substitution helpers where practical.
- Forward and inverse substitution maps must be checked as inverse polynomial-ring automorphisms on the ring generators used by the context.
- If bounded search cannot find a verified candidate, return a staged diagnostic record instead of guessing.
- A selected monic entry outside coordinate 1 must produce non-empty elementary coordinate-move factors.
- Replay must prove coordinate-move and inverse-substituted factors act over the original ring and send the original column exactly to `e_n`.
- Focused verification command is `julia --project=. -e 'include("test/expert/ecp_monicity_normalization.jl")'`.
- Package verification command is `julia --project=. -e 'using Pkg; Pkg.test()'`.

---

## File Structure

- Modify `src/algorithm/column_reduction.jl`: add `ECPMonicityNormalization`, `ECPMonicityNormalizationFailure`, context-level constructor, coordinate-move helpers, replay verifier, and updated monicity stage replay.
- Create `test/expert/ecp_monicity_normalization.jl`: focused positive and negative tests for #244.
- Modify `test/expert/ecp_variable_change_replay.jl`: update existing expectations from `:not_moved` to concrete coordinate-move replay.
- Modify `test/runtests.jl`: register the new expert test file.
- Add workflow docs in `docs/superpowers/specs/2026-07-01-issue-244-ecp-monicity-normalization-design.md` and this plan.

### Task 1: Add Red Monicity Normalization Tests

**Files:**
- Create: `test/expert/ecp_monicity_normalization.jl`
- Modify: `test/runtests.jl`

**Interfaces:**
- Consumes future `Suslin.ECPMonicityNormalization`, `Suslin.ecp_monicity_normalization`, and `Suslin.verify_ecp_monicity_normalization`.
- Produces focused tests for identity normalization, bounded variable-change normalization, coordinate moves, replay verification, and tamper rejection.

- [ ] **Step 1: Write the failing test file**

Create `test/expert/ecp_monicity_normalization.jl` with helpers that load `test/fixtures/ecp_column_cases.jl` and `test/fixtures/ecp_mainline_cases.jl`, construct #243 contexts, apply factors to columns, verify substitution-map inverse identities on ring generators, and rebuild a tampered `ECPMonicityNormalization` by replacing one struct field.

The positive tests must construct:

```julia
mainline_cases["ecp-mainline-gf2-hard-slice"]
```

as an already-first-monic #242 context with selected variable `x`.

They must construct a bounded variable-change case from the inverse substitution of `column_cases["ecp-variable-change-monic-gf2"]`, with `selected_monic_index = 1`, selected variable `y`, and `max_shift_power = 2`.

They must construct:

```julia
column_cases["ecp-variable-change-permuted-gf2"]
```

as a non-first monic-entry case with selected variable `y`, and assert `selected_monic_index != 1`, `!isempty(coordinate_move_factors)`, and `normalized_column[1] == selected_monic_entry`.

- [ ] **Step 2: Add negative controls**

In the same test file, corrupt:

```julia
:inverse_substitution
:selected_variable
:selected_monic_entry
:coordinate_move_factors
```

For the non-first case, also replace `coordinate_move_factors` with an empty vector and assert verification rejects the record.

- [ ] **Step 3: Register the expert test**

Add `"expert/ecp_monicity_normalization.jl"` to `TEST_GROUP_FILES["expert"]` near the other ECP expert tests.

- [ ] **Step 4: Run RED verification**

Run:

```bash
julia --project=. -e 'include("test/expert/ecp_monicity_normalization.jl")'
```

Expected: FAIL because `Suslin.ecp_monicity_normalization` and `Suslin.verify_ecp_monicity_normalization` do not exist yet.

### Task 2: Implement Record, Coordinate Moves, And Replay

**Files:**
- Modify: `src/algorithm/column_reduction.jl`
- Modify: `test/expert/ecp_variable_change_replay.jl`
- Test: `test/expert/ecp_monicity_normalization.jl`

**Interfaces:**
- Produces `ECPMonicityNormalization`, `ECPMonicityNormalizationFailure`, `ecp_monicity_normalization(context::ECPInputContext; ...)`, `_ecp_monicity_normalization_replay_summary(record)`, and `verify_ecp_monicity_normalization(record)::Bool`.
- Updates `:monicity_normalization` stage construction so `first_coordinate_move_factors` and `first_coordinate_column` are replayed concrete data.

- [ ] **Step 1: Add structs and public-internal entry points**

Add the two unexported structs near the existing ECP structs. Implement `ecp_monicity_normalization(context::ECPInputContext; ...)` so it verifies the context, resolves the selected variable, normalizes the target-last variable order, and returns either a verified success record or an `ECPMonicityNormalizationFailure`.

- [ ] **Step 2: Add coordinate-move helpers**

Implement helpers for first-coordinate elementary moves, substitution-map inverse checks, selected-monic-index hint validation, and target-last variable ordering.

- [ ] **Step 3: Update monicity candidate construction**

Change `_ecp_monicity_candidate_stage` so it applies coordinate moves before reducing the transformed column, stores inverse-substituted coordinate moves, and makes `stage.factors` the full inverse-substituted pipeline.

- [ ] **Step 4: Update replay**

Update `_ecp_replay_stage` for `:monicity_normalization` to recompute coordinate moves, normalized columns, inverse-substituted coordinate moves, inverse-substituted reduction factors, the full original-ring factor sequence, and the expanded verification tuple.

- [ ] **Step 5: Add standalone record replay**

Implement `_ecp_monicity_normalization_replay_summary(record)` and `verify_ecp_monicity_normalization(record)::Bool` by replaying the stored stage, checking record fields match the stage, checking the context, and requiring exact reduction to `e_n`.

- [ ] **Step 6: Update existing variable-change replay expectations**

In `test/expert/ecp_variable_change_replay.jl`, replace the old `:not_moved`/empty-factor assertion with the new rule: empty only when `selected_monic_index == 1`, otherwise `:moved_to_first` with non-empty elementary factors and `first_coordinate_column[1] == selected_monic_entry`.

- [ ] **Step 7: Run GREEN focused verification**

Run:

```bash
julia --project=. -e 'include("test/expert/ecp_monicity_normalization.jl")'
julia --project=. -e 'include("test/expert/ecp_variable_change_replay.jl")'
```

Expected: both PASS.

### Task 3: Package Verification, Review, And Commit

**Files:**
- Modify: `src/algorithm/column_reduction.jl`
- Modify: `test/expert/ecp_variable_change_replay.jl`
- Create: `test/expert/ecp_monicity_normalization.jl`
- Modify: `test/runtests.jl`
- Create: `docs/superpowers/specs/2026-07-01-issue-244-ecp-monicity-normalization-design.md`
- Create: `docs/superpowers/plans/2026-07-01-issue-244-ecp-monicity-normalization.md`

**Interfaces:**
- Consumes the finished implementation from Task 2.
- Produces a verified branch ready for PR.

- [ ] **Step 1: Run focused verification**

Run:

```bash
julia --project=. -e 'include("test/expert/ecp_monicity_normalization.jl")'
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
git add docs/superpowers/specs/2026-07-01-issue-244-ecp-monicity-normalization-design.md docs/superpowers/plans/2026-07-01-issue-244-ecp-monicity-normalization.md src/algorithm/column_reduction.jl test/expert/ecp_monicity_normalization.jl test/expert/ecp_variable_change_replay.jl test/runtests.jl
git commit -m "Add ECP monicity normalization replay"
```

Expected: commit succeeds with only the intended files.

## Plan Self-Review

- Every #244 hard requirement maps to a focused assertion or replay check.
- The plan preserves TDD by adding the failing focused test before implementation.
- The new API is internal and context-level.
- Search failure is diagnostic and non-guessing.
- No placeholder markers remain.
