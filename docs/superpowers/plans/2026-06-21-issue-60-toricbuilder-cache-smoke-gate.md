# Issue 60 ToricBuilder Cache Smoke Gate Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an optional local ToricBuilder cache smoke gate for selected cache cases and blocks.

**Architecture:** Keep the smoke command under `example/toric_decoupling` and load ToricBuilder dynamically from a local checkout. The command emits stable `TORIC_SMOKE` key-value rows and a local optional test file verifies cache-error handling plus live cache summaries when the local environment is present.

**Tech Stack:** Julia, Suslin, optional local ToricBuilder checkout, Oscar matrices, Test stdlib.

## Global Constraints

- Do not add ToricBuilder to Suslin's main `Project.toml`.
- Do not register the optional smoke gate in `test/runtests.jl` or GitHub Actions.
- Default selected cases are `case_001` and `case_004`.
- Negative cache control `case_013` must report `CACHE_ERROR` or another stable unavailable-case status and must not report a false successful factorization.
- Output rows must include case id, block role, size, determinant classification, normalization status, `SL` core status, optional `GL` certificate status, factor count, exact verification result, and stable failure code.
- Stable status tokens must include `SL_CORE_PASS`, `GL_CERT_PASS`, `UNSUPPORTED_STAGED`, and `CACHE_ERROR`.
- Determinant-one blocks use the recursive column-peel route.
- Laurent monomial-unit blocks use determinant normalization and the Laurent `GL_n` certificate path.
- Do not commit `Manifest.toml`.
- Focused optional verification commands are `julia --project=example/toric_decoupling example/toric_decoupling/runtests.jl`, `julia --project=example/toric_decoupling example/toric_decoupling/try_column_block_decoupling.jl --case=case_001,case_004`, and `julia --project=example/toric_decoupling example/toric_decoupling/try_column_block_decoupling.jl --case=case_013`.
- Package verification command is `julia --project=. -e 'using Pkg; Pkg.test()'`.

---

## File Structure

- Create `example/toric_decoupling/Project.toml`: tiny optional example environment containing `Test` only; Suslin and ToricBuilder are loaded by local `LOAD_PATH` entries inside the script.
- Create `example/toric_decoupling/try_column_block_decoupling.jl`: optional smoke command, dynamic dependency loading, cache case loading, block extraction, route execution, stable summary formatting, and guarded CLI entry point.
- Create `example/toric_decoupling/runtests.jl`: optional tests for stable cache errors and live cache summaries when the local ToricBuilder cache is available.

---

### Task 1: Add the Optional ToricBuilder Cache Smoke Gate

**Files:**
- Create: `example/toric_decoupling/Project.toml`
- Create: `example/toric_decoupling/try_column_block_decoupling.jl`
- Create: `example/toric_decoupling/runtests.jl`

**Interfaces:**
- Consumes: optional local ToricBuilder checkout, `load_cached_toric_case`, `ToricBuilder.laurent_conjugate`, `Suslin.classify_laurent_determinant`, `Suslin.normalize_laurent_gl_matrix`, `Suslin._factor_laurent_sl_column_peel`, and `Suslin.laurent_gl_factorization_certificate`.
- Produces: guarded command `main(args=ARGS; io=stdout)` and reusable `run_smoke(cases; toricbuilder_dir, cache_dir, io)` returning row named tuples.

- [ ] **Step 1: Write the failing optional tests**

Create `example/toric_decoupling/runtests.jl` with tests that include
`try_column_block_decoupling.jl`, call `run_smoke(["case_013"]; toricbuilder_dir="/definitely/missing", cache_dir="/definitely/missing", io=io)`, and assert:

```julia
@test only(rows).case == "case_013"
@test only(rows).status == "FAIL"
@test only(rows).failure == "CACHE_ERROR"
@test occursin("TORIC_SMOKE case=case_013", String(take!(io)))
```

Also add a live-cache test guarded by `local_toricbuilder_available()` that runs
`case_001` and `case_004`, asserts rows exist for `column_Q` and
`pair_mix_2_1`, and checks that every `pair_mix_2_1` row has
`sl_core == "SL_CORE_PASS"` and `verified == "true"`.

- [ ] **Step 2: Run the focused test and confirm RED**

Run:

```bash
julia --project=example/toric_decoupling example/toric_decoupling/runtests.jl
```

Expected: FAIL because `try_column_block_decoupling.jl` is missing or does not
yet define `run_smoke`.

- [ ] **Step 3: Add the optional example project**

Create `example/toric_decoupling/Project.toml`:

```toml
name = "SuslinToricDecouplingExample"
uuid = "9e38cbf1-3fa0-481f-865b-fd27b9f0e7d5"
version = "0.1.0"

[deps]
Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
```

- [ ] **Step 4: Implement the smoke command**

Create `example/toric_decoupling/try_column_block_decoupling.jl` implementing:

```julia
const DEFAULT_CASES = ["case_001", "case_004"]
const DEFAULT_BLOCKS = ["column_Q", "pair_mix_2_1"]
const DEFAULT_TORICBUILDER_DIR = joinpath(homedir(), "jcode", "topological-code-decoupling", "julia_code", "ToricBuilder")

main(args=ARGS; io=stdout)
run_smoke(cases=DEFAULT_CASES; toricbuilder_dir=default_toricbuilder_dir(), cache_dir=default_cache_dir(toricbuilder_dir), io=stdout)
local_toricbuilder_available(; toricbuilder_dir=default_toricbuilder_dir(), cache_dir=default_cache_dir(toricbuilder_dir))
```

The implementation must:

- Parse `--case=case_001,case_004`, `--toricbuilder-dir=...`, and
  `--cache-dir=...`.
- Push the Suslin repository root and ToricBuilder checkout onto `LOAD_PATH`
  before dynamically importing packages.
- Return `CACHE_ERROR` rows for missing ToricBuilder, missing cache files,
  failed cache status, missing transfer fields, or block extraction errors.
- Extract `column_Q` from the upper-left decoder block and `pair_mix_2_1` from
  the deterministic toric pair mix.
- Emit deterministic `TORIC_SMOKE` key-value rows with no raw stack traces.
- Use `Suslin._factor_laurent_sl_column_peel(A)` for `det=one` blocks.
- Use `Suslin.laurent_gl_factorization_certificate(A)` for
  `det=laurent_monomial_unit` blocks.
- Map staged `ArgumentError` route failures to `UNSUPPORTED_STAGED`.

- [ ] **Step 5: Run focused optional verification**

Run:

```bash
julia --project=example/toric_decoupling example/toric_decoupling/runtests.jl
julia --project=example/toric_decoupling example/toric_decoupling/try_column_block_decoupling.jl --case=case_001,case_004
julia --project=example/toric_decoupling example/toric_decoupling/try_column_block_decoupling.jl --case=case_013
```

Expected: tests pass; `case_001` and `case_004` print deterministic rows with
stable route tokens; `case_013` reports `CACHE_ERROR` and no successful factors.

- [ ] **Step 6: Run package verification**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: package tests pass and remain independent of ToricBuilder.

- [ ] **Step 7: Commit**

Run:

```bash
git add docs/superpowers/specs/2026-06-21-issue-60-toricbuilder-cache-smoke-gate-design.md docs/superpowers/plans/2026-06-21-issue-60-toricbuilder-cache-smoke-gate.md example/toric_decoupling/Project.toml example/toric_decoupling/try_column_block_decoupling.jl example/toric_decoupling/runtests.jl
git commit -m "Add optional ToricBuilder cache smoke gate"
```

Expected: commit succeeds.

## Self-Review

- Spec coverage: the plan covers the optional command, local test command,
  stable summary tokens, determinant-one and monomial-unit routes, unavailable
  case handling, and package-test independence.
- Placeholder scan: no unresolved placeholders remain.
- Type consistency: planned function names and token strings match the design.
