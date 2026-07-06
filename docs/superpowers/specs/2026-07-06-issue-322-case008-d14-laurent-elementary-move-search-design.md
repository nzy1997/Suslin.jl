# Issue 322 Case008 D14 Laurent Elementary Move Search Design

## Goal

Add an expert-only, deterministic bounded search report for `case_008 d=14`
Laurent elementary entry additions. The report evaluates the explicit default
window from issue #322 against the #321 descent measure contract and records
only replay-verified decreasing candidates, or an exhausted boundary if none
exist.

## Context

No repository `AGENTS.md` file is present in this worktree. Issue #322 has no
comments. Relevant merged PR context:

- #318 added the checked-in `case_008 d=14` Laurent boundary fixture.
- #319 exposed the Laurent-native ECP boundary diagnostic.
- #320 added the d14 Laurent descent profile.
- #324 added the d14 Laurent descent measure contract.

The baseline focused command
`include("test/expert/case008_d14_laurent_descent_measure_contract.jl")`
passes in this sandbox when run with the direct Julia 1.12.5 binary and a
writable temporary depot.

## Approach Options

Recommended: create `test/expert/case008_d14_laurent_elementary_move_search.jl`
as a test-local report module. Reuse the d14 profile and measure contract,
scan the bounded operation window using exact Laurent support-set arithmetic
over `GF(2)`, and replay only operations that appear to strictly decrease the
measure before recording them as candidates.

Alternative: recompute the full profile from polynomial columns for every one
of the 1638 operations. That is simpler conceptually but too slow: the current
d14 profile computation takes about 95 seconds on its own in this environment.

Alternative: move the search into `src/` as a production reducer step. That is
out of scope for this issue and would imply public reducer behavior before the
future Laurent link witness and endpoint reduction issues exist.

Chosen approach: expert-only support scan plus replay verification for recorded
candidates. This preserves exactness for the declared `GF(2)` coefficient set,
keeps the report deterministic and offline, and avoids production API churn.

## Report Shape

Expose:

`case008_d14_laurent_elementary_move_search_report(fixture = boundary_fixture())`

The returned named tuple has stable fields:

- `case_id = "case_008"`;
- `dimension = 14`;
- `input_measure` equal to `case008_d14_laurent_descent_measure(...)`;
- `operation_families = (:entry_addition,)`;
- `search_bounds`, including all ordered `(target_index, source_index)` pairs
  with distinct indices, exponent radius `1`, exponent vectors
  `(-1:1) x (-1:1)`, coefficient family `(1,)`, and `GF(2)`;
- `status in (:candidate_found, :exhausted)`;
- `candidate_count`;
- `checked_operation_count`;
- `replay_verified_count`;
- `candidates`, each containing the elementary operation and before/after
  measures.

The default window checks exactly `14 * 13 * 9 == 1638` operations. The scan
does not stop early when it finds a candidate; it records all replay-verified
decreasing candidates inside the declared default bounds.

## Operation Semantics

The only operation family is:

`target_entry <- target_entry + coefficient * u^a * v^b * source_entry`

with ordered `target_index != source_index`, exponent vector `(a, b)`, and
`coefficient == 1` in `GF(2)`.

The report stores operations as named tuples with:

- `family = :entry_addition`;
- `target_index`;
- `source_index`;
- `coefficient = 1`;
- `exponent = (a, b)`;
- `ring_generators = ("u", "v")`.

Replay validates the operation shape and indices before applying it to a fresh
copy of the original column.

## Measure Evaluation

The input measure comes directly from the #321 contract. For transformed
columns, the search uses the same component order and value definitions:

- `whole_support_count`;
- `max_entry_terms`;
- `valuation_span`;
- `leading_exponent`;
- `leading_entry_index`.

For speed, non-candidate operations are evaluated through exact support-set
updates. Over `GF(2)`, adding a monomial multiple of the source entry changes
the target support by symmetric difference with the shifted source support.
This is equivalent to polynomial addition for the declared coefficient set.

Any operation that appears to strictly decrease the measure is replayed against
the original polynomial column. Its after-measure is recomputed from the
replayed column before it can be recorded. Candidate verification rejects:

- malformed operations;
- operations whose replayed column does not match the recorded operation;
- candidates whose before or after measure was not recomputed from replay;
- non-decreasing candidates.

## Controls

Add synthetic calibration controls over `GF(2)[u^+/-1, v^+/-1]`:

- known-good: a two-entry column where adding the first entry to the second
  strictly decreases `max_entry_terms`;
- known-bad: a non-decreasing operation is rejected;
- malformed operation: wrong source or target index is rejected;
- stale after-measure: a candidate with an after-measure not recomputed from
  replay is rejected.

## Runner Registration

Register the new expert report file in `test/runtests.jl` immediately after
`expert/case008_d14_laurent_descent_measure_contract.jl`.

## Verification

Required focused command:

```bash
julia --project=. -e 'include("test/expert/case008_d14_laurent_elementary_move_search.jl")'
```

Required package command:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

In this Agent Desk sandbox, the PATH `julia` launcher cannot create its juliaup
lockfile. Use the installed direct Julia 1.12.5 binary with
`JULIA_DEPOT_PATH=/private/tmp/suslin-julia-depot` for local verification; the
logical command remains the issue's required Julia invocation.
