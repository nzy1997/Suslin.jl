# Issue 245 ECP Link Witness Extraction Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Extract verified Park-Woodburn ECP link witness records from #244-normalized ordinary-polynomial contexts.

**Architecture:** Extend the existing internal ECP link-witness code in `src/algorithm/column_reduction.jl`. The extractor consumes a verified `ECPMonicityNormalization`, generates bounded tail-combination candidates, computes resultants, obtains Bezout and cover coordinates through Oscar ideals, stores the equations in `ECPLinkWitnessRecord`, and returns a staged diagnostic when no bounded cover is proved.

**Tech Stack:** Julia, Oscar polynomial rings and ideals, existing Suslin ECP normalization/link helpers, Julia `Test`.

## Global Constraints

- Do not realize link elementary factors, consume the #184 `SL_3` route, or change the public reducer.
- Keep new APIs internal and unexported, accessed in expert tests as `Suslin.<name>`.
- Input support is ordinary polynomial rings only; Laurent columns stay rejected.
- The extractor input must be a #244 `ECPMonicityNormalization` whose normalized first entry is monic in the selected variable.
- Supplied link evidence must still be accepted first and exactly verified.
- Missing automatic coverage must return `ECPLinkWitnessExtractionFailure`, not fixture ids and not guessed witness data.
- Default automatic search bounds are total tail-coefficient monomial degree `1`, at most `2` nonzero tail-coefficient terms per candidate, and at most `3` link candidates in a cover subset.
- Use Oscar `resultant`, `ideal`, `in`, and `coordinates` for exact resultants, Bezout coordinates, and cover multipliers.
- Focused verification commands are `julia --project=. -e 'include("test/expert/ecp_link_witness_general.jl")'` and `julia --project=. -e 'include("test/internal/ecp_mainline_fixtures.jl")'`.
- Package verification command is `julia --project=. -e 'using Pkg; Pkg.test()'`.

---

## File Structure

- Modify `src/algorithm/column_reduction.jl`: add the extraction failure record, normalization overload, bounded tail-combination search helpers, ideal-coordinate helpers, extracted-record builder, and metadata verifier update.
- Create `test/expert/ecp_link_witness_general.jl`: focused #245 extraction tests and negative controls.
- Modify `test/expert/ecp_link_witnesses.jl`: update the old "missing metadata throws" assertion to the new extraction/diagnostic behavior.
- Modify `test/runtests.jl`: register the new expert test file.
- Add this plan and `docs/superpowers/specs/2026-07-01-issue-245-ecp-link-witness-extraction-design.md`.

### Task 1: Add Red Link-Extraction Tests

**Files:**
- Create: `test/expert/ecp_link_witness_general.jl`
- Modify: `test/expert/ecp_link_witnesses.jl`
- Modify: `test/runtests.jl`

**Interfaces:**
- Consumes future `Suslin.ECPLinkWitnessExtractionFailure` and `Suslin.ecp_link_witness(normalization::Suslin.ECPMonicityNormalization; ...)`.
- Produces tests for #242 QQ extraction, #242 length-4 multivariate ordinary extraction, exact equation replay, staged diagnostics, and tamper rejection.

- [ ] **Step 1: Write the new failing expert test**

Create `test/expert/ecp_link_witness_general.jl` with helpers that:

```julia
using Test
using Oscar
using Suslin

include(joinpath(@__DIR__, "..", "fixtures", "ecp_mainline_cases.jl"))

function _lw_column(entry)
    return [getproperty(entry.column_entries, name) for name in entry.column_order]
end

function _lw_context(entry)
    return Suslin.ecp_input_context(
        _lw_column(entry),
        entry.ring.object;
        variable_order = entry.ring.generators,
        selected_variable = entry.selected_variable.generator,
        unimodularity_witness = entry.unimodularity.coefficients,
    )
end

function _lw_normalization(entry)
    return Suslin.ecp_monicity_normalization(
        _lw_context(entry);
        selected_variable = entry.selected_variable.generator,
        max_shift_power = 2,
    )
end

function _lw_recompute_tail(record, idx::Int)
    total = zero(record.ring)
    tail_entries = record.original_column[2:end]
    for tail_idx in eachindex(tail_entries)
        total += record.tail_reductions[idx].lifted_tail_coefficients[tail_idx] * tail_entries[tail_idx]
    end
    return total
end
```

The main testset must extract from:

```julia
cases["ecp-mainline-qq-link-bezout"]
cases["ecp-mainline-length4-coupled-qq"]
```

For each record, assert:

```julia
@test record isa Suslin.ECPLinkWitnessRecord
@test record.metadata.source == :extracted_link_witness
@test Suslin.verify_ecp_link_witness(record)
@test record.verification.tail_reduction_ok
@test record.verification.resultants_ok
@test record.verification.bezout_ok
@test record.verification.coverage_ok
@test record.verification.path_ok
```

Then recompute every stored equation:

```julia
for idx in eachindex(record.tail_reductions)
    G = _lw_recompute_tail(record, idx)
    @test record.tail_reductions[idx].G == G
    @test record.tail_reductions[idx].tilde_G == G
    @test record.resultants[idx] == resultant(record.selected_monic_entry, G, record.selected_variable_index)
    bezout = record.bezout_coefficients[idx]
    @test bezout.f * record.selected_monic_entry + bezout.h * G == record.resultants[idx]
end
coverage_total = sum(record.coverage_multipliers[idx] * record.resultants[idx] for idx in eachindex(record.resultants); init = zero(record.ring))
@test coverage_total == one(record.ring)
@test first(record.path_points) == zero(record.ring)
@test last(record.path_points) == record.selected_variable
for idx in eachindex(record.resultants)
    @test record.path_points[idx + 1] - record.path_points[idx] ==
          record.coverage_multipliers[idx] * record.resultants[idx] * record.selected_variable
end
```

- [ ] **Step 2: Add negative controls and diagnostic checks**

In the same test file, define a record-field replacement helper and verify:

```julia
@test !Suslin.verify_ecp_link_witness(_lw_replace_record_field(record, :resultants, _lw_replace_tuple_entry(record.resultants, 1, record.resultants[1] + one(record.ring))))
@test !Suslin.verify_ecp_link_witness(_lw_replace_record_field(record, :bezout_coefficients, _lw_replace_tuple_entry(record.bezout_coefficients, 1, merge(record.bezout_coefficients[1], (; f = record.bezout_coefficients[1].f + one(record.ring))))))
@test !Suslin.verify_ecp_link_witness(_lw_replace_record_field(record, :path_points, _lw_replace_tuple_entry(record.path_points, 2, record.path_points[2] + one(record.ring))))
@test !Suslin.verify_ecp_link_witness(_lw_replace_record_field(record, :tail_reductions, _lw_replace_tuple_entry(record.tail_reductions, 1, merge(record.tail_reductions[1], (; lifted_tail_coefficients = _lw_replace_tuple_entry(record.tail_reductions[1].lifted_tail_coefficients, 1, record.tail_reductions[1].lifted_tail_coefficients[1] + one(record.ring)))))))
```

For the length-4 case, call extraction with `max_tail_terms = 1` and assert it returns `Suslin.ECPLinkWitnessExtractionFailure` with `kind == :link_witness_cover_not_proved`.

- [ ] **Step 3: Update the older supplied-witness expert test**

In `test/expert/ecp_link_witnesses.jl`, replace the final `@test err isa ArgumentError` missing-metadata expectation with a bounded diagnostic expectation on a case that cannot prove coverage at `max_tail_terms = 1`, or with an extraction success assertion for `ecp-monic-first-entry-qq`.

- [ ] **Step 4: Register the new expert test**

Add `"expert/ecp_link_witness_general.jl"` in `TEST_GROUP_FILES["expert"]` next to the existing ECP link witness tests.

- [ ] **Step 5: Run RED verification**

Run:

```bash
julia --project=. -e 'include("test/expert/ecp_link_witness_general.jl")'
```

Expected: FAIL because `ECPLinkWitnessExtractionFailure` and the normalization overload do not exist yet, and missing supplied metadata still throws.

### Task 2: Implement Bounded Exact Link-Witness Extraction

**Files:**
- Modify: `src/algorithm/column_reduction.jl`
- Test: `test/expert/ecp_link_witness_general.jl`
- Test: `test/expert/ecp_link_witnesses.jl`

**Interfaces:**
- Consumes: `ECPMonicityNormalization`, `verify_ecp_monicity_normalization`, `_ecp_normalize_variable_order`, `_is_monic_in_variable`, `_ecp_link_witness_replay_summary`.
- Produces: `ECPLinkWitnessExtractionFailure`, `ecp_link_witness(normalization::ECPMonicityNormalization; ...)`, automatic extraction for the vector entry point when no supplied witness is present, and helper routines for bounded tail candidates and Oscar ideal coordinates.

- [ ] **Step 1: Add the failure record**

Add near `ECPLinkWitnessRecord`:

```julia
struct ECPLinkWitnessExtractionFailure
    kind::Symbol
    original_column
    ring
    variable_order
    selected_variable_index::Int
    selected_variable
    selected_monic_index::Int
    selected_monic_entry
    max_tail_coefficient_degree::Int
    max_tail_terms::Int
    max_cover_witnesses::Int
    attempted_tail_reductions::Int
    valid_resultants
    message::String
end
```

- [ ] **Step 2: Refactor stored-record construction**

Extract the common record verification path into:

```julia
function _ecp_link_witness_record_from_data(
    column,
    R,
    normalized_order,
    selected_variable_index::Int,
    selected_monic_index::Int,
    residue_probes,
    tail_reductions,
    resultants,
    bezout_coefficients,
    coverage_multipliers,
    path_points,
    metadata;
    failure_message = "Park-Woodburn ECP link witness data failed exact replay verification",
)
    replay_record = ECPLinkWitnessRecord(
        tuple(column...),
        R,
        tuple(normalized_order...),
        selected_variable_index,
        gens(R)[selected_variable_index],
        selected_monic_index,
        column[selected_monic_index],
        tuple(residue_probes...),
        tuple(tail_reductions...),
        tuple(resultants...),
        tuple(bezout_coefficients...),
        tuple(coverage_multipliers...),
        tuple(path_points...),
        metadata,
        nothing,
    )
    verification = _ecp_link_witness_replay_summary(replay_record)
    verification.overall_ok || throw(ArgumentError(failure_message))
    stored = ECPLinkWitnessRecord(
        replay_record.original_column,
        replay_record.ring,
        replay_record.variable_order,
        replay_record.selected_variable_index,
        replay_record.selected_variable,
        replay_record.selected_monic_index,
        replay_record.selected_monic_entry,
        replay_record.residue_probes,
        replay_record.tail_reductions,
        replay_record.resultants,
        replay_record.bezout_coefficients,
        replay_record.coverage_multipliers,
        replay_record.path_points,
        replay_record.metadata,
        verification,
    )
    verify_ecp_link_witness(stored) || throw(ArgumentError("stored Park-Woodburn ECP link witness data failed exact replay verification"))
    return stored
end
```

Update supplied-witness construction to call this helper with `metadata = (; source = :supplied_link_witness)`.

- [ ] **Step 3: Add extraction entry points**

Add a normalization overload:

```julia
function ecp_link_witness(
    normalization::ECPMonicityNormalization;
    supplied_link_witness = nothing,
    max_tail_coefficient_degree::Integer = 1,
    max_tail_terms::Integer = 2,
    max_cover_witnesses::Integer = 3,
)
    verify_ecp_monicity_normalization(normalization) ||
        throw(ArgumentError("ECP link witness extraction requires a verified monicity normalization record"))
    return ecp_link_witness(
        collect(normalization.normalized_column),
        normalization.ring;
        variable_order = normalization.variable_order,
        selected_variable = normalization.selected_variable,
        selected_monic_index = 1,
        supplied_link_witness,
        max_tail_coefficient_degree,
        max_tail_terms,
        max_cover_witnesses,
    )
end
```

Extend the vector entry point keyword list with the three bounds. When
`supplied_link_witness === nothing`, call `_ecp_extract_link_witness(...)`
instead of throwing.

- [ ] **Step 4: Add bounded candidate and ideal-coordinate helpers**

Implement these helpers in `column_reduction.jl`:

```julia
_ecp_link_small_scalars(R)
_ecp_link_exponent_tuples(var_count::Int, max_degree::Int)
_ecp_link_monomial_basis(R, max_degree::Int)
_ecp_link_tail_reduction_candidates(tail_entries, R; max_tail_coefficient_degree::Int, max_tail_terms::Int)
_ecp_link_combinations(indices, width::Int)
_ecp_link_coordinates_tuple(coordinates_value, R, expected_length::Int, label::AbstractString)
_ecp_link_bezout_for_resultant(v1, G, resultant_value, R)
_ecp_link_cover_multipliers(resultants, R)
```

`_ecp_link_tail_reduction_candidates` must generate deterministic coefficient
tuples from monomial/scalar atoms, skip zero combinations, deduplicate
coefficient tuples, and return named tuples containing `lifted_tail_coefficients`
and `G`.

`_ecp_link_bezout_for_resultant` must call `ideal(R, [v1, G])`, require
`resultant_value in ideal`, extract two coordinates, and verify
`f * v1 + h * G == resultant_value`.

`_ecp_link_cover_multipliers` must call `ideal(R, collect(resultants))`,
require `one(R) in ideal`, extract one multiplier per resultant, and verify the
coverage sum equals `one(R)`.

- [ ] **Step 5: Build extracted records or diagnostics**

Implement `_ecp_extract_link_witness` so it:

1. validates nonnegative bounds and selected first-entry monicity;
2. computes valid candidate records with nonzero resultants and Bezout data;
3. searches cover subsets of width `1:min(max_cover_witnesses, length(valid))`;
4. builds residue probes with unique ids and `kind = :bounded_tail_combination`;
5. builds path points from cumulative cover terms;
6. calls `_ecp_link_witness_record_from_data` with extracted metadata.

If no cover is found, return `ECPLinkWitnessExtractionFailure` with
`kind = :link_witness_cover_not_proved` and the valid resultants.

- [ ] **Step 6: Update verifier metadata acceptance**

Change `_ecp_link_witness_replay_summary` so `metadata_ok` accepts:

```julia
metadata_ok = hasproperty(record.metadata, :source) &&
    record.metadata.source in (:supplied_link_witness, :extracted_link_witness)
```

Keep all existing tail/resultant/Bezout/coverage/path checks exact.

- [ ] **Step 7: Run GREEN focused verification**

Run:

```bash
julia --project=. -e 'include("test/expert/ecp_link_witness_general.jl")'
julia --project=. -e 'include("test/expert/ecp_link_witnesses.jl")'
```

Expected: both PASS.

### Task 3: Package Verification, Review, And Commit

**Files:**
- Modify: `src/algorithm/column_reduction.jl`
- Create: `test/expert/ecp_link_witness_general.jl`
- Modify: `test/expert/ecp_link_witnesses.jl`
- Modify: `test/runtests.jl`
- Create: `docs/superpowers/specs/2026-07-01-issue-245-ecp-link-witness-extraction-design.md`
- Create: `docs/superpowers/plans/2026-07-01-issue-245-ecp-link-witness-extraction.md`

**Interfaces:**
- Consumes the finished extractor and tests.
- Produces a verified branch ready for PR.

- [ ] **Step 1: Run issue verification**

Run:

```bash
julia --project=. -e 'include("test/expert/ecp_link_witness_general.jl")'
julia --project=. -e 'include("test/internal/ecp_mainline_fixtures.jl")'
```

Expected: both PASS.

- [ ] **Step 2: Run package verification**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: PASS.

- [ ] **Step 3: Run diff hygiene**

Run:

```bash
git diff --check
git status --short
```

Expected: no whitespace errors and only intended files changed.

- [ ] **Step 4: Commit implementation**

Run:

```bash
git add docs/superpowers/plans/2026-07-01-issue-245-ecp-link-witness-extraction.md src/algorithm/column_reduction.jl test/expert/ecp_link_witness_general.jl test/expert/ecp_link_witnesses.jl test/runtests.jl
git commit -m "Extract ECP link witness records"
```

Expected: commit succeeds with no generated dependency artifacts.

## Plan Self-Review

- Every issue requirement maps to a task or verification command.
- The plan preserves TDD by adding failing expert tests before implementation.
- The first extractor is exact and bounded, and failure is diagnostic.
- No public reducer or link realization work is included.
- No incomplete markers remain.
