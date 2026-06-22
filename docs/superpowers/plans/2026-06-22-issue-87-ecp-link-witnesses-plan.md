# Issue 87 ECP Link Witnesses Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add validated Park-Woodburn ECP link-theorem witness records for supported ordinary polynomial columns with a selected monic first entry.

**Architecture:** Extend the existing ECP reducer module with a non-exported expert witness record and replay verifier. The first implementation stage supports supplied link witness data only; extraction remains staged with a clear `ArgumentError`.

**Tech Stack:** Julia, Oscar multivariate polynomial rings, existing Suslin ECP helpers in `src/algorithm/column_reduction.jl`, and Julia `Test`.

## Global Constraints

- Extend `src/algorithm/column_reduction.jl`; do not add a parallel reducer.
- Keep `reduce_unimodular_column(v, R)` returning only factors.
- Names added for this issue remain non-exported expert/internal APIs accessed as `Suslin.<name>` in expert tests.
- Support ordinary polynomial rings only; Laurent columns stay out of scope.
- Require selected monic first entry: `selected_monic_index == 1` and `v[1]` is monic in the selected variable.
- Support `:supplied_link_witness` metadata now; missing extraction must throw `ArgumentError`.
- Do not implement the link-step operation mapping `v(b_{i-1})` to `v(b_i)`.
- Do not implement Quillen patching or the final Park-Woodburn matrix driver.
- Verification command required by issue #87: `julia --project=. -e 'include("test/expert/ecp_link_witnesses.jl")'`.
- Verification command required by Agent Desk: `julia --project=. -e 'using Pkg; Pkg.test()'`.

---

## File Structure

- Modify `src/algorithm/column_reduction.jl`: define `ECPLinkWitnessRecord`, constructor/replay helpers, and `verify_ecp_link_witness(record)::Bool`.
- Create `test/expert/ecp_link_witnesses.jl`: fixture-backed positive and negative tests for supplied link theorem witnesses.
- Modify `test/runtests.jl`: register `expert/ecp_link_witnesses.jl` in the expert group.

### Task 1: Add Link Witness Record And Verifier

**Files:**
- Modify: `src/algorithm/column_reduction.jl`

**Interfaces:**
- Consumes: `_validated_unimodular_column(v, R)`, `_is_laurent_polynomial_ring(R)`, `_ecp_normalize_variable_order(R, variable_order)`, `_ecp_generator_index(R, variable)`, `_is_monic_in_variable(p, R, variable_index)`, `_coerce_into_ring(R, value, label)`.
- Produces: `ECPLinkWitnessRecord`, `ecp_link_witness(v, R; variable_order, selected_variable, selected_monic_index, supplied_link_witness)`, `verify_ecp_link_witness(record)::Bool`.

- [ ] **Step 1: Write a failing smoke test in scratch**

Run this temporary probe before production edits to prove the API is absent:

```bash
julia --project=. -e 'using Suslin, Oscar; R,(x,y)=Oscar.polynomial_ring(QQ, ["x", "y"]); v=[x^2 + y + one(R), x, y]; Suslin.ecp_link_witness(v, R; variable_order=(x, y), selected_variable=x, supplied_link_witness=(; source=:supplied_link_witness))'
```

Expected: FAIL with `UndefVarError: ecp_link_witness not defined`.

- [ ] **Step 2: Add the record type**

Add this near the existing ECP structs:

```julia
struct ECPLinkWitnessRecord
    original_column
    ring
    variable_order
    selected_variable_index::Int
    selected_variable
    selected_monic_index::Int
    selected_monic_entry
    residue_probes
    tail_reductions
    resultants
    bezout_coefficients
    coverage_multipliers
    path_points
    metadata
    verification
end
```

- [ ] **Step 3: Add the constructor**

Add `ecp_link_witness` after `_deterministic_ecp_monicity_search` helpers or before replay helpers:

```julia
function ecp_link_witness(
    v::AbstractVector,
    R;
    variable_order = tuple(gens(R)...),
    selected_variable = nothing,
    selected_monic_index::Integer = 1,
    supplied_link_witness = nothing,
)
    _is_laurent_polynomial_ring(R) &&
        throw(ArgumentError("ECP link witnesses currently support ordinary polynomial columns only"))
    supplied_link_witness === nothing &&
        throw(ArgumentError("Park-Woodburn ECP link witness extraction is not implemented; pass supplied_link_witness metadata with source = :supplied_link_witness"))

    column = _validated_unimodular_column(v, R)
    normalized_order = _ecp_normalize_variable_order(R, variable_order)
    selected_variable = selected_variable === nothing ? first(normalized_order) : selected_variable
    selected_variable_index = _ecp_selected_variable_index(R, selected_variable)
    selected_monic_index = Int(selected_monic_index)
    selected_monic_index == 1 ||
        throw(ArgumentError("Park-Woodburn ECP link witnesses require the selected monic entry to be first"))
    _is_monic_in_variable(column[selected_monic_index], R, selected_variable_index) ||
        throw(ArgumentError("selected first entry must be monic in the selected variable"))

    metadata = (; source = _ecp_link_field(supplied_link_witness, :source))
    metadata.source == :supplied_link_witness ||
        throw(ArgumentError("supplied ECP link witness metadata must use source = :supplied_link_witness"))
    record = ECPLinkWitnessRecord(
        tuple(column...),
        R,
        tuple(normalized_order...),
        selected_variable_index,
        gens(R)[selected_variable_index],
        selected_monic_index,
        column[selected_monic_index],
        tuple(_ecp_link_field(supplied_link_witness, :residue_probes)...),
        tuple(_ecp_link_field(supplied_link_witness, :tail_reductions)...),
        tuple(_ecp_link_field(supplied_link_witness, :resultants)...),
        tuple(_ecp_link_field(supplied_link_witness, :bezout_coefficients)...),
        tuple(_ecp_link_field(supplied_link_witness, :coverage_multipliers)...),
        tuple(_ecp_link_field(supplied_link_witness, :path_points)...),
        metadata,
        nothing,
    )
    verification = _ecp_link_witness_replay_summary(record)
    verification.overall_ok ||
        throw(ArgumentError("supplied Park-Woodburn ECP link witness data failed exact replay verification"))
    stored = ECPLinkWitnessRecord(
        record.original_column,
        record.ring,
        record.variable_order,
        record.selected_variable_index,
        record.selected_variable,
        record.selected_monic_index,
        record.selected_monic_entry,
        record.residue_probes,
        record.tail_reductions,
        record.resultants,
        record.bezout_coefficients,
        record.coverage_multipliers,
        record.path_points,
        record.metadata,
        verification,
    )
    verify_ecp_link_witness(stored) ||
        throw(ArgumentError("stored Park-Woodburn ECP link witness data failed exact replay verification"))
    return stored
end
```

- [ ] **Step 4: Add replay helpers**

Add helpers with these exact behaviors:

```julia
function verify_ecp_link_witness(record)::Bool
    try
        replay = _ecp_link_witness_replay_summary(record)
        return replay.overall_ok && record.verification == replay
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end

function _ecp_link_field(source, field::Symbol)
    hasproperty(source, field) || throw(ArgumentError("supplied ECP link witness missing field $(field)"))
    return getproperty(source, field)
end

function _ecp_selected_variable_index(R, variable)
    ring_gens = collect(gens(R))
    idx = _ecp_variable_order_match_index(ring_gens, variable)
    idx === nothing && throw(ArgumentError("selected_variable must be a generator of R"))
    return idx
end
```

`_ecp_link_witness_replay_summary(record)` must compute and return a `NamedTuple`
containing at least:

```julia
(;
    overall_ok,
    metadata_ok,
    selected_monic_ok,
    lengths_ok,
    tail_reduction_ok,
    resultants_ok,
    bezout_ok,
    coverage_ok,
    path_ok,
    recomputed_tail_tilde_Gs,
    recomputed_resultants,
    coverage_total,
)
```

For each tail reduction, use `tail.lifted_tail_coefficients` and `record.original_column[2:end]` to recompute `tilde_G`; check `tail.G == tail.tilde_G == recomputed_tilde_G`.

- [ ] **Step 5: Run the scratch probe**

Re-run the Step 1 probe with a complete supplied witness from Task 2's test helper if Task 2 already exists, or skip to Task 2 and let the test drive the remaining details.

### Task 2: Add Fixture-Backed Expert Tests

**Files:**
- Create: `test/expert/ecp_link_witnesses.jl`

**Interfaces:**
- Consumes: `test/fixtures/ecp_column_cases.jl`, `Suslin.ecp_link_witness`, `Suslin.verify_ecp_link_witness`.
- Produces: Positive tests for `ecp-variable-change-monic-gf2` and `ecp-monic-first-entry-qq`, and negative controls for corrupted resultant, Bezout coefficient, coverage multiplier, and path point.

- [ ] **Step 1: Write the failing test file**

Create `test/expert/ecp_link_witnesses.jl` with:

```julia
using Test
using Oscar
using Suslin

const ECP_COLUMN_CATALOG_PATH = joinpath(@__DIR__, "..", "fixtures", "ecp_column_cases.jl")
include(ECP_COLUMN_CATALOG_PATH)

function _case_by_id(id::AbstractString)
    return ECPColumnFixtureCatalog.cases_by_id()[id]
end

function _column(entry)
    return [getproperty(entry.entries, name) for name in entry.column_order]
end

function _gf2_link_witness(entry)
    R = entry.ring.object
    x, y = entry.ring.generators
    v = _column(entry)
    G = y * v[2] + v[3]
    return (;
        source = :supplied_link_witness,
        residue_probes = ((; id = :gf2_fixture_probe, kind = :deterministic_fixture, maximal_ideal_generators = (y,)),),
        tail_reductions = ((;
            probe_id = :gf2_fixture_probe,
            G,
            lifted_tail_coefficients = (y, one(R)),
            tilde_G = G,
        ),),
        resultants = (one(R),),
        bezout_coefficients = ((; f = x, h = one(R)),),
        coverage_multipliers = (one(R),),
        path_points = (zero(R), x),
    )
end

function _qq_link_witness(entry)
    R = entry.ring.object
    x, y = entry.ring.generators
    return (;
        source = :supplied_link_witness,
        residue_probes = (
            (; id = :qq_y_probe, kind = :deterministic_fixture, maximal_ideal_generators = (y,)),
            (; id = :qq_x_probe, kind = :deterministic_fixture, maximal_ideal_generators = (x,)),
        ),
        tail_reductions = (
            (; probe_id = :qq_y_probe, G = y, lifted_tail_coefficients = (zero(R), one(R)), tilde_G = y),
            (; probe_id = :qq_x_probe, G = x, lifted_tail_coefficients = (one(R), zero(R)), tilde_G = x),
        ),
        resultants = (y^2, y + one(R)),
        bezout_coefficients = (
            (; f = zero(R), h = y),
            (; f = one(R), h = -x),
        ),
        coverage_multipliers = (one(R), one(R) - y),
        path_points = (zero(R), y^2 * x, x),
    )
end
```

Add tests asserting no unit entries, source metadata, replay fields, and `Suslin.verify_ecp_link_witness(record) == true`.

- [ ] **Step 2: Verify RED**

Run:

```bash
julia --project=. -e 'include("test/expert/ecp_link_witnesses.jl")'
```

Expected: FAIL before Task 1 is complete because `Suslin.ecp_link_witness` is not defined or verifier details are missing.

- [ ] **Step 3: Add negative controls**

Add helpers:

```julia
function _replace_tuple_entry(values::Tuple, idx::Int, value)
    return ntuple(j -> j == idx ? value : values[j], length(values))
end

function _mutate_witness(witness; resultants = witness.resultants, bezout_coefficients = witness.bezout_coefficients, coverage_multipliers = witness.coverage_multipliers, path_points = witness.path_points)
    return merge(witness, (;
        resultants,
        bezout_coefficients,
        coverage_multipliers,
        path_points,
    ))
end
```

For the QQ witness, assert each mutation is rejected by construction:

```julia
@test_throws ArgumentError Suslin.ecp_link_witness(v, R; variable_order = entry.ring.generators, selected_variable = x, supplied_link_witness = _mutate_witness(witness; resultants = (y^2 + one(R), y + one(R))))
@test_throws ArgumentError Suslin.ecp_link_witness(v, R; variable_order = entry.ring.generators, selected_variable = x, supplied_link_witness = _mutate_witness(witness; bezout_coefficients = ((; f = one(R), h = y), witness.bezout_coefficients[2])))
@test_throws ArgumentError Suslin.ecp_link_witness(v, R; variable_order = entry.ring.generators, selected_variable = x, supplied_link_witness = _mutate_witness(witness; coverage_multipliers = (one(R), one(R))))
@test_throws ArgumentError Suslin.ecp_link_witness(v, R; variable_order = entry.ring.generators, selected_variable = x, supplied_link_witness = _mutate_witness(witness; path_points = (zero(R), y^2 * x + one(R), x)))
```

Also tamper a stored valid record and assert `!Suslin.verify_ecp_link_witness(tampered)`.

- [ ] **Step 4: Verify GREEN**

Run:

```bash
julia --project=. -e 'include("test/expert/ecp_link_witnesses.jl")'
```

Expected: PASS.

### Task 3: Register Expert Test And Final Verification

**Files:**
- Modify: `test/runtests.jl`

**Interfaces:**
- Consumes: `test/expert/ecp_link_witnesses.jl`.
- Produces: expert test registration.

- [ ] **Step 1: Register the expert test**

Add `"expert/ecp_link_witnesses.jl"` to the `expert_tests` list in `test/runtests.jl`, near the other ECP expert files.

- [ ] **Step 2: Run issue verification**

Run:

```bash
julia --project=. -e 'include("test/expert/ecp_link_witnesses.jl")'
```

Expected: PASS.

- [ ] **Step 3: Run default package verification**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: PASS.

- [ ] **Step 4: Run expert harness**

Run:

```bash
julia --project=. test/runtests.jl expert
```

Expected: PASS.

- [ ] **Step 5: Check whitespace**

Run:

```bash
git diff --check
```

Expected: no output and exit 0.

---

## Plan Self-Review

- Spec coverage: The plan covers the record, supplied witness boundary, monicity, tail reductions, resultants, Bezout identities, coverage, paths, negative controls, and registration.
- Placeholder scan: No placeholder steps remain.
- Type consistency: Public names are `ECPLinkWitnessRecord`, `ecp_link_witness`, and `verify_ecp_link_witness`; tests use the same names through `Suslin.<name>`.
