# Issue 208 Murthy Local Input Context Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an internal, replayable Murthy local input context for special-form `SL_3` matrices.

**Architecture:** Add a non-exported context record, constructor, and verifier beside the existing `SL3Local...` records in `src/algorithm/sl3_local.jl`. The context validates ordinary `QQ[X]` cases without local evidence and validates `QQ[u, X]` local-contract cases only when supplied local-unit or Bezout evidence proves the needed local units.

**Tech Stack:** Julia, Oscar polynomial rings, Suslin internal local `SL_3` helpers, Test stdlib.

## Global Constraints

- Do not produce elementary factors in this issue.
- Do not implement q(0)-unit, q(0)-nonunit, Quillen denominator covers, or global patching.
- Keep the new context internal: do not add exports in `src/Suslin.jl`.
- Reuse the #206 local-unit witness schema exactly.
- The constructor and verifier must recompute target shape, determinant one, selected variable, monic `p`, degrees, constant terms, and local/global unit classification.
- Existing ordinary `QQ[X]` Murthy cases must validate without local-unit witnesses.
- `QQ[u, X]` local-contract examples that need local units must validate only when local-unit evidence is supplied.
- Negative controls must reject missing local-unit evidence, wrong selected variable, non-monic `p`, corrupted stored `q0`, and corrupted supplied local-unit or Bezout witnesses.
- Focused verification command is `julia --project=. -e 'include("test/expert/sl3_local_murthy_context.jl")'`.
- Package verification command is `julia --project=. -e 'using Pkg; Pkg.test()'`.

---

## File Structure

- Modify `src/algorithm/sl3_local.jl`: define `SL3LocalMurthyInputContext`, matrix and entry constructors, verifier, witness normalization, local-unit witness verification, and optional split/Bezout witness verification helpers.
- Create `test/expert/sl3_local_murthy_context.jl`: focused context tests for ordinary and local-contract examples plus tamper cases.
- Modify `test/runtests.jl`: register the new expert test file.
- Keep `src/Suslin.jl` unchanged except for natural includes already present.

---

### Task 1: Context Tests

**Files:**
- Create: `test/expert/sl3_local_murthy_context.jl`
- Modify: `test/runtests.jl`

**Interfaces:**
- Consumes existing fixture catalog `test/fixtures/sl3_murthy_gupta_cases.jl`.
- Produces tests against `Suslin.SL3LocalMurthyInputContext`, `Suslin.sl3_local_murthy_input_context`, and `Suslin.verify_sl3_local_murthy_input_context`.

- [ ] **Step 1: Write the failing context test**

Create `test/expert/sl3_local_murthy_context.jl` with these test helpers and assertions:

```julia
using Test
using Suslin
using Oscar

const SL3_CONTEXT_FIXTURE_PATH =
    joinpath(@__DIR__, "..", "fixtures", "sl3_murthy_gupta_cases.jl")

function _context_as_namedtuple(context)
    names = propertynames(context)
    return NamedTuple{names}(Tuple(getproperty(context, name) for name in names))
end

function _captured_error(f)
    try
        f()
        return nothing
    catch err
        return err
    end
end

@testset "Murthy local input context" begin
    include(SL3_CONTEXT_FIXTURE_PATH)
    catalog = SL3MurthyGuptaFixtureCatalog.catalog()
    by_id = Dict(entry.id => entry for entry in catalog.cases)

    ordinary = by_id["mg-q0-unit-recursion"]
    ordinary_context = Suslin.sl3_local_murthy_input_context(
        ordinary.target,
        ordinary.variable,
    )
    @test ordinary_context.R === ordinary.ring.object
    @test ordinary_context.X == ordinary.variable
    @test ordinary_context.entries == ordinary.entries
    @test ordinary_context.target == ordinary.target
    @test ordinary_context.determinant == one(ordinary.ring.object)
    @test ordinary_context.degree_p == 1
    @test ordinary_context.degree_q == 0
    @test ordinary_context.p0 == one(ordinary.ring.object)
    @test ordinary_context.q0 == one(ordinary.ring.object)
    @test ordinary_context.p_monic == true
    @test ordinary_context.global_units.q0 == true
    @test ordinary_context.local_units.q0 == true
    @test Suslin.verify_sl3_local_murthy_input_context(ordinary_context)

    local_q0_unit = by_id["mg-local-q0-unit-at-u"]
    local_witness = first(local_q0_unit.witnesses)
    local_context = Suslin.sl3_local_murthy_input_context(
        local_q0_unit.entries.p,
        local_q0_unit.entries.q,
        local_q0_unit.entries.r,
        local_q0_unit.entries.s,
        local_q0_unit.variable;
        witness = local_witness,
    )
    @test local_context.target == local_q0_unit.target
    @test local_context.p0 == local_witness.p0
    @test local_context.q0 == local_witness.q0
    @test local_context.global_units.q0 == false
    @test local_context.local_units.q0 == true
    @test local_context.local_unit_witnesses.q0 == local_witness.local_unit_witness
    @test Suslin.verify_sl3_local_murthy_input_context(local_context)

    missing_local_witness_err = _captured_error(() -> Suslin.sl3_local_murthy_input_context(
        local_q0_unit.target,
        local_q0_unit.variable,
    ))
    @test missing_local_witness_err isa ArgumentError
    @test occursin("local-unit witness", sprint(showerror, missing_local_witness_err))

    wrong_variable_err = _captured_error(() -> Suslin.sl3_local_murthy_input_context(
        local_q0_unit.target,
        first(gens(local_q0_unit.ring.object));
        witness = local_witness,
    ))
    @test wrong_variable_err isa ArgumentError
    @test occursin("p must be monic", sprint(showerror, wrong_variable_err))

    R, (X,) = Oscar.polynomial_ring(QQ, ["X"])
    nonmonic_err = _captured_error(() -> Suslin.sl3_local_murthy_input_context(
        2 * X + one(R),
        X,
        R(2),
        one(R),
        X,
    ))
    @test nonmonic_err isa ArgumentError
    @test occursin("p must be monic", sprint(showerror, nonmonic_err))

    corrupted_q0_context = Suslin.SL3LocalMurthyInputContext(
        values(merge(
            _context_as_namedtuple(local_context),
            (; q0 = local_context.q0 + one(local_context.R)),
        ))...,
    )
    @test !Suslin.verify_sl3_local_murthy_input_context(corrupted_q0_context)

    corrupted_local_witness = merge(
        local_witness,
        (;
            local_unit_witness = merge(
                local_witness.local_unit_witness,
                (; residue_difference_coefficients = (zero(local_context.R),)),
            ),
        ),
    )
    @test_throws ArgumentError Suslin.sl3_local_murthy_input_context(
        local_q0_unit.target,
        local_q0_unit.variable;
        witness = corrupted_local_witness,
    )

    local_nonunit = by_id["mg-local-q0-nonunit-bezout-at-u"]
    bezout_witness = first(local_nonunit.witnesses)
    bezout_context = Suslin.sl3_local_murthy_input_context(
        local_nonunit.target,
        local_nonunit.variable;
        witness = bezout_witness,
    )
    @test bezout_context.global_units.branch_unit == false
    @test bezout_context.local_units.branch_unit == true
    @test bezout_context.bezout_witness == bezout_witness
    @test Suslin.verify_sl3_local_murthy_input_context(bezout_context)

    corrupted_bezout_witness = merge(
        bezout_witness,
        (; p_prime = bezout_witness.p_prime + one(local_nonunit.ring.object)),
    )
    @test_throws ArgumentError Suslin.sl3_local_murthy_input_context(
        local_nonunit.target,
        local_nonunit.variable;
        witness = corrupted_bezout_witness,
    )
end
```

Add `"expert/sl3_local_murthy_context.jl"` to the expert file list in
`test/runtests.jl` immediately after `"expert/sl3_local_murthy_resultant.jl"`.

- [ ] **Step 2: Run focused test to verify RED**

Run:

```bash
julia --project=. -e 'include("test/expert/sl3_local_murthy_context.jl")'
```

Expected: FAIL with `UndefVarError` for `sl3_local_murthy_input_context` or
`SL3LocalMurthyInputContext`.

- [ ] **Step 3: Commit failing test**

```bash
git add test/expert/sl3_local_murthy_context.jl test/runtests.jl
git commit -m "test: cover Murthy local input context"
```

---

### Task 2: Context Constructor And Verifier

**Files:**
- Modify: `src/algorithm/sl3_local.jl`

**Interfaces:**
- Produces `SL3LocalMurthyInputContext`.
- Produces `sl3_local_murthy_input_context(A, X; witness, local_unit_witnesses, split_witness, bezout_witness)`.
- Produces `sl3_local_murthy_input_context(p, q, r, s, X; witness, local_unit_witnesses, split_witness, bezout_witness)`.
- Produces `verify_sl3_local_murthy_input_context(context)::Bool`.

- [ ] **Step 1: Add the context struct and public internal wrappers**

In `src/algorithm/sl3_local.jl`, place this after `SL3LocalMurthyQ0NonunitReduction`:

```julia
struct SL3LocalMurthyInputContext
    R
    X
    var_idx::Int
    entries::NamedTuple
    target
    determinant
    degree_p::Int
    degree_q::Int
    p0
    q0
    p_monic::Bool
    global_units::NamedTuple
    local_units::NamedTuple
    local_unit_witnesses::NamedTuple
    split_witness
    bezout_witness
end
```

Add constructor wrappers before `realize_sl3_local`:

```julia
function sl3_local_murthy_input_context(A, X; witness=nothing, local_unit_witnesses=(;),
        split_witness=nothing, bezout_witness=nothing)
    nrows(A) == 3 || throw(ArgumentError("Murthy local input context requires a 3x3 special-form matrix"))
    ncols(A) == 3 || throw(ArgumentError("Murthy local input context requires a 3x3 special-form matrix"))
    entries = _sl3_local_target_entries(A)
    entries === nothing && throw(ArgumentError("Murthy local input context requires a special-form SL_3 target"))
    return sl3_local_murthy_input_context(
        entries.p,
        entries.q,
        entries.r,
        entries.s,
        X;
        witness,
        local_unit_witnesses,
        split_witness,
        bezout_witness,
    )
end

function sl3_local_murthy_input_context(p, q, r, s, X; witness=nothing,
        local_unit_witnesses=(;), split_witness=nothing, bezout_witness=nothing)
    return _sl3_local_murthy_input_context(
        p,
        q,
        r,
        s,
        X;
        witness,
        local_unit_witnesses,
        split_witness,
        bezout_witness,
    )
end

function verify_sl3_local_murthy_input_context(context)::Bool
    try
        return _sl3_local_murthy_input_context_verification(context).overall_ok
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end
```

- [ ] **Step 2: Implement context construction helpers**

Add helpers near the other `_sl3_local_*` helpers:

```julia
function _sl3_local_murthy_input_context(p, q, r, s, X; witness, local_unit_witnesses,
        split_witness, bezout_witness)
    R = parent(p)
    parent(q) === R || throw(ArgumentError("Murthy local input context entries and X must lie in the same polynomial ring"))
    parent(r) === R || throw(ArgumentError("Murthy local input context entries and X must lie in the same polynomial ring"))
    parent(s) === R || throw(ArgumentError("Murthy local input context entries and X must lie in the same polynomial ring"))
    parent(X) === R || throw(ArgumentError("Murthy local input context entries and X must lie in the same polynomial ring"))

    ring_gens = collect(gens(R))
    var_idx = findfirst(isequal(X), ring_gens)
    var_idx === nothing && throw(ArgumentError("Murthy local input context variable must be a ring generator"))
    _is_laurent_polynomial_ring(R) &&
        throw(ArgumentError("Murthy local input context is only supported for ordinary polynomial rings"))

    target = _sl3_local_special_form_target(R, p, q, r, s)
    determinant = det(target)
    determinant == one(R) || throw(ArgumentError("Murthy local input context target must have determinant 1"))
    p_monic = _is_monic_in_variable(p, var_idx, R)
    p_monic || throw(ArgumentError("Murthy local input context p must be monic in X"))

    normalized_witnesses = _sl3_local_murthy_normalize_witness_data(
        witness,
        local_unit_witnesses,
        split_witness,
        bezout_witness,
    )
    p0 = _sl3_local_constant_coefficient(p, var_idx, R)
    q0 = _sl3_local_constant_coefficient(q, var_idx, R)
    degree_p = degree(p, var_idx)
    degree_q = degree(q, var_idx)
    global_units = _sl3_local_murthy_global_units(p, q, r, s, p0, q0)
    local_units = _sl3_local_murthy_local_units(
        R,
        X,
        p,
        q,
        p0,
        q0,
        degree_p,
        degree_q,
        normalized_witnesses.local_unit_witnesses,
        normalized_witnesses.bezout_witness,
    )
    _sl3_local_murthy_validate_required_local_evidence(
        R,
        degree_p,
        degree_q,
        q0,
        global_units,
        local_units,
        normalized_witnesses.bezout_witness,
    )
    _sl3_local_murthy_verify_split_witness(R, normalized_witnesses.split_witness)

    context = SL3LocalMurthyInputContext(
        R,
        X,
        var_idx,
        (; p, q, r, s),
        target,
        determinant,
        degree_p,
        degree_q,
        p0,
        q0,
        p_monic,
        global_units,
        local_units,
        normalized_witnesses.local_unit_witnesses,
        normalized_witnesses.split_witness,
        normalized_witnesses.bezout_witness,
    )
    verify_sl3_local_murthy_input_context(context) ||
        error("internal Murthy local input context verification failed")
    return context
end
```

Implement `_sl3_local_murthy_normalize_witness_data`, `_sl3_local_murthy_global_units`,
`_sl3_local_murthy_local_units`, `_sl3_local_murthy_validate_required_local_evidence`,
`_sl3_local_murthy_verify_local_unit_witness`, `_sl3_local_murthy_bezout_data`,
and `_sl3_local_murthy_verify_split_witness` following the spec. Use `merge` on
named tuples to combine explicit witnesses with extracted #206 witness fields.

- [ ] **Step 3: Implement verifier**

Add `_sl3_local_murthy_input_context_verification(context)` that:

1. checks target size and special-form entries;
2. checks `context.R === base_ring(context.target)`;
3. checks `parent(context.X) === context.R` and `context.var_idx` matches
   `findfirst(isequal(context.X), collect(gens(context.R)))`;
4. recomputes determinant, degrees, `p0`, `q0`, monicity, global-unit named
   tuple, local-unit named tuple, split witness validity, and Bezout witness
   validity;
5. returns a named tuple with `overall_ok` and individual boolean fields.

The verifier must return `false` through `verify_sl3_local_murthy_input_context`
when any field is corrupted.

- [ ] **Step 4: Run focused test to verify GREEN**

Run:

```bash
julia --project=. -e 'include("test/expert/sl3_local_murthy_context.jl")'
```

Expected: PASS.

- [ ] **Step 5: Commit implementation**

```bash
git add src/algorithm/sl3_local.jl
git commit -m "feat: add Murthy local input context"
```

---

### Task 3: Verification Sweep

**Files:**
- Modify only if Task 2 reveals a registration or verifier gap: `test/runtests.jl`, `src/algorithm/sl3_local.jl`, `test/expert/sl3_local_murthy_context.jl`.

**Interfaces:**
- Consumes the context constructor and verifier from Task 2.
- Produces verified package state and any needed polish fixes.

- [ ] **Step 1: Run focused context verification**

Run:

```bash
julia --project=. -e 'include("test/expert/sl3_local_murthy_context.jl")'
```

Expected: PASS.

- [ ] **Step 2: Run package verification**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: PASS.

- [ ] **Step 3: Inspect git diff for unrelated edits**

Run:

```bash
git status --short
git diff --stat origin/main...HEAD
```

Expected: only the design, plan, local context implementation, context test, and test registration are present.

- [ ] **Step 4: Commit verification fixes if needed**

If Step 1 or Step 2 required changes, commit them:

```bash
git add src/algorithm/sl3_local.jl test/expert/sl3_local_murthy_context.jl test/runtests.jl
git commit -m "test: verify Murthy local input context"
```

If no changes were needed, do not create an empty commit.
