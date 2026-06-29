# Issue 206 Murthy Local Contract Catalog Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Extend the Murthy-Gupta fixture catalog with Park-Woodburn Section 5 local-contract cases and negative controls without changing solver behavior.

**Architecture:** Keep all new data in `test/fixtures/sl3_murthy_gupta_cases.jl` and all checks in `test/internal/sl3_murthy_gupta_fixtures.jl`. Existing ordinary `QQ[X]` cases remain stable pass cases, while new `QQ[u, X]` local-contract cases are marked as staged current-solver failures and validated through explicit q-degree, local-unit, and Bezout witness equations.

**Tech Stack:** Julia, Oscar polynomial rings, Suslin local `SL_3` helpers, Test stdlib.

## Global Constraints

- Do not implement new Murthy solving, denominator-aware replay, Quillen patching, ECP, recursive `SL_n`, or public documentation claims.
- Keep existing fixture ids stable.
- Add new local-contract ids rather than reclassifying old pass cases.
- New local-contract cases must set `local_contract = true`.
- New local-contract cases must set `expected_current_solver = (; status = :staged_fail, message_substring = "staged local SL_3 solver failure")`.
- New local-contract cases must set `consumer_issue_ids = ("#182", "#208", "#207", "#209", "#210")`.
- Existing ordinary `QQ[X]` pass cases must still pass through the current solver.
- The validator must reconstruct every positive target from `(p, q, r, s)`, check determinant one, check `p` monic in the selected variable, and validate q-degree, split-lemma, q(0)-unit, and q(0)-nonunit witness equations exactly.
- The validator must reject negative controls for non-monic `p`, determinant not one, corrupted split witness, corrupted local-unit witness, and corrupted Bezout equality.
- Focused validator command is `julia --project=. -e 'include("test/internal/sl3_murthy_gupta_fixtures.jl")'`.
- Package verification command is `julia --project=. -e 'using Pkg; Pkg.test()'`.

---

## File Structure

- Modify `test/internal/sl3_murthy_gupta_fixtures.jl`: selected-variable constant helper, local-context/local-unit witness validation, local-contract metadata checks, q(0)-unit and q(0)-nonunit witness extensions, required id checks, and negative-control tests.
- Modify `test/fixtures/sl3_murthy_gupta_cases.jl`: generalized ring metadata, three `QQ[u, X]` local-contract cases, local witness constructors, and catalog negative controls.

---

### Task 1: Local Contract Catalog And Validator

**Files:**
- Modify: `test/internal/sl3_murthy_gupta_fixtures.jl`
- Modify: `test/fixtures/sl3_murthy_gupta_cases.jl`

**Interfaces:**
- Consumes: `SL3MurthyGuptaFixtureCatalog.catalog()` returning `(; ring, cases)` today.
- Produces: `SL3MurthyGuptaFixtureCatalog.catalog()` returning `(; ring, cases, negative_controls)` while preserving the existing `cases` interface.
- Produces: internal validator functions that continue to expose `validate_sl3_murthy_gupta_fixture(entry)` and `validate_sl3_murthy_gupta_fixture_catalog(catalog)`.

- [ ] **Step 1: Write failing validator expectations**

In `test/internal/sl3_murthy_gupta_fixtures.jl`, add required-id constants near the existing branch constants:

```julia
const REQUIRED_SL3_MURTHY_GUPTA_LOCAL_CONTRACT_IDS = Set([
    "mg-local-q-degree-qq-u-x",
    "mg-local-q0-unit-at-u",
    "mg-local-q0-nonunit-bezout-at-u",
])

const REQUIRED_SL3_MURTHY_GUPTA_NEGATIVE_IDS = Set([
    "mg-negative-nonmonic-p",
    "mg-negative-determinant-not-one",
    "mg-negative-corrupted-split-witness",
    "mg-negative-corrupted-local-unit-witness",
    "mg-negative-corrupted-bezout-equality",
])
```

In the `"Murthy-Gupta local SL3 fixture catalog"` testset, after `by_id` is built, add:

```julia
@test REQUIRED_SL3_MURTHY_GUPTA_LOCAL_CONTRACT_IDS ⊆ Set(keys(by_id))

local_contract_entries = [by_id[id] for id in REQUIRED_SL3_MURTHY_GUPTA_LOCAL_CONTRACT_IDS]
@test length(local_contract_entries) == 3
for entry in local_contract_entries
    @test entry.local_contract == true
    @test entry.expected_current_solver.status == :staged_fail
    @test ("#182", "#208", "#207", "#209", "#210") == entry.consumer_issue_ids
end

@test hasproperty(catalog, :negative_controls)
negative_by_id = Dict(entry.id => entry for entry in catalog.negative_controls)
@test REQUIRED_SL3_MURTHY_GUPTA_NEGATIVE_IDS ⊆ Set(keys(negative_by_id))
for entry in values(negative_by_id)
    @test_throws ArgumentError validate_sl3_murthy_gupta_fixture(entry)
end
```

- [ ] **Step 2: Run focused test to verify RED**

Run:

```bash
julia --project=. -e 'include("test/internal/sl3_murthy_gupta_fixtures.jl")'
```

Expected: FAIL because the new local-contract ids and `negative_controls` field have not been added to the catalog.

- [ ] **Step 3: Implement selected-variable constants and local witnesses**

In `test/internal/sl3_murthy_gupta_fixtures.jl`, add a helper that evaluates the `X = 0` coefficient while preserving other variables:

```julia
function _sl3_mg_constant_in_variable(value, X)
    R = parent(value)
    parent(X) == R || throw(ArgumentError("value and variable must lie in same ring"))
    var_idx = findfirst(isequal(X), collect(gens(R)))
    var_idx === nothing && throw(ArgumentError("variable must be a generator of the ambient ring"))

    vars = collect(gens(R))
    constant = zero(R)
    for (coeff, exponents) in zip(AbstractAlgebra.coefficients(value), AbstractAlgebra.exponent_vectors(value))
        exponents[var_idx] == 0 || continue
        term = R(coeff)
        for idx in eachindex(vars)
            idx == var_idx && continue
            exponent = exponents[idx]
            exponent == 0 || (term *= vars[idx]^exponent)
        end
        constant += term
    end
    return constant
end
```

Update q(0)-unit and q(0)-nonunit checks to use `_sl3_mg_constant_in_variable(value, entry.variable)` instead of the old all-variables constant term.

Add exact local witness validation:

```julia
function _sl3_mg_assert_local_unit_witness(entry, witness, expected_unit; label)
    R = entry.ring.object
    context = _sl3_mg_field(witness, :context)
    unit = _sl3_mg_field(witness, :unit)
    residue_unit = _sl3_mg_field(witness, :residue_unit)
    residue_inverse = _sl3_mg_field(witness, :residue_inverse)
    generators = _sl3_mg_field(witness, :maximal_ideal_generators)
    coefficients = _sl3_mg_field(witness, :residue_difference_coefficients)
    global_unit = hasproperty(witness, :global_unit) ? witness.global_unit : is_unit(unit)

    context.kind == :localization_at_maximal_ideal ||
        throw(ArgumentError("fixture $(entry.id) $(label) local context kind is unsupported"))
    context.selected_variable == entry.variable ||
        throw(ArgumentError("fixture $(entry.id) $(label) local context variable mismatch"))
    unit == expected_unit ||
        throw(ArgumentError("fixture $(entry.id) $(label) local unit does not match expected value"))
    length(generators) == length(coefficients) ||
        throw(ArgumentError("fixture $(entry.id) $(label) local witness generator/coefficient length mismatch"))
    all(parent(generator) === R for generator in generators) ||
        throw(ArgumentError("fixture $(entry.id) $(label) local witness generator ring mismatch"))
    all(parent(coefficient) === R for coefficient in coefficients) ||
        throw(ArgumentError("fixture $(entry.id) $(label) local witness coefficient ring mismatch"))

    difference = zero(R)
    for (coefficient, generator) in zip(coefficients, generators)
        difference += coefficient * generator
    end
    unit - residue_unit == difference ||
        throw(ArgumentError("fixture $(entry.id) $(label) local residue equation failed"))
    residue_unit * residue_inverse == one(R) ||
        throw(ArgumentError("fixture $(entry.id) $(label) local residue inverse equation failed"))
    is_unit(unit) == global_unit ||
        throw(ArgumentError("fixture $(entry.id) $(label) global unit flag is incorrect"))
    return true
end
```

Extend `_sl3_mg_assert_metadata(entry)` so entries with `local_contract = true` must have staged current-solver status and the exact consumer issue tuple from the global constraints.

Extend `_sl3_mg_assert_q0_unit_witness` with two schemas:

```julia
if hasproperty(witness, :local_unit_witness)
    p0 == _sl3_mg_constant_in_variable(p, entry.variable) || throw(...)
    q0 == _sl3_mg_constant_in_variable(q, entry.variable) || throw(...)
    _sl3_mg_assert_local_unit_witness(entry, witness.local_unit_witness, q0; label = "q0-unit")
    _sl3_mg_degree_in_variable(q, entry.variable) < _sl3_mg_degree_in_variable(p, entry.variable) ||
        throw(ArgumentError("fixture $(entry.id) q0-unit local witness must already satisfy q-degree guard"))
    return true
end
```

Keep the existing global-unit schema for old `QQ[X]` pass cases.

Extend `_sl3_mg_assert_q0_nonunit_bezout_witness` so `resultant` and `branch_unit` may be validated by `resultant_unit_witness` or `branch_unit_witness` when they are not global units.

- [ ] **Step 4: Implement catalog cases**

In `test/fixtures/sl3_murthy_gupta_cases.jl`, generalize ring metadata helpers to accept descriptions and variables without changing existing call sites:

```julia
function _ring_metadata(description, R, generators)
    return (;
        description = description,
        object = R,
        generator_names = tuple((string(generator) for generator in generators)...),
        generators = generators,
    )
end

function _ring_metadata(R, X)
    return _ring_metadata("QQ[X]", R, (X,))
end

function _ring_constructor_metadata(coefficient = "QQ", variables = ("X",))
    return (;
        function_name = :polynomial_ring,
        coefficient = coefficient,
        variables = variables,
    )
end
```

Allow `_case(...)` to accept keyword overrides:

```julia
function _case(id, branch, variable, entries, target, witness, expected_current_solver, consumer_issue_ids;
        ring_constructor = _ring_constructor_metadata(),
        ring = _ring_metadata(parent(variable), variable),
        extra = (;))
    base = (;
        id,
        branch,
        ring_constructor,
        ring,
        variable,
        entries,
        target,
        murthy_path = true,
        expected_current_solver,
        witnesses = witness,
        source_refs = ("Park-Woodburn arXiv:alg-geom/9405003 section 5",),
        consumer_issue_ids,
    )
    return merge(base, extra)
end
```

Inside `catalog()`, create `RU, (u, UX) = Oscar.polynomial_ring(QQ, ["u", "X"])`, local metadata, and this helper:

```julia
local_context = (;
    kind = :localization_at_maximal_ideal,
    description = "QQ[u] localized at (u), with X as the Section 5 variable",
    selected_variable = UX,
    maximal_ideal_generators = (u,),
    residue_description = "u => 0",
)

function _local_unit(unit, residue_unit, residue_inverse, generators, coefficients)
    return (;
        context = merge(local_context, (; maximal_ideal_generators = generators)),
        unit,
        residue_unit,
        residue_inverse,
        maximal_ideal_generators = generators,
        residue_difference_coefficients = coefficients,
        global_unit = false,
    )
end
```

Add the multivariate q-degree case:

```julia
p = UX^2 + u * UX + one(RU)
q = UX * p + one(RU)
r = -one(RU)
s = -UX
witness = (;
    quotient = UX,
    remainder = one(RU),
    normalized_s = zero(RU),
)
```

Add the local q(0)-unit case:

```julia
p = UX^2 + (u + one(RU)) * UX + one(RU)
q = UX + u + one(RU)
r = UX + p * UX
s = p
witness = (;
    p0 = one(RU),
    q0 = u + one(RU),
    local_unit_witness = _local_unit(u + one(RU), one(RU), one(RU), (u,), (one(RU),)),
    formal_right_e21_coefficient = "-1/(1 + u)",
)
```

Add the q(0)-nonunit Bezout/resultant case:

```julia
p = UX^2 + u * UX + one(RU)
q = UX + u
r = UX + p * UX
s = p
witness = (;
    p0 = one(RU),
    q0 = u,
    p_prime = one(RU),
    q_prime = UX,
    resultant = one(RU),
    p_prime_degree = 0,
    q_prime_degree = 1,
    branch_unit = u + one(RU),
    branch_unit_witness = _local_unit(u + one(RU), one(RU), one(RU), (u,), (one(RU),)),
    case1_entries = (;
        p = p + UX,
        q = q + one(RU),
        r = UX,
        s = one(RU),
    ),
)
```

Each new case uses:

```julia
expected_local_failure = (; status = :staged_fail, message_substring = "staged local SL_3 solver failure")
local_consumer_issues = ("#182", "#208", "#207", "#209", "#210")
extra = (;
    local_contract = true,
    requires_local_units = true_or_false,
    requires_bezout_witness = true_or_false,
)
```

- [ ] **Step 5: Implement catalog negative controls**

In `test/fixtures/sl3_murthy_gupta_cases.jl`, add:

```julia
function _negative_control(id, base_case_id, reason, entry)
    return merge(entry, (; id, base_case_id, reason))
end
```

Return `negative_controls` from `catalog()` with these entries:

```julia
nonmonic_entries = (; p = 2 * X + one(R), q = X, r = R(2), s = one(R))
nonmonic = _case(
    "mg-negative-nonmonic-p",
    :open_slice_control,
    X,
    nonmonic_entries,
    _target(R, nonmonic_entries.p, nonmonic_entries.q, nonmonic_entries.r, nonmonic_entries.s),
    (),
    now_supported,
    ("#206",),
)

det_bad_entries = (; p = X + one(R), q = zero(R), r = zero(R), s = one(R))
det_bad = _case(
    "mg-negative-determinant-not-one",
    :open_slice_control,
    X,
    det_bad_entries,
    _target(R, det_bad_entries.p, det_bad_entries.q, det_bad_entries.r, det_bad_entries.s),
    (),
    now_supported,
    ("#206",),
)
```

For the remaining controls, merge existing positive cases:

```julia
split_negative = _negative_control(
    "mg-negative-corrupted-split-witness",
    split_lemma_case.id,
    "split witness a no longer reconstructs the target",
    merge(split_lemma_case, (; witnesses = (merge(first(split_lemma_case.witnesses), (; a = first(split_lemma_case.witnesses).a + one(R))),))),
)

local_unit_negative = _negative_control(
    "mg-negative-corrupted-local-unit-witness",
    local_q0_unit_case.id,
    "local unit residue equation is corrupted",
    merge(local_q0_unit_case, (; witnesses = (merge(first(local_q0_unit_case.witnesses), (;
        local_unit_witness = merge(first(local_q0_unit_case.witnesses).local_unit_witness, (;
            residue_difference_coefficients = (zero(RU),),
        )),
    )),))),
)

bezout_negative = _negative_control(
    "mg-negative-corrupted-bezout-equality",
    local_q0_nonunit_case.id,
    "p_prime*p - q_prime*q no longer equals the resultant",
    merge(local_q0_nonunit_case, (; witnesses = (merge(first(local_q0_nonunit_case.witnesses), (; p_prime = one(RU) + u)),))),
)
```

Return:

```julia
return (;
    ring = _ring_metadata(R, X),
    cases = [
        existing_cases...,
        local_q_degree_case,
        local_q0_unit_case,
        local_q0_nonunit_case,
    ],
    negative_controls = [
        nonmonic,
        det_bad,
        split_negative,
        local_unit_negative,
        bezout_negative,
    ],
)
```

- [ ] **Step 6: Run focused test to verify GREEN**

Run:

```bash
julia --project=. -e 'include("test/internal/sl3_murthy_gupta_fixtures.jl")'
```

Expected: PASS with the Murthy-Gupta catalog testset.

- [ ] **Step 7: Run package verification**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: PASS.

- [ ] **Step 8: Commit**

Run:

```bash
git add test/fixtures/sl3_murthy_gupta_cases.jl test/internal/sl3_murthy_gupta_fixtures.jl docs/superpowers/plans/2026-06-29-issue-206-murthy-local-contract-catalog.md
git commit -m "test: extend murthy local contract catalog"
```

Expected: commit succeeds.

## Plan Self-Review

- Spec coverage: the task covers additive local-contract cases, local-unit and Bezout witness validation, negative controls, focused verification, package verification, and no solver behavior changes.
- Marker scan: no unresolved markers remain.
- Type consistency: `local_contract`, `requires_local_units`, `requires_bezout_witness`, `local_unit_witness`, `branch_unit_witness`, `negative_controls`, and the required ids are named consistently across the plan.
