# Issue 207 SL3 Local Denominator Factor Replay Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an internal denominator-aware local elementary factor record and replay verifier for SL3 local certificates without changing public ordinary factorization behavior.

**Architecture:** Keep ordinary certificate factors as matrices for all currently supported `realize_sl3_local(...)` cases. Add a separate local factor record plus replay summary in `src/algorithm/sl3_local.jl`; denominator-one records materialize to ordinary matrices, while non-one denominators verify through an explicit denominator-cleared equality and #208 local-unit witness checks.

**Tech Stack:** Julia, Oscar.jl matrices and polynomial rings, existing Suslin SL3 local certificate helpers, `Test`.

## Global Constraints

- Do not change public `elementary_matrix`, `verify_factorization`, or ordinary `realize_sl3_local(...)` semantics.
- Do not route q(0)-unit or q(0)-nonunit Murthy branches through this representation yet.
- Keep the new API internal and unexported.
- Reuse the #208 local-unit witness verifier and schema.
- Avoid storing an ambiguous mixture of raw matrices and local factor records without an explicit adapter layer.
- Verification commands required by issue #207: `julia --project=. -e 'include("test/expert/sl3_local_local_factors.jl")'` and `julia --project=. -e 'include("test/expert/sl3_local_certificate.jl")'`.
- Package verification required by Agent Desk: `julia --project=. -e 'using Pkg; Pkg.test()'`.

---

### Task 1: Local Factor Records And Denominator-Cleared Replay

**Files:**
- Create: `test/expert/sl3_local_local_factors.jl`
- Modify: `src/algorithm/sl3_local.jl`

**Interfaces:**
- Consumes: `_sl3_local_murthy_verify_local_unit_witness(R, X, witness, denominator; label = "...")`, `_same_base_ring`, `_coerce_into_ring`, `_sl3_local_factor_product`.
- Produces:
  - `SL3LocalElementaryFactor`
  - `SL3LocalElementaryFactorReplay`
  - `sl3_local_elementary_factor(row, col, numerator, denominator, X; local_unit_witness = nothing, n::Int = 3)`
  - `sl3_local_materialize_elementary_factor(record)`
  - `sl3_local_denominator_one_records_from_matrices(factors, X)`
  - `sl3_local_elementary_factor_replay(target, records, X)`
  - `verify_sl3_local_elementary_factor(record)::Bool`
  - `verify_sl3_local_elementary_factor_replay(replay)::Bool`

- [ ] **Step 1: Write the failing focused test**

Add `test/expert/sl3_local_local_factors.jl` with this content:

```julia
using Test
using Suslin
using Oscar

function _local_factor_product(factors, R)
    product = identity_matrix(R, 3)
    for factor in factors
        product *= factor
    end
    return product
end

function _local_factor_error(f)
    try
        f()
        return nothing
    catch err
        return err
    end
end

function _local_factor_witness(unit, selected_variable, generator)
    R = parent(unit)
    return (;
        context = (;
            kind = :localization_at_maximal_ideal,
            selected_variable,
            maximal_ideal_generators = (generator,),
        ),
        unit,
        residue_unit = one(R),
        residue_inverse = one(R),
        maximal_ideal_generators = (generator,),
        residue_difference_coefficients = (one(R),),
        global_unit = is_unit(unit),
    )
end

@testset "SL3 local elementary factor replay" begin
    R, (X,) = Oscar.polynomial_ring(QQ, ["X"])
    coefficient = X + one(R)
    ordinary_factor = Suslin.elementary_matrix(3, 1, 2, coefficient, R)
    local_factor = Suslin.sl3_local_elementary_factor(1, 2, coefficient, one(R), X)

    @test Suslin.verify_sl3_local_elementary_factor(local_factor)
    @test Suslin.sl3_local_materialize_elementary_factor(local_factor) == ordinary_factor

    ordinary_replay = Suslin.sl3_local_elementary_factor_replay(
        ordinary_factor,
        [local_factor],
        X,
    )
    @test ordinary_replay.mode == :ordinary
    @test ordinary_replay.denominator_product == one(R)
    @test ordinary_replay.materialized_factors == [ordinary_factor]
    @test ordinary_replay.cleared_product == ordinary_factor
    @test Suslin.verify_sl3_local_elementary_factor_replay(ordinary_replay)

    adapted = Suslin.sl3_local_denominator_one_records_from_matrices([ordinary_factor], X)
    @test length(adapted) == 1
    @test Suslin.sl3_local_materialize_elementary_factor(first(adapted)) == ordinary_factor
    adapted_replay = Suslin.sl3_local_elementary_factor_replay(ordinary_factor, adapted, X)
    @test adapted_replay.mode == :ordinary
    @test Suslin.verify_sl3_local_elementary_factor_replay(adapted_replay)

    RU, (u, UX) = Oscar.polynomial_ring(QQ, ["u", "X"])
    denominator = one(RU) + u
    local_unit_witness = _local_factor_witness(denominator, UX, u)
    first_local = Suslin.sl3_local_elementary_factor(
        1,
        2,
        one(RU),
        denominator,
        UX;
        local_unit_witness,
    )
    second_local = Suslin.sl3_local_elementary_factor(
        1,
        2,
        -one(RU),
        denominator,
        UX;
        local_unit_witness,
    )
    identity_target = identity_matrix(RU, 3)
    cleared_replay = Suslin.sl3_local_elementary_factor_replay(
        identity_target,
        [first_local, second_local],
        UX,
    )
    @test cleared_replay.mode == :denominator_cleared
    @test cleared_replay.materialized_factors === nothing
    @test cleared_replay.denominator_product == denominator^2
    @test cleared_replay.cleared_product == denominator^2 * identity_target
    @test Suslin.verify_sl3_local_elementary_factor_replay(cleared_replay)

    nonmaterialized_err =
        _local_factor_error(() -> Suslin.sl3_local_materialize_elementary_factor(first_local))
    @test nonmaterialized_err isa ArgumentError
    @test occursin("cannot materialize", sprint(showerror, nonmaterialized_err))

    corrupted_numerator = Suslin.SL3LocalElementaryFactor(
        first_local.R,
        first_local.n,
        first_local.row,
        first_local.col,
        first_local.numerator + one(RU),
        first_local.denominator,
        first_local.selected_variable,
        first_local.local_unit_witness,
    )
    @test Suslin.verify_sl3_local_elementary_factor(corrupted_numerator)
    bad_numerator_replay = Suslin.SL3LocalElementaryFactorReplay(
        cleared_replay.target,
        [corrupted_numerator, second_local],
        cleared_replay.selected_variable,
        cleared_replay.mode,
        cleared_replay.denominator_product,
        cleared_replay.cleared_product,
        cleared_replay.materialized_factors,
    )
    @test !Suslin.verify_sl3_local_elementary_factor_replay(bad_numerator_replay)

    bad_row = Suslin.SL3LocalElementaryFactor(
        first_local.R,
        first_local.n,
        4,
        first_local.col,
        first_local.numerator,
        first_local.denominator,
        first_local.selected_variable,
        first_local.local_unit_witness,
    )
    @test !Suslin.verify_sl3_local_elementary_factor(bad_row)
    @test_throws ArgumentError Suslin.sl3_local_elementary_factor(
        1,
        1,
        one(RU),
        denominator,
        UX;
        local_unit_witness,
    )

    corrupted_witness = merge(
        local_unit_witness,
        (; residue_difference_coefficients = (zero(RU),)),
    )
    @test_throws ArgumentError Suslin.sl3_local_elementary_factor(
        1,
        2,
        one(RU),
        denominator,
        UX;
        local_unit_witness = corrupted_witness,
    )
    bad_witness_record = Suslin.SL3LocalElementaryFactor(
        first_local.R,
        first_local.n,
        first_local.row,
        first_local.col,
        first_local.numerator,
        first_local.denominator,
        first_local.selected_variable,
        corrupted_witness,
    )
    @test !Suslin.verify_sl3_local_elementary_factor(bad_witness_record)

    wrong_denominator_err = _local_factor_error(() -> Suslin.sl3_local_elementary_factor(
        1,
        2,
        one(RU),
        denominator + u,
        UX;
        local_unit_witness,
    ))
    @test wrong_denominator_err isa ArgumentError
    @test occursin("unit does not match expected value", sprint(showerror, wrong_denominator_err))
end
```

- [ ] **Step 2: Run the focused test to verify it fails for missing implementation**

Run:

```bash
julia --project=. -e 'include("test/expert/sl3_local_local_factors.jl")'
```

Expected: FAIL with `UndefVarError` for `sl3_local_elementary_factor` or another missing local factor helper.

- [ ] **Step 3: Implement the local factor records and verifier**

In `src/algorithm/sl3_local.jl`, add the two structs near the existing SL3 local structs:

```julia
struct SL3LocalElementaryFactor
    R
    n::Int
    row::Int
    col::Int
    numerator
    denominator
    selected_variable
    local_unit_witness
end

struct SL3LocalElementaryFactorReplay
    target
    factors::Vector{SL3LocalElementaryFactor}
    selected_variable
    mode::Symbol
    denominator_product
    cleared_product
    materialized_factors
end
```

Add helper functions after `_sl3_local_factor_product` so they can reuse local product helpers:

```julia
function sl3_local_elementary_factor(row, col, numerator, denominator, X;
        local_unit_witness=nothing, n::Int=3)
    R = parent(X)
    record = SL3LocalElementaryFactor(
        R,
        n,
        Int(row),
        Int(col),
        _coerce_into_ring(R, numerator, "local elementary factor numerator"),
        _coerce_into_ring(R, denominator, "local elementary factor denominator"),
        X,
        local_unit_witness,
    )
    _require_sl3_local_elementary_factor(record)
    return record
end

function verify_sl3_local_elementary_factor(record)::Bool
    try
        _require_sl3_local_elementary_factor(record)
        return true
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end

function sl3_local_materialize_elementary_factor(record)
    _require_sl3_local_elementary_factor(record)
    record.denominator == one(record.R) ||
        throw(ArgumentError("local elementary factor cannot materialize over the ordinary base ring"))
    return elementary_matrix(record.n, record.row, record.col, record.numerator, record.R)
end

function sl3_local_denominator_one_records_from_matrices(factors, X)
    return [
        _sl3_local_denominator_one_record_from_matrix(factor, X)
        for factor in factors
    ]
end

function sl3_local_elementary_factor_replay(target, records, X)
    collected = SL3LocalElementaryFactor[records...]
    replay = _sl3_local_elementary_factor_replay(target, collected, X)
    verify_sl3_local_elementary_factor_replay(replay) ||
        error("internal local elementary factor replay verification failed")
    return replay
end

function verify_sl3_local_elementary_factor_replay(replay)::Bool
    try
        return _sl3_local_elementary_factor_replay_verification(replay).overall_ok
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end
```

Implement the private helpers with these exact rules:

```julia
function _require_sl3_local_elementary_factor(record)
    record.n == 3 || throw(ArgumentError("local elementary factor size must be 3"))
    1 <= record.row <= record.n || throw(ArgumentError("local elementary factor row is out of bounds"))
    1 <= record.col <= record.n || throw(ArgumentError("local elementary factor column is out of bounds"))
    record.row != record.col || throw(ArgumentError("local elementary factor row and column must differ"))
    parent(record.selected_variable) === record.R ||
        throw(ArgumentError("local elementary factor variable must lie in the factor ring"))
    parent(record.numerator) === record.R ||
        throw(ArgumentError("local elementary factor numerator ring mismatch"))
    parent(record.denominator) === record.R ||
        throw(ArgumentError("local elementary factor denominator ring mismatch"))
    iszero(record.denominator) &&
        throw(ArgumentError("local elementary factor denominator must be nonzero"))
    if record.denominator == one(record.R)
        record.local_unit_witness === nothing ||
            _sl3_local_murthy_verify_local_unit_witness(
                record.R,
                record.selected_variable,
                record.local_unit_witness,
                record.denominator;
                label = "local elementary factor denominator witness",
            )
    else
        record.local_unit_witness === nothing &&
            throw(ArgumentError("local elementary factor denominator requires a local-unit witness"))
        _sl3_local_murthy_verify_local_unit_witness(
            record.R,
            record.selected_variable,
            record.local_unit_witness,
            record.denominator;
            label = "local elementary factor denominator witness",
        )
    end
    return record
end

function _sl3_local_cleared_elementary_factor(record)
    _require_sl3_local_elementary_factor(record)
    cleared = record.denominator * identity_matrix(record.R, record.n)
    cleared[record.row, record.col] += record.numerator
    return cleared
end

function _sl3_local_denominator_one_record_from_matrix(factor, X)
    nrows(factor) == 3 || throw(ArgumentError("ordinary elementary factor must be 3x3"))
    ncols(factor) == 3 || throw(ArgumentError("ordinary elementary factor must be 3x3"))
    R = base_ring(factor)
    parent(X) === R || throw(ArgumentError("ordinary elementary factor variable ring mismatch"))
    identity = identity_matrix(R, 3)
    row = 0
    col = 0
    coefficient = zero(R)
    for i in 1:3, j in 1:3
        if i == j
            factor[i, j] == one(R) || throw(ArgumentError("ordinary factor diagonal is not identity"))
        elseif factor[i, j] != zero(R)
            row == 0 || throw(ArgumentError("ordinary factor has more than one off-diagonal entry"))
            row = i
            col = j
            coefficient = factor[i, j]
        end
    end
    row != 0 || throw(ArgumentError("ordinary identity factor has no elementary row/column data"))
    return sl3_local_elementary_factor(row, col, coefficient, one(R), X)
end

function _sl3_local_elementary_factor_replay(target, records::Vector{SL3LocalElementaryFactor}, X)
    nrows(target) == 3 || throw(ArgumentError("local elementary factor replay target must be 3x3"))
    ncols(target) == 3 || throw(ArgumentError("local elementary factor replay target must be 3x3"))
    R = base_ring(target)
    parent(X) === R || throw(ArgumentError("local elementary factor replay variable ring mismatch"))
    denominator_product = one(R)
    cleared_product = identity_matrix(R, 3)
    all_materializable = true
    materialized = Any[]
    for record in records
        _require_sl3_local_elementary_factor(record)
        record.R === R || throw(ArgumentError("local elementary factor replay ring mismatch"))
        record.selected_variable == X ||
            throw(ArgumentError("local elementary factor replay variable mismatch"))
        denominator_product *= record.denominator
        cleared_product *= _sl3_local_cleared_elementary_factor(record)
        if record.denominator == one(R)
            push!(materialized, sl3_local_materialize_elementary_factor(record))
        else
            all_materializable = false
        end
    end
    mode = all_materializable ? :ordinary : :denominator_cleared
    return SL3LocalElementaryFactorReplay(
        target,
        records,
        X,
        mode,
        denominator_product,
        cleared_product,
        all_materializable ? collect(materialized) : nothing,
    )
end
```

Implement `_sl3_local_elementary_factor_replay_verification(replay)` so it recomputes `denominator_product`, `cleared_product`, mode, and materialized factors from `replay.factors`, then requires:

```julia
denominator_cleared_ok = expected_cleared_product == replay.denominator_product * replay.target
ordinary_ok =
    replay.mode == :ordinary &&
    replay.materialized_factors !== nothing &&
    _sl3_local_factor_product(replay.materialized_factors, R) == replay.target
```

The final `overall_ok` is true when shape/ring/variable checks pass, factor checks pass, stored fields match recomputed values, and either the ordinary replay or denominator-cleared replay equality holds.

- [ ] **Step 4: Run the focused local factor test**

Run:

```bash
julia --project=. -e 'include("test/expert/sl3_local_local_factors.jl")'
```

Expected: PASS with one `SL3 local elementary factor replay` testset.

- [ ] **Step 5: Commit Task 1**

Run:

```bash
git add src/algorithm/sl3_local.jl test/expert/sl3_local_local_factors.jl
git commit -m "feat: add SL3 local elementary factor replay"
```

Expected: commit succeeds.

---

### Task 2: Certificate Adapter Coverage And Test Registration

**Files:**
- Modify: `test/expert/sl3_local_certificate.jl`
- Modify: `test/runtests.jl`

**Interfaces:**
- Consumes: Task 1 `sl3_local_denominator_one_records_from_matrices`, `sl3_local_elementary_factor_replay`, and `verify_sl3_local_elementary_factor_replay`.
- Produces: Expert-suite registration for `expert/sl3_local_local_factors.jl` and certificate coverage proving current ordinary factors replay through denominator-one records.

- [ ] **Step 1: Write the failing certificate adapter assertion and register the new test**

In `test/expert/sl3_local_certificate.jl`, inside `_assert_sl3_certificate_replays(cert)` after the existing three assertions, add:

```julia
    local_records = Suslin.sl3_local_denominator_one_records_from_matrices(
        cert.factors,
        cert.selected_variable,
    )
    local_replay = Suslin.sl3_local_elementary_factor_replay(
        cert.target,
        local_records,
        cert.selected_variable,
    )
    @test local_replay.mode == :ordinary
    @test local_replay.materialized_factors == cert.factors
    @test Suslin.verify_sl3_local_elementary_factor_replay(local_replay)
```

In `test/runtests.jl`, add the new expert file immediately after `expert/sl3_local_certificate.jl`:

```julia
        "expert/sl3_local_local_factors.jl",
```

- [ ] **Step 2: Run the focused certificate test**

Run:

```bash
julia --project=. -e 'include("test/expert/sl3_local_certificate.jl")'
```

Expected: PASS. If it fails because an existing ordinary factor is an identity matrix and cannot be adapted, update Task 1's adapter to support an explicit future row/column adapter rather than guessing identity row/column data; keep the certificate assertion scoped to factors that have unambiguous matrix form.

- [ ] **Step 3: Run the expert test group enough to exercise registration**

Run:

```bash
julia --project=. test/runtests.jl expert
```

Expected: PASS for the expert group.

- [ ] **Step 4: Commit Task 2**

Run:

```bash
git add test/expert/sl3_local_certificate.jl test/runtests.jl
git commit -m "test: cover SL3 local factor replay adapters"
```

Expected: commit succeeds.

---

### Task 3: Final Verification And Branch Review

**Files:**
- No planned code changes unless verification or review finds a defect.

**Interfaces:**
- Consumes: Task 1 and Task 2 implementation.
- Produces: Verified branch ready for pull request creation.

- [ ] **Step 1: Run issue-required focused verification**

Run:

```bash
julia --project=. -e 'include("test/expert/sl3_local_local_factors.jl")'
julia --project=. -e 'include("test/expert/sl3_local_certificate.jl")'
```

Expected: both commands exit 0.

- [ ] **Step 2: Run Agent Desk package verification**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: exit 0 with public and internal test groups passing.

- [ ] **Step 3: Inspect git status and branch diff**

Run:

```bash
git status --short
git log --oneline origin/main..HEAD
git diff --stat origin/main..HEAD
```

Expected: no uncommitted tracked changes, and branch commits limited to the issue #207 docs, implementation, and tests.

- [ ] **Step 4: Run final code review**

Use the superpowers:requesting-code-review reviewer on the full branch diff from `git merge-base origin/main HEAD` to `HEAD`. Fix Critical and Important findings, then repeat the focused verification covering the fix.

## Plan Self-Review

- Every issue #207 requirement maps to Task 1 or Task 2.
- The plan keeps public ordinary APIs unchanged.
- The local-unit evidence path explicitly reuses the #208 verifier.
- Negative controls cover numerator, denominator, row/column, and witness corruption.
- Verification commands include the two issue-required commands and the Agent Desk package command.
