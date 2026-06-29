# Issue 209 Local q-Degree and q(0)-Unit Replay Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Route q-degree normalization and the q(0)-unit Murthy branch through checked local contexts and denominator-aware local replay without changing ordinary `QQ[X]` behavior.

**Architecture:** Add internal context methods in `src/algorithm/sl3_local.jl`. Q-degree contexts reuse the existing reduction certificate. Nontrivial local q(0)-unit contexts translate an exact fraction-field Murthy certificate into #207 `SL3LocalElementaryFactor` records and verify the final certificate by denominator-cleared replay.

**Tech Stack:** Julia, Oscar/AbstractAlgebra polynomial rings, existing `SL3LocalMurthyInputContext`, `SL3LocalRealizationCertificate`, and `SL3LocalElementaryFactorReplay`.

## Global Constraints

- Existing ordinary `QQ[X]` behavior must stay unchanged.
- Nontrivial local-unit q(0)-unit acceptance must use localized replay/certificate verification, not ordinary `verify_factorization(A, factors)` over the original base ring.
- Do not implement the q(0)-nonunit local Bezout/resultant branch.
- Do not integrate with Quillen patching or public Park-Woodburn `SL_3` routing.
- Do not export new public APIs from `src/Suslin.jl`.
- Focused verification commands are `julia --project=. -e 'include("test/expert/sl3_local_q_degree_normalization.jl")'` and `julia --project=. -e 'include("test/expert/sl3_local_murthy_q_unit.jl")'`.
- Package verification command is `julia --project=. -e 'using Pkg; Pkg.test()'`.

---

## File Structure

- Modify `src/algorithm/sl3_local.jl`: add context method overloads, fraction-field translation helpers, `SL3LocalMurthyQUnitLocalReduction`, local q(0)-unit certificate construction, and verifier integration.
- Modify `test/expert/sl3_local_q_degree_normalization.jl`: add #206 local q-degree context replay coverage and context guard checks.
- Modify `test/expert/sl3_local_murthy_q_unit.jl`: add #206 local q(0)-unit certificate coverage and negative controls for local replay corruption.
- Modify `test/expert/sl3_local_murthy_context.jl`: add small route smoke checks for context-consuming certificate dispatch.

### Task 1: Context q-Degree Normalization

**Files:**
- Modify: `src/algorithm/sl3_local.jl`
- Modify: `test/expert/sl3_local_q_degree_normalization.jl`

**Interfaces:**
- Consumes: `SL3LocalMurthyInputContext`, `verify_sl3_local_murthy_input_context`, existing `sl3_local_q_degree_normalization(p, q, r, s, X)`.
- Produces: `sl3_local_q_degree_normalization(context::SL3LocalMurthyInputContext)` and `sl3_local_q_degree_normalization_certificate(context::SL3LocalMurthyInputContext)`.

- [x] **Step 1: Write the failing local q-degree context test**

Add this block to `test/expert/sl3_local_q_degree_normalization.jl` inside the existing testset after the ordinary fixture certificate checks:

```julia
    local_fixture = by_id["mg-local-q-degree-qq-u-x"]
    local_context = Suslin.sl3_local_murthy_input_context(
        local_fixture.target,
        local_fixture.variable;
        witness = first(local_fixture.witnesses),
    )
    local_record = Suslin.sl3_local_q_degree_normalization(local_context)
    local_witness = first(local_fixture.witnesses)
    _assert_qdegree_record(
        local_record,
        local_fixture.entries.p,
        local_fixture.entries.q,
        local_fixture.entries.r,
        local_fixture.entries.s,
        local_fixture.variable,
        local_witness.quotient,
        local_witness.remainder,
    )
    local_cert = Suslin.sl3_local_q_degree_normalization_certificate(local_context)
    @test local_cert.branch == :murthy_q_degree_normalization
    @test local_cert.target == local_fixture.target
    @test Suslin.verify_sl3_local_realization(local_cert)
```

Add this guard check near the q-degree negative controls:

```julia
    local_q0_unit = by_id["mg-local-q0-unit-at-u"]
    local_q0_context = Suslin.sl3_local_murthy_input_context(
        local_q0_unit.target,
        local_q0_unit.variable;
        witness = first(local_q0_unit.witnesses),
    )
    @test_throws ArgumentError Suslin.sl3_local_q_degree_normalization(local_q0_context)
```

- [x] **Step 2: Run the focused q-degree test and verify RED**

Run:

```bash
julia --project=. -e 'include("test/expert/sl3_local_q_degree_normalization.jl")'
```

Expected before implementation: failure with no method matching `sl3_local_q_degree_normalization(::SL3LocalMurthyInputContext)`.

- [x] **Step 3: Add the q-degree context implementation**

Add these methods near the existing q-degree helper methods in `src/algorithm/sl3_local.jl`:

```julia
function sl3_local_q_degree_normalization(context::SL3LocalMurthyInputContext)
    verify_sl3_local_murthy_input_context(context) ||
        throw(ArgumentError("Murthy q-degree normalization requires a verified local input context"))
    context.degree_q >= context.degree_p ||
        throw(ArgumentError("Murthy q-degree normalization context requires deg(q) >= deg(p)"))
    return sl3_local_q_degree_normalization(
        context.entries.p,
        context.entries.q,
        context.entries.r,
        context.entries.s,
        context.X,
    )
end

function sl3_local_q_degree_normalization_certificate(context::SL3LocalMurthyInputContext)
    return sl3_local_q_degree_normalization_certificate(
        sl3_local_q_degree_normalization(context),
    )
end
```

- [x] **Step 4: Run the focused q-degree test and verify GREEN**

Run:

```bash
julia --project=. -e 'include("test/expert/sl3_local_q_degree_normalization.jl")'
```

Expected after implementation: all tests pass.

- [x] **Step 5: Commit Task 1**

Run:

```bash
git add src/algorithm/sl3_local.jl test/expert/sl3_local_q_degree_normalization.jl
git commit -m "feat: route q-degree normalization through Murthy context"
```

### Task 2: Local q(0)-Unit Certificate Replay

**Files:**
- Modify: `src/algorithm/sl3_local.jl`
- Modify: `test/expert/sl3_local_murthy_q_unit.jl`

**Interfaces:**
- Consumes: Task 1 context dispatch, #208 local-unit witnesses, #207 local elementary factor replay, existing ordinary q(0)-unit certificate.
- Produces: `SL3LocalMurthyQUnitLocalReduction`, `realize_sl3_local_certificate(context::SL3LocalMurthyInputContext)`, and local q(0)-unit `:murthy_q0_unit` certificates whose factors are `SL3LocalElementaryFactor` records.

- [x] **Step 1: Write the failing local q(0)-unit replay test**

Add helper assertions to `test/expert/sl3_local_murthy_q_unit.jl`:

```julia
function _assert_local_q0_certificate(cert)
    @test cert.branch == :murthy_q0_unit
    @test Suslin.verify_sl3_local_realization(cert)
    @test !isempty(cert.factors)
    @test all(factor -> factor isa Suslin.SL3LocalElementaryFactor, cert.factors)
    reduction = cert.witness.reduction
    @test reduction isa Suslin.SL3LocalMurthyQUnitLocalReduction
    @test reduction.local_factor_replay.mode == :denominator_cleared
    @test reduction.local_factor_replay.target == cert.target
    @test reduction.local_factor_replay.factors == cert.factors
    @test Suslin.verify_sl3_local_elementary_factor_replay(reduction.local_factor_replay)
    @test Suslin.verify_sl3_local_realization(reduction.source_certificate)
    @test reduction.split_certificate.branch == :murthy_split_lemma
    @test Suslin.verify_sl3_local_realization(reduction.split_certificate)
    @test Suslin.verify_sl3_local_realization(reduction.split_certificate.witness.first_child_certificate)
    @test Suslin.verify_sl3_local_realization(reduction.split_certificate.witness.second_child_certificate)
end
```

Add this local-contract case after the ordinary q0-unit fixture assertions:

```julia
    local_fixture = by_id["mg-local-q0-unit-at-u"]
    local_context = Suslin.sl3_local_murthy_input_context(
        local_fixture.target,
        local_fixture.variable;
        witness = first(local_fixture.witnesses),
    )
    local_cert = Suslin.realize_sl3_local_certificate(local_context)
    _assert_local_q0_certificate(local_cert)
    @test local_cert.target == local_fixture.target
```

- [x] **Step 2: Run the focused q0-unit test and verify RED**

Run:

```bash
julia --project=. -e 'include("test/expert/sl3_local_murthy_q_unit.jl")'
```

Expected before implementation: failure because `realize_sl3_local_certificate(::SL3LocalMurthyInputContext)` and `SL3LocalMurthyQUnitLocalReduction` do not exist.

- [x] **Step 3: Add the local reduction type and context route**

Add `SL3LocalMurthyQUnitLocalReduction` near `SL3LocalMurthyQUnitReduction`:

```julia
struct SL3LocalMurthyQUnitLocalReduction
    target
    context::SL3LocalMurthyInputContext
    q0
    q0_inverse
    p0
    right_e21_coefficient
    elimination_factor::SL3LocalElementaryFactor
    inverse_elimination_factor::SL3LocalElementaryFactor
    source_certificate
    split_certificate
    local_factor_replay::SL3LocalElementaryFactorReplay
    selected_variable
    degree_p::Int
    degree_q::Int
end
```

Add `realize_sl3_local_certificate(context::SL3LocalMurthyInputContext)` after the existing certificate entry points:

```julia
function realize_sl3_local_certificate(context::SL3LocalMurthyInputContext)
    verify_sl3_local_murthy_input_context(context) ||
        throw(ArgumentError("local Murthy realization requires a verified input context"))
    if context.degree_q >= context.degree_p
        return sl3_local_q_degree_normalization_certificate(context)
    end
    context.global_units.q0 || context.local_units.q0 ||
        throw(ArgumentError("local Murthy q(0)-unit realization requires q(0) to be a global or local unit"))
    if context.global_units.q0 && _sl3_local_supports_murthy_q0_unit_branch(context.R, context.var_idx)
        return realize_sl3_local_certificate(context.target, context.X)
    end
    return _realize_sl3_local_murthy_q0_unit_local_certificate(context)
end
```

- [x] **Step 4: Add fraction-field translation helpers and local certificate construction**

Implement helpers in `src/algorithm/sl3_local.jl`:

```julia
_sl3_local_fraction_model(context)
_sl3_local_to_fraction_polynomial(value, model)
_sl3_local_fraction_polynomial_to_ratio(value, context)
_sl3_local_coefficient_model_to_ring(value, context, coefficient_names)
_sl3_local_derive_local_unit_witness(context, unit)
_sl3_local_local_factor_from_fraction_matrix(factor, context)
_sl3_local_local_factors_from_fraction_matrices(factors, context)
_sl3_local_source_q0_reduction(source_certificate)
_realize_sl3_local_murthy_q0_unit_local_certificate(context)
_sl3_local_murthy_q0_unit_local_reduction(context)
```

The construction must:

```julia
source_certificate = realize_sl3_local_certificate(fraction_target, fraction_variable)
records = _sl3_local_local_factors_from_fraction_matrices(source_certificate.factors, context)
replay = sl3_local_elementary_factor_replay(context.target, records, context.X)
source_reduction = _sl3_local_source_q0_reduction(source_certificate)
elimination_factor = first(_sl3_local_local_factors_from_fraction_matrices(
    [source_reduction.elimination_factor],
    context,
))
inverse_elimination_factor = first(_sl3_local_local_factors_from_fraction_matrices(
    [source_reduction.inverse_elimination_factor],
    context,
))
```

For denominator witness derivation, support exactly the checked one-generator
localization used by #206 first. If the context has no q0 local-unit witness,
or has multiple maximal-ideal generators, throw `ArgumentError` with a message
containing `unsupported local-unit denominator witness`.

- [x] **Step 5: Teach certificate verification about local q0 reductions**

Add a specific method:

```julia
function verify_sl3_local_murthy_q_unit_reduction(
        reduction::SL3LocalMurthyQUnitLocalReduction,
)::Bool
    try
        return _sl3_local_murthy_q0_unit_local_reduction_verification(reduction).overall_ok
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end
```

Update `_sl3_local_certificate_expected_factors` so a `:murthy_q0_unit`
certificate with `reduction isa SL3LocalMurthyQUnitLocalReduction` returns
`reduction.local_factor_replay.factors`.

Update `_sl3_local_realization_verification` so the final factor check uses
`verify_sl3_local_elementary_factor_replay(reduction.local_factor_replay)` for
local q0 reductions and keeps `verify_factorization` for every existing branch.

- [x] **Step 6: Run the focused q0-unit test and verify GREEN**

Run:

```bash
julia --project=. -e 'include("test/expert/sl3_local_murthy_q_unit.jl")'
```

Expected after implementation: all tests pass.

- [x] **Step 7: Commit Task 2**

Run:

```bash
git add src/algorithm/sl3_local.jl test/expert/sl3_local_murthy_q_unit.jl
git commit -m "feat: replay local Murthy q0-unit certificates"
```

### Task 3: Negative Controls and Context Route Smoke Tests

**Files:**
- Modify: `src/algorithm/sl3_local.jl`
- Modify: `test/expert/sl3_local_murthy_q_unit.jl`
- Modify: `test/expert/sl3_local_murthy_context.jl`

**Interfaces:**
- Consumes: Task 2 local reduction verifier.
- Produces: focused corruption coverage for `q0_inverse`, right `E21` coefficient, source split child, and local factor replay.

- [x] **Step 1: Add q0-unit local negative controls**

In `test/expert/sl3_local_murthy_q_unit.jl`, after the local positive case,
construct corrupted reductions and certificates:

```julia
    local_reduction = local_cert.witness.reduction
    bad_q0_inverse = Suslin.SL3LocalMurthyQUnitLocalReduction(
        local_reduction.target,
        local_reduction.context,
        local_reduction.q0,
        merge(local_reduction.q0_inverse, (; denominator = local_reduction.q0 + one(local_reduction.context.R))),
        local_reduction.p0,
        local_reduction.right_e21_coefficient,
        local_reduction.elimination_factor,
        local_reduction.inverse_elimination_factor,
        local_reduction.source_certificate,
        local_reduction.split_certificate,
        local_reduction.local_factor_replay,
        local_reduction.selected_variable,
        local_reduction.degree_p,
        local_reduction.degree_q,
    )
    @test !Suslin.verify_sl3_local_murthy_q_unit_reduction(bad_q0_inverse)

    bad_e21 = Suslin.SL3LocalMurthyQUnitLocalReduction(
        local_reduction.target,
        local_reduction.context,
        local_reduction.q0,
        local_reduction.q0_inverse,
        local_reduction.p0,
        merge(local_reduction.right_e21_coefficient, (; numerator = local_reduction.right_e21_coefficient.numerator + one(local_reduction.context.R))),
        local_reduction.elimination_factor,
        local_reduction.inverse_elimination_factor,
        local_reduction.source_certificate,
        local_reduction.split_certificate,
        local_reduction.local_factor_replay,
        local_reduction.selected_variable,
        local_reduction.degree_p,
        local_reduction.degree_q,
    )
    @test !Suslin.verify_sl3_local_murthy_q_unit_reduction(bad_e21)
```

Also corrupt one local factor in `local_factor_replay.factors` and assert the
reduction verifier and certificate verifier both return `false`.

- [x] **Step 2: Add source split child corruption coverage**

Build a copy of `local_reduction.split_certificate.witness.first_child_certificate`
with a mismatched target and rewrap it in a copied split certificate. Rebuild
the local reduction with that split certificate and assert:

```julia
@test !Suslin.verify_sl3_local_murthy_q_unit_reduction(bad_split_reduction)
```

- [x] **Step 3: Add context route smoke tests**

In `test/expert/sl3_local_murthy_context.jl`, after the local q0 context
positive assertions, add:

```julia
    local_context_cert = Suslin.realize_sl3_local_certificate(local_context)
    @test local_context_cert.branch == :murthy_q0_unit
    @test Suslin.verify_sl3_local_realization(local_context_cert)
```

For the local q0-nonunit context already built in that file, add:

```julia
    @test_throws ArgumentError Suslin.realize_sl3_local_certificate(bezout_context)
```

- [x] **Step 4: Run focused tests**

Run:

```bash
julia --project=. -e 'include("test/expert/sl3_local_q_degree_normalization.jl")'
julia --project=. -e 'include("test/expert/sl3_local_murthy_q_unit.jl")'
julia --project=. -e 'include("test/expert/sl3_local_murthy_context.jl")'
```

Expected: all commands pass.

- [x] **Step 5: Run package tests**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: default package tests pass.

- [x] **Step 6: Commit Task 3**

Run:

```bash
git add src/algorithm/sl3_local.jl test/expert/sl3_local_murthy_q_unit.jl test/expert/sl3_local_murthy_context.jl
git commit -m "test: harden local Murthy q0-unit replay"
```

## Self-Review

- Spec coverage: Task 1 covers q-degree context consumption; Task 2 covers local q0-unit replay with #207 records; Task 3 covers negative controls and context route smoke tests.
- Placeholder scan: no task contains placeholder markers or deferred implementation language.
- Type consistency: `SL3LocalMurthyQUnitLocalReduction`, `local_factor_replay`, and `realize_sl3_local_certificate(context::SL3LocalMurthyInputContext)` are introduced before later tasks consume them.
