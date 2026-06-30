# Issue 210 Local q(0)-Nonunit Resultant Replay Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Route supplied q(0)-nonunit Bezout/resultant Murthy local contexts through replayable local certificates while preserving ordinary `QQ[X]` extraction.

**Architecture:** Extend the existing q(0)-nonunit reduction instead of adding a second branch type. Ordinary univariate inputs keep materialized matrix factors; nontrivial local contexts use supplied checked witness data, construct a local q(0)-unit child context, and verify the final sequence through `SL3LocalElementaryFactorReplay`.

**Tech Stack:** Julia, Oscar/AbstractAlgebra polynomial rings, existing `SL3LocalMurthyInputContext`, `SL3LocalRealizationCertificate`, `SL3LocalMurthyQUnitLocalReduction`, and `SL3LocalElementaryFactorReplay`.

## Global Constraints

- Existing ordinary `QQ[X]` supplied-witness and `gcdx` extraction behavior must stay unchanged.
- Broader local contexts must prefer supplied witness data and stage-fail unsupported automatic extraction.
- Nontrivial local-witness acceptance must use localized replay/certificate verification, not ordinary `verify_factorization(A, factors)` over the original base ring.
- The branch verifier must recompute Bezout equality, degree bounds, child q(0)-unit condition, child replay, elementary identities, and final replay metadata.
- Do not implement Quillen local-to-global patching or the general `SL_3` driver.
- Do not export new public APIs from `src/Suslin.jl`.
- Focused verification command is `julia --project=. -e 'include("test/expert/sl3_local_murthy_resultant.jl")'`.
- Package verification command is `julia --project=. -e 'using Pkg; Pkg.test()'`.

---

## File Structure

- Modify `src/algorithm/sl3_local.jl`: add local q(0)-nonunit context routing, local replay metadata on the reduction record, helper functions that adapt prefix matrices plus child local factors into local elementary replay, and verifier updates.
- Modify `test/expert/sl3_local_murthy_resultant.jl`: add the #206 supplied-witness local-contract case and negative controls.
- Modify `test/expert/sl3_local_murthy_context.jl`: update the old staged-failure assertion for the #206 q(0)-nonunit context to expect a verifying certificate.

### Task 1: RED Coverage for Local q(0)-Nonunit Replay

**Files:**
- Modify: `test/expert/sl3_local_murthy_resultant.jl`
- Modify: `test/expert/sl3_local_murthy_context.jl`

**Interfaces:**
- Consumes: existing fixture id `mg-local-q0-nonunit-bezout-at-u`, `sl3_local_murthy_input_context`, `realize_sl3_local_certificate`.
- Produces: failing tests that describe local q(0)-nonunit replay acceptance and verifier rejection cases.

- [ ] **Step 1: Extend the resultant assertion helper**

In `test/expert/sl3_local_murthy_resultant.jl`, update `_assert_resultant_certificate` so it accepts local replay certificates:

```julia
    if getproperty(reduction, :local_factor_replay, nothing) === nothing
        @test Suslin.verify_factorization(cert.target, cert.factors)
        @test _resultant_product(cert.factors, base_ring(cert.target)) == cert.target
        @test reduction.branch_unit * reduction.branch_unit_inverse == one(base_ring(cert.target))
    else
        @test all(factor -> factor isa Suslin.SL3LocalElementaryFactor, cert.factors)
        @test reduction.local_factor_replay.target == cert.target
        @test reduction.local_factor_replay.factors == cert.factors
        @test Suslin.verify_sl3_local_elementary_factor_replay(reduction.local_factor_replay)
        child_reduction = reduction.child_certificate.witness.reduction
        @test child_reduction isa Suslin.SL3LocalMurthyQUnitLocalReduction
        @test child_reduction.context.q0 == reduction.branch_unit
        @test child_reduction.context.local_units.q0
    end
```

- [ ] **Step 2: Add the local positive case**

Add this block after the extracted ordinary case:

```julia
    local_fixture = by_id["mg-local-q0-nonunit-bezout-at-u"]
    local_context = Suslin.sl3_local_murthy_input_context(
        local_fixture.target,
        local_fixture.variable;
        witness = first(local_fixture.witnesses),
    )
    local_cert = Suslin.realize_sl3_local_certificate(local_context)
    _assert_resultant_certificate(local_cert; expected_source = :supplied_bezout_witness)
    @test local_cert.target == local_fixture.target
```

- [ ] **Step 3: Add local negative controls**

Add tests that merge the local reduction named tuple and verify rejection for
corrupted `q_prime`, corrupted `resultant`, corrupted `degree_p_prime`, a child
certificate with corrupted local q0 witness evidence, and a child certificate
with a wrong target.

Use direct struct reconstruction for the corrupted child q0-unit local
reduction so construction does not normalize away the bad witness:

```julia
function _local_q0_reduction_copy(reduction; context = reduction.context)
    return Suslin.SL3LocalMurthyQUnitLocalReduction(
        reduction.target,
        context,
        reduction.q0,
        reduction.q0_inverse,
        reduction.p0,
        reduction.right_e21_coefficient,
        reduction.elimination_factor,
        reduction.inverse_elimination_factor,
        reduction.source_certificate,
        reduction.split_certificate,
        reduction.local_factor_replay,
        reduction.selected_variable,
        reduction.degree_p,
        reduction.degree_q,
    )
end
```

- [ ] **Step 4: Update the context route smoke test**

In `test/expert/sl3_local_murthy_context.jl`, replace the staged-failure
assertion for `bezout_context` with:

```julia
    bezout_context_cert = Suslin.realize_sl3_local_certificate(bezout_context)
    @test bezout_context_cert.branch == :murthy_q0_nonunit_bezout_resultant
    @test Suslin.verify_sl3_local_realization(bezout_context_cert)
```

- [ ] **Step 5: Run focused tests and verify RED**

Run:

```bash
julia --project=. -e 'include("test/expert/sl3_local_murthy_resultant.jl")'
julia --project=. -e 'include("test/expert/sl3_local_murthy_context.jl")'
```

Expected before implementation: failure because q(0)-nonunit local contexts are still rejected or because the reduction record has no local replay metadata.

### Task 2: Implement Local q(0)-Nonunit Context Routing and Replay

**Files:**
- Modify: `src/algorithm/sl3_local.jl`

**Interfaces:**
- Consumes: `SL3LocalMurthyInputContext`, `_sl3_local_murthy_bezout_data`, `realize_sl3_local_certificate(context::SL3LocalMurthyInputContext)`, `sl3_local_denominator_one_records_from_matrices`, `SL3LocalMurthyQUnitLocalReduction`.
- Produces: local q(0)-nonunit certificate construction and verifier-owned local replay.

- [ ] **Step 1: Add replay metadata to the reduction record**

Add `local_factor_replay` to `SL3LocalMurthyQ0NonunitReduction` before
`witness_source::Symbol`. Ordinary branch construction must pass `nothing`.

- [ ] **Step 2: Route q(0)-nonunit contexts**

Change `realize_sl3_local_certificate(context::SL3LocalMurthyInputContext)` so
the q0-unit guard becomes an explicit branch:

```julia
    if context.global_units.q0 || context.local_units.q0
        if context.global_units.q0 &&
                _sl3_local_supports_murthy_q0_unit_branch(context.R, context.var_idx)
            return realize_sl3_local_certificate(context.target, context.X)
        end
        return _realize_sl3_local_murthy_q0_unit_local_certificate(context)
    end
    if _sl3_local_supports_murthy_q0_unit_branch(context.R, context.var_idx)
        return realize_sl3_local_certificate(
            context.target,
            context.X;
            murthy_q0_nonunit_witness = context.bezout_witness,
        )
    end
    context.bezout_witness === nothing &&
        _throw_staged_sl3_local_failure("Murthy q(0)-nonunit local Bezout/resultant extraction is unsupported")
    return _realize_sl3_local_murthy_q0_nonunit_local_certificate(context)
```

- [ ] **Step 3: Add local q(0)-nonunit construction helpers**

Implement helpers that:

- verify the context;
- call `_sl3_local_murthy_bezout_data` and require supplied witness data for
  non-univariate local contexts;
- build `bezout_target`, `child_link_target`, `left_factor`, and
  `first_elementary_factor`;
- build the child context with `context.local_unit_witnesses.branch_unit` mapped
  to `q0` when the branch unit is local;
- call `realize_sl3_local_certificate(child_context)`;
- build `local_factor_replay` from denominator-one prefix records plus child
  local elementary records.

The main helper signatures are:

```julia
function _realize_sl3_local_murthy_q0_nonunit_local_certificate(
        context::SL3LocalMurthyInputContext,
)

function _sl3_local_murthy_q0_nonunit_local_reduction(
        context::SL3LocalMurthyInputContext,
)

function _sl3_local_murthy_q0_nonunit_child_context(
        context::SL3LocalMurthyInputContext,
        child_link_target,
)

function _sl3_local_q0_nonunit_local_factor_records(reduction)
```

- [ ] **Step 4: Update expected-factor and certificate verification paths**

When `reduction.local_factor_replay !== nothing`, return
`reduction.local_factor_replay.factors` from
`_sl3_local_certificate_expected_factors` and make
`_sl3_local_realization_verification` validate the local replay instead of
ordinary `verify_factorization`.

- [ ] **Step 5: Relax and harden q(0)-nonunit reduction verification**

Update `_sl3_local_murthy_q0_nonunit_reduction_verification` to verify ordinary
polynomial rings with any selected generator, allow `branch_unit_inverse ===
nothing` only when the child q0-unit certificate carries checked local q0
evidence, recompute `local_factor_replay` when present, and keep the ordinary
matrix-product check when replay metadata is absent.

- [ ] **Step 6: Run focused tests and verify GREEN**

Run:

```bash
julia --project=. -e 'include("test/expert/sl3_local_murthy_resultant.jl")'
julia --project=. -e 'include("test/expert/sl3_local_murthy_context.jl")'
```

Expected after implementation: both commands exit 0.

- [ ] **Step 7: Commit implementation**

Run:

```bash
git add src/algorithm/sl3_local.jl test/expert/sl3_local_murthy_resultant.jl test/expert/sl3_local_murthy_context.jl docs/superpowers/plans/2026-06-30-issue-210-local-q0-nonunit-resultant-replay.md
git commit -m "feat: route local q0 nonunit resultant replay"
```

### Task 3: Final Verification and PR Preparation

**Files:**
- Modify only if verification uncovers a defect in the Task 2 files.

**Interfaces:**
- Consumes: all Task 1 and Task 2 outputs.
- Produces: verified branch ready for PR.

- [ ] **Step 1: Run issue-required focused verification**

Run:

```bash
julia --project=. -e 'include("test/expert/sl3_local_murthy_resultant.jl")'
```

Expected: command exits 0.

- [ ] **Step 2: Run package verification**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: command exits 0.

- [ ] **Step 3: Review final diff**

Run:

```bash
git status --short
git diff --stat origin/main...HEAD
git diff --check
```

Expected: only issue #210 source, test, and Superpowers docs changes; no whitespace errors.

- [ ] **Step 4: Use finishing workflow**

Choose "Push and create a Pull Request" when the finishing workflow presents
integration options.

## Plan Self-Review

- The plan covers the spec requirements for supplied local witness routing,
  ordinary extraction preservation, exact elementary identities, child q0-unit
  replay, and negative controls.
- No placeholder markers remain.
- The task interfaces use existing internal types and do not require exported
  APIs.
