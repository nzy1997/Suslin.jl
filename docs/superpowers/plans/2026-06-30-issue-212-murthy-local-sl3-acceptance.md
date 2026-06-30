# Issue 212 Murthy Local SL3 Acceptance Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add #182 closeout acceptance coverage and staged-boundary documentation for the Murthy local `SL_3` solver.

**Architecture:** Keep production behavior unchanged unless the new acceptance checks expose a real gap. Add direct non-fixture ordinary/local witness tests to the existing Murthy expert acceptance file, add one public materializable ordinary-factor smoke, and update README/docs wording plus the docs-boundary smoke test.

**Tech Stack:** Julia, Oscar/AbstractAlgebra polynomial rings, existing `SL3LocalRealizationCertificate`, `SL3LocalMurthyInputContext`, `SL3LocalElementaryFactorReplay`, `elementary_factorization`, and Documenter markdown docs.

## Global Constraints

- Murthy local `SL_3` certificates are supported for the proven ordinary/local-witness contract.
- Ordinary factor vectors are exposed only when the certificate can materialize them over the base ring.
- Nontrivial local-witness acceptance must use localized/denominator-cleared replay, not ordinary `verify_factorization(A, factors)` over the original base ring.
- Quillen automatic patching (#183), general `SL_3` (#184), ECP (#185), recursive `SL_n` (#186), and full public Park-Woodburn acceptance (#187) remain open.
- Laurent/ToricBuilder mainline acceptance and Steinberg factor-count optimization remain staged.
- Do not implement #183-#187, Steinberg factor-count optimization, or broader Laurent/ToricBuilder support.
- Skip the optional #211 smoke check unless a callable Murthy-to-Quillen adapter already exists.
- Focused verification command is `julia --project=. -e 'include("test/expert/sl3_local_murthy_gupta.jl")'`.
- Full verification commands are `julia --project=. test/runtests.jl all` and `julia --project=. -e 'using Pkg; Pkg.test()'`.

---

## File Structure

- Modify `test/expert/sl3_local_murthy_gupta.jl`: add #182 acceptance helpers, non-fixture ordinary branch cases, local witness replay cases, and negative controls.
- Modify `test/public/factorization_driver_shell.jl`: add a public ordinary/materializable non-fixture Murthy case.
- Modify `README.md`: update current-scope boundary for #182 support and #183-#187 staged work.
- Modify `docs/src/index.md`: mirror README scope boundary.
- Modify `test/expert/polynomial_normality_support_boundary.jl`: update exact required phrases for the new docs wording.

### Task 1: Expert Murthy Closeout Acceptance

**Files:**
- Modify: `test/expert/sl3_local_murthy_gupta.jl`

**Interfaces:**
- Consumes: `realize_sl3_local_certificate`, `realize_sl3_local`, `sl3_local_murthy_input_context`, `verify_sl3_local_realization`, `verify_sl3_local_elementary_factor_replay`.
- Produces: an Issue #182 closeout testset that proves non-fixture ordinary branch coverage, local witness replay coverage, and required negative controls.

- [ ] **Step 1: Add expert acceptance helpers**

Add these helpers after `_mg_acceptance_special_form_target`:

```julia
function _mg_acceptance_captured_error(f)
    try
        f()
        return nothing
    catch err
        return err
    end
end

function _mg_acceptance_as_namedtuple(value)
    names = propertynames(value)
    return NamedTuple{names}(Tuple(getproperty(value, name) for name in names))
end

function _mg_acceptance_local_unit_witness(unit, residue_unit, residue_inverse, generator, coefficient, X)
    return (;
        context = (;
            kind = :localization_at_maximal_ideal,
            selected_variable = X,
            maximal_ideal_generators = (generator,),
        ),
        unit,
        residue_unit,
        residue_inverse,
        maximal_ideal_generators = (generator,),
        residue_difference_coefficients = (coefficient,),
        global_unit = is_unit(unit),
    )
end

function _mg_acceptance_assert_local_replay(certificate)
    @test Suslin.verify_sl3_local_realization(certificate)
    reduction = certificate.witness.reduction
    replay = reduction.local_factor_replay
    @test replay.target == certificate.target
    @test replay.factors == certificate.factors
    @test replay.mode == :denominator_cleared
    @test replay.materialized_factors === nothing
    @test all(factor -> factor isa Suslin.SL3LocalElementaryFactor, replay.factors)
    @test Suslin.verify_sl3_local_elementary_factor_replay(replay)

    denominator_index = findfirst(factor -> factor.denominator != one(factor.R), replay.factors)
    @test denominator_index !== nothing
    if denominator_index !== nothing
        @test_throws ArgumentError Suslin.sl3_local_materialize_elementary_factor(
            replay.factors[denominator_index],
        )
    end
    return replay
end
```

- [ ] **Step 2: Add the #182 closeout positive cases**

Append this testset after the existing "Issue 61 Murthy-Gupta local SL3 acceptance" testset:

```julia
@testset "Issue 182 Murthy local SL3 closeout acceptance" begin
    R, (X,) = Oscar.polynomial_ring(QQ, ["X"])

    p_qdegree = X^2 + X + 1
    q_qdegree = X * p_qdegree + one(R)
    qdegree_target = _mg_acceptance_special_form_target(
        R,
        p_qdegree,
        q_qdegree,
        -one(R),
        -X,
    )
    qdegree_cert = Suslin.realize_sl3_local_certificate(qdegree_target, X)
    @test qdegree_cert.branch == :murthy_q0_unit
    @test qdegree_cert.witness.normalization !== nothing
    _mg_acceptance_assert_q0_unit_certificate(qdegree_cert)
    _mg_acceptance_assert_elementary_sequence(
        qdegree_target,
        Suslin.realize_sl3_local(qdegree_target, X),
    )

    p_q0_unit = X^2 + X + 1
    q0_unit_target = _mg_acceptance_special_form_target(
        R,
        p_q0_unit,
        one(R),
        p_q0_unit^2 - one(R),
        p_q0_unit,
    )
    q0_unit_cert = Suslin.realize_sl3_local_certificate(q0_unit_target, X)
    @test q0_unit_cert.branch == :murthy_q0_unit
    @test q0_unit_cert.witness.normalization === nothing
    _mg_acceptance_assert_q0_unit_certificate(q0_unit_cert)
    _mg_acceptance_assert_elementary_sequence(
        q0_unit_target,
        Suslin.realize_sl3_local(q0_unit_target, X),
    )

    p_nonunit = X^2 + X + 1
    q_nonunit = X
    nonunit_witness = (;
        p0 = one(R),
        q0 = zero(R),
        p_prime = one(R),
        q_prime = X + 1,
        resultant = one(R),
        p_prime_degree = 0,
        q_prime_degree = 1,
        branch_unit = one(R),
        case1_entries = (;
            p = p_nonunit + X + 1,
            q = q_nonunit + one(R),
            r = X + 1,
            s = one(R),
        ),
    )
    nonunit_target = _mg_acceptance_special_form_target(
        R,
        p_nonunit,
        q_nonunit,
        X + 1 + p_nonunit,
        X + 1,
    )
    supplied_nonunit_cert = Suslin.realize_sl3_local_certificate(
        nonunit_target,
        X;
        murthy_q0_nonunit_witness = nonunit_witness,
    )
    _mg_acceptance_assert_resultant_certificate(
        supplied_nonunit_cert;
        expected_source = :supplied_bezout_witness,
    )
    _mg_acceptance_assert_elementary_sequence(
        nonunit_target,
        Suslin.realize_sl3_local(
            nonunit_target,
            X;
            murthy_q0_nonunit_witness = nonunit_witness,
        ),
    )

    RU, (u, Y) = Oscar.polynomial_ring(QQ, ["u", "X"])

    local_q0_unit = Y + u + 2
    local_q0_unit_p = Y * local_q0_unit + one(RU)
    local_q0_unit_target = _mg_acceptance_special_form_target(
        RU,
        local_q0_unit_p,
        local_q0_unit,
        Y + local_q0_unit_p * Y,
        local_q0_unit_p,
    )
    local_q0_unit_witness = (;
        p0 = one(RU),
        q0 = u + 2,
        local_unit_witness = _mg_acceptance_local_unit_witness(
            u + 2,
            RU(2),
            RU(QQ(1) // QQ(2)),
            u,
            one(RU),
            Y,
        ),
        formal_right_e21_coefficient = "-1/(u + 2)",
    )
    local_q0_unit_context = Suslin.sl3_local_murthy_input_context(
        local_q0_unit_target,
        Y;
        witness = local_q0_unit_witness,
    )
    local_q0_unit_cert = Suslin.realize_sl3_local_certificate(local_q0_unit_context)
    @test local_q0_unit_cert.branch == :murthy_q0_unit
    _mg_acceptance_assert_local_replay(local_q0_unit_cert)

    local_nonunit_q = Y + 2 * u
    local_nonunit_p = Y * local_nonunit_q + one(RU)
    local_nonunit_target = _mg_acceptance_special_form_target(
        RU,
        local_nonunit_p,
        local_nonunit_q,
        Y + local_nonunit_p * Y,
        local_nonunit_p,
    )
    local_nonunit_witness = (;
        p0 = one(RU),
        q0 = 2 * u,
        p_prime = one(RU),
        q_prime = Y,
        resultant = one(RU),
        p_prime_degree = 0,
        q_prime_degree = 1,
        branch_unit = one(RU) + 2 * u,
        branch_unit_witness = _mg_acceptance_local_unit_witness(
            one(RU) + 2 * u,
            one(RU),
            one(RU),
            u,
            RU(2),
            Y,
        ),
        case1_entries = (;
            p = local_nonunit_p + Y,
            q = local_nonunit_q + one(RU),
            r = Y,
            s = one(RU),
        ),
    )
    local_nonunit_context = Suslin.sl3_local_murthy_input_context(
        local_nonunit_target,
        Y;
        witness = local_nonunit_witness,
    )
    local_nonunit_cert = Suslin.realize_sl3_local_certificate(local_nonunit_context)
    @test local_nonunit_cert.branch == :murthy_q0_nonunit_bezout_resultant
    _mg_acceptance_assert_resultant_certificate(
        local_nonunit_cert;
        expected_source = :supplied_bezout_witness,
    )
    _mg_acceptance_assert_local_replay(local_nonunit_cert)
end
```

- [ ] **Step 3: Add closeout negative controls**

Append this testset after the positive closeout testset:

```julia
@testset "Issue 182 Murthy local SL3 closeout negative controls" begin
    R, (X,) = Oscar.polynomial_ring(QQ, ["X"])
    determinant_bad = _mg_acceptance_special_form_target(
        R,
        X + one(R),
        zero(R),
        zero(R),
        one(R),
    )
    determinant_err = _mg_acceptance_captured_error(
        () -> Suslin.realize_sl3_local_certificate(determinant_bad, X),
    )
    @test determinant_err isa ArgumentError
    @test occursin("determinant", sprint(showerror, determinant_err))

    nonmonic = _mg_acceptance_special_form_target(R, 2 * X + one(R), X, R(2), one(R))
    nonmonic_err = _mg_acceptance_captured_error(
        () -> Suslin.realize_sl3_local_certificate(nonmonic, X),
    )
    @test nonmonic_err isa ArgumentError
    @test occursin("p must be monic", sprint(showerror, nonmonic_err))

    RU, (u, Y) = Oscar.polynomial_ring(QQ, ["u", "X"])
    q_local = Y + u + 2
    p_local = Y * q_local + one(RU)
    missing_witness_target = _mg_acceptance_special_form_target(
        RU,
        p_local,
        q_local,
        Y + p_local * Y,
        p_local,
    )
    missing_witness_err = _mg_acceptance_captured_error(
        () -> Suslin.sl3_local_murthy_input_context(missing_witness_target, Y),
    )
    @test missing_witness_err isa ArgumentError
    @test occursin("local-unit witness", sprint(showerror, missing_witness_err))

    q_nonunit = Y + 2 * u
    p_nonunit = Y * q_nonunit + one(RU)
    unsupported_extraction_target = _mg_acceptance_special_form_target(
        RU,
        p_nonunit,
        q_nonunit,
        Y + p_nonunit * Y,
        p_nonunit,
    )
    unsupported_extraction_err = _mg_acceptance_captured_error(
        () -> Suslin.sl3_local_murthy_input_context(unsupported_extraction_target, Y),
    )
    @test unsupported_extraction_err isa ArgumentError
    @test occursin("unsupported local Bezout/resultant extraction", sprint(showerror, unsupported_extraction_err)) ||
        occursin("staged local SL_3 solver failure", sprint(showerror, unsupported_extraction_err))

    p_supported = X^2 + X + 1
    q_supported = one(R)
    supported_target = _mg_acceptance_special_form_target(
        R,
        p_supported,
        q_supported,
        p_supported^2 - one(R),
        p_supported,
    )
    supported_cert = Suslin.realize_sl3_local_certificate(supported_target, X)
    corrupted_factors = copy(supported_cert.factors)
    corrupted_factors[1] =
        corrupted_factors[1] * elementary_matrix(3, 1, 3, one(R), R)
    corrupted_cert = Suslin.SL3LocalRealizationCertificate(
        supported_cert.target,
        supported_cert.branch,
        corrupted_factors,
        supported_cert.selected_variable,
        supported_cert.witness,
    )
    @test !Suslin.verify_sl3_local_realization(corrupted_cert)
end
```

- [ ] **Step 4: Run focused acceptance**

Run:

```bash
julia --project=. -e 'include("test/expert/sl3_local_murthy_gupta.jl")'
```

Expected: pass. This issue is a closeout gate over already-landed dependency behavior, so no production-code RED is expected.

### Task 2: Public Ordinary Boundary and Documentation

**Files:**
- Modify: `test/public/factorization_driver_shell.jl`
- Modify: `README.md`
- Modify: `docs/src/index.md`
- Modify: `test/expert/polynomial_normality_support_boundary.jl`

**Interfaces:**
- Consumes: public `elementary_factorization`, `verify_factorization`, docs wording.
- Produces: public smoke for materializable non-fixture Murthy factors and docs that state the #182 support boundary.

- [ ] **Step 1: Add the public ordinary/materializable smoke**

In `test/public/factorization_driver_shell.jl`, after the first supported
`fast_local_sl3` case, add:

```julia
    nonfixture_p = X^2 + X + one(R)
    nonfixture_q = X
    nonfixture_local_sl3 = matrix(R, [
        nonfixture_p       nonfixture_q zero(R);
        X + one(R) + nonfixture_p X + one(R) zero(R);
        zero(R)            zero(R)      one(R)
    ])
    nonfixture_factors = elementary_factorization(nonfixture_local_sl3)
    @test verify_factorization(nonfixture_local_sl3, nonfixture_factors)
    nonfixture_cert = Suslin._polynomial_factorization_route_certificate(nonfixture_local_sl3)
    @test nonfixture_cert.route == :fast_local_sl3
    @test nonfixture_cert.evidence isa Suslin.SL3LocalRealizationCertificate
    @test nonfixture_cert.evidence.branch == :murthy_q0_nonunit_bezout_resultant
    @test nonfixture_factors == nonfixture_cert.factors
```

- [ ] **Step 2: Update README current scope**

In `README.md`, replace the staged-boundary bullet that currently says
`#182 through #187 are complete. Murthy local solving...` with wording that
states:

```markdown
- The Murthy local `SL_3` solver (#182) is supported for the proven
  ordinary/local-witness contract: ordinary factor vectors are exposed only
  when the certificate can materialize them over the base ring, while
  nontrivial local-witness cases are verified through localized
  denominator-cleared certificate replay.
- The implementation is not yet the full Park-Woodburn algorithm for arbitrary
  `SL_n(k[x_1, ..., x_m])`, `n >= 3`: Quillen automatic patching (#183),
  general `SL_3` (#184), the general ECP reducer (#185), recursive `SL_n`
  (#186), full public Park-Woodburn acceptance (#187), coefficient-ring
  support beyond exact field-backed ordinary polynomial rings, arbitrary
  Laurent `GL_n` determinant correction, Laurent/ToricBuilder mainline
  acceptance, and Steinberg factor-count optimization remain staged
  boundaries.
```

- [ ] **Step 3: Mirror docs index wording**

Apply the same two bullets to `docs/src/index.md` under `## Scope`.

- [ ] **Step 4: Update documentation smoke phrases**

In `test/expert/polynomial_normality_support_boundary.jl`, replace the two old
phrases:

```julia
        "Murthy local solving, general Quillen local realizability, the general ECP reducer",
        "Laurent/ToricBuilder mainline acceptance, and Steinberg factor-count optimization remain staged boundaries",
```

with:

```julia
        "Murthy local `SL_3` solver (#182) is supported for the proven ordinary/local-witness contract",
        "ordinary factor vectors are exposed only when the certificate can materialize them over the base ring",
        "Quillen automatic patching (#183), general `SL_3` (#184), the general ECP reducer (#185)",
        "recursive `SL_n` (#186), full public Park-Woodburn acceptance (#187)",
        "Laurent/ToricBuilder mainline acceptance, and Steinberg factor-count optimization remain staged boundaries",
```

- [ ] **Step 5: Run focused public/docs checks**

Run:

```bash
julia --project=. -e 'include("test/public/factorization_driver_shell.jl")'
julia --project=. -e 'include("test/expert/polynomial_normality_support_boundary.jl")'
```

Expected: both commands exit 0.

### Task 3: Full Verification and Cleanup

**Files:**
- Verify all changed files.

**Interfaces:**
- Consumes: Tasks 1 and 2.
- Produces: final verified branch ready for PR.

- [ ] **Step 1: Run issue focused acceptance**

Run:

```bash
julia --project=. -e 'include("test/expert/sl3_local_murthy_gupta.jl")'
```

Expected: exits 0.

- [ ] **Step 2: Run full test runner**

Run:

```bash
julia --project=. test/runtests.jl all
```

Expected: exits 0.

- [ ] **Step 3: Run package test entry**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: exits 0.

- [ ] **Step 4: Check diff hygiene**

Run:

```bash
git diff --check
git status --short
```

Expected: no whitespace errors; only intentional source/docs/test changes are present before commit.
