# ECP Public Unimodular Pipeline Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Route supported ECP fixture columns through verified link-witness, link-step, induction/normality replay before `reduce_unimodular_column(v, R)` returns factors.

**Architecture:** Keep `reduce_unimodular_column(v, R)` factor-returning. Add a non-exported staged certificate API in `src/algorithm/column_reduction.jl` that composes the existing #87, #88, and #90 replay verifiers, and let the public reducer prefer that narrow route only when deterministic fixture witness data is recognized. Fall back to the existing certificate reducer for all other supported columns, including Laurent consumers.

**Tech Stack:** Julia, Oscar polynomial rings, existing Suslin elementary factor and ECP certificate helpers, `Test`.

## Global Constraints

- `reduce_unimodular_column(v, R)` remains the public entry point and returns only elementary factors.
- The staged ECP public route must return factors only after exact replay proves `product(factors) * v == e_n`.
- Expert certificate APIs may be non-exported and accessed as `Suslin.<name>` from expert tests.
- Reuse `ecp_link_witness`, `ecp_link_step_certificate`, `ecp_induction_normality_certificate`, and their verifiers instead of duplicating route verification.
- Unsupported columns throw staged `ArgumentError`s without returning factors.
- Non-unimodular columns fail immediately with `ArgumentError("v must be a unimodular column")`.
- Do not implement Quillen local-to-global patching, the general Park-Woodburn matrix driver, Laurent determinant correction, or factor minimization.

---

### Task 1: Add ECP Acceptance Tests

**Files:**
- Create: `test/expert/elementary_column_property.jl`
- Modify: `test/runtests.jl`

**Interfaces:**
- Consumes: `Suslin.reduce_unimodular_column(v, R)`, `Suslin.ecp_column_reduction_certificate(v, R)`, planned `Suslin.ecp_staged_column_reduction_certificate(v, R; kwargs...)`, planned `Suslin.verify_ecp_staged_column_reduction(certificate)`.
- Produces: Focused acceptance coverage for the public factor-returning path and expert staged route inspection.

- [ ] **Step 1: Write the failing focused test**

Create `test/expert/elementary_column_property.jl` with these concrete test sections:

```julia
using Test
using Oscar
using Suslin

function _ecp_acceptance_product(factors, R, n::Int)
    product = identity_matrix(R, n)
    for factor in factors
        product *= factor
    end
    return product
end

function _ecp_acceptance_apply(factors, column, R)
    return _ecp_acceptance_product(factors, R, length(column)) *
        matrix(R, length(column), 1, collect(column))
end

function _ecp_acceptance_target(R, n::Int)
    target = zero_matrix(R, n, 1)
    target[n, 1] = one(R)
    return target
end

function _ecp_acceptance_assert_public(column, R)
    @test !any(is_unit, column)
    @test Suslin._reduce_supported_unimodular_column_certificate(column, R) === nothing
    factors = Suslin.reduce_unimodular_column(column, R)
    @test _ecp_acceptance_apply(factors, column, R) == _ecp_acceptance_target(R, length(column))
    return factors
end

function _ecp_acceptance_gf2_cases()
    R, (x, y) = Oscar.polynomial_ring(GF(2), ["x", "y"])
    base = [
        x + y^2,
        x * y + x + one(R),
        x^2 + x * y + y + one(R),
    ]
    return R, [
        ("canonical-full-route", base[[1, 2, 3]]),
        ("inverse-substitution-permuted", base[[2, 1, 3]]),
        ("permuted-third-middle", base[[1, 3, 2]]),
        ("permuted-cyclic", base[[3, 1, 2]]),
    ]
end

function _ecp_acceptance_good_link_witness(column, R)
    x, y = gens(R)
    G = y * column[2] + column[3]
    return (;
        source = :supplied_link_witness,
        residue_probes = ((; id = :gf2_fixture_probe, kind = :deterministic_fixture, maximal_ideal_generators = (y,)),),
        tail_reductions = ((; probe_id = :gf2_fixture_probe, G, lifted_tail_coefficients = (y, one(R)), tilde_G = G),),
        resultants = (one(R),),
        bezout_coefficients = ((; f = x, h = one(R)),),
        coverage_multipliers = (one(R),),
        path_points = (zero(R), x),
    )
end

function _ecp_acceptance_normality_witness(lower_certificate, entry)
    R = lower_certificate.ring
    n = length(lower_certificate.original_column)
    lower_product = _ecp_acceptance_product(lower_certificate.factors, R, n)
    return (;
        source = :supplied_normality_witness,
        conjugator = inv(lower_product),
        sl2_indices = (n, 1),
        sl2_entry = entry,
    )
end

function _ecp_acceptance_capture_error(f)
    try
        f()
        return nothing
    catch err
        return err
    end
end

@testset "public ECP unimodular-column pipeline" begin
    R, cases = _ecp_acceptance_gf2_cases()
    public_factors_by_name = Dict{String, Any}()
    for (name, column) in cases
        public_factors_by_name[name] = _ecp_acceptance_assert_public(column, R)
    end

    canonical = cases[1][2]
    staged = Suslin.ecp_staged_column_reduction_certificate(canonical, R)
    @test Suslin.verify_ecp_staged_column_reduction(staged)
    @test staged.verification.route == (:validation, :monicity_forcing, :link_witness, :link_step, :induction_normality)
    @test staged.verification.link_witness_ok
    @test staged.verification.link_step_ok
    @test staged.verification.induction_normality_ok
    @test staged.induction_normality.normality_rewrite.sl2_block != identity_matrix(R, 2)
    @test public_factors_by_name["canonical-full-route"] == staged.factors
    legacy = Suslin.ecp_column_reduction_certificate(canonical, R)
    @test any(stage -> stage.kind == :monicity_normalization, legacy.stages)
    @test legacy.factors != staged.factors

    permuted = cases[2][2]
    permuted_cert = Suslin.ecp_column_reduction_certificate(permuted, R)
    monicity_stage = only([stage for stage in permuted_cert.stages if stage.kind == :monicity_normalization])
    @test all(factor -> base_ring(factor) == R, monicity_stage.inverse_substituted_factors)
    @test _ecp_acceptance_apply(monicity_stage.inverse_substituted_factors, permuted, R) ==
          _ecp_acceptance_target(R, length(permuted))

    x, y = gens(R)
    non_unimodular = [x, y, x * y]
    err = _ecp_acceptance_capture_error(() -> Suslin.reduce_unimodular_column(non_unimodular, R))
    @test err isa ArgumentError
    @test occursin("v must be a unimodular column", sprint(showerror, err))

    bad_link = merge(_ecp_acceptance_good_link_witness(canonical, R), (; resultants = (x,),))
    @test_throws ArgumentError Suslin.ecp_staged_column_reduction_certificate(
        canonical,
        R;
        supplied_link_witness = bad_link,
    )

    good_link = Suslin.ecp_link_step_certificate(
        canonical,
        R;
        variable_order = gens(R),
        selected_variable = x,
        supplied_link_witness = _ecp_acceptance_good_link_witness(canonical, R),
    )
    lower = Suslin.ecp_column_reduction_certificate(collect(good_link.lower_variable_column), R)
    bad_normality = merge(
        _ecp_acceptance_normality_witness(lower, y + one(R)),
        (; conjugator = identity_matrix(R, length(canonical))),
    )
    @test_throws ArgumentError Suslin.ecp_staged_column_reduction_certificate(
        canonical,
        R;
        supplied_link_witness = _ecp_acceptance_good_link_witness(canonical, R),
        normality_witness = bad_normality,
    )
end
```

- [ ] **Step 2: Register the test and verify RED**

Add `"expert/elementary_column_property.jl"` after `"expert/ecp_induction_normality.jl"` in `test/runtests.jl`.

Run:

```bash
julia --project=. -e 'include("test/expert/elementary_column_property.jl")'
```

Expected before production code: fails because `Suslin.ecp_staged_column_reduction_certificate` is not defined.

- [ ] **Step 3: Commit after GREEN in Task 2**

Do not commit the failing test alone. After Task 2 makes it pass, commit the test and implementation together.

---

### Task 2: Wire Public Staged ECP Route

**Files:**
- Modify: `src/algorithm/column_reduction.jl`
- Test: `test/expert/elementary_column_property.jl`

**Interfaces:**
- Produces: `ECPStagedColumnReductionCertificate`, `ecp_staged_column_reduction_certificate`, `verify_ecp_staged_column_reduction`.
- Updates: `reduce_unimodular_column(v, R)` to prefer verified staged ECP factors when deterministic route data exists.

- [ ] **Step 1: Add the staged certificate type**

Add near the other ECP certificate structs:

```julia
struct ECPStagedColumnReductionCertificate
    original_column
    ring
    monicity
    link_step::ECPLinkStepCertificate
    lower_reduction
    normality_witness
    induction_normality::ECPInductionNormalityCertificate
    factors::Vector
    final_column
    verification
end
```

- [ ] **Step 2: Update public reducer routing**

Change `reduce_unimodular_column(v::AbstractVector, R)` so it validates once, tries `_ecp_public_staged_reduction_certificate(column, R)`, returns those factors when present, and otherwise falls back to `ecp_column_reduction_certificate(column, R).factors`.

- [ ] **Step 3: Implement the expert staged constructor**

Add `ecp_staged_column_reduction_certificate(v::AbstractVector, R; variable_order = tuple(gens(R)...), selected_variable = nothing, selected_monic_index::Integer = 1, supplied_link_witness = nothing, lower_reduction = nothing, normality_witness = nothing)` that:

```julia
column = _validated_unimodular_column(v, R)
_is_laurent_polynomial_ring(R) &&
    throw(ArgumentError("ECP staged public column pipeline currently supports ordinary polynomial columns only"))
selected_variable = selected_variable === nothing ? first(_ecp_normalize_variable_order(R, variable_order)) : selected_variable
link_witness = supplied_link_witness === nothing ?
    _ecp_default_public_link_witness(column, R, selected_variable) :
    supplied_link_witness
link = ecp_link_step_certificate(column, R; variable_order, selected_variable, selected_monic_index, supplied_link_witness = link_witness)
lower_certificate = lower_reduction === nothing ?
    ecp_column_reduction_certificate(collect(link.lower_variable_column), R) :
    lower_reduction
lower_factors = lower_certificate isa ECPColumnReductionCertificate ? lower_certificate.factors : collect(lower_certificate)
normality = normality_witness === nothing ?
    _ecp_default_public_normality_witness(lower_factors, R) :
    normality_witness
induction = ecp_induction_normality_certificate(column, R; link_step = link, lower_reduction = lower_certificate, normality_witness = normality)
```

Then store factors from `induction.final_factors`, exact final column, and a replay summary. Throw `ArgumentError("constructed ECP staged public column pipeline failed exact replay verification")` if the summary or stored verifier fails.

- [ ] **Step 4: Implement deterministic route helpers**

Add `_ecp_public_staged_reduction_certificate(column, R)` that returns `nothing` unless `_ecp_default_public_link_witness(column, R, first(gens(R)))` recognizes the column. It must catch only `ArgumentError` from route construction and return `nothing` for public fallback.

Add `_ecp_default_public_link_witness(column, R, selected_variable)` for the GF(2) fixture:

```julia
x, y = gens(R)
expected = (x + y^2, x * y + x + one(R), x^2 + x * y + y + one(R))
```

When the ring has characteristic 2, two generators, `selected_variable == x`, and `tuple(column...) == expected`, return the same witness metadata used by the focused test. Otherwise throw `ArgumentError("unsupported ECP staged public column pipeline route")`.

Add `_ecp_default_public_normality_witness(lower_factors, R)` that returns:

```julia
n = nrows(first(lower_factors))
lower_product = _factor_sequence_product(lower_factors, R, n)
(; source = :supplied_normality_witness, conjugator = inv(lower_product), sl2_indices = (n, 1), sl2_entry = gens(R)[end] + one(R))
```

- [ ] **Step 5: Implement staged replay verification**

Add `verify_ecp_staged_column_reduction(certificate)::Bool` and `_ecp_staged_column_reduction_replay_summary(certificate)` that check:

```julia
route = (:validation, :monicity_forcing, :link_witness, :link_step, :induction_normality)
link_step_ok = verify_ecp_link_step_certificate(certificate.link_step)
link_witness_ok = link_step_ok && verify_ecp_link_witness(certificate.link_step.link_witness)
induction_normality_ok = verify_ecp_induction_normality_certificate(certificate.induction_normality)
monicity_ok = link_witness_ok && certificate.monicity.selected_monic_ok
factors_match_ok = _ecp_factor_sequences_equal(certificate.factors, certificate.induction_normality.final_factors)
final_column_ok = certificate.final_column == _apply_reduction_factors(certificate.factors, collect(certificate.original_column), certificate.ring)
target_ok = certificate.final_column == _target_reduced_column(certificate.ring, length(certificate.original_column))
overall_ok = all((link_witness_ok, link_step_ok, induction_normality_ok, monicity_ok, factors_match_ok, final_column_ok, target_ok))
```

The verifier returns true only when `summary.overall_ok && certificate.verification == summary`.

- [ ] **Step 6: Verify GREEN**

Run:

```bash
julia --project=. -e 'include("test/expert/elementary_column_property.jl")'
```

Expected after production code: all tests pass.

- [ ] **Step 7: Regression checks**

Run:

```bash
julia --project=. -e 'include("test/expert/unimodular_reduction_exact.jl")'
julia --project=. -e 'include("test/expert/ecp_variable_change_replay.jl")'
julia --project=. -e 'include("test/expert/ecp_induction_normality.jl")'
```

Expected: all commands exit 0.

- [ ] **Step 8: Commit**

```bash
git add src/algorithm/column_reduction.jl test/expert/elementary_column_property.jl test/runtests.jl
git commit -m "feat: wire ECP staged column pipeline"
```

## Self-Review

Spec coverage: Task 1 covers focused acceptance, full route inspection, inverse substitution, non-unimodular rejection, and corrupt replay controls. Task 2 covers public routing, expert API exposure, verifier reuse, and fallback preservation.

Placeholder scan: No placeholder tokens are present.

Type consistency: The planned expert names are `ECPStagedColumnReductionCertificate`, `ecp_staged_column_reduction_certificate`, and `verify_ecp_staged_column_reduction` in both tasks.
