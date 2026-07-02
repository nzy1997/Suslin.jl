# Issue 248 Public ECP Column Route Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Route public ordinary-polynomial unimodular column reduction through the checked general ECP pipeline while preserving legacy simple and Laurent paths.

**Architecture:** Keep `ECPColumnReductionCertificate` as the public certificate type and add a replayable `:ecp_pipeline` stage containing #243-#247 evidence plus route metadata. Preserve unit-entry, witness-unit, embedded three-block, and Laurent routes before attempting the general ECP route. For length greater than three link endpoints, record and verify direct elementary endpoint transport metadata instead of requiring unsupported #186 matrix factorization.

**Tech Stack:** Julia, Oscar exact ordinary polynomial rings, existing Suslin ECP and polynomial route certificate internals, Test stdlib.

## Global Constraints

- Repository has no `AGENTS.md`; follow README test commands.
- Preserve `reduce_unimodular_column(v, R)` returning a factor sequence.
- Preserve `ecp_column_reduction_certificate(v, R)` returning `ECPColumnReductionCertificate`.
- Preserve unit-entry, witness-unit, embedded three-block, and Laurent fallback behavior.
- Non-unimodular columns must fail before route work with `ArgumentError("v must be a unimodular column")`.
- Unsupported but unimodular ordinary-polynomial columns must fail with a staged `ArgumentError`.
- Do not implement #186 recursive matrix factorization.
- Do not implement #187 final public Park-Woodburn matrix acceptance.

---

## File Structure

- Modify `src/algorithm/column_reduction.jl`: add the general public ECP stage builder, replay support, staged diagnostics, and direct elementary endpoint transport metadata.
- Modify `src/algorithm/polynomial_column_peel.jl`: store an ECP column certificate on peel steps while keeping the existing positional constructor compatible.
- Modify `test/expert/elementary_column_property.jl`: add public ECP route tests, route metadata checks, and tamper rejection.
- Modify `test/expert/unimodular_reduction_exact.jl`: preserve old exact/Laurent behavior and unsupported diagnostics.
- `test/runtests.jl` already includes both expert files; modify only if test file placement changes.

### Task 1: Public Route Tests

**Files:**
- Modify: `test/expert/elementary_column_property.jl`
- Modify: `test/expert/unimodular_reduction_exact.jl`

**Interfaces:**
- Consumes: existing `Suslin.reduce_unimodular_column`, `Suslin.ecp_column_reduction_certificate`, `Suslin.verify_ecp_column_reduction`, `Suslin._reduce_supported_unimodular_column_certificate`, and `Suslin._reduce_via_supported_three_block_certificate`.
- Produces: failing acceptance tests for `:ecp_pipeline` route metadata and staged diagnostics.

- [ ] **Step 1: Add the failing non-fixture length-four acceptance helper**

Add this helper near the other helpers in `test/expert/elementary_column_property.jl`:

```julia
function _ecp_acceptance_length4_general_case()
    R, (x, y) = Oscar.polynomial_ring(QQ, ["x", "y"])
    column = [
        x * y,
        x * (one(R) - y),
        (one(R) - x) * y,
        (one(R) - x) * (one(R) - y),
    ]
    for (target, source, coeff) in (
        (2, 4, -x * y - y),
        (3, 2, -x * y - y),
        (4, 3, x^2),
        (2, 3, x * y + y),
        (2, 3, -y - one(R)),
    )
        column[target] += coeff * column[source]
    end
    return R, column
end
```

- [ ] **Step 2: Add the failing public ECP route test**

In `test/expert/elementary_column_property.jl`, add this block inside the `"public ECP unimodular-column pipeline"` testset after the canonical staged checks:

```julia
    general_R, general_column = _ecp_acceptance_length4_general_case()
    @test length(general_column) > 3
    @test Suslin.is_unimodular_column(general_column, general_R)
    @test !any(is_unit, general_column)
    @test Suslin._reduce_supported_unimodular_column_certificate(general_column, general_R) === nothing
    @test Suslin._reduce_via_supported_three_block_certificate(general_column, general_R) === nothing

    general_cert = Suslin.ecp_column_reduction_certificate(general_column, general_R)
    @test Suslin.verify_ecp_column_reduction(general_cert)
    @test general_cert.stages[end].kind == :ecp_pipeline
    @test general_cert.stages[end].route_metadata.route == :general_ecp_pipeline
    @test general_cert.stages[end].route_metadata.link_route_mode == :direct_elementary
    @test general_cert.stages[end].route_metadata.normalized_column_length == length(general_column)
    @test _ecp_acceptance_apply(general_cert.factors, general_column, general_R) ==
          _ecp_acceptance_target(general_R, length(general_column))
    @test Suslin.reduce_unimodular_column(general_column, general_R) == general_cert.factors

    tampered_factors = copy(general_cert.factors)
    tampered_factors[1] = identity_matrix(general_R, length(general_column))
    tampered_cert = Suslin.ECPColumnReductionCertificate(
        general_cert.original_column,
        general_cert.ring,
        general_cert.stages,
        tampered_factors,
        general_cert.final_column,
        general_cert.verification,
    )
    @test !Suslin.verify_ecp_column_reduction(tampered_cert)
```

- [ ] **Step 3: Strengthen staged diagnostic coverage**

In `test/expert/unimodular_reduction_exact.jl`, change the unsupported ordinary-column assertion to require the staged ECP diagnostic substring while preserving the existing exact-reduction substring:

```julia
    @test occursin("unsupported exact unimodular column reduction", sprint(showerror, unsupported_err))
    @test occursin("general ECP pipeline", sprint(showerror, unsupported_err))
    @test !occursin("not unimodular", sprint(showerror, unsupported_err))
```

- [ ] **Step 4: Run RED tests**

Run:

```bash
julia --project=. -e 'include("test/expert/elementary_column_property.jl")'
julia --project=. -e 'include("test/expert/unimodular_reduction_exact.jl")'
```

Expected: the first command fails because `:ecp_pipeline` and direct elementary route metadata are not implemented; the second command fails until the unsupported diagnostic mentions the general ECP pipeline.

- [ ] **Step 5: Commit**

```bash
git add test/expert/elementary_column_property.jl test/expert/unimodular_reduction_exact.jl
git commit -m "test: cover public general ECP column route"
```

### Task 2: General ECP Public Certificate Route

**Files:**
- Modify: `src/algorithm/column_reduction.jl`
- Test: `test/expert/elementary_column_property.jl`
- Test: `test/expert/unimodular_reduction_exact.jl`

**Interfaces:**
- Consumes: #243-#247 helpers `ecp_input_context`, `ecp_monicity_normalization`, `ecp_link_witness`, `ecp_link_step_certificate`, and `ecp_induction_normality_certificate`.
- Produces: `_reduce_via_general_ecp_pipeline_certificate`, `:ecp_pipeline` stage replay, direct elementary link transport route metadata, and staged ordinary-polynomial failure diagnostics.

- [ ] **Step 1: Add direct endpoint-transport route support**

In `src/algorithm/column_reduction.jl`, update `_ecp_link_step_resolve_route_mode`, `_ecp_link_step_supported_family`, `_ecp_link_step_endpoint_transport`, and `_ecp_link_step_segment_verification` so:

```julia
route_mode in (:auto, :legacy_fixture, :polynomial_sl3, :direct_elementary)
```

`route_mode = :auto` keeps legacy fixture matching first, uses `:direct_elementary` for columns of length greater than three, and uses `:polynomial_sl3` otherwise. The direct route returns `raw_factors` as the endpoint transport factors, stores empty `sl3_route_*` collections, and records direct route metadata in the segment `support_family` plus the outer `:ecp_pipeline` stage `route_metadata`. Do not store a standalone NamedTuple in `sl3_route_metadata`; that field remains aligned with `sl3_route_matrices`.

```julia
route_metadata = (;
    source = :direct_elementary_endpoint_transport,
    route = :direct_elementary,
    factor_count = length(raw_factors),
)
```

In verification, treat `:direct_elementary_endpoint_transport` as supported when the forward factors are elementary over the route ring and the endpoint map checks pass.

- [ ] **Step 2: Add ECP pipeline stage builder**

Add `_reduce_via_general_ecp_pipeline_certificate(column, R)` that:

1. rejects Laurent rings by returning `nothing`;
2. creates `context = ecp_input_context(column, R)`;
3. chooses `selected_variable = context.selected_variable === nothing ? first(context.variable_order) : context.selected_variable`;
4. creates `normalization = ecp_monicity_normalization(context; selected_variable)`;
5. throws a staged `ArgumentError` if normalization returns `ECPMonicityNormalizationFailure`;
6. creates `link_witness = ecp_link_witness(normalization)` and throws staged `ArgumentError` if extraction returns `ECPLinkWitnessExtractionFailure`;
7. creates `link_step = ecp_link_step_certificate(collect(normalization.normalized_column), R; link_witness)`;
8. creates `induction = ecp_induction_normality_certificate(collect(normalization.normalized_column), R; link_step)`;
9. inverse-substitutes `induction.final_factors` through `normalization.inverse_substitution`;
10. returns a certificate-stage tuple with `kind = :ecp_pipeline`, `route_metadata`, nested evidence, `factors`, and `output_column`.

The final factors must be:

```julia
vcat(inverse_substituted_induction_factors,
     normalization.inverse_substituted_coordinate_move_factors)
```

- [ ] **Step 3: Route polynomial certificate construction**

Refactor `_reduce_polynomial_unimodular_column_exact_certificate(column, R)` to:

1. return unit-entry and witness-unit certificates first;
2. return embedded three-block certificates first for length greater than three;
3. attempt `_reduce_via_general_ecp_pipeline_certificate(column, R)`;
4. fall back to existing monicity-normalization exact support if the general ECP route is unsupported;
5. return `nothing` only after all supported ordinary-polynomial stages fail.

Keep Laurent callers on `_reduce_laurent_unimodular_column_certificate`.

- [ ] **Step 4: Add `:ecp_pipeline` replay**

Extend `_ecp_replay_stage(stage, input_column, R)` with a `stage.kind == :ecp_pipeline` branch that recomputes and verifies:

- stage keys,
- `verify_ecp_monicity_normalization(stage.normalization)`,
- `verify_ecp_link_witness(stage.link_witness)`,
- `verify_ecp_link_step_certificate(stage.link_step)`,
- `verify_ecp_induction_normality_certificate(stage.induction_normality)`,
- inverse-substituted induction factors,
- final factor sequence,
- exact final output.

- [ ] **Step 5: Update unsupported diagnostics**

Update ordinary-polynomial diagnostics and `_unsupported_unimodular_column_reduction_message` so unsupported but unimodular ordinary columns mention the general ECP pipeline while retaining the prefix `"unsupported exact unimodular column reduction"`.

- [ ] **Step 6: Run GREEN tests**

Run:

```bash
julia --project=. -e 'include("test/expert/elementary_column_property.jl")'
julia --project=. -e 'include("test/expert/unimodular_reduction_exact.jl")'
```

Expected: both commands exit 0.

- [ ] **Step 7: Commit**

```bash
git add src/algorithm/column_reduction.jl test/expert/elementary_column_property.jl test/expert/unimodular_reduction_exact.jl
git commit -m "feat: route public columns through general ECP"
```

### Task 3: Column Peel Metadata Compatibility

**Files:**
- Modify: `src/algorithm/polynomial_column_peel.jl`
- Modify: `test/expert/park_woodburn_polynomial_column_peel.jl` if required by constructor shape

**Interfaces:**
- Consumes: `ecp_column_reduction_certificate(last_column, R)`.
- Produces: `PolynomialColumnPeelStep.left_certificate` while keeping the old constructor usable.

- [ ] **Step 1: Store the column certificate**

Add a `left_certificate` field to `PolynomialColumnPeelStep`. Add an outer constructor with the old eight-argument positional signature that sets `left_certificate = nothing` so existing tests and internal helper constructors keep working.

- [ ] **Step 2: Use the certificate in peel steps**

In `_polynomial_column_peel_step`, replace:

```julia
left_factors = reduce_unimodular_column(last_column, R)
```

with:

```julia
left_certificate = ecp_column_reduction_certificate(last_column, R)
left_factors = left_certificate.factors
```

and pass `left_certificate` into the stored step.

- [ ] **Step 3: Verify optional certificate metadata**

Update `_is_valid_polynomial_column_peel_step_data` only if needed to accept the new field through the step object. In `_polynomial_column_peel_core_verification`, add a check that any non-`nothing` `left_certificate` verifies, matches the step last column, and has factors equal to `left_factors`.

- [ ] **Step 4: Run peel and required route tests**

Run:

```bash
julia --project=. -e 'include("test/expert/park_woodburn_polynomial_column_peel.jl")'
julia --project=. -e 'include("test/expert/elementary_column_property.jl")'
```

Expected: both commands exit 0.

- [ ] **Step 5: Commit**

```bash
git add src/algorithm/polynomial_column_peel.jl test/expert/park_woodburn_polynomial_column_peel.jl
git commit -m "feat: expose ECP route metadata to column peel"
```

### Task 4: Final Verification

**Files:**
- No planned source edits unless verification reveals a defect.

**Interfaces:**
- Consumes: all changes from Tasks 1-3.
- Produces: final evidence for the PR.

- [ ] **Step 1: Run issue verification**

```bash
julia --project=. -e 'include("test/expert/elementary_column_property.jl")'
julia --project=. -e 'include("test/expert/unimodular_reduction_exact.jl")'
```

Expected: both commands exit 0.

- [ ] **Step 2: Run package verification**

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: command exits 0.

- [ ] **Step 3: Check diff hygiene**

```bash
git diff --check
git status --short
```

Expected: `git diff --check` exits 0 and `git status --short` shows only intentional changes before final commit.

- [ ] **Step 4: Commit if verification fixes were needed**

If Step 1-3 required fixes, commit them:

```bash
git add src/algorithm/column_reduction.jl src/algorithm/polynomial_column_peel.jl test/expert/elementary_column_property.jl test/expert/unimodular_reduction_exact.jl test/expert/park_woodburn_polynomial_column_peel.jl
git commit -m "fix: stabilize public ECP column route verification"
```
