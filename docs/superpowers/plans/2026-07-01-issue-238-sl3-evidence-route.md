# Issue 238 SL3 Evidence Route Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Route supported multivariate `3 x 3` ordinary-polynomial `SL_3` inputs through the checked #235 context, #236 witness, #237 Murthy/Quillen provider, #219 adapter consumption, and #220 global patch before returning public factors.

**Architecture:** Keep the public route tag `:quillen_patch`, but add a new internal evidence wrapper that proves the patch came from the SL3 Murthy/Quillen driver path. The route builds checked local evidence first, consumes it into a verified Quillen patch, adapts the patch through the existing polynomial patch adapter, and then lets the standard route certificate verify factor products.

**Tech Stack:** Julia, Oscar, existing Suslin factorization, SL3 local, Murthy/Quillen adapter, and Quillen patch APIs.

## Global Constraints

- Preserve the existing cheap `:fast_local_sl3` route.
- Return public factors only after all #235, #236, #237, #219, and #220 certificates replay.
- Do not bypass #236 by matching raw matrices.
- Do not return Murthy factors or local sequence factors before global patch verification.
- Do not implement general ECP, recursive `SL_n`, Laurent/ToricBuilder support, broad coordinate-change search, or Steinberg optimization.
- Unsupported determinant-one inputs must fail with a staged diagnostic naming the missing boundary.
- Verification commands required by the issue:
  `julia --project=. -e 'include("test/expert/park_woodburn_route_certificate.jl")'`
  `julia --project=. -e 'include("test/public/factorization_driver_shell.jl")'`
  `julia --project=. -e 'include("test/public/park_woodburn_polynomial_factorization.jl")'`
- Agent Desk package gate:
  `julia --project=. -e 'using Pkg; Pkg.test()'`

---

## File Structure

- Modify `src/algorithm/factorization.jl`: add the route evidence wrapper, automatic SL3 evidence route construction, `:already_handled` base-term policy verification, staged diagnostics, and route-certificate verification support.
- Modify `test/expert/park_woodburn_route_certificate.jl`: add positive SL3 evidence route checks and tamper rejection for provider, consumption, and patch metadata.
- Modify `test/public/factorization_driver_shell.jl`: add public acceptance and staged negative controls.
- Modify `test/public/park_woodburn_polynomial_factorization.jl`: add acceptance for the same multivariate SL3 evidence route.
- Keep `test/runtests.jl` unchanged unless a new test file is created.

### Task 1: Red Tests For SL3 Evidence Route Acceptance And Rejection

**Files:**
- Modify: `test/expert/park_woodburn_route_certificate.jl`
- Modify: `test/public/factorization_driver_shell.jl`
- Modify: `test/public/park_woodburn_polynomial_factorization.jl`

**Interfaces:**
- Consumes: future `Suslin.PolynomialSL3QuillenMurthyRouteEvidence` and future automatic `_polynomial_factorization_route_certificate(A)` SL3 route.
- Produces: failing tests that require public factors to come from verified provider and global patch evidence.

- [ ] **Step 1: Add shared issue-238 test helpers to the expert route test**

Insert these helpers after `_pw_corrupt_route_peel_evidence` in `test/expert/park_woodburn_route_certificate.jl`:

```julia
function _issue238_sl3_route_case()
    R, (X, r, g) = Oscar.polynomial_ring(QQ, ["X", "r", "g"])
    p = X + r * g + one(R)
    q = one(R)
    s = one(R)
    lower = X + r * g
    A = matrix(R, [
        p q zero(R);
        lower s zero(R);
        zero(R) zero(R) one(R)
    ])
    @assert det(A) == one(R)
    return (; R, X, r, g, p, q, s, lower, A)
end

function _issue238_assert_sl3_route_evidence(cert, A)
    @test cert.route == :quillen_patch
    @test cert.status == :supported
    @test cert.evidence isa Suslin.PolynomialSL3QuillenMurthyRouteEvidence
    @test Suslin._verify_polynomial_factorization_route_certificate(cert)
    @test verify_factorization(A, cert.factors)

    evidence = cert.evidence
    @test evidence.target == A
    @test evidence.route == :quillen_patch
    @test Suslin._verify_sl3_realization_input_context(evidence.context)
    @test Suslin._verify_sl3_local_form_witness_selection(evidence.witness_selection)
    @test Suslin._verify_sl3_murthy_quillen_local_evidence_provider(
        evidence.local_evidence_provider,
    )
    @test Suslin.verify_quillen_murthy_adapter_consumption(evidence.quillen_consumption)
    @test Suslin._verify_polynomial_quillen_patch_route_adapter(
        evidence.quillen_route_adapter,
    )
    @test evidence.local_evidence_provider.staged_diagnostic.status == :supported
    @test evidence.local_evidence_provider.murthy_adapter.mode ==
          :ordinary_quillen_factor_sequence
    @test evidence.quillen_consumption.patch == evidence.quillen_route_adapter.quillen_patch
    @test evidence.quillen_route_adapter.quillen_patch.replay_metadata.metadata.source ==
          :sl3_quillen_murthy_polynomial_route
    @test evidence.quillen_route_adapter.quillen_patch.replay_metadata.metadata.context_issue_id ==
          "#235"
    @test evidence.quillen_route_adapter.quillen_patch.replay_metadata.metadata.witness_issue_id ==
          "#236"
    @test evidence.quillen_route_adapter.quillen_patch.replay_metadata.metadata.provider_issue_id ==
          "#237"
    @test evidence.quillen_route_adapter.quillen_patch.replay_metadata.metadata.patch_issue_id ==
          "#220"
end
```

- [ ] **Step 2: Add the expert positive and tamper tests**

Append this block inside the existing `@testset "Park-Woodburn polynomial route certificates"` after the nonfixture Quillen assertions:

```julia
    issue238 = _issue238_sl3_route_case()
    sl3_cert = Suslin._polynomial_factorization_route_certificate(
        issue238.A;
        allow_recursive_column_peel = false,
    )
    _issue238_assert_sl3_route_evidence(sl3_cert, issue238.A)
    @test sl3_cert.factors == sl3_cert.evidence.quillen_route_adapter.global_elementary_factors
    @test sl3_cert.evidence.base_term_policy == :already_handled
    @test isempty(sl3_cert.evidence.base_term_factors)
    @test sl3_cert.evidence.quillen_route_adapter.quillen_patch.base_term_policy ==
          :already_handled

    bad_provider = _pw_rebuild(
        sl3_cert.evidence.local_evidence_provider;
        staged_diagnostic = merge(
            sl3_cert.evidence.local_evidence_provider.staged_diagnostic,
            (; status = :staged),
        ),
    )
    bad_provider_evidence = _pw_rebuild(
        sl3_cert.evidence;
        local_evidence_provider = bad_provider,
    )
    bad_provider_cert = _pw_replace_certificate(sl3_cert; evidence = bad_provider_evidence)
    @test !Suslin._verify_polynomial_factorization_route_certificate(bad_provider_cert)

    tampered_consumption = _pw_rebuild(
        sl3_cert.evidence.quillen_consumption;
        replay_metadata = (; source = :tampered_consumption),
    )
    bad_consumption_evidence = _pw_rebuild(
        sl3_cert.evidence;
        quillen_consumption = tampered_consumption,
    )
    bad_consumption_cert = _pw_replace_certificate(sl3_cert; evidence = bad_consumption_evidence)
    @test !Suslin._verify_polynomial_factorization_route_certificate(bad_consumption_cert)

    tampered_patch = _pw_rebuild(
        sl3_cert.evidence.quillen_route_adapter.quillen_patch;
        replay_metadata = (; source = :tampered_patch),
    )
    tampered_adapter = _pw_rebuild(
        sl3_cert.evidence.quillen_route_adapter;
        quillen_patch = tampered_patch,
    )
    bad_patch_evidence = _pw_rebuild(
        sl3_cert.evidence;
        quillen_route_adapter = tampered_adapter,
    )
    bad_patch_cert = _pw_replace_certificate(sl3_cert; evidence = bad_patch_evidence)
    @test !Suslin._verify_polynomial_factorization_route_certificate(bad_patch_cert)
```

- [ ] **Step 3: Add public driver acceptance and staged diagnostics**

Append this block inside `test/public/factorization_driver_shell.jl` before the existing `quillen_unsupported` check:

```julia
    issue238_R, (issue238_X, issue238_r, issue238_g) =
        Oscar.polynomial_ring(QQ, ["X", "r", "g"])
    issue238_p = issue238_X + issue238_r * issue238_g + one(issue238_R)
    issue238_A = matrix(issue238_R, [
        issue238_p one(issue238_R) zero(issue238_R);
        issue238_X + issue238_r * issue238_g one(issue238_R) zero(issue238_R);
        zero(issue238_R) zero(issue238_R) one(issue238_R)
    ])
    issue238_factors = elementary_factorization(issue238_A)
    @test verify_factorization(issue238_A, issue238_factors)
    issue238_cert = Suslin._polynomial_factorization_route_certificate(issue238_A)
    @test issue238_cert.route == :quillen_patch
    @test issue238_cert.evidence isa Suslin.PolynomialSL3QuillenMurthyRouteEvidence
    @test issue238_factors == issue238_cert.factors
    @test Suslin._verify_polynomial_factorization_route_certificate(issue238_cert)

    missing_witness = identity_matrix(issue238_R, 3)
    missing_witness_err = _captured_error(() -> elementary_factorization(missing_witness))
    @test missing_witness_err isa ArgumentError
    @test occursin("missing Quillen/local realizability witness", sprint(showerror, missing_witness_err))
    @test occursin("#236 local-form witness", sprint(showerror, missing_witness_err))

```

- [ ] **Step 4: Add public Park-Woodburn acceptance coverage**

Append this block inside `test/public/park_woodburn_polynomial_factorization.jl` before the negative controls:

```julia
    issue238_R, (issue238_X, issue238_r, issue238_g) =
        Oscar.polynomial_ring(QQ, ["X", "r", "g"])
    issue238_p = issue238_X + issue238_r * issue238_g + one(issue238_R)
    issue238_A = matrix(issue238_R, [
        issue238_p one(issue238_R) zero(issue238_R);
        issue238_X + issue238_r * issue238_g one(issue238_R) zero(issue238_R);
        zero(issue238_R) zero(issue238_R) one(issue238_R)
    ])
    issue238_factors, issue238_err = _pw_acceptance_result_or_error(issue238_A)
    @test issue238_err === nothing
    @test issue238_factors !== nothing
    @test verify_factorization(issue238_A, issue238_factors)
    issue238_cert = Suslin._polynomial_factorization_route_certificate(issue238_A)
    @test issue238_cert.route == :quillen_patch
    @test issue238_cert.evidence isa Suslin.PolynomialSL3QuillenMurthyRouteEvidence
    @test issue238_factors == issue238_cert.factors
    @test Suslin._verify_polynomial_factorization_route_certificate(issue238_cert)
```

- [ ] **Step 5: Run red tests**

If Julia reports missing packages, run `julia --project=. -e 'using Pkg; Pkg.instantiate()'` once, then rerun the focused commands.

Run:

```bash
julia --project=. -e 'include("test/expert/park_woodburn_route_certificate.jl")'
julia --project=. -e 'include("test/public/factorization_driver_shell.jl")'
julia --project=. -e 'include("test/public/park_woodburn_polynomial_factorization.jl")'
```

Expected: FAIL because `PolynomialSL3QuillenMurthyRouteEvidence` and the automatic SL3 route are not implemented.

### Task 2: Implement The SL3 Murthy/Quillen Route Evidence Wrapper

**Files:**
- Modify: `src/algorithm/factorization.jl`

**Interfaces:**
- Consumes: `SL3RealizationInputContext`, `SL3LocalFormWitnessSelection`, `SL3MurthyQuillenLocalEvidenceProvider`, `consume_murthy_quillen_adapters_for_patch`, `PolynomialQuillenPatchRouteAdapter`.
- Produces: `PolynomialSL3QuillenMurthyRouteEvidence`, `_verify_polynomial_sl3_quillen_murthy_route_evidence`, `_polynomial_sl3_quillen_murthy_route_certificate`.

- [ ] **Step 1: Add the evidence struct**

Add this after `PolynomialQuillenPatchRouteAdapter`:

```julia
struct PolynomialSL3QuillenMurthyRouteEvidence
    target
    route::Symbol
    context::SL3RealizationInputContext
    witness_selection::SL3LocalFormWitnessSelection
    local_evidence_provider::SL3MurthyQuillenLocalEvidenceProvider
    quillen_consumption
    quillen_route_adapter::PolynomialQuillenPatchRouteAdapter
    base_term_policy::Symbol
    base_term_factors::Vector
    replay_metadata
    verification
end
```

- [ ] **Step 2: Add route witness helpers**

Add these helpers near `_polynomial_quillen_supplied_evidence_route_certificate`:

```julia
function _polynomial_sl3_quillen_murthy_selected_variable(A)
    R = base_ring(A)
    entries = _sl3_local_target_entries(A)
    entries === nothing &&
        throw(ArgumentError("SL_3 Quillen/Murthy route is missing a #236 local-form witness: input is not an already-special-form SL_3 context"))
    for (index, generator) in enumerate(gens(R))
        if _is_monic_in_variable(entries.p, R, index)
            return generator, index, entries
        end
    end
    throw(ArgumentError("SL_3 Quillen/Murthy route is missing a #236 local-form witness: special-form p is not monic in any selected generator"))
end

function _polynomial_sl3_quillen_murthy_local_form_witness(A, selected_variable, entries)
    return (;
        entries,
        source_matrix = A,
        selected_variable,
        replay_steps = ((; kind = :issue238_already_special_form_replay),),
        witness_issue_id = "#236",
        route_issue_id = "#238",
    )
end

```

- [ ] **Step 3: Add construction fields and metadata**

Add:

```julia
function _polynomial_sl3_quillen_murthy_route_metadata(
    context,
    selection,
    provider,
    base_term_policy::Symbol,
    base_term_factors,
    metadata,
)
    return (;
        source = :sl3_quillen_murthy_polynomial_route,
        route_issue_id = "#238",
        context_issue_id = "#235",
        witness_issue_id = "#236",
        provider_issue_id = "#237",
        adapter_issue_id = "#219",
        patch_issue_id = "#220",
        context_metadata = context.catalog_metadata,
        witness_metadata = selection.local_form_witness,
        provider_replay_metadata = provider.replay_metadata,
        base_term_policy,
        base_term_factor_count = length(base_term_factors),
        metadata,
    )
end

function _polynomial_sl3_quillen_murthy_route_fields(A; metadata = (;))
    nrows(A) == 3 && ncols(A) == 3 ||
        throw(ArgumentError("SL_3 Quillen/Murthy route requires a 3 x 3 input"))
    R = base_ring(A)
    _factorization_ring_profile(R) == :polynomial ||
        throw(ArgumentError("SL_3 Quillen/Murthy route requires an ordinary polynomial ring"))
    _polynomial_exact_field_backed_ring(R) ||
        throw(ArgumentError(_polynomial_unsupported_coefficient_ring_message()))
    _require_polynomial_sl_determinant(A)
    selected_variable, selected_variable_index, entries =
        _polynomial_sl3_quillen_murthy_selected_variable(A)
    local_form_witness =
        _polynomial_sl3_quillen_murthy_local_form_witness(A, selected_variable, entries)
    context = _sl3_realization_input_context(
        A;
        selected_variable = (;
            name = string(selected_variable),
            generator = selected_variable,
            index = selected_variable_index,
            status = :passes,
        ),
        catalog_metadata = (;
            source = :sl3_quillen_murthy_polynomial_route,
            route_issue_id = "#238",
            context_issue_id = "#235",
            driver_issue_id = "#184",
        ),
        local_form_witness,
    )
    selection = _select_sl3_local_form_witness(context)
    provider = _sl3_murthy_quillen_local_evidence_provider(
        selection;
        metadata = (;
            source = :sl3_quillen_murthy_polynomial_route,
            route_issue_id = "#238",
            provider_issue_id = "#237",
        ),
    )
    provider.staged_diagnostic.status == :supported ||
        throw(ArgumentError("SL_3 Quillen/Murthy route is missing #237 ordinary Quillen local evidence: $(provider.staged_diagnostic.message)"))
    provider.murthy_adapter.mode == :ordinary_quillen_factor_sequence ||
        throw(ArgumentError("SL_3 Quillen/Murthy route requires ordinary Quillen local sequence evidence before #220 patching"))
    base_term_factors = typeof(A)[]
    base_term_policy = :already_handled
    route_metadata = _polynomial_sl3_quillen_murthy_route_metadata(
        context,
        selection,
        provider,
        base_term_policy,
        base_term_factors,
        metadata,
    )
    consumption = consume_murthy_quillen_adapters_for_patch(
        A,
        selected_variable,
        [provider.murthy_adapter];
        base_term_policy,
        base_term_factors,
        metadata = route_metadata,
    )
    verify_quillen_murthy_adapter_consumption(consumption) ||
        throw(ArgumentError("SL_3 Quillen/Murthy route #219/#220 adapter consumption and global patch evidence does not replay"))
    adapter = _polynomial_quillen_patch_route_adapter(A, consumption.patch)
    replay_metadata = (;
        source = :sl3_quillen_murthy_polynomial_route,
        route_issue_id = "#238",
        route_metadata,
        consumption_replay_metadata = consumption.replay_metadata,
        patch_replay_metadata = adapter.quillen_patch.replay_metadata,
    )
    return (;
        target = adapter.target_matrix,
        route = :quillen_patch,
        context,
        witness_selection = selection,
        local_evidence_provider = provider,
        quillen_consumption = consumption,
        quillen_route_adapter = adapter,
        base_term_policy,
        base_term_factors,
        replay_metadata,
    )
end
```

- [ ] **Step 4: Add evidence verification and constructor**

Add:

```julia
function _polynomial_sl3_quillen_murthy_route_core_verification(evidence)
    route_ok = evidence.route == :quillen_patch
    context_ok =
        _verify_sl3_realization_input_context(evidence.context) &&
        evidence.context.matrix == evidence.target
    selection_ok =
        context_ok &&
        _verify_sl3_local_form_witness_selection(evidence.witness_selection) &&
        evidence.witness_selection.context == evidence.context
    provider_ok =
        selection_ok &&
        _verify_sl3_murthy_quillen_local_evidence_provider(evidence.local_evidence_provider) &&
        evidence.local_evidence_provider.context == evidence.context &&
        evidence.local_evidence_provider.witness_selection == evidence.witness_selection &&
        evidence.local_evidence_provider.staged_diagnostic.status == :supported &&
        evidence.local_evidence_provider.murthy_adapter.mode == :ordinary_quillen_factor_sequence
    consumption_ok =
        provider_ok &&
        verify_quillen_murthy_adapter_consumption(evidence.quillen_consumption) &&
        evidence.quillen_consumption.original_input == evidence.target &&
        length(evidence.quillen_consumption.murthy_adapters) == 1 &&
        _same_murthy_quillen_local_adapter_data(
            first(evidence.quillen_consumption.murthy_adapters),
            evidence.local_evidence_provider.murthy_adapter,
        )
    adapter_ok =
        consumption_ok &&
        _verify_polynomial_quillen_patch_route_adapter(evidence.quillen_route_adapter) &&
        evidence.quillen_route_adapter.target_matrix == evidence.target &&
        evidence.quillen_route_adapter.quillen_patch == evidence.quillen_consumption.patch
    base_term_ok =
        adapter_ok &&
        evidence.quillen_route_adapter.quillen_patch.base_term_policy == evidence.base_term_policy &&
        evidence.base_term_policy == :already_handled &&
        isempty(evidence.base_term_factors)
    expected_metadata =
        adapter_ok ?
        (;
            source = :sl3_quillen_murthy_polynomial_route,
            route_issue_id = "#238",
            route_metadata = evidence.quillen_route_adapter.quillen_patch.replay_metadata.metadata,
            consumption_replay_metadata = evidence.quillen_consumption.replay_metadata,
            patch_replay_metadata = evidence.quillen_route_adapter.quillen_patch.replay_metadata,
        ) :
        nothing
    replay_metadata_ok = evidence.replay_metadata == expected_metadata
    overall_core_ok =
        route_ok &&
        context_ok &&
        selection_ok &&
        provider_ok &&
        consumption_ok &&
        adapter_ok &&
        base_term_ok &&
        replay_metadata_ok
    return (;
        route_ok,
        context_ok,
        selection_ok,
        provider_ok,
        consumption_ok,
        adapter_ok,
        base_term_ok,
        replay_metadata_ok,
        overall_core_ok,
    )
end

function _polynomial_sl3_quillen_murthy_route_verification(evidence)
    core = _polynomial_sl3_quillen_murthy_route_core_verification(evidence)
    stored_verification_ok = evidence.verification == core
    return merge(core, (; stored_verification_ok, overall_ok = core.overall_core_ok && stored_verification_ok))
end

function _verify_polynomial_sl3_quillen_murthy_route_evidence(evidence)::Bool
    try
        return _polynomial_sl3_quillen_murthy_route_verification(evidence).overall_ok
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end

function _polynomial_sl3_quillen_murthy_route_evidence(A; metadata = (;))
    fields = _polynomial_sl3_quillen_murthy_route_fields(A; metadata)
    raw = PolynomialSL3QuillenMurthyRouteEvidence(values(merge(fields, (; verification = nothing,)))...)
    verification = _polynomial_sl3_quillen_murthy_route_core_verification(raw)
    evidence = PolynomialSL3QuillenMurthyRouteEvidence(values(merge(fields, (; verification,)))...)
    _verify_polynomial_sl3_quillen_murthy_route_evidence(evidence) ||
        error("internal SL_3 Quillen/Murthy route evidence verification failed")
    return evidence
end

function _polynomial_sl3_quillen_murthy_route_certificate(A)
    evidence = _polynomial_sl3_quillen_murthy_route_evidence(A)
    adapter = evidence.quillen_route_adapter
    factors = copy(adapter.global_elementary_factors)
    return _polynomial_route_certificate(
        adapter.target_matrix,
        :quillen_patch,
        factors,
        adapter.product,
        evidence,
        :supported,
    )
end
```

- [ ] **Step 5: Route automatic certificate construction through the new helper**

In `_polynomial_factorization_route_certificate`, after the existing supplied-evidence Quillen route attempt and before recursive column peel, add:

```julia
        if n == 3
            try
                return _polynomial_sl3_quillen_murthy_route_certificate(A)
            catch err
                err isa InterruptException && rethrow()
                err isa ArgumentError || rethrow()
            end
        end
```

In `_polynomial_route_evidence_ok`, replace the `:quillen_patch` branch with:

```julia
        elseif cert.route == :quillen_patch
            if cert.evidence isa PolynomialQuillenPatchRouteAdapter
                return cert.evidence.target_matrix == cert.matrix &&
                    _verify_polynomial_quillen_patch_route_adapter(cert.evidence) &&
                    _polynomial_route_factor_sequences_equal(
                        cert.factors,
                        cert.evidence.global_elementary_factors,
                    )
            elseif cert.evidence isa PolynomialSL3QuillenMurthyRouteEvidence
                return cert.evidence.target == cert.matrix &&
                    _verify_polynomial_sl3_quillen_murthy_route_evidence(cert.evidence) &&
                    _polynomial_route_factor_sequences_equal(
                        cert.factors,
                        cert.evidence.quillen_route_adapter.global_elementary_factors,
                    )
            end
            return false
```

- [ ] **Step 6: Add staged diagnostics for failed SL3 route attempts**

Add:

```julia
function _polynomial_sl3_quillen_murthy_route_error(A)
    try
        _polynomial_sl3_quillen_murthy_route_certificate(A)
        return nothing
    catch err
        err isa InterruptException && rethrow()
        err isa ArgumentError || rethrow()
        return err
    end
end
```

In `_throw_staged_factorization_failure`, replace the multivariate `length(gens) > 1` throw with:

```julia
    if length(collect(gens(base_ring(A)))) > 1
        sl3_error = n == 3 ? _polynomial_sl3_quillen_murthy_route_error(A) : nothing
        detail = sl3_error === nothing ? "" : " ($(sprint(showerror, sl3_error)))"
        throw(ArgumentError("determinant-one polynomial input is outside the implemented Quillen/local evidence route: missing Quillen/local realizability witness$(detail)"))
    end
```

In `_polynomial_staged_failure_evidence`, before the `nrows(A) > 3` block, add:

```julia
    if nrows(A) == 3 && length(collect(gens(R))) > 1
        sl3_error = _polynomial_sl3_quillen_murthy_route_error(A)
        if sl3_error === nothing
            return (; error_type = :none, message = "")
        end
        return (;
            error_type = :ArgumentError,
            message = "determinant-one polynomial input is outside the implemented Quillen/local evidence route: missing Quillen/local realizability witness ($(sprint(showerror, sl3_error)))",
        )
    end
```

- [ ] **Step 7: Run focused green tests**

Run:

```bash
julia --project=. -e 'include("test/expert/park_woodburn_route_certificate.jl")'
julia --project=. -e 'include("test/public/factorization_driver_shell.jl")'
julia --project=. -e 'include("test/public/park_woodburn_polynomial_factorization.jl")'
```

Expected: all three exit 0.

### Task 3: Full Verification, Hygiene, Commit, And PR

**Files:**
- Modify: `src/algorithm/factorization.jl`
- Modify: `test/expert/park_woodburn_route_certificate.jl`
- Modify: `test/public/factorization_driver_shell.jl`
- Modify: `test/public/park_woodburn_polynomial_factorization.jl`
- Create: `docs/superpowers/plans/2026-07-01-issue-238-sl3-evidence-route.md`

**Interfaces:**
- Consumes: completed route helper and tests.
- Produces: verified branch and pull request.

- [ ] **Step 1: Run required issue verification**

Run:

```bash
julia --project=. -e 'include("test/expert/park_woodburn_route_certificate.jl")'
julia --project=. -e 'include("test/public/factorization_driver_shell.jl")'
julia --project=. -e 'include("test/public/park_woodburn_polynomial_factorization.jl")'
```

Expected: all commands exit 0.

- [ ] **Step 2: Run package verification**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: exit 0.

- [ ] **Step 3: Run diff hygiene**

Run:

```bash
git diff --check
git status --short
```

Expected: no whitespace errors. Only the intended source, test, and plan files are modified beyond the already committed design doc.

- [ ] **Step 4: Commit implementation**

Run:

```bash
git add docs/superpowers/plans/2026-07-01-issue-238-sl3-evidence-route.md src/algorithm/factorization.jl test/expert/park_woodburn_route_certificate.jl test/public/factorization_driver_shell.jl test/public/park_woodburn_polynomial_factorization.jl
git commit -m "Route SL3 evidence through polynomial factorization"
```

Expected: commit succeeds.

- [ ] **Step 5: Push and create the pull request**

Run:

```bash
git push -u origin agent/issue-238-route-evidence-backed-sl3-contexts-through-the-p-run-1
```

Then create a PR targeting `main` with title:

```text
Route SL3 evidence through polynomial factorization
```

Use this PR body:

```markdown
## Summary

- Route supported multivariate `SL_3` ordinary-polynomial inputs through checked #235 context, #236 witness, #237 Murthy/Quillen provider evidence, #219 adapter consumption, and #220 global patch verification before returning factors.
- Reuse the `:quillen_patch` route tag with a new internal evidence wrapper that verifies provider, consumption, patch adapter, base-term, and replay metadata.
- Add expert and public coverage for the evidence-backed route plus staged and tamper controls.

Closes #238.

## Verification

- `julia --project=. -e 'include("test/expert/park_woodburn_route_certificate.jl")'`
- `julia --project=. -e 'include("test/public/factorization_driver_shell.jl")'`
- `julia --project=. -e 'include("test/public/park_woodburn_polynomial_factorization.jl")'`
- `julia --project=. -e 'using Pkg; Pkg.test()'`
- `git diff --check`
```

Expected: PR URL is returned. Stop after PR creation.

## Plan Self-Review

- The plan maps every issue requirement to tests or implementation steps.
- The route returns factors only after provider, adapter consumption, patch, and route certificate verification.
- The existing fast-local and supplied-evidence Quillen paths remain intact.
- Negative controls cover determinant, coefficient ring, missing #236 witness, missing #237 ordinary local evidence, missing #220 global patch evidence through tampered patch/consumption records, tampered provider, and corrupted route certificate verification.
- No broad coordinate-change search, ECP, recursive `SL_n`, Laurent/ToricBuilder, or direct Murthy-factor return path is introduced.
