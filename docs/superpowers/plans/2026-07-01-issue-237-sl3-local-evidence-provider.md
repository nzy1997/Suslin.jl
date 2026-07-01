# Issue 237 SL3 Local Evidence Provider Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build an internal `SL_3` driver evidence provider that turns verified #236 witness selections into Murthy-to-Quillen local evidence or a staged adapter record.

**Architecture:** Add a thin provider in `src/algorithm/factorization.jl` beside the existing #235/#236 driver context records. The provider validates the #235 context and #236 selection, calls the existing Murthy input/certificate APIs, adapts through `_murthy_quillen_local_adapter`, and calls `quillen_local_sequences_from_murthy_adapters` only for ordinary materializable adapter mode.

**Tech Stack:** Julia, Oscar, existing Suslin internal Murthy and Quillen APIs, `Test`.

## Global Constraints

- The provider is internal and must not change public `elementary_factorization`.
- Do not reimplement Murthy branch logic, q-degree normalization, q(0)-unit checks, Bezout/resultant logic, denominator extraction, denominator-cover solving, #218 patch assembly, #219 adapter consumption, or Quillen patching.
- The provider input is a verified `SL3LocalFormWitnessSelection`; the output is verified Quillen-consumable local sequence evidence for ordinary replay or a verified staged provider/adapter record for localized replay.
- The provider must carry selected variable, local product, denominator/local-unit metadata, Murthy certificate provenance, #236 witness provenance, and #235 context metadata.
- Localized denominator-cleared Murthy replay must stage with a diagnostic instead of being handed to ordinary Quillen factor-sequence conversion.
- Verification commands required by the issue:
  `julia --project=. -e 'include("test/expert/park_woodburn_sl3_local_evidence_provider.jl")'`
  `julia --project=. -e 'include("test/expert/park_woodburn_quillen_route_adapter.jl")'`

---

### Task 1: Red Tests For Provider Evidence And Tamper Rejection

**Files:**
- Create: `test/expert/park_woodburn_sl3_local_evidence_provider.jl`
- Modify: `test/expert/park_woodburn_quillen_route_adapter.jl`
- Modify: `test/runtests.jl`

**Interfaces:**
- Consumes: existing `Suslin._sl3_realization_input_context`, `Suslin._select_sl3_local_form_witness`, `Suslin.quillen_local_sequences_from_murthy_adapters`.
- Produces: failing expectations for `Suslin._sl3_murthy_quillen_local_evidence_provider` and `Suslin._verify_sl3_murthy_quillen_local_evidence_provider`.

- [ ] **Step 1: Write the failing provider test**

Create `test/expert/park_woodburn_sl3_local_evidence_provider.jl` with a non-catalog multivariate `SL_3` special-form case:

```julia
using Test
using Oscar
using Suslin

function _issue237_rebuild(record; kwargs...)
    overrides = Dict{Symbol,Any}(pair.first => pair.second for pair in kwargs)
    values = map(fieldnames(typeof(record))) do name
        get(overrides, name, getproperty(record, name))
    end
    return typeof(record)(values...)
end

function _issue237_context_case()
    R, (X, u, v) = Oscar.polynomial_ring(QQ, ["X", "u", "v"])
    p = X + u * v + one(R)
    q = one(R)
    r = X + u * v
    s = one(R)
    A = matrix(R, [
        p q zero(R);
        r s zero(R);
        zero(R) zero(R) one(R)
    ])
    context_metadata = (;
        fixture_id = "issue-237-non-fixture-sl3-context",
        context_issue_id = "#235",
        driver_issue_id = "#184",
        original_matrix_label = :issue237_non_fixture,
    )
    witness_metadata = (;
        entries = (; p, q, r, s),
        source_matrix = A,
        selected_variable = X,
        replay_steps = ((; kind = :issue236_supplied_local_form),),
        witness_issue_id = "#236",
    )
    context = Suslin._sl3_realization_input_context(
        A;
        selected_variable = (; name = "X", generator = X, index = 1, status = :passes),
        catalog_metadata = context_metadata,
        local_form_witness = witness_metadata,
    )
    selection = Suslin._select_sl3_local_form_witness(context)
    return (; R, X, u, v, p, q, r, s, A, context, selection, context_metadata, witness_metadata)
end

@testset "Park-Woodburn SL3 Murthy-to-Quillen local evidence provider" begin
    case = _issue237_context_case()
    provider = Suslin._sl3_murthy_quillen_local_evidence_provider(
        case.selection;
        metadata = (; provider_test = :non_fixture, route_issue_id = "#184"),
    )

    @test provider.context == case.context
    @test provider.witness_selection == case.selection
    @test provider.selected_variable == case.X
    @test provider.selected_variable_index == 1
    @test provider.local_product == case.A
    @test provider.murthy_context.target == case.A
    @test Suslin.verify_sl3_local_murthy_input_context(provider.murthy_context)
    @test Suslin.verify_sl3_local_realization(provider.murthy_certificate)
    @test Suslin._verify_murthy_quillen_local_adapter(provider.murthy_adapter)
    @test provider.murthy_adapter.mode == :ordinary_quillen_factor_sequence
    @test provider.staged_diagnostic.status == :supported
    @test provider.denominator_metadata.denominator_product == one(case.R)
    @test provider.denominator_metadata.factor_denominators ==
          [one(case.R) for _ in provider.murthy_adapter.local_factor_replay.factors]
    @test provider.witness_metadata.context_metadata == case.context_metadata
    @test provider.witness_metadata.local_form_witness == case.witness_metadata
    @test provider.witness_metadata.witness_source == :already_special_form
    @test provider.replay_metadata.original_matrix == case.A
    @test provider.replay_metadata.local_product == case.A
    @test provider.replay_metadata.context_metadata == case.context_metadata
    @test length(provider.quillen_local_sequences) == 1
    @test all(Suslin.verify_quillen_local_factor_sequence_certificate, provider.quillen_local_sequences)
    @test provider.quillen_local_sequences[1].local_product == case.A
    @test provider.quillen_local_sequences[1].witness_metadata == provider.witness_metadata
    @test Suslin._verify_sl3_murthy_quillen_local_evidence_provider(provider)

    bad_selection = _issue237_rebuild(case.selection; selected_variable_index = 2)
    @test !Suslin._verify_sl3_local_form_witness_selection(bad_selection)
    @test_throws ArgumentError Suslin._sl3_murthy_quillen_local_evidence_provider(bad_selection)
    @test !Suslin._verify_sl3_murthy_quillen_local_evidence_provider(
        _issue237_rebuild(provider; witness_selection = bad_selection),
    )

    bad_certificate = Suslin.SL3LocalRealizationCertificate(
        provider.murthy_certificate.target,
        provider.murthy_certificate.branch,
        reverse(provider.murthy_certificate.factors),
        provider.murthy_certificate.selected_variable,
        provider.murthy_certificate.witness,
    )
    @test !Suslin._verify_sl3_murthy_quillen_local_evidence_provider(
        _issue237_rebuild(provider; murthy_certificate = bad_certificate),
    )

    @test !Suslin._verify_sl3_murthy_quillen_local_evidence_provider(
        _issue237_rebuild(provider; selected_variable = case.u),
    )

    bad_denominator_metadata = merge(
        provider.denominator_metadata,
        (; denominator_product = provider.denominator_metadata.denominator_product + case.X),
    )
    @test !Suslin._verify_sl3_murthy_quillen_local_evidence_provider(
        _issue237_rebuild(provider; denominator_metadata = bad_denominator_metadata),
    )

    bad_context = Suslin.SL3RealizationInputContext(
        values(merge(
            NamedTuple{propertynames(case.context)}(
                Tuple(getproperty(case.context, name) for name in propertynames(case.context)),
            ),
            (; catalog_metadata = (; fixture_id = "tampered-context"),),
        ))...,
    )
    @test !Suslin._verify_sl3_realization_input_context(bad_context)
    @test !Suslin._verify_sl3_murthy_quillen_local_evidence_provider(
        _issue237_rebuild(provider; context = bad_context),
    )
end
```

- [ ] **Step 2: Extend the route adapter expert test**

Append a new testset to `test/expert/park_woodburn_quillen_route_adapter.jl`:

```julia
@testset "SL3 provider handoff to Murthy adapter consumption" begin
    R, (X, u, v) = Oscar.polynomial_ring(QQ, ["X", "u", "v"])
    p = X + u * v + one(R)
    q = one(R)
    r = X + u * v
    s = one(R)
    A = matrix(R, [p q zero(R); r s zero(R); zero(R) zero(R) one(R)])
    context = Suslin._sl3_realization_input_context(
        A;
        selected_variable = (; name = "X", generator = X, index = 1, status = :passes),
        catalog_metadata = (; fixture_id = "issue-237-route-adapter-provider"),
        local_form_witness = (;
            entries = (; p, q, r, s),
            source_matrix = A,
            selected_variable = X,
            replay_steps = ((; kind = :issue236_route_adapter_replay),),
        ),
    )
    selection = Suslin._select_sl3_local_form_witness(context)
    provider = Suslin._sl3_murthy_quillen_local_evidence_provider(selection)
    sequences = Suslin.quillen_local_sequences_from_murthy_adapters(
        provider.local_product,
        provider.selected_variable,
        [provider.murthy_adapter],
    )
    @test provider.staged_diagnostic.status == :supported
    @test Suslin._same_quillen_local_factor_sequence_certificates(
        provider.quillen_local_sequences,
        sequences,
    )

    localized_fixture = Dict(entry.id => entry for entry in SL3MurthyGuptaFixtureCatalog.catalog().cases)[
        "mg-local-q0-unit-at-u"
    ]
    localized_R = base_ring(localized_fixture.target)
    localized_index = findfirst(isequal(localized_fixture.variable), collect(gens(localized_R)))
    localized_context = Suslin._sl3_realization_input_context(
        localized_fixture.target;
        selected_variable = (;
            name = string(localized_fixture.variable),
            generator = localized_fixture.variable,
            index = localized_index,
            status = :passes,
        ),
        catalog_metadata = (; fixture_id = "issue-237-localized-provider"),
        local_form_witness = (;
            entries = localized_fixture.entries,
            source_matrix = localized_fixture.target,
            selected_variable = localized_fixture.variable,
            replay_steps = ((; kind = :localized_fixture_replay),),
        ),
    )
    localized_selection = Suslin._select_sl3_local_form_witness(localized_context)
    localized_provider = Suslin._sl3_murthy_quillen_local_evidence_provider(
        localized_selection;
        witness = first(localized_fixture.witnesses),
    )
    @test localized_provider.murthy_adapter.mode == :localized_replay_handoff
    @test isempty(localized_provider.quillen_local_sequences)
    @test localized_provider.staged_diagnostic.status == :staged
    @test occursin("localized denominator-cleared", localized_provider.staged_diagnostic.message)
    @test Suslin._verify_sl3_murthy_quillen_local_evidence_provider(localized_provider)
    @test_throws ArgumentError Suslin.quillen_local_sequences_from_murthy_adapters(
        localized_provider.local_product,
        localized_provider.selected_variable,
        [localized_provider.murthy_adapter],
    )
end
```

- [ ] **Step 3: Register the new expert test**

Add `"expert/park_woodburn_sl3_local_evidence_provider.jl"` immediately after `"expert/park_woodburn_sl3_witness_selection.jl"` in `test/runtests.jl`.

- [ ] **Step 4: Run red tests**

Run:

```bash
julia --project=. -e 'include("test/expert/park_woodburn_sl3_local_evidence_provider.jl")'
julia --project=. -e 'include("test/expert/park_woodburn_quillen_route_adapter.jl")'
```

Expected: both fail because `_sl3_murthy_quillen_local_evidence_provider` is not defined yet.

### Task 2: Thin Provider Implementation

**Files:**
- Modify: `src/algorithm/factorization.jl`

**Interfaces:**
- Consumes: `SL3LocalFormWitnessSelection`, `sl3_local_murthy_input_context`, `realize_sl3_local_certificate`, `_murthy_quillen_local_adapter`, `quillen_local_sequences_from_murthy_adapters`.
- Produces: `SL3MurthyQuillenLocalEvidenceProvider`, `_sl3_murthy_quillen_local_evidence_provider`, `_verify_sl3_murthy_quillen_local_evidence_provider`.

- [ ] **Step 1: Add the provider struct**

Add this struct after `SL3LocalFormWitnessSelection`:

```julia
struct SL3MurthyQuillenLocalEvidenceProvider
    context::SL3RealizationInputContext
    witness_selection::SL3LocalFormWitnessSelection
    murthy_context
    murthy_certificate
    murthy_adapter
    quillen_local_sequences::Vector
    selected_variable
    selected_variable_index
    local_product
    denominator_metadata
    witness_metadata
    replay_metadata
    staged_diagnostic::NamedTuple
    metadata
    verification
end
```

- [ ] **Step 2: Add provider construction helpers**

Add helper functions in `factorization.jl` near the #236 witness-selection helpers. They must compute witness metadata, denominator metadata, replay metadata, staged diagnostic, and provider fields from the verified selection. The main fields helper must call:

```julia
murthy_context = sl3_local_murthy_input_context(
    selection.entries.p,
    selection.entries.q,
    selection.entries.r,
    selection.entries.s,
    selection.selected_variable;
    witness,
    local_unit_witnesses,
    split_witness,
    bezout_witness,
)
murthy_certificate = realize_sl3_local_certificate(murthy_context)
murthy_adapter = _murthy_quillen_local_adapter(
    murthy_certificate,
    selection.local_form_matrix,
    selection.selected_variable;
    witness_metadata,
)
quillen_local_sequences = murthy_adapter.mode == :ordinary_quillen_factor_sequence ?
    quillen_local_sequences_from_murthy_adapters(
        selection.local_form_matrix,
        selection.selected_variable,
        [murthy_adapter],
    ) :
    Any[]
```

- [ ] **Step 3: Add provider verification**

Add `_sl3_murthy_quillen_local_evidence_provider_core_verification(provider)` that recomputes the expected certificate, adapter, local sequences, metadata, diagnostic, and denominator data from the stored provider fields. The final `_verify_sl3_murthy_quillen_local_evidence_provider(provider)::Bool` must return `false` on malformed or tampered records and rethrow `InterruptException`.

- [ ] **Step 4: Add provider constructor**

Add:

```julia
function _sl3_murthy_quillen_local_evidence_provider(
    selection::SL3LocalFormWitnessSelection;
    witness = nothing,
    local_unit_witnesses = (;),
    split_witness = nothing,
    bezout_witness = nothing,
    metadata = (;),
)
```

The constructor must build an unchecked provider with `verification = nothing`, compute core verification, rebuild the checked provider, and throw an internal error if checked verification fails.

- [ ] **Step 5: Run green tests**

Run:

```bash
julia --project=. -e 'include("test/expert/park_woodburn_sl3_local_evidence_provider.jl")'
julia --project=. -e 'include("test/expert/park_woodburn_quillen_route_adapter.jl")'
```

Expected: both exit 0.

### Task 3: Full Verification And PR Preparation

**Files:**
- Modify: no additional source files unless verification exposes a defect.

**Interfaces:**
- Consumes: completed provider and tests.
- Produces: verified branch and PR.

- [ ] **Step 1: Run the required package test gate**

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: exit 0.

- [ ] **Step 2: Run expert file through the test runner**

```bash
julia --project=. test/runtests.jl expert
```

Expected: exit 0.

- [ ] **Step 3: Review git status**

```bash
git status --short
```

Expected: only intended source, test, and Superpowers doc files changed. `Manifest.toml` must remain ignored/uncommitted.
