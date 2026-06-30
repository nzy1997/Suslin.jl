# Issue 211 Murthy-Quillen Adapter Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an internal Quillen-facing adapter for verified Murthy local `SL_3` certificates.

**Architecture:** The adapter lives in `src/algorithm/quillen_induction.jl`, where both Murthy local certificate types and Quillen local certificate types are available. It accepts only verified `SL3LocalRealizationCertificate` values, extracts or rebuilds the Murthy local factor replay, materializes denominator-one ordinary factors into a `QuillenLocalFactorSequenceCertificate`, and keeps nontrivial localized replays as guarded handoff records for #183.

**Tech Stack:** Julia, Oscar, existing `SL3LocalRealizationCertificate`, `SL3LocalElementaryFactorReplay`, `QuillenLocalFactorSequenceCertificate`, and `QuillenLocalRealizationCertificate` helpers.

## Global Constraints

- Keep the adapter internal; do not export new names from `src/Suslin.jl`.
- Accept only `SL3LocalRealizationCertificate`; raw local factor vectors must not construct adapter data.
- Verify Murthy certificates before constructing Quillen-facing data.
- Materialize ordinary Quillen factors only when `SL3LocalElementaryFactorReplay.mode == :ordinary`.
- Localized `:denominator_cleared` Murthy replays must not pass through ordinary-factor-only Quillen APIs.
- Do not discover denominator covers, assemble global Quillen patches, implement a general `SL_3` driver, or modify public `elementary_factorization` route selection.
- Preserve existing #100 Quillen local certificate behavior.
- Required verification commands:
  - `julia --project=. -e 'include("test/expert/park_woodburn_quillen_route_adapter.jl")'`
  - `julia --project=. -e 'include("test/expert/quillen_local_certificate.jl")'`
  - `julia --project=. -e 'using Pkg; Pkg.test()'`

---

### Task 1: Red Tests For Murthy-Quillen Adapter

**Files:**
- Modify: `test/expert/park_woodburn_quillen_route_adapter.jl`
- Read: `test/fixtures/sl3_murthy_gupta_cases.jl`

**Interfaces:**
- Consumes: existing `SL3MurthyGuptaFixtureCatalog`, `realize_sl3_local_certificate`, `sl3_local_murthy_input_context`, and `verify_sl3_local_realization`.
- Produces: failing tests that define `_murthy_quillen_local_adapter`, `_verify_murthy_quillen_local_adapter`, `_murthy_quillen_local_factor_sequence_certificate`, and `_murthy_quillen_local_realization_certificate`.

- [ ] **Step 1: Add Murthy fixture include and helpers**

Add near the top of `test/expert/park_woodburn_quillen_route_adapter.jl`:

```julia
const SL3_MURTHY_QUILLEN_ADAPTER_FIXTURE_PATH =
    joinpath(@__DIR__, "..", "fixtures", "sl3_murthy_gupta_cases.jl")

if !isdefined(Main, :SL3MurthyGuptaFixtureCatalog)
    include(SL3_MURTHY_QUILLEN_ADAPTER_FIXTURE_PATH)
end

function _pwq_rebuild_sl3_certificate(
        cert;
        target = cert.target,
        branch = cert.branch,
        factors = cert.factors,
        selected_variable = cert.selected_variable,
        witness = cert.witness)
    return Suslin.SL3LocalRealizationCertificate(
        target,
        branch,
        factors,
        selected_variable,
        witness,
    )
end
```

- [ ] **Step 2: Add ordinary and localized adapter tests**

Append a new `@testset` to `test/expert/park_woodburn_quillen_route_adapter.jl`.
Leave the testset open; Step 3 adds the negative controls and closes it.

```julia
@testset "Murthy local SL3 Quillen handoff adapter" begin
    catalog = SL3MurthyGuptaFixtureCatalog.catalog()
    by_id = Dict(entry.id => entry for entry in catalog.cases)

    ordinary_fixture = by_id["mg-q0-unit-recursion"]
    ordinary_cert = Suslin.realize_sl3_local_certificate(
        ordinary_fixture.entries.p,
        ordinary_fixture.entries.q,
        ordinary_fixture.entries.r,
        ordinary_fixture.entries.s,
        ordinary_fixture.variable,
    )
    @test ordinary_cert isa Suslin.SL3LocalRealizationCertificate
    @test ordinary_cert.branch == :murthy_q0_unit
    @test Suslin.verify_sl3_local_realization(ordinary_cert)

    ordinary_adapter = Suslin._murthy_quillen_local_adapter(
        ordinary_cert,
        ordinary_fixture.target,
        ordinary_fixture.variable;
        witness_metadata = (;
            fixture_id = ordinary_fixture.id,
            consumer_issue = 211,
        ),
    )
    @test ordinary_adapter isa Suslin.MurthyQuillenLocalAdapter
    @test ordinary_adapter.mode == :ordinary_quillen_factor_sequence
    @test ordinary_adapter.selected_variable == ordinary_fixture.variable
    @test ordinary_adapter.local_product == ordinary_fixture.target
    @test ordinary_adapter.local_correction == ordinary_fixture.target
    @test ordinary_adapter.local_factor_replay.mode == :ordinary
    @test ordinary_adapter.local_factor_replay.materialized_factors == ordinary_cert.factors
    @test ordinary_adapter.witness_metadata.fixture_id == ordinary_fixture.id
    @test ordinary_adapter.replay_metadata.murthy_branch == ordinary_cert.branch
    @test ordinary_adapter.replay_metadata.denominator_product == one(base_ring(ordinary_fixture.target))
    @test Suslin._verify_murthy_quillen_local_adapter(ordinary_adapter)

    ordinary_sequence = Suslin._murthy_quillen_local_factor_sequence_certificate(ordinary_adapter)
    @test ordinary_sequence isa Suslin.QuillenLocalFactorSequenceCertificate
    @test Suslin.verify_quillen_local_factor_sequence_certificate(ordinary_sequence)
    @test ordinary_sequence.selected_variable == ordinary_fixture.variable
    @test ordinary_sequence.local_product == ordinary_fixture.target
    @test ordinary_sequence.local_correction == ordinary_fixture.target
    @test ordinary_sequence.witness_metadata.fixture_id == ordinary_fixture.id
    @test ordinary_sequence.verification.local_product == ordinary_fixture.target

    local_fixture = by_id["mg-local-q0-nonunit-bezout-at-u"]
    local_context = Suslin.sl3_local_murthy_input_context(
        local_fixture.target,
        local_fixture.variable;
        witness = first(local_fixture.witnesses),
    )
    localized_cert = Suslin.realize_sl3_local_certificate(local_context)
    @test localized_cert.branch == :murthy_q0_nonunit_bezout_resultant
    @test Suslin.verify_sl3_local_realization(localized_cert)

    localized_adapter = Suslin._murthy_quillen_local_adapter(
        localized_cert,
        local_fixture.target,
        local_fixture.variable;
        witness_metadata = (;
            fixture_id = local_fixture.id,
            consumer_issue = 211,
        ),
    )
    @test localized_adapter.mode == :localized_replay_handoff
    @test localized_adapter.local_factor_replay.mode == :denominator_cleared
    @test localized_adapter.quillen_factor_sequence === nothing
    @test localized_adapter.quillen_local_certificate === nothing
    @test localized_adapter.local_product === nothing
    @test localized_adapter.local_correction == local_fixture.target
    @test localized_adapter.replay_metadata.denominator_product != one(base_ring(local_fixture.target))
    @test localized_adapter.replay_metadata.cleared_product ==
          localized_adapter.local_factor_replay.cleared_product
    @test Suslin._verify_murthy_quillen_local_adapter(localized_adapter)
    @test_throws ArgumentError Suslin._murthy_quillen_local_factor_sequence_certificate(localized_adapter)
    @test_throws ArgumentError Suslin._murthy_quillen_local_realization_certificate(localized_adapter)
```

- [ ] **Step 3: Add negative controls**

Continue the same testset with:

```julia
    R = base_ring(ordinary_fixture.target)
    tampered_factors = copy(ordinary_cert.factors)
    tampered_factors[1] =
        tampered_factors[1] * elementary_matrix(3, 1, 3, one(R), R)
    tampered_cert = _pwq_rebuild_sl3_certificate(
        ordinary_cert;
        factors = tampered_factors,
    )
    @test !Suslin.verify_sl3_local_realization(tampered_cert)
    @test_throws ArgumentError Suslin._murthy_quillen_local_adapter(
        tampered_cert,
        ordinary_fixture.target,
        ordinary_fixture.variable,
    )

    local_generators = collect(gens(base_ring(local_fixture.target)))
    mismatched_variable = first(filter(gen -> gen != local_fixture.variable, local_generators))
    @test_throws ArgumentError Suslin._murthy_quillen_local_adapter(
        localized_cert,
        local_fixture.target,
        mismatched_variable,
    )

    @test_throws MethodError Suslin._murthy_quillen_local_adapter(
        ordinary_cert.factors,
        ordinary_fixture.target,
        ordinary_fixture.variable,
    )
end
```

- [ ] **Step 4: Run RED verification**

Run:

```bash
julia --project=. -e 'include("test/expert/park_woodburn_quillen_route_adapter.jl")'
```

Expected: FAIL because `MurthyQuillenLocalAdapter` or `_murthy_quillen_local_adapter` is not defined.

- [ ] **Step 5: Commit tests**

After RED failure is confirmed:

```bash
git add test/expert/park_woodburn_quillen_route_adapter.jl
git commit -m "test: cover murthy quillen adapter"
```

### Task 2: Implement Internal Murthy-Quillen Adapter

**Files:**
- Modify: `src/algorithm/quillen_induction.jl`

**Interfaces:**
- Consumes: tests from Task 1.
- Produces:
  - `MurthyQuillenLocalAdapter`
  - `_murthy_quillen_local_adapter(::SL3LocalRealizationCertificate, original_input, selected_variable; witness_metadata=(;), provenance=(; source = :murthy_quillen_local_adapter))`
  - `_verify_murthy_quillen_local_adapter(adapter)::Bool`
  - `_murthy_quillen_local_factor_sequence_certificate(adapter)`
  - `_murthy_quillen_local_realization_certificate(adapter)`

- [ ] **Step 1: Add adapter struct after Quillen local sequence structs**

Insert after `struct QuillenLocalDenominatorSupport`:

```julia
struct MurthyQuillenLocalAdapter
    original_input
    ring
    size::Int
    selected_variable
    murthy_certificate::SL3LocalRealizationCertificate
    local_factor_replay::SL3LocalElementaryFactorReplay
    mode::Symbol
    materialized_factors
    local_product
    local_correction
    quillen_factor_sequence
    quillen_local_certificate
    witness_metadata
    replay_metadata
    verification
end
```

- [ ] **Step 2: Add extraction and metadata helpers**

Add near the other Quillen local helper functions:

```julia
function _murthy_quillen_local_replay(certificate::SL3LocalRealizationCertificate)
    if certificate.branch == :murthy_q0_unit &&
            hasproperty(certificate.witness, :reduction) &&
            certificate.witness.reduction isa SL3LocalMurthyQUnitLocalReduction
        return certificate.witness.reduction.local_factor_replay
    elseif certificate.branch == :murthy_q0_nonunit_bezout_resultant &&
            hasproperty(certificate.witness, :reduction) &&
            hasproperty(certificate.witness.reduction, :local_factor_replay) &&
            certificate.witness.reduction.local_factor_replay !== nothing
        return certificate.witness.reduction.local_factor_replay
    elseif all(factor -> factor isa SL3LocalElementaryFactor, certificate.factors)
        return sl3_local_elementary_factor_replay(
            certificate.target,
            SL3LocalElementaryFactor[certificate.factors...],
            certificate.selected_variable,
        )
    else
        records = sl3_local_denominator_one_records_from_matrices(
            certificate.factors,
            certificate.selected_variable,
        )
        return sl3_local_elementary_factor_replay(
            certificate.target,
            records,
            certificate.selected_variable,
        )
    end
end

function _murthy_quillen_local_record_certificate(record::SL3LocalElementaryFactor)
    return LocalCertificate([record.row, record.col], [record.denominator, record.denominator])
end

function _murthy_quillen_local_sequence_factor(
        record::SL3LocalElementaryFactor,
        index::Int,
        witness_metadata,
)
    return QuillenLocalElementaryFactor(
        record.row,
        record.col,
        record.numerator,
        record.denominator,
        one(record.R),
        _murthy_quillen_local_record_certificate(record),
        (;
            source = :murthy_local_sl3,
            factor_index = index,
            murthy_denominator = record.denominator,
            murthy_local_unit_witness = record.local_unit_witness,
        ),
        (;
            source = :murthy_quillen_local_adapter,
            witness_metadata,
        ),
    )
end

function _murthy_quillen_local_replay_metadata(certificate, replay, mode, witness_metadata)
    return (;
        source = :murthy_quillen_local_adapter,
        murthy_branch = certificate.branch,
        replay_mode = replay.mode,
        adapter_mode = mode,
        selected_variable = replay.selected_variable,
        factor_count = length(replay.factors),
        denominator_product = replay.denominator_product,
        cleared_product = replay.cleared_product,
        ordinary_materialized = replay.materialized_factors !== nothing,
        witness_metadata,
    )
end
```

- [ ] **Step 3: Add ordinary Quillen construction helpers**

Add:

```julia
function _murthy_quillen_local_factor_sequence(
        original_input,
        selected_variable,
        replay::SL3LocalElementaryFactorReplay,
        witness_metadata,
        provenance,
)
    replay.mode == :ordinary ||
        throw(ArgumentError("Murthy local replay is not materializable over the ordinary base ring"))
    factors = [
        _murthy_quillen_local_sequence_factor(record, index, witness_metadata)
        for (index, record) in enumerate(replay.factors)
    ]
    return quillen_local_factor_sequence_certificate(
        original_input,
        selected_variable;
        factors,
        local_correction = replay.target,
        witness_metadata,
        local_evidence = (;
            source = :murthy_local_sl3,
            replay_mode = replay.mode,
            denominator_product = replay.denominator_product,
            expected_product = replay.target,
        ),
        provenance,
    )
end

function _murthy_quillen_local_single_realization(
        original_input,
        selected_variable,
        replay::SL3LocalElementaryFactorReplay,
        witness_metadata,
)
    replay.mode == :ordinary || return nothing
    length(replay.factors) == 1 || return nothing
    record = only(replay.factors)
    return quillen_local_realization_certificate(
        original_input,
        selected_variable;
        local_certificate = _murthy_quillen_local_record_certificate(record),
        denominator = record.denominator,
        coverage_multiplier = one(record.R),
        correction = QuillenElementaryCorrection(record.row, record.col, record.numerator),
        factors = replay.materialized_factors,
        local_correction = replay.target,
        witness_metadata,
    )
end
```

- [ ] **Step 4: Add constructor and verifier**

Add:

```julia
function _murthy_quillen_local_adapter(
        certificate::SL3LocalRealizationCertificate,
        original_input,
        selected_variable;
        witness_metadata = (;),
        provenance = (; source = :murthy_quillen_local_adapter),
)
    verify_sl3_local_realization(certificate) ||
        throw(ArgumentError("Murthy local certificate does not replay"))
    R, n = _quillen_local_input_ring_size(original_input)
    n == 3 || throw(ArgumentError("Murthy local Quillen adapter requires a 3x3 input"))
    _same_base_ring(R, base_ring(certificate.target)) ||
        throw(ArgumentError("Murthy local Quillen adapter ring mismatch"))
    original_input == certificate.target ||
        throw(ArgumentError("Murthy local Quillen adapter requires the original input to match the Murthy target"))
    selected = _require_substitution_generator(R, selected_variable)
    selected == certificate.selected_variable ||
        throw(ArgumentError("Murthy local Quillen adapter selected variable mismatch"))

    replay = _murthy_quillen_local_replay(certificate)
    verify_sl3_local_elementary_factor_replay(replay) ||
        throw(ArgumentError("Murthy local factor replay does not verify"))
    replay.target == certificate.target ||
        throw(ArgumentError("Murthy local replay target mismatch"))
    replay.selected_variable == selected ||
        throw(ArgumentError("Murthy local replay selected variable mismatch"))
    replay.factors == certificate.factors ||
        throw(ArgumentError("Murthy local replay factor mismatch"))

    mode = replay.mode == :ordinary ?
        :ordinary_quillen_factor_sequence :
        :localized_replay_handoff
    quillen_sequence = mode == :ordinary_quillen_factor_sequence ?
        _murthy_quillen_local_factor_sequence(
            original_input,
            selected,
            replay,
            witness_metadata,
            provenance,
        ) :
        nothing
    quillen_local = mode == :ordinary_quillen_factor_sequence ?
        _murthy_quillen_local_single_realization(
            original_input,
            selected,
            replay,
            witness_metadata,
        ) :
        nothing
    local_product = quillen_sequence === nothing ? nothing : quillen_sequence.local_product
    local_correction = quillen_sequence === nothing ? replay.target : quillen_sequence.local_correction
    replay_metadata = _murthy_quillen_local_replay_metadata(
        certificate,
        replay,
        mode,
        witness_metadata,
    )
    provisional = MurthyQuillenLocalAdapter(
        original_input,
        R,
        n,
        selected,
        certificate,
        replay,
        mode,
        replay.materialized_factors,
        local_product,
        local_correction,
        quillen_sequence,
        quillen_local,
        witness_metadata,
        replay_metadata,
        nothing,
    )
    verification = _murthy_quillen_local_adapter_summary(provisional)
    verification.overall_ok ||
        throw(ArgumentError("Murthy local Quillen adapter data does not replay"))
    return MurthyQuillenLocalAdapter(
        provisional.original_input,
        provisional.ring,
        provisional.size,
        provisional.selected_variable,
        provisional.murthy_certificate,
        provisional.local_factor_replay,
        provisional.mode,
        provisional.materialized_factors,
        provisional.local_product,
        provisional.local_correction,
        provisional.quillen_factor_sequence,
        provisional.quillen_local_certificate,
        provisional.witness_metadata,
        provisional.replay_metadata,
        verification,
    )
end

function _murthy_quillen_local_adapter_summary(adapter::MurthyQuillenLocalAdapter)
    certificate_ok = verify_sl3_local_realization(adapter.murthy_certificate)
    replay_ok = verify_sl3_local_elementary_factor_replay(adapter.local_factor_replay)
    input_ok = adapter.original_input == adapter.murthy_certificate.target
    selected_variable_ok = adapter.selected_variable == adapter.murthy_certificate.selected_variable
    replay_alignment_ok =
        adapter.local_factor_replay.target == adapter.murthy_certificate.target &&
        adapter.local_factor_replay.selected_variable == adapter.selected_variable &&
        adapter.local_factor_replay.factors == adapter.murthy_certificate.factors
    expected_mode = adapter.local_factor_replay.mode == :ordinary ?
        :ordinary_quillen_factor_sequence :
        :localized_replay_handoff
    mode_ok = adapter.mode == expected_mode
    materialized_ok = adapter.materialized_factors == adapter.local_factor_replay.materialized_factors
    sequence_ok =
        adapter.mode == :ordinary_quillen_factor_sequence ?
        adapter.quillen_factor_sequence isa QuillenLocalFactorSequenceCertificate &&
            verify_quillen_local_factor_sequence_certificate(adapter.quillen_factor_sequence) &&
            adapter.quillen_factor_sequence.original_input == adapter.original_input &&
            adapter.quillen_factor_sequence.selected_variable == adapter.selected_variable &&
            adapter.quillen_factor_sequence.local_product == adapter.local_factor_replay.target &&
            adapter.quillen_factor_sequence.local_correction == adapter.local_factor_replay.target :
        adapter.quillen_factor_sequence === nothing
    local_certificate_ok =
        adapter.quillen_local_certificate === nothing ||
        (
            adapter.quillen_local_certificate isa QuillenLocalRealizationCertificate &&
            verify_quillen_local_certificate(adapter.quillen_local_certificate)
        )
    product_ok =
        adapter.mode == :ordinary_quillen_factor_sequence ?
        adapter.local_product == adapter.local_factor_replay.target :
        adapter.local_product === nothing
    correction_ok = adapter.local_correction == adapter.local_factor_replay.target
    expected_metadata = _murthy_quillen_local_replay_metadata(
        adapter.murthy_certificate,
        adapter.local_factor_replay,
        adapter.mode,
        adapter.witness_metadata,
    )
    replay_metadata_ok = adapter.replay_metadata == expected_metadata
    overall_ok =
        certificate_ok &&
        replay_ok &&
        input_ok &&
        selected_variable_ok &&
        replay_alignment_ok &&
        mode_ok &&
        materialized_ok &&
        sequence_ok &&
        local_certificate_ok &&
        product_ok &&
        correction_ok &&
        replay_metadata_ok
    return (;
        certificate_ok,
        replay_ok,
        input_ok,
        selected_variable_ok,
        replay_alignment_ok,
        mode_ok,
        materialized_ok,
        sequence_ok,
        local_certificate_ok,
        product_ok,
        correction_ok,
        replay_metadata_ok,
        overall_ok,
    )
end

function _verify_murthy_quillen_local_adapter(adapter)::Bool
    try
        replay = _murthy_quillen_local_adapter_summary(adapter)
        return replay.overall_ok && adapter.verification == replay
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end
```

- [ ] **Step 5: Add ordinary-only accessors**

Add:

```julia
function _murthy_quillen_local_factor_sequence_certificate(
        adapter::MurthyQuillenLocalAdapter,
)
    _verify_murthy_quillen_local_adapter(adapter) ||
        throw(ArgumentError("Murthy local Quillen adapter does not replay"))
    adapter.quillen_factor_sequence !== nothing ||
        throw(ArgumentError("Murthy local adapter contains localized denominator-cleared replay; #183 must define the localized Quillen local certificate shape before ordinary factor conversion"))
    return adapter.quillen_factor_sequence
end

function _murthy_quillen_local_realization_certificate(
        adapter::MurthyQuillenLocalAdapter,
)
    _verify_murthy_quillen_local_adapter(adapter) ||
        throw(ArgumentError("Murthy local Quillen adapter does not replay"))
    adapter.quillen_local_certificate !== nothing ||
        throw(ArgumentError("Murthy local adapter does not contain a length-one ordinary Quillen local realization certificate"))
    return adapter.quillen_local_certificate
end
```

- [ ] **Step 6: Run GREEN verification**

Run:

```bash
julia --project=. -e 'include("test/expert/park_woodburn_quillen_route_adapter.jl")'
julia --project=. -e 'include("test/expert/quillen_local_certificate.jl")'
```

Expected: both pass.

- [ ] **Step 7: Commit implementation**

```bash
git add src/algorithm/quillen_induction.jl
git commit -m "feat: add murthy quillen local adapter"
```

### Task 3: Final Verification And Branch Review

**Files:**
- No planned edits unless verification finds a defect.

**Interfaces:**
- Consumes: committed tests and implementation from Tasks 1 and 2.
- Produces: verified branch ready for PR.

- [ ] **Step 1: Run required expert verification**

Run:

```bash
julia --project=. -e 'include("test/expert/park_woodburn_quillen_route_adapter.jl")'
julia --project=. -e 'include("test/expert/quillen_local_certificate.jl")'
```

Expected: both exit 0.

- [ ] **Step 2: Run package tests**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: exits 0.

- [ ] **Step 3: Run diff hygiene**

Run:

```bash
git diff --check
```

Expected: exits 0.

- [ ] **Step 4: Commit any verification fixes**

If verification required code or test fixes, commit them with a targeted message:

```bash
git add src/algorithm/quillen_induction.jl test/expert/park_woodburn_quillen_route_adapter.jl
git commit -m "fix: harden murthy quillen adapter replay"
```
