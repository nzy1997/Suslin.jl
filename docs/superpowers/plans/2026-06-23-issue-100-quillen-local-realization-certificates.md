# Issue 100 Quillen Local Realization Certificates Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add replayable Quillen local realization certificates that consume supplied local fixture data and reject tampering.

**Architecture:** Keep the new API expert/internal by defining unexported certificate names in `src/algorithm/quillen_induction.jl`. The constructor normalizes existing `LocalCertificate` and `QuillenElementaryCorrection` data through current Quillen helpers, records supplied factors and optional patched-substitution witnesses, and stores a replay summary that the verifier recomputes exactly.

**Tech Stack:** Julia, Oscar exact polynomial rings, Suslin Quillen helpers, Test stdlib, Issue 99 fixture catalog.

## Global Constraints

- Preserve the existing `LocalCertificate` constructor.
- Do not choose denominator covers.
- Do not assemble global factors.
- Do not call the local `SL_3` solver automatically.
- Consume supplied local certificate data from deterministic fixtures.
- Keep backend-sensitive shared algebra helpers in `src/core/groebner_tools.jl` only when shared outside Quillen patching.
- New certificate names remain unexported; do not modify `src/Suslin.jl` exports or `test/public/api_surface.jl`.
- Tests use qualified `Suslin.<name>` access for new expert/internal names.
- Carry forward patched-substitution witness fields from #99: `matrix`, `variable`, `denominator`, `exponent`, `shift`, and `expected_matrix`.
- Focused command is `julia --project=. -e 'include("test/expert/quillen_local_certificate.jl")'`.
- Expert group command is `julia --project=. test/runtests.jl expert`.
- Full package command is `julia --project=. -e 'using Pkg; Pkg.test()'`.

---

## File Structure

- Modify `src/algorithm/quillen_induction.jl`: owns `QuillenLocalRealizationCertificate`, constructor, replay helper, patched-substitution witness replay, and `verify_quillen_local_certificate`.
- Create `test/expert/quillen_local_certificate.jl`: owns fixture-backed positive tests and tamper controls.
- Modify `test/runtests.jl`: registers `expert/quillen_local_certificate.jl` in the expert group near existing Quillen expert tests.

---

### Task 1: Expert Test And RED Verification

**Files:**
- Create: `test/expert/quillen_local_certificate.jl`

**Interfaces:**
- Consumes: `test/fixtures/quillen_patch_cases.jl`.
- Produces: a failing test that specifies `Suslin.QuillenLocalRealizationCertificate`, `Suslin.quillen_local_realization_certificate`, and `Suslin.verify_quillen_local_certificate`.

- [ ] **Step 1: Create fixture-backed test helpers**

Create `test/expert/quillen_local_certificate.jl` with imports, fixture include, and helper functions:

```julia
using Test
using Suslin
using Oscar

include(joinpath(@__DIR__, "..", "fixtures", "quillen_patch_cases.jl"))

function quillen_local_test_product(factors, R, n)
    product = identity_matrix(R, n)
    for factor in factors
        product *= factor
    end
    return product
end

function local_certificate_from_fixture(local_factor)
    return Suslin.LocalCertificate(
        local_factor.certificate.indices,
        local_factor.certificate.denominators,
    )
end

function local_correction_from_fixture(local_factor)
    correction = local_factor.correction
    return Suslin.QuillenElementaryCorrection(correction.row, correction.col, correction.entry)
end

function realization_certificate_from_fixture(entry; local_index::Int=1)
    local_factor = entry.local_factors[local_index]
    return Suslin.quillen_local_realization_certificate(
        entry.target_matrix,
        entry.substitution_variable;
        local_certificate = local_certificate_from_fixture(local_factor),
        denominator = local_factor.denominator,
        coverage_multiplier = local_factor.coverage_multiplier,
        correction = local_correction_from_fixture(local_factor),
        factors = [local_factor.factor],
        patched_substitution_witness = entry.patched_substitution_witness,
        witness_metadata = (;
            fixture_id = entry.id,
            local_index = local_index,
            source_refs = entry.source_refs,
            consumer_issue_ids = entry.consumer_issue_ids,
        ),
    )
end

function rebuild_local_certificate(cert; kwargs...)
    fields = merge((;
        original_input = cert.original_input,
        ring = cert.ring,
        size = cert.size,
        selected_variable = cert.selected_variable,
        local_certificate = cert.local_certificate,
        denominator = cert.denominator,
        coverage_multiplier = cert.coverage_multiplier,
        correction = cert.correction,
        factors = cert.factors,
        local_product = cert.local_product,
        local_correction = cert.local_correction,
        patched_substitution_witness = cert.patched_substitution_witness,
        witness_metadata = cert.witness_metadata,
        verification = cert.verification,
    ), kwargs)
    return Suslin.QuillenLocalRealizationCertificate(
        fields.original_input,
        fields.ring,
        fields.size,
        fields.selected_variable,
        fields.local_certificate,
        fields.denominator,
        fields.coverage_multiplier,
        fields.correction,
        fields.factors,
        fields.local_product,
        fields.local_correction,
        fields.patched_substitution_witness,
        fields.witness_metadata,
        fields.verification,
    )
end
```

- [ ] **Step 2: Add positive fixture tests**

Add a testset that builds certificates for at least these Issue 99 fixture ids:

```julia
@testset "Quillen local realization certificates" begin
    entries = QuillenPatchFixtureCatalog.cases_by_id()
    fixture_ids = [
        "quillen-two-open-cover-qq",
        "quillen-supplied-local-certificate-gf2",
        "quillen-patched-substitution-witness-qq",
    ]

    certs = [realization_certificate_from_fixture(entries[id]) for id in fixture_ids]
    @test all(cert -> cert isa Suslin.QuillenLocalRealizationCertificate, certs)

    for cert in certs
        @test Suslin.verify_quillen_local_certificate(cert)
        @test quillen_local_test_product(cert.factors, cert.ring, cert.size) == cert.local_correction
        @test cert.local_product == cert.local_correction
        @test cert.verification.local_product == cert.local_correction
        @test cert.verification.local_correction_ok
        @test cert.verification.factors_ok
        @test cert.verification.overall_ok
        @test Suslin._coerce_into_ring(cert.ring, cert.denominator, "test denominator") == cert.denominator
        for denominator in cert.local_certificate.denominators
            @test Suslin._coerce_into_ring(cert.ring, denominator, "test certificate denominator") == denominator
        end
    end

    patched_cert = certs[3]
    @test patched_cert.patched_substitution_witness !== nothing
    @test patched_cert.verification.patched_substitution_ok
    witness = patched_cert.patched_substitution_witness
    @test Suslin.patched_substitution(
        witness.matrix,
        witness.variable,
        witness.denominator,
        witness.exponent,
        witness.shift,
    ) == witness.expected_matrix
end
```

- [ ] **Step 3: Add tamper controls**

In the same testset, add separate negative checks:

```julia
base_cert = certs[1]
R = base_cert.ring
bad_factor = base_cert.factors[1] * elementary_matrix(base_cert.size, 1, 3, one(R), R)
@test !Suslin.verify_quillen_local_certificate(rebuild_local_certificate(base_cert; factors = [bad_factor]))

@test !Suslin.verify_quillen_local_certificate(rebuild_local_certificate(
    base_cert;
    denominator = base_cert.denominator + one(R),
))

generators = collect(gens(R))
replacement_variable = first(filter(gen -> gen != base_cert.selected_variable, generators))
@test !Suslin.verify_quillen_local_certificate(rebuild_local_certificate(
    base_cert;
    selected_variable = replacement_variable,
))

patched_cert = certs[3]
bad_witness = merge(patched_cert.patched_substitution_witness, (;
    exponent = patched_cert.patched_substitution_witness.exponent + 1,
))
@test !Suslin.verify_quillen_local_certificate(rebuild_local_certificate(
    patched_cert;
    patched_substitution_witness = bad_witness,
))
```

- [ ] **Step 4: Run RED verification**

Run:

```bash
julia --project=. -e 'include("test/expert/quillen_local_certificate.jl")'
```

Expected: FAIL with `UndefVarError` or `is not defined` for the new certificate constructor/type.

### Task 2: Quillen Local Certificate Replay Implementation

**Files:**
- Modify: `src/algorithm/quillen_induction.jl`

**Interfaces:**
- Produces: `QuillenLocalRealizationCertificate`.
- Produces: `quillen_local_realization_certificate(original_input, selected_variable; local_certificate, denominator, coverage_multiplier, correction, factors=nothing, local_correction=nothing, patched_substitution_witness=nothing, witness_metadata=(;), ring=nothing, size=nothing)`.
- Produces: `verify_quillen_local_certificate(certificate)::Bool`.

- [ ] **Step 1: Add certificate struct**

Append these definitions after `patched_substitution`:

```julia
struct QuillenLocalRealizationCertificate
    original_input
    ring
    size::Int
    selected_variable
    local_certificate::LocalCertificate
    denominator
    coverage_multiplier
    correction::QuillenElementaryCorrection
    factors::Vector
    local_product
    local_correction
    patched_substitution_witness
    witness_metadata
    verification
end
```

- [ ] **Step 2: Add input normalization helpers**

Add helpers that derive a target ring and size from a matrix input or from an
elementary correction plus explicit `ring` and `size` keywords:

```julia
function _quillen_local_input_ring_size(original_input; ring = nothing, size = nothing)
    if original_input isa QuillenElementaryCorrection
        ring === nothing && throw(ArgumentError("ring is required when original input is an elementary correction"))
        size === nothing && throw(ArgumentError("size is required when original input is an elementary correction"))
        R = _require_supported_quillen_ring(ring)
        size isa Integer && size >= 2 || throw(ArgumentError("certificate size must be at least 2"))
        return R, Int(size)
    end

    matrix_size = _require_square_matrix(original_input, "original input")
    R = _require_supported_quillen_ring(base_ring(original_input))
    return R, matrix_size
end

function _quillen_local_require_factor_matrix(factor, R, n::Int, label::AbstractString)
    nrows(factor) == n && ncols(factor) == n ||
        throw(ArgumentError("$(label) must be a square matrix of certificate size"))
    base_ring(factor) == R ||
        throw(ArgumentError("$(label) must be defined over the certificate ring"))
    return factor
end
```

- [ ] **Step 3: Add factor and witness replay helpers**

Add helpers:

```julia
function _quillen_local_factor_vector(factors, local_correction, R, n::Int)
    if factors === nothing
        local_correction === nothing &&
            throw(ArgumentError("either factors or local_correction must be supplied"))
        return [_quillen_local_require_factor_matrix(local_correction, R, n, "local correction")]
    end
    collected = collect(factors)
    isempty(collected) && throw(ArgumentError("local factors must be nonempty"))
    return [
        _quillen_local_require_factor_matrix(factor, R, n, "local factor")
        for factor in collected
    ]
end

function _quillen_local_patched_witness_summary(witness, R, n::Int, selected_variable)
    witness === nothing && return (; present = false, ok = true)
    for field in (:matrix, :variable, :denominator, :exponent, :shift, :expected_matrix)
        hasproperty(witness, field) ||
            throw(ArgumentError("patched-substitution witness missing field $(field)"))
    end
    _quillen_local_require_factor_matrix(witness.matrix, R, n, "patched-substitution witness matrix")
    _quillen_local_require_factor_matrix(witness.expected_matrix, R, n, "patched-substitution witness expected matrix")
    witness_variable = _require_substitution_generator(R, witness.variable)
    denominator = _coerce_into_ring(R, witness.denominator, "patched-substitution witness denominator")
    shift = _coerce_into_ring(R, witness.shift, "patched-substitution witness shift")
    actual = patched_substitution(witness.matrix, witness_variable, denominator, witness.exponent, shift)
    ok = witness_variable == selected_variable && actual == witness.expected_matrix
    return (;
        present = true,
        variable = witness_variable,
        denominator = denominator,
        exponent = witness.exponent,
        shift = shift,
        actual_matrix = actual,
        expected_matrix = witness.expected_matrix,
        ok = ok,
    )
end
```

- [ ] **Step 4: Add replay summary and constructor**

Implement replay by reusing existing Quillen contribution normalization:

```julia
function _quillen_local_certificate_replay_summary(certificate::QuillenLocalRealizationCertificate)
    R = _require_supported_quillen_ring(certificate.ring)
    n = certificate.size
    n >= 2 || throw(ArgumentError("certificate size must be at least 2"))
    selected_variable = _require_substitution_generator(R, certificate.selected_variable)

    if certificate.original_input isa QuillenElementaryCorrection
        certificate.original_input == certificate.correction ||
            throw(ArgumentError("original elementary correction must match recorded correction"))
    else
        _quillen_local_require_factor_matrix(certificate.original_input, R, n, "original input")
    end

    normalized = _normalize_quillen_contribution(
        QuillenLocalContribution(
            certificate.local_certificate,
            certificate.denominator,
            certificate.coverage_multiplier,
            certificate.correction,
        ),
        R,
        n,
    )
    expected_factor = _quillen_factors(R, n, [normalized])[1]
    factors = [
        _quillen_local_require_factor_matrix(factor, R, n, "local factor")
        for factor in certificate.factors
    ]
    local_product = _quillen_product(R, n, factors)
    local_correction = _quillen_local_require_factor_matrix(
        certificate.local_correction,
        R,
        n,
        "local correction",
    )
    witness = _quillen_local_patched_witness_summary(
        certificate.patched_substitution_witness,
        R,
        n,
        selected_variable,
    )
    denominator_ok = normalized.denominator == certificate.denominator
    correction_ok = expected_factor == local_correction
    factors_ok = local_product == local_correction
    stored_product_ok = certificate.local_product == local_product
    overall_ok = denominator_ok && correction_ok && factors_ok && stored_product_ok && witness.ok
    return (;
        selected_variable = selected_variable,
        denominator = normalized.denominator,
        coverage_multiplier = normalized.coverage_multiplier,
        expected_local_correction = expected_factor,
        local_product = local_product,
        local_correction = local_correction,
        denominator_ok = denominator_ok,
        local_correction_ok = correction_ok,
        factors_ok = factors_ok,
        stored_product_ok = stored_product_ok,
        patched_substitution = witness,
        patched_substitution_ok = witness.ok,
        witness_metadata = certificate.witness_metadata,
        overall_ok = overall_ok,
    )
end

function quillen_local_realization_certificate(
    original_input,
    selected_variable;
    local_certificate::LocalCertificate,
    denominator,
    coverage_multiplier,
    correction::QuillenElementaryCorrection,
    factors = nothing,
    local_correction = nothing,
    patched_substitution_witness = nothing,
    witness_metadata = (;),
    ring = nothing,
    size = nothing,
)
    R, n = _quillen_local_input_ring_size(original_input; ring, size)
    selected = _require_substitution_generator(R, selected_variable)
    normalized = _normalize_quillen_contribution(
        QuillenLocalContribution(local_certificate, denominator, coverage_multiplier, correction),
        R,
        n,
    )
    normalized_factors = _quillen_local_factor_vector(factors, local_correction, R, n)
    product = _quillen_product(R, n, normalized_factors)
    recorded_correction = local_correction === nothing ?
        product :
        _quillen_local_require_factor_matrix(local_correction, R, n, "local correction")
    provisional = QuillenLocalRealizationCertificate(
        original_input,
        R,
        n,
        selected,
        normalized.certificate,
        normalized.denominator,
        normalized.coverage_multiplier,
        normalized.correction,
        normalized_factors,
        product,
        recorded_correction,
        patched_substitution_witness,
        witness_metadata,
        nothing,
    )
    verification = _quillen_local_certificate_replay_summary(provisional)
    verification.overall_ok ||
        throw(ArgumentError("Quillen local realization certificate data does not replay"))
    return QuillenLocalRealizationCertificate(
        provisional.original_input,
        provisional.ring,
        provisional.size,
        provisional.selected_variable,
        provisional.local_certificate,
        provisional.denominator,
        provisional.coverage_multiplier,
        provisional.correction,
        provisional.factors,
        provisional.local_product,
        provisional.local_correction,
        provisional.patched_substitution_witness,
        provisional.witness_metadata,
        verification,
    )
end
```

- [ ] **Step 5: Add verifier**

Add:

```julia
function verify_quillen_local_certificate(certificate)::Bool
    try
        replay = _quillen_local_certificate_replay_summary(certificate)
        return replay.overall_ok && certificate.verification == replay
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end
```

- [ ] **Step 6: Run focused GREEN verification**

Run:

```bash
julia --project=. -e 'include("test/expert/quillen_local_certificate.jl")'
```

Expected: exits 0 with all Quillen local certificate tests passing.

### Task 3: Expert Registration And Branch Verification

**Files:**
- Modify: `test/runtests.jl`

**Interfaces:**
- Consumes: focused test from Task 1.
- Produces: expert group coverage for the new replay certificate.

- [ ] **Step 1: Register the expert test**

In `test/runtests.jl`, add the new file after `"expert/quillen_patching_exact.jl",`:

```julia
"expert/quillen_local_certificate.jl",
```

- [ ] **Step 2: Run required focused command**

Run:

```bash
julia --project=. -e 'include("test/expert/quillen_local_certificate.jl")'
```

Expected: exits 0.

- [ ] **Step 3: Run expert group**

Run:

```bash
julia --project=. test/runtests.jl expert
```

Expected: exits 0 and includes the new expert test.

- [ ] **Step 4: Run full package tests**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: exits 0.

- [ ] **Step 5: Inspect diff and public API**

Run:

```bash
git status --short
git diff --stat origin/main...HEAD
git diff --name-only origin/main...HEAD
rg -n "QuillenLocalRealizationCertificate|verify_quillen_local_certificate|quillen_local_realization_certificate" src/Suslin.jl test/public/api_surface.jl
```

Expected: source, expert test, test registration, and Superpowers docs are changed; no public export or API surface lines are added for the new expert/internal names.

- [ ] **Step 6: Commit implementation**

Commit the implementation and tests:

```bash
git add src/algorithm/quillen_induction.jl test/expert/quillen_local_certificate.jl test/runtests.jl docs/superpowers/plans/2026-06-23-issue-100-quillen-local-realization-certificates.md
git commit -m "feat: add quillen local realization certificates"
```

Expected: commit succeeds on branch `agent/issue-100-add-replayable-local-realization-certificates-fo-run-1`.
