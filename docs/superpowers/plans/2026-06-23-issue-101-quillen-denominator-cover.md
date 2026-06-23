# Issue 101 Quillen Denominator Cover Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add replayable exact denominator cover certificates for deterministic Quillen patching fixtures.

**Architecture:** Keep the new cover certificate expert/internal and define it in `src/algorithm/quillen_induction.jl`, next to existing Quillen induction helpers. The constructor normalizes supplied denominators and multipliers into an exact ordinary polynomial ring, stores the exact sum, and exposes replay/verification helpers without exporting new public API.

**Tech Stack:** Julia, Oscar ordinary polynomial rings, existing Suslin coercion helpers, Test stdlib.

## Global Constraints

- No `AGENTS.md`, `CLAUDE.md`, or `GEMINI.md` file is present in this checkout.
- Dependency #99 is closed and merged; reuse `test/fixtures/quillen_patch_cases.jl`.
- Do not implement a general Nullstellensatz or cover-solving engine.
- Do not assemble global Quillen factors.
- Do not solve local realizability.
- Keep new cover certificate names expert/internal and do not export them.
- Tests must use qualified `Suslin.<name>` access.
- Include a positive cover with denominators `r` and `1-r`.
- Include a positive cover where a multiplier is not `1`.
- Make uncovered and inexact-ring failures explicit `ArgumentError`s.
- Focused command is `julia --project=. -e 'include("test/expert/quillen_denominator_cover.jl")'`.
- Expert group command is `julia --project=. test/runtests.jl expert`.
- Full package command is `julia --project=. -e 'using Pkg; Pkg.test()'`.

---

## File Structure

- Modify `src/algorithm/quillen_induction.jl`: add `QuillenDenominatorCoverVerification`, `QuillenDenominatorCoverCertificate`, `quillen_denominator_cover_certificate`, `replay_quillen_denominator_cover`, and `verify_quillen_denominator_cover`.
- Create `test/expert/quillen_denominator_cover.jl`: build certificates from #99 fixture catalog entries and test replay plus negative controls.
- Modify `test/runtests.jl`: add the focused test to the `expert` group after `expert/quillen_patching_exact.jl`.
- Leave `src/Suslin.jl` and `test/public/api_surface.jl` unchanged.

---

### Task 1: Expert Cover Certificate Test

**Files:**
- Create: `test/expert/quillen_denominator_cover.jl`
- Modify: `test/runtests.jl`

**Interfaces:**
- Consumes: `QuillenPatchFixtureCatalog.cases_by_id()`.
- Produces: RED tests for `Suslin.quillen_denominator_cover_certificate`, `Suslin.replay_quillen_denominator_cover`, and `Suslin.verify_quillen_denominator_cover`.

- [ ] **Step 1: Write the failing test**

Create `test/expert/quillen_denominator_cover.jl` with tests that:

```julia
using Test
using Suslin
using Oscar

const QUILLEN_PATCH_CATALOG_PATH = joinpath(@__DIR__, "..", "fixtures", "quillen_patch_cases.jl")
if !isdefined(Main, :QuillenPatchFixtureCatalog)
    include(QUILLEN_PATCH_CATALOG_PATH)
end

function _cover_inputs(entry)
    return (
        [data.denominator for data in entry.denominator_data],
        [data.coverage_multiplier for data in entry.denominator_data],
    )
end

@testset "Quillen denominator cover certificates" begin
    entries = Main.QuillenPatchFixtureCatalog.cases_by_id()
    for id in ("quillen-two-open-cover-qq", "quillen-nontrivial-multipliers-qq")
        entry = entries[id]
        R = entry.ring.object
        denominators, multipliers = _cover_inputs(entry)
        certificate = Suslin.quillen_denominator_cover_certificate(R, denominators, multipliers)
        replay = Suslin.replay_quillen_denominator_cover(certificate)

        @test Suslin.verify_quillen_denominator_cover(certificate)
        @test certificate.coverage_sum == one(R)
        @test replay.coverage_sum == certificate.coverage_sum
        @test replay.coverage_ok
        @test certificate.denominators == [R(denominator) for denominator in denominators]
        @test certificate.coverage_multipliers == [R(multiplier) for multiplier in multipliers]
        @test all(denominator -> parent(denominator) == R, certificate.denominators)
        @test all(multiplier -> parent(multiplier) == R, certificate.coverage_multipliers)
    end

    two_open = entries["quillen-two-open-cover-qq"]
    R_two = two_open.ring.object
    r = two_open.denominator_data[1].denominator
    @test [data.denominator for data in two_open.denominator_data] == [r, one(R_two) - r]

    nontrivial = entries["quillen-nontrivial-multipliers-qq"]
    R_nontrivial = nontrivial.ring.object
    @test any(data -> data.coverage_multiplier != one(R_nontrivial), nontrivial.denominator_data)

    denominators, multipliers = _cover_inputs(two_open)
    @test_throws ArgumentError Suslin.quillen_denominator_cover_certificate(
        R_two,
        denominators[1:1],
        multipliers[1:1],
    )
    @test_throws ArgumentError Suslin.quillen_denominator_cover_certificate(
        R_two,
        denominators,
        [one(R_two), zero(R_two)],
    )

    RR, (Y, s) = Oscar.polynomial_ring(RealField(), ["Y", "s"])
    @test_throws ArgumentError Suslin.quillen_denominator_cover_certificate(
        RR,
        [s, one(RR) - s],
        [one(RR), one(RR)],
    )

    certificate = Suslin.quillen_denominator_cover_certificate(R_two, denominators, multipliers)
    tampered = Suslin.QuillenDenominatorCoverCertificate(
        certificate.ring,
        certificate.denominators,
        [one(R_two), zero(R_two)],
        certificate.coverage_sum,
        certificate.verification,
    )
    @test !Suslin.verify_quillen_denominator_cover(tampered)
end
```

- [ ] **Step 2: Register the expert test**

In `test/runtests.jl`, add:

```julia
        "expert/quillen_denominator_cover.jl",
```

immediately after `expert/quillen_patching_exact.jl`.

- [ ] **Step 3: Run RED verification**

Run:

```bash
julia --project=. -e 'include("test/expert/quillen_denominator_cover.jl")'
```

Expected: FAIL with `UndefVarError` for the new `Suslin.quillen_denominator_cover_certificate` helper.

### Task 2: Cover Certificate Implementation

**Files:**
- Modify: `src/algorithm/quillen_induction.jl`

**Interfaces:**
- Consumes: `_coerce_into_ring(R, value, label)`.
- Produces:
  - `QuillenDenominatorCoverVerification`
  - `QuillenDenominatorCoverCertificate`
  - `quillen_denominator_cover_certificate(R, denominators, coverage_multipliers)`
  - `replay_quillen_denominator_cover(certificate)`
  - `verify_quillen_denominator_cover(certificate)`

- [ ] **Step 1: Add internal data structs**

Append to `src/algorithm/quillen_induction.jl`:

```julia
struct QuillenDenominatorCoverVerification
    denominator_count::Int
    multiplier_count::Int
    parent_ring_ok::Bool
    exact_ring_ok::Bool
    coverage_terms::Vector
    coverage_sum
    coverage_ok::Bool
end

struct QuillenDenominatorCoverCertificate
    ring
    denominators::Vector
    coverage_multipliers::Vector
    coverage_sum
    verification::QuillenDenominatorCoverVerification
end
```

- [ ] **Step 2: Add ring and normalization helpers**

Add helpers that accept only exact ordinary polynomial rings:

```julia
function _require_quillen_denominator_cover_ring(R)
    R isa MPolyRing || R isa PolyRing ||
        throw(ArgumentError("staged Quillen denominator cover certificates require a supported exact ordinary polynomial ring"))
    Oscar.is_exact_type(typeof(zero(coefficient_ring(R)))) ||
        throw(ArgumentError("staged Quillen denominator cover certificates require a supported exact ordinary polynomial ring"))
    return R
end

function _normalize_quillen_cover_elements(R, values, label::AbstractString)
    collected = collect(values)
    return [_coerce_into_ring(R, value, label) for value in collected]
end
```

- [ ] **Step 3: Add replay and constructor**

Add exact replay logic:

```julia
function _quillen_denominator_cover_verification(R, denominators, coverage_multipliers)
    denominator_count = length(denominators)
    multiplier_count = length(coverage_multipliers)
    exact_ring_ok = Oscar.is_exact_type(typeof(zero(coefficient_ring(R))))
    parent_ring_ok =
        denominator_count == multiplier_count &&
        all(denominator -> parent(denominator) == R, denominators) &&
        all(multiplier -> parent(multiplier) == R, coverage_multipliers)
    coverage_terms = parent_ring_ok ?
        [coverage_multipliers[idx] * denominators[idx] for idx in eachindex(denominators)] :
        Any[]
    coverage_sum = parent_ring_ok ? sum(coverage_terms; init = zero(R)) : zero(R)
    coverage_ok = parent_ring_ok && exact_ring_ok && coverage_sum == one(R)
    return QuillenDenominatorCoverVerification(
        denominator_count,
        multiplier_count,
        parent_ring_ok,
        exact_ring_ok,
        coverage_terms,
        coverage_sum,
        coverage_ok,
    )
end

function quillen_denominator_cover_certificate(R, denominators, coverage_multipliers)
    _require_quillen_denominator_cover_ring(R)
    normalized_denominators = _normalize_quillen_cover_elements(R, denominators, "cover denominator")
    normalized_multipliers = _normalize_quillen_cover_elements(R, coverage_multipliers, "cover coverage multiplier")
    isempty(normalized_denominators) &&
        throw(ArgumentError("staged Quillen denominator cover certificates require at least one denominator"))
    length(normalized_denominators) == length(normalized_multipliers) ||
        throw(ArgumentError("staged Quillen denominator cover certificates require matching denominator and multiplier counts"))
    verification = _quillen_denominator_cover_verification(R, normalized_denominators, normalized_multipliers)
    verification.coverage_ok ||
        throw(ArgumentError("staged Quillen denominator cover certificate requires coverage sum to equal one"))
    return QuillenDenominatorCoverCertificate(
        R,
        normalized_denominators,
        normalized_multipliers,
        verification.coverage_sum,
        verification,
    )
end
```

- [ ] **Step 4: Add replay equality and verifier**

Add:

```julia
function replay_quillen_denominator_cover(certificate::QuillenDenominatorCoverCertificate)
    _require_quillen_denominator_cover_ring(certificate.ring)
    return _quillen_denominator_cover_verification(
        certificate.ring,
        certificate.denominators,
        certificate.coverage_multipliers,
    )
end

function _same_quillen_denominator_cover_verification(left::QuillenDenominatorCoverVerification, right::QuillenDenominatorCoverVerification)::Bool
    return left.denominator_count == right.denominator_count &&
           left.multiplier_count == right.multiplier_count &&
           left.parent_ring_ok == right.parent_ring_ok &&
           left.exact_ring_ok == right.exact_ring_ok &&
           left.coverage_terms == right.coverage_terms &&
           left.coverage_sum == right.coverage_sum &&
           left.coverage_ok == right.coverage_ok
end

function verify_quillen_denominator_cover(certificate::QuillenDenominatorCoverCertificate)::Bool
    try
        replay = replay_quillen_denominator_cover(certificate)
        return certificate.coverage_sum == replay.coverage_sum &&
               _same_quillen_denominator_cover_verification(certificate.verification, replay) &&
               replay.parent_ring_ok &&
               replay.exact_ring_ok &&
               replay.coverage_ok
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end
```

- [ ] **Step 5: Run GREEN verification**

Run:

```bash
julia --project=. -e 'include("test/expert/quillen_denominator_cover.jl")'
```

Expected: exits 0.

### Task 3: Group and Package Verification

**Files:**
- Verify all changed files.

**Interfaces:**
- Produces: branch ready for review and PR.

- [ ] **Step 1: Run expert group**

Run:

```bash
julia --project=. test/runtests.jl expert
```

Expected: exits 0 and includes `expert/quillen_denominator_cover.jl`.

- [ ] **Step 2: Run package tests**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: exits 0.

- [ ] **Step 3: Inspect diff**

Run:

```bash
git status --short
git diff --stat origin/main
git diff --name-only origin/main
```

Expected: only the issue #101 spec, plan, implementation file, new expert test,
and `test/runtests.jl` changed.

## Plan Self-Review

- Spec coverage: all issue-required positive and negative controls map to Task
  1 tests and Task 2 implementation.
- Placeholder scan: no incomplete or deferred implementation markers are
  present.
- Type consistency: function and type names match between tests and
  implementation tasks.

## Execution Choice

Plan complete and saved. Under the Standing Answer Policy, choose
Subagent-Driven (recommended).
