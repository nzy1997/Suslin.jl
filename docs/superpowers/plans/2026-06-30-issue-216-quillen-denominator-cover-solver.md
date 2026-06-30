# Issue 216 Quillen Denominator Cover Solver Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a bounded exact solver that turns #215 raw denominator-cover candidates into replayable ordinary-polynomial denominator cover certificates.

**Architecture:** Extend `src/algorithm/quillen_induction.jl` beside the existing denominator-cover certificate and #215 candidate replay code. The solver searches bounded powers with Oscar ideal coordinates, stores raw denominators separately from powered cover terms, and adapts proven powered terms through the existing `QuillenDenominatorCoverCertificate` verifier.

**Tech Stack:** Julia, Oscar ordinary polynomial rings and ideals, existing Suslin Quillen certificate structs, `Test`.

## Global Constraints

- Input is a #215 denominator-cover candidate over an exact ordinary polynomial ring, plus bounded solver options such as maximum exponent and optional supplied multipliers.
- Output records raw extracted denominators `r_i`, chosen exponent `l`, powered cover terms `r_i^l`, coverage multipliers `g_i`, exact coverage terms `g_i * r_i^l`, and a replayable `QuillenDenominatorCoverCertificate`.
- Keep raw denominators and powered cover terms distinct.
- If automatic membership cannot prove a cover within the bound, throw `ArgumentError` with a message containing `coverage not proven`.
- If supplied multipliers are provided, replay the exact identity and reject stale or tampered multipliers.
- Do not produce patched-substitution brackets, global elementary factors, base-term factors, or public route selection.
- Keep new names expert/internal and do not export them from `src/Suslin.jl`.
- Focused solver command is `julia --project=. -e 'include("test/expert/quillen_denominator_cover_solver.jl")'`.
- Existing verifier command is `julia --project=. -e 'include("test/expert/quillen_denominator_cover.jl")'`.
- Full package command is `julia --project=. -e 'using Pkg; Pkg.test()'`.

---

## File Structure

- Modify `src/algorithm/quillen_induction.jl`: add solver result structs, bounded solver, Oscar coordinate helper, replay, and verification.
- Create `test/expert/quillen_denominator_cover_solver.jl`: positive automatic and supplied-multiplier solver tests, plus uncovered and tampered negative controls.
- Modify `test/expert/quillen_denominator_cover.jl`: add a small powered-denominator certificate replay check.
- Modify `test/runtests.jl`: add the new expert solver test after `expert/quillen_denominator_cover.jl`.

### Task 1: Add Red Solver Tests

**Files:**
- Create: `test/expert/quillen_denominator_cover_solver.jl`
- Modify: `test/expert/quillen_denominator_cover.jl`
- Modify: `test/runtests.jl`

**Interfaces:**
- Consumes: future `Suslin.solve_quillen_denominator_cover`, `Suslin.verify_quillen_denominator_cover_solver_result`, `Suslin.replay_quillen_denominator_cover_solver_result`, and `Suslin.QuillenDenominatorCoverSolverResult`.
- Produces: failing behavioral coverage for the solver implementation.

- [ ] **Step 1: Write the failing solver test**

Create `test/expert/quillen_denominator_cover_solver.jl`:

```julia
using Test
using Suslin
using Oscar

const QUILLEN_PATCH_CATALOG_PATH = joinpath(@__DIR__, "..", "fixtures", "quillen_patch_cases.jl")
if !isdefined(Main, :QuillenPatchFixtureCatalog)
    include(QUILLEN_PATCH_CATALOG_PATH)
end

function _solver_tamper_multiplier(result)
    multipliers = copy(result.coverage_multipliers)
    multipliers[1] += one(result.ring)
    return Suslin.QuillenDenominatorCoverSolverResult(
        result.source_candidate,
        result.ring,
        result.raw_denominators,
        result.exponent,
        result.powered_denominators,
        multipliers,
        result.coverage_terms,
        result.coverage_sum,
        result.cover_certificate,
        result.verification,
    )
end

function _solver_tamper_exponent(result)
    return Suslin.QuillenDenominatorCoverSolverResult(
        result.source_candidate,
        result.ring,
        result.raw_denominators,
        result.exponent + 1,
        result.powered_denominators,
        result.coverage_multipliers,
        result.coverage_terms,
        result.coverage_sum,
        result.cover_certificate,
        result.verification,
    )
end

function _coverage_error_message(err)
    return err isa ArgumentError && occursin("coverage not proven", sprint(showerror, err))
end

@testset "Quillen denominator cover solver" begin
    entries = Main.QuillenPatchFixtureCatalog.cases_by_id()

    two_open = entries["quillen-two-open-cover-qq"]
    R = two_open.ring.object
    raw = [data.denominator for data in two_open.denominator_data]
    result = Suslin.solve_quillen_denominator_cover(R, raw; max_exponent = 2)

    @test result isa Suslin.QuillenDenominatorCoverSolverResult
    @test Suslin.verify_quillen_denominator_cover_solver_result(result)
    @test Suslin.verify_quillen_denominator_cover(result.cover_certificate)
    @test result.raw_denominators == raw
    @test result.exponent == 1
    @test result.powered_denominators == raw
    @test result.coverage_multipliers == [one(R), one(R)]
    @test result.coverage_terms == [result.coverage_multipliers[i] * result.powered_denominators[i] for i in eachindex(raw)]
    @test result.coverage_sum == one(R)
    @test result.cover_certificate.denominators == result.powered_denominators

    squared = Suslin.solve_quillen_denominator_cover(
        R,
        raw;
        max_exponent = 2,
        exponent = 2,
    )
    @test squared.exponent == 2
    @test squared.powered_denominators == [denominator^2 for denominator in raw]
    @test squared.raw_denominators == raw
    @test squared.cover_certificate.denominators == squared.powered_denominators
    @test Suslin.verify_quillen_denominator_cover_solver_result(squared)

    nontrivial = entries["quillen-nontrivial-multipliers-qq"]
    R_nontrivial = nontrivial.ring.object
    raw_nontrivial = [data.denominator for data in nontrivial.denominator_data]
    supplied = [data.coverage_multiplier for data in nontrivial.denominator_data]
    supplied_result = Suslin.solve_quillen_denominator_cover(
        R_nontrivial,
        raw_nontrivial;
        max_exponent = 1,
        coverage_multipliers = supplied,
    )
    @test supplied_result.coverage_multipliers == supplied
    @test any(multiplier -> multiplier != one(R_nontrivial), supplied_result.coverage_multipliers)
    @test supplied_result.coverage_sum == one(R_nontrivial)
    @test Suslin.verify_quillen_denominator_cover_solver_result(supplied_result)

    @test !Suslin.verify_quillen_denominator_cover_solver_result(_solver_tamper_multiplier(result))
    @test !Suslin.verify_quillen_denominator_cover_solver_result(_solver_tamper_exponent(result))

    R_bad, (X_bad, r_bad, s_bad) = Oscar.polynomial_ring(QQ, ["X", "r", "s"])
    try
        Suslin.solve_quillen_denominator_cover(R_bad, [r_bad, s_bad]; max_exponent = 2)
        @test false
    catch err
        @test _coverage_error_message(err)
    end

    try
        Suslin.solve_quillen_denominator_cover(
            R_nontrivial,
            raw_nontrivial;
            max_exponent = 1,
            coverage_multipliers = [supplied[1] + one(R_nontrivial), supplied[2]],
        )
        @test false
    catch err
        @test _coverage_error_message(err)
    end
end
```

- [ ] **Step 2: Add a powered-term replay check to the existing verifier test**

Append inside the existing `@testset` in `test/expert/quillen_denominator_cover.jl`:

```julia
    powered_certificate = Suslin.quillen_denominator_cover_certificate(
        R_two,
        [r^2, (one(R_two) - r)^2],
        [-2 * r + 3 * one(R_two), 2 * r + one(R_two)],
    )
    @test Suslin.verify_quillen_denominator_cover(powered_certificate)
    @test powered_certificate.coverage_sum == one(R_two)
    tampered_powered = Suslin.QuillenDenominatorCoverCertificate(
        powered_certificate.ring,
        powered_certificate.denominators,
        powered_certificate.coverage_multipliers,
        powered_certificate.coverage_sum + r,
        powered_certificate.verification,
    )
    @test !Suslin.verify_quillen_denominator_cover(tampered_powered)
```

- [ ] **Step 3: Register the red test**

In `test/runtests.jl`, add:

```julia
        "expert/quillen_denominator_cover_solver.jl",
```

immediately after:

```julia
        "expert/quillen_denominator_cover.jl",
```

- [ ] **Step 4: Run the focused red test**

Run:

```bash
julia --project=. -e 'include("test/expert/quillen_denominator_cover_solver.jl")'
```

Expected: FAIL with `UndefVarError` for `QuillenDenominatorCoverSolverResult` or `solve_quillen_denominator_cover`.

### Task 2: Implement Bounded Exact Solver

**Files:**
- Modify: `src/algorithm/quillen_induction.jl`

**Interfaces:**
- Consumes: `_require_quillen_denominator_cover_ring`, `_normalize_quillen_cover_elements`, `quillen_denominator_cover_certificate`, `verify_quillen_denominator_cover`, `QuillenDenominatorCoverCandidate`, and `verify_quillen_denominator_cover_candidate`.
- Produces: `QuillenDenominatorCoverSolverVerification`, `QuillenDenominatorCoverSolverResult`, `solve_quillen_denominator_cover`, `replay_quillen_denominator_cover_solver_result`, and `verify_quillen_denominator_cover_solver_result`.

- [ ] **Step 1: Add solver replay structs**

Add after `QuillenDenominatorCoverCandidate`:

```julia
struct QuillenDenominatorCoverSolverVerification
    raw_denominator_count::Int
    multiplier_count::Int
    raw_denominators::Vector
    exponent::Int
    powered_denominators::Vector
    coverage_multipliers::Vector
    parent_ring_ok::Bool
    exact_ring_ok::Bool
    exponent_ok::Bool
    source_candidate_ok::Bool
    coverage_terms::Vector
    coverage_sum
    coverage_ok::Bool
    cover_certificate_ok::Bool
    cover_certificate_matches::Bool
    overall_ok::Bool
end

struct QuillenDenominatorCoverSolverResult
    source_candidate
    ring
    raw_denominators::Vector
    exponent::Int
    powered_denominators::Vector
    coverage_multipliers::Vector
    coverage_terms::Vector
    coverage_sum
    cover_certificate::QuillenDenominatorCoverCertificate
    verification::QuillenDenominatorCoverSolverVerification
end
```

- [ ] **Step 2: Add solver verification helpers**

Implement helpers with these signatures:

```julia
function _quillen_denominator_cover_solver_verification(source_candidate, R, raw_denominators, exponent::Int, coverage_multipliers, cover_certificate)
    # compute parent/exact/exponent/candidate checks, powered terms, coverage terms,
    # certificate replay, and overall status
end

function _same_quillen_denominator_cover_solver_verification(
    left::QuillenDenominatorCoverSolverVerification,
    right::QuillenDenominatorCoverSolverVerification,
)::Bool
    # compare every stored field
end
```

The verification must compute:

```julia
powered_denominators = [denominator^exponent for denominator in raw_denominators]
coverage_terms = [coverage_multipliers[idx] * powered_denominators[idx] for idx in eachindex(powered_denominators)]
coverage_sum = sum(coverage_terms; init = zero(R))
coverage_ok = parent_ring_ok && exact_ring_ok && exponent_ok && coverage_sum == one(R)
cover_certificate_ok = cover_certificate isa QuillenDenominatorCoverCertificate &&
    verify_quillen_denominator_cover(cover_certificate)
cover_certificate_matches = cover_certificate_ok &&
    cover_certificate.ring == R &&
    cover_certificate.denominators == powered_denominators &&
    cover_certificate.coverage_multipliers == coverage_multipliers &&
    cover_certificate.coverage_sum == coverage_sum
overall_ok = coverage_ok && source_candidate_ok && cover_certificate_matches
```

- [ ] **Step 3: Add replay and result verification**

Implement:

```julia
function replay_quillen_denominator_cover_solver_result(result::QuillenDenominatorCoverSolverResult)
    _require_quillen_denominator_cover_ring(result.ring)
    return _quillen_denominator_cover_solver_verification(
        result.source_candidate,
        result.ring,
        result.raw_denominators,
        result.exponent,
        result.coverage_multipliers,
        result.cover_certificate,
    )
end

function verify_quillen_denominator_cover_solver_result(result)::Bool
    try
        replay = replay_quillen_denominator_cover_solver_result(result)
        return replay.overall_ok &&
               result.raw_denominators == replay.raw_denominators &&
               result.exponent == replay.exponent &&
               result.powered_denominators == replay.powered_denominators &&
               result.coverage_multipliers == replay.coverage_multipliers &&
               result.coverage_terms == replay.coverage_terms &&
               result.coverage_sum == replay.coverage_sum &&
               _same_quillen_denominator_cover_solver_verification(result.verification, replay)
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end
```

- [ ] **Step 4: Add bounded solver helpers**

Implement helpers that normalize options, build results, and call Oscar:

```julia
function _quillen_cover_exponent_range(max_exponent::Integer, exponent)
    bound = Int(max_exponent)
    bound >= 1 || throw(ArgumentError("coverage not proven: max_exponent must be positive"))
    if exponent === nothing
        return 1:bound
    end
    exponent isa Integer || throw(ArgumentError("coverage not proven: exponent must be an integer"))
    chosen = Int(exponent)
    1 <= chosen <= bound ||
        throw(ArgumentError("coverage not proven: requested exponent is outside the configured bound"))
    return chosen:chosen
end

function _quillen_supplied_cover_multipliers(coverage_multipliers, supplied_multipliers)
    coverage_multipliers !== nothing && supplied_multipliers !== nothing &&
        throw(ArgumentError("coverage not proven: provide only one supplied multiplier collection"))
    return coverage_multipliers === nothing ? supplied_multipliers : coverage_multipliers
end
```

Use `Oscar.coordinates(one(R), ideal(R, powered_denominators))` for automatic
membership and catch the non-membership exception by returning `nothing`.

- [ ] **Step 5: Add solver overloads**

Implement:

```julia
function solve_quillen_denominator_cover(
    R,
    raw_denominators;
    max_exponent::Integer = 4,
    exponent = nothing,
    coverage_multipliers = nothing,
    supplied_multipliers = nothing,
    source_candidate = nothing,
)
    # normalize raw denominators, search requested exponents, return first replaying result
end

function solve_quillen_denominator_cover(
    candidate::QuillenDenominatorCoverCandidate;
    max_exponent::Integer = 4,
    exponent = nothing,
    coverage_multipliers = nothing,
    supplied_multipliers = nothing,
)
    verify_quillen_denominator_cover_candidate(candidate) ||
        throw(ArgumentError("coverage not proven: denominator-cover candidate does not replay"))
    return solve_quillen_denominator_cover(
        candidate.ring,
        candidate.raw_denominators;
        max_exponent,
        exponent,
        coverage_multipliers,
        supplied_multipliers,
        source_candidate = candidate,
    )
end
```

All failure exits after bounded search must throw an `ArgumentError` whose
message includes `coverage not proven`.

- [ ] **Step 6: Run the focused solver test**

Run:

```bash
julia --project=. -e 'include("test/expert/quillen_denominator_cover_solver.jl")'
```

Expected: PASS.

### Task 3: Final Verification, Review, and Commit

**Files:**
- Modify: all files changed by Tasks 1-2.

**Interfaces:**
- Consumes: completed solver implementation.
- Produces: committed worker-branch change ready for pull request creation.

- [ ] **Step 1: Run focused issue verification**

Run:

```bash
julia --project=. -e 'include("test/expert/quillen_denominator_cover_solver.jl")'
julia --project=. -e 'include("test/expert/quillen_denominator_cover.jl")'
```

Expected: both commands exit 0.

- [ ] **Step 2: Run package entry-point verification**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: exit 0.

- [ ] **Step 3: Run diff hygiene**

Run:

```bash
git diff --check
```

Expected: exit 0 with no whitespace errors.

- [ ] **Step 4: Review the changed files**

Run:

```bash
git status --short
git diff --stat
git diff -- src/algorithm/quillen_induction.jl test/expert/quillen_denominator_cover_solver.jl test/expert/quillen_denominator_cover.jl test/runtests.jl
```

Expected: only the solver implementation, solver tests, existing verifier test addition, runtests registration, and Superpowers spec/plan are changed.

- [ ] **Step 5: Commit**

Run:

```bash
git add docs/superpowers/specs/2026-06-30-issue-216-quillen-denominator-cover-solver-design.md \
    docs/superpowers/plans/2026-06-30-issue-216-quillen-denominator-cover-solver.md \
    src/algorithm/quillen_induction.jl \
    test/expert/quillen_denominator_cover_solver.jl \
    test/expert/quillen_denominator_cover.jl \
    test/runtests.jl
git commit -m "Implement exact Quillen denominator cover solver"
```

Expected: commit succeeds on the worker branch.
