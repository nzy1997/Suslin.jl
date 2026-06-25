# Issue 135 Case 010 Laurent GL Certificate Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Promote original-input ToricBuilder `case_010` from a Q-block route error to a verified Laurent GL certificate pass.

**Architecture:** Reuse the existing certified `:laurent_unit_creation` column-reduction stage added for #134. The remaining full-matrix boundary is a later length-4 column with the same exact unit-creation relation, so the implementation broadens that stage's length guard, deterministically prefers downstream-compatible unit-creation candidates, and updates acceptance/status coverage.

**Tech Stack:** Julia, Oscar, Suslin internal Laurent GL certificates, Test stdlib.

## Global Constraints

- Repository has no `AGENTS.md`, `CLAUDE.md`, `GEMINI.md`, or `CONVENTIONS.md`; follow `README.md` test commands and existing Suslin style.
- Keep `public_elementary_status == :staged_boundary` for original `case_010`.
- Do not make `elementary_factorization(case_010)` return factors for the original Laurent `GL_n` input.
- Do not add a separate broad Laurent column algorithm.
- Reuse the existing `:laurent_unit_creation` certificate/replay path.
- When several exact unit-creation candidates exist, prefer the larger `pivot_index`, then the smaller `source_index` for ties, so the original `case_010` recursive peel reaches the existing supported final block.
- The issue focused command must be `julia --project=. -e 'include("test/internal/toricbuilder_cache_case010_certificate.jl")'`.
- The report verification command must be `julia --project=. scripts/report_toricbuilder_cache_q_blocks.jl --exercise=case_010 --output=/tmp/case010-q-block-status.md`.
- The package verification command must be `julia --project=. -e 'using Pkg; Pkg.test()'`.

---

## File Structure

- Create `test/internal/toricbuilder_cache_case010_certificate.jl`: focused full-matrix certificate acceptance test and corrupted-matrix negative control.
- Modify `test/runtests.jl`: register the new internal test in the default package test set.
- Modify `src/algorithm/column_reduction_case010.jl`: allow the existing exact Laurent unit-creation stage for the observed length-4 boundary as well as the prior length-5 boundary, and prefer downstream-compatible candidates when several exact row pairs work.
- Modify `test/internal/toricbuilder_cache_status_report.jl`: promote `case_010` expectations from route error to GL certificate pass.
- Modify `docs/audits/2026-06-24-toricbuilder-cache-q-block-status.md`: refresh the generated markdown report after the route passes.
- Keep `docs/superpowers/specs/2026-06-25-issue-135-case010-laurent-gl-certificate-design.md` and this plan in the PR.

---

### Task 1: Route Case 010 Through The Laurent GL Certificate

**Files:**
- Create: `test/internal/toricbuilder_cache_case010_certificate.jl`
- Modify: `test/runtests.jl`
- Modify: `src/algorithm/column_reduction_case010.jl`
- Modify: `test/internal/toricbuilder_cache_status_report.jl`
- Modify: `docs/audits/2026-06-24-toricbuilder-cache-q-block-status.md`

**Interfaces:**
- Consumes: `ToricBuilderCacheQBlockStatusReport.ToricBuilderCacheQBlocks.materialize_matrix(entry)`, `laurent_gl_factorization_certificate(A)`, `verify_laurent_gl_factorization_certificate(certificate)`, `elementary_factorization(A)`, and `ToricBuilderCacheQBlockStatusReport._exercised_row(entry)`.
- Produces: `test/internal/toricbuilder_cache_case010_certificate.jl`, a widened internal `:laurent_unit_creation` length guard with deterministic candidate preference, updated default internal status-report expectations, and refreshed generated audit markdown.

- [ ] **Step 1: Write the failing focused certificate test**

Create `test/internal/toricbuilder_cache_case010_certificate.jl` with this content:

```julia
using Test
using Suslin
using Oscar

const TORICBUILDER_CACHE_CASE010_REPORT_SCRIPT =
    joinpath(@__DIR__, "..", "..", "scripts", "report_toricbuilder_cache_q_blocks.jl")

include(TORICBUILDER_CACHE_CASE010_REPORT_SCRIPT)

function _case010_entry()
    return only(filter(
        entry -> entry.id == "case_010",
        ToricBuilderCacheQBlockStatusReport.ToricBuilderCacheQBlocks.catalog().cases,
    ))
end

function _case010_product(factors, R, n::Int)
    product = identity_matrix(R, n)
    for factor in factors
        product *= factor
    end
    return product
end

function _case010_certificate_with_original(certificate, original_matrix)
    return Suslin.LaurentGLFactorizationCertificate(
        original_matrix,
        certificate.determinant_profile,
        certificate.normalization,
        certificate.correction,
        certificate.inverse_correction,
        certificate.normalized_core,
        certificate.core_factorization,
        certificate.core_factors,
        certificate.reconstructed_product,
        certificate.verification,
    )
end

function _case010_corrupt_matrix(A)
    corrupted = copy(A)
    corrupted[1, 1] += one(base_ring(A))
    return corrupted
end

function _case010_corrupt_entry(entry)
    corrupted_entries = collect(entry.sparse_entries)
    row, col, _ = corrupted_entries[1]
    corrupted_entries[1] = (row, col, "0")
    return merge(entry, (; sparse_entries = corrupted_entries, sparse_entry_count = length(corrupted_entries)))
end

function _case010_verify_or_staged_false(certificate)
    try
        return verify_laurent_gl_factorization_certificate(certificate)
    catch err
        err isa InterruptException && rethrow()
        err isa ArgumentError || rethrow()
        message = sprint(showerror, err)
        occursin("staged", message) || occursin("unsupported", message) || rethrow()
        return false
    end
end

@testset "ToricBuilder case_010 Laurent GL certificate" begin
    entry = _case010_entry()
    A = ToricBuilderCacheQBlockStatusReport.ToricBuilderCacheQBlocks.materialize_matrix(entry)
    R = base_ring(A)

    profile = classify_laurent_determinant(A)
    @test profile.classification == :laurent_monomial_unit
    @test profile.determinant == only(gens(R)) * collect(gens(R))[2]

    certificate = laurent_gl_factorization_certificate(A)
    @test certificate.original_matrix == A
    @test certificate.normalized_core == certificate.normalization.normalized_matrix
    @test det(certificate.normalized_core) == one(R)
    @test length(certificate.core_factors) > 0
    @test length(certificate.core_factorization.peel_steps) == nrows(A) - 2
    @test any(step -> step.dimension == 4, certificate.core_factorization.peel_steps)
    @test _case010_product(certificate.core_factors, R, nrows(A)) == certificate.normalized_core
    @test certificate.reconstructed_product == A
    @test verify_laurent_gl_factorization_certificate(certificate)

    public_error = try
        elementary_factorization(A)
        nothing
    catch err
        err
    end
    @test public_error isa ArgumentError
    @test occursin("Laurent GL_n normalization boundary", sprint(showerror, public_error))

    corrupted_matrix = _case010_corrupt_matrix(A)
    corrupted_certificate = _case010_certificate_with_original(certificate, corrupted_matrix)
    @test !_case010_verify_or_staged_false(corrupted_certificate)

    corrupted_entry = _case010_corrupt_entry(entry)
    corrupted_row = ToricBuilderCacheQBlockStatusReport._exercised_row(corrupted_entry)
    @test !(corrupted_row.route_status == :gl_certificate_pass && corrupted_row.verified)
end
```

- [ ] **Step 2: Register the focused test in the internal group**

In `test/runtests.jl`, add the new internal test immediately after
`"internal/toricbuilder_cache_status_report.jl",`:

```julia
        "internal/toricbuilder_cache_case010_certificate.jl",
```

- [ ] **Step 3: Run the focused test to verify RED**

Run:

```bash
julia --project=. -e 'include("test/internal/toricbuilder_cache_case010_certificate.jl")'
```

Expected: FAIL before implementation with `unsupported exact unimodular column reduction for Laurent-normalized column of length 4`. If the failure is a syntax or import error, fix the test file and rerun until the failure proves the current route boundary.

- [ ] **Step 4: Implement the minimal production change**

In `src/algorithm/column_reduction_case010.jl`, add the candidate preference
helper:

```julia
function _preferred_laurent_unit_creation_candidate(left, right)
    left === nothing && return right
    # Prefer lower rows so the recursive column-peel keeps the exposed case_010
    # trailing blocks in the existing supported Laurent families.
    right.pivot_index > left.pivot_index && return right
    right.pivot_index == left.pivot_index && right.source_index < left.source_index && return right
    return left
end
```

Then change:

```julia
    length(column) == 5 || return nothing
```

to:

```julia
    length(column) in (4, 5) || return nothing
    target_unit = one(R)
    candidate = nothing
```

and replace the current immediate `return (; pivot_index = ..., ...)` inside
the nested loop with:

```julia
        candidate = _preferred_laurent_unit_creation_candidate(candidate, (;
            pivot_index = pivot_idx,
            source_index = source_idx,
            target_unit,
            creation_coefficient = coeff,
        ))
```

Return `candidate` after the loop. Do not change the stage kind, replay method,
factor ordering, or public API.

- [ ] **Step 5: Run the focused test to verify GREEN**

Run:

```bash
julia --project=. -e 'include("test/internal/toricbuilder_cache_case010_certificate.jl")'
```

Expected: PASS with all assertions in `ToricBuilder case_010 Laurent GL certificate`.

- [ ] **Step 6: Update status-report expectations**

In `test/internal/toricbuilder_cache_status_report.jl`, replace the `case_010`
route-error assertions:

```julia
    @test by_id["case_010"].route_status == :route_error
    @test by_id["case_010"].public_elementary_status == :staged_boundary
    @test by_id["case_010"].determinant_class == :laurent_monomial_unit
    @test by_id["case_010"].verified == false
    @test by_id["case_010"].runtime_seconds > 0
    @test occursin("unsupported exact unimodular column reduction", by_id["case_010"].error_details)
```

with:

```julia
    @test by_id["case_010"].route_status == :gl_certificate_pass
    @test by_id["case_010"].public_elementary_status == :staged_boundary
    @test by_id["case_010"].determinant_class == :laurent_monomial_unit
    @test by_id["case_010"].verified == true
    @test by_id["case_010"].decomposed_base_matrix_count > 0
    @test by_id["case_010"].runtime_seconds > 0
    @test by_id["case_010"].error_details == "none"
```

Replace the `case_010` markdown regex:

```julia
    @test occursin(r"\| case_010 \| 6x6 \| 34 \| default_contract \| route_error \| staged_boundary \| laurent_monomial_unit \| 0 \| [0-9]+\.[0-9]{3} \|", markdown)
```

with:

```julia
    @test occursin(r"\| case_010 \| 6x6 \| 34 \| default_contract \| gl_certificate_pass \| staged_boundary \| laurent_monomial_unit \| [1-9][0-9]* \| [0-9]+\.[0-9]{3} \|", markdown)
```

Replace the route-error-details assertions:

```julia
    @test occursin("## Route Error Details", markdown)
    @test occursin("unsupported exact unimodular column reduction", markdown)
```

with:

```julia
    @test !occursin("## Route Error Details", markdown)
    @test !occursin("unsupported exact unimodular column reduction", markdown)
```

- [ ] **Step 7: Run status-report test**

Run:

```bash
julia --project=. -e 'include("test/internal/toricbuilder_cache_status_report.jl")'
```

Expected: PASS with the updated `case_010` status expectations.

- [ ] **Step 8: Refresh the generated markdown report**

Run:

```bash
julia --project=. scripts/report_toricbuilder_cache_q_blocks.jl
```

Expected: `docs/audits/2026-06-24-toricbuilder-cache-q-block-status.md` is rewritten and its `case_010` row contains `gl_certificate_pass`, `staged_boundary`, `laurent_monomial_unit`, and a positive decomposed base-matrix count. The markdown must not contain a `Route Error Details` section.

- [ ] **Step 9: Run issue report verification**

Run:

```bash
julia --project=. scripts/report_toricbuilder_cache_q_blocks.jl --exercise=case_010 --output=/tmp/case010-q-block-status.md
```

Expected: PASS, and `/tmp/case010-q-block-status.md` contains `case_010`, `gl_certificate_pass`, `verified `true` in the exercised evidence line, a positive decomposed base-matrix count, and no `Route Error Details` section.

- [ ] **Step 10: Run package verification**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: PASS for public and internal test groups.

- [ ] **Step 11: Run diff hygiene check**

Run:

```bash
git diff --check origin/main
```

Expected: PASS with no whitespace errors.

- [ ] **Step 12: Commit**

Run:

```bash
git add src/algorithm/column_reduction_case010.jl test/internal/toricbuilder_cache_case010_certificate.jl test/internal/toricbuilder_cache_status_report.jl test/runtests.jl docs/audits/2026-06-24-toricbuilder-cache-q-block-status.md docs/superpowers/plans/2026-06-25-issue-135-case010-laurent-gl-certificate.md
git commit -m "feat: route case010 Laurent GL certificate"
```

Expected: commit succeeds and contains the implementation, focused test, status-report expectation update, refreshed generated report, and this plan.

---

## Self-Review

- Spec coverage: Task 1 covers the full-matrix certificate route, verifier result, exact reconstruction, determinant normalization, positive decomposed count, public staged boundary, corrupted-matrix negative control, status-report promotion, markdown refresh, issue commands, package verification, and diff hygiene.
- Placeholder scan: no `TBD`, `TODO`, or open-ended implementation placeholders are present.
- Type consistency: the plan consistently uses `LaurentGLFactorizationCertificate`, `:laurent_unit_creation`, `ToricBuilderCacheQBlockStatusReport`, `_preferred_laurent_unit_creation_candidate`, and the existing report row fields.
