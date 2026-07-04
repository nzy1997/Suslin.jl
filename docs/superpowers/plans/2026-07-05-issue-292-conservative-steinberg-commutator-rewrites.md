# Conservative Steinberg Commutator Rewrites Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an internal certificate-backed optimizer for exact four-factor Steinberg commutator rewrite windows from issue #292.

**Architecture:** Extend `src/algorithm/redundancy.jl` beside the existing adjacent optimizer. The new private helper scans canonical elementary factor records left to right, accepts only exact four-factor commutator windows, checks each accepted window by exact matrix multiplication, and returns the existing `SteinbergOptimizationCertificate`. Expert tests cover the three positive commutator forms, catalog fixtures, and negative controls that must remain unchanged.

**Tech Stack:** Julia, Oscar exact polynomial rings and matrices, existing Suslin canonical factor records, existing Steinberg certificate helpers, `Test`.

## Global Constraints

- Keep the commutator optimizer internal; do not export a public API.
- Accept only ordinary-polynomial elementary factor sequences already accepted by `_steinberg_sequence_context`.
- Implement only exact windows of the forms from Section 6:
  `E_ij(a) E_jl(b) E_ij(-a) E_jl(-b) -> E_il(a*b)`,
  `E_ij(a) E_li(b) E_ij(-a) E_li(-b) -> E_lj(-a*b)`, and the disjoint commutator identity.
- Skip overlapping or ambiguous windows unless they match exactly.
- Before accepting a commutator rewrite, verify the local original and replacement products by exact matrix multiplication.
- Do not search globally for optimal factor count.
- Do not change `elementary_factorization(A)` or optimize public routes by default.

---

## File Structure

- Modify `test/expert/steinberg_factor_count_optimization.jl`: add failing tests for the private commutator optimizer, including fixture-backed positives and negative controls.
- Modify `src/algorithm/redundancy.jl`: add private commutator matching helpers and `_steinberg_commutator_rewrite_optimization_certificate(factors)`.
- Existing files only; no new source module or public export is needed.

### Task 1: Expert Tests For Conservative Commutator Rewrites

**Files:**
- Modify: `test/expert/steinberg_factor_count_optimization.jl`

**Interfaces:**
- Consumes: `elementary_matrix`, existing `SteinbergOptimizationFixtureCatalog`, existing certificate verifier `Suslin._verify_steinberg_optimization_certificate`.
- Produces: failing tests for `Suslin._steinberg_commutator_rewrite_optimization_certificate(factors)`.

- [ ] **Step 1: Add local test helpers**

Append this helper block after the adjacent optimizer testset in `test/expert/steinberg_factor_count_optimization.jl`:

```julia
function _assert_valid_commutator_certificate(
    certificate,
    original_factors,
    expected_optimized_factors,
    expected_rule_names,
)
    @test certificate isa Suslin.SteinbergOptimizationCertificate
    @test Suslin._verify_steinberg_optimization_certificate(certificate)
    @test certificate.original_factors == original_factors
    @test certificate.optimized_factors == expected_optimized_factors
    @test certificate.original_product == certificate.optimized_product
    @test certificate.comparison_summary.products_equal
    @test certificate.comparison_summary.original_factor_count == length(original_factors)
    @test certificate.comparison_summary.optimized_factor_count == length(expected_optimized_factors)
    @test certificate.comparison_summary.factor_count_delta ==
          length(expected_optimized_factors) - length(original_factors)
    @test [rewrite.rule_name for rewrite in certificate.applied_rewrites] == expected_rule_names
    @test all(rewrite -> rewrite.metadata.local_products_equal, certificate.applied_rewrites)
    return certificate
end
```

- [ ] **Step 2: Add positive commutator tests**

Append this testset after the helper:

```julia
@testset "Steinberg conservative commutator optimizer positives" begin
    R, (x, y) = Oscar.polynomial_ring(QQ, ["x", "y"])
    a = x + one(R)
    b = y

    forward_factors = [
        elementary_matrix(3, 1, 2, a, R),
        elementary_matrix(3, 2, 3, b, R),
        elementary_matrix(3, 1, 2, -a, R),
        elementary_matrix(3, 2, 3, -b, R),
    ]
    forward_expected = [
        elementary_matrix(3, 1, 3, a * b, R),
    ]
    forward_certificate =
        Suslin._steinberg_commutator_rewrite_optimization_certificate(forward_factors)
    _assert_valid_commutator_certificate(
        forward_certificate,
        forward_factors,
        forward_expected,
        [:commutator_forward],
    )
    @test forward_certificate.applied_rewrites[1].original_span == (start = 1, stop = 4)
    @test forward_certificate.applied_rewrites[1].optimized_span == (start = 1, stop = 1)
    @test forward_certificate.applied_rewrites[1].metadata.indices == (i = 1, j = 2, l = 3)

    reverse_factors = [
        elementary_matrix(3, 2, 3, x, R),
        elementary_matrix(3, 1, 2, y + one(R), R),
        elementary_matrix(3, 2, 3, -x, R),
        elementary_matrix(3, 1, 2, -(y + one(R)), R),
    ]
    reverse_expected = [
        elementary_matrix(3, 1, 3, -(x * (y + one(R))), R),
    ]
    reverse_certificate =
        Suslin._steinberg_commutator_rewrite_optimization_certificate(reverse_factors)
    _assert_valid_commutator_certificate(
        reverse_certificate,
        reverse_factors,
        reverse_expected,
        [:commutator_reverse],
    )
    @test reverse_certificate.applied_rewrites[1].metadata.indices == (l = 1, i = 2, j = 3)

    disjoint_factors = [
        elementary_matrix(4, 1, 2, x, R),
        elementary_matrix(4, 3, 4, y + one(R), R),
        elementary_matrix(4, 1, 2, -x, R),
        elementary_matrix(4, 3, 4, -(y + one(R)), R),
    ]
    disjoint_certificate =
        Suslin._steinberg_commutator_rewrite_optimization_certificate(disjoint_factors)
    _assert_valid_commutator_certificate(
        disjoint_certificate,
        disjoint_factors,
        typeof(first(disjoint_factors))[],
        [:disjoint_commutator_identity],
    )
    @test disjoint_certificate.applied_rewrites[1].optimized_span == (start = 1, stop = 0)
    @test disjoint_certificate.applied_rewrites[1].metadata.indices == (i = 1, j = 2, l = 3, p = 4)
end
```

- [ ] **Step 3: Add fixture-backed positive checks**

Append this testset after the positive tests:

```julia
@testset "Steinberg commutator optimizer fixture catalog positives" begin
    entries = SteinbergOptimizationFixtureCatalog.cases_by_id()
    for id in (
        "steinberg-commutator-forward-qq",
        "steinberg-commutator-reverse-qq",
        "steinberg-disjoint-commutator-identity-qq",
    )
        entry = entries[id]
        original_factors = collect(entry.factors)
        expected_factors = collect(entry.expected_rewrite_factors)
        certificate =
            Suslin._steinberg_commutator_rewrite_optimization_certificate(original_factors)

        _assert_valid_commutator_certificate(
            certificate,
            original_factors,
            expected_factors,
            [entry.rule_name],
        )
        @test certificate.optimized_product == entry.rewritten_product
        @test certificate.original_product == entry.original_product
    end
end
```

- [ ] **Step 4: Add negative controls**

Append this testset after the fixture checks:

```julia
@testset "Steinberg commutator optimizer negative controls" begin
    R, (x, y) = Oscar.polynomial_ring(QQ, ["x", "y"])
    a = x + one(R)
    b = y + one(R)

    reordered_factors = [
        elementary_matrix(3, 1, 2, a, R),
        elementary_matrix(3, 1, 2, -a, R),
        elementary_matrix(3, 2, 3, b, R),
        elementary_matrix(3, 2, 3, -b, R),
    ]
    reordered_certificate =
        Suslin._steinberg_commutator_rewrite_optimization_certificate(reordered_factors)
    _assert_valid_commutator_certificate(reordered_certificate, reordered_factors, reordered_factors, Symbol[])

    wrong_inverse_factors = [
        elementary_matrix(3, 1, 2, a, R),
        elementary_matrix(3, 2, 3, b, R),
        elementary_matrix(3, 1, 2, -(a + one(R)), R),
        elementary_matrix(3, 2, 3, -b, R),
    ]
    wrong_inverse_certificate =
        Suslin._steinberg_commutator_rewrite_optimization_certificate(wrong_inverse_factors)
    _assert_valid_commutator_certificate(
        wrong_inverse_certificate,
        wrong_inverse_factors,
        wrong_inverse_factors,
        Symbol[],
    )

    invalid_forward_indices = [
        elementary_matrix(3, 1, 2, a, R),
        elementary_matrix(3, 2, 1, b, R),
        elementary_matrix(3, 1, 2, -a, R),
        elementary_matrix(3, 2, 1, -b, R),
    ]
    invalid_forward_certificate =
        Suslin._steinberg_commutator_rewrite_optimization_certificate(invalid_forward_indices)
    _assert_valid_commutator_certificate(
        invalid_forward_certificate,
        invalid_forward_indices,
        invalid_forward_indices,
        Symbol[],
    )

    invalid_disjoint_indices = [
        elementary_matrix(3, 1, 2, a, R),
        elementary_matrix(3, 3, 1, b, R),
        elementary_matrix(3, 1, 2, -a, R),
        elementary_matrix(3, 3, 1, -b, R),
    ]
    invalid_disjoint_certificate =
        Suslin._steinberg_commutator_rewrite_optimization_certificate(invalid_disjoint_indices)
    _assert_valid_commutator_certificate(
        invalid_disjoint_certificate,
        invalid_disjoint_indices,
        invalid_disjoint_indices,
        Symbol[],
    )
end
```

- [ ] **Step 5: Run the focused test to verify RED**

Run:

```bash
julia --project=. -e 'include("test/expert/steinberg_factor_count_optimization.jl")'
```

Expected: FAIL with `UndefVarError: _steinberg_commutator_rewrite_optimization_certificate not defined`, proving the new tests cover the missing helper.

- [ ] **Step 6: Commit the failing tests**

Run:

```bash
git add test/expert/steinberg_factor_count_optimization.jl
git commit -m "test: cover Steinberg commutator rewrites"
```

### Task 2: Private Commutator Optimizer

**Files:**
- Modify: `src/algorithm/redundancy.jl`

**Interfaces:**
- Consumes: `_steinberg_sequence_context(factors, label)`, `_elementary_factor_record_matrix(record)`, `_steinberg_factor_product(factors, R, n)`, `_steinberg_adjacent_rewrite_record(...)`, `_steinberg_optimization_certificate(original_factors, optimized_factors, applied_rewrite_metadata)`.
- Produces: `_steinberg_commutator_rewrite_optimization_certificate(factors)`.

- [ ] **Step 1: Add commutator matching helpers**

Insert this block after `_steinberg_adjacent_rewrite_record` in `src/algorithm/redundancy.jl`:

```julia
function _steinberg_commutator_inverse_tail_matches(first, second, third, fourth)::Bool
    return third.row == first.row &&
           third.col == first.col &&
           third.coefficient == -first.coefficient &&
           fourth.row == second.row &&
           fourth.col == second.col &&
           fourth.coefficient == -second.coefficient
end

function _steinberg_commutator_forward_candidate(first, second)
    i = first.row
    j = first.col
    l = second.col
    first.col == second.row || return nothing
    i != l || return nothing

    return (;
        rule_name = :commutator_forward,
        replacement_records = (;
            kind = :elementary,
            n = first.n,
            ring = first.ring,
            row = i,
            col = l,
            coefficient = first.coefficient * second.coefficient,
        ),
        metadata = (;
            indices = (; i, j, l),
            a = first.coefficient,
            b = second.coefficient,
        ),
    )
end

function _steinberg_commutator_reverse_candidate(first, second)
    l = second.row
    i = first.row
    j = first.col
    second.col == i || return nothing
    j != l || return nothing

    return (;
        rule_name = :commutator_reverse,
        replacement_records = (;
            kind = :elementary,
            n = first.n,
            ring = first.ring,
            row = l,
            col = j,
            coefficient = -(first.coefficient * second.coefficient),
        ),
        metadata = (;
            indices = (; l, i, j),
            a = first.coefficient,
            b = second.coefficient,
        ),
    )
end

function _steinberg_commutator_disjoint_candidate(first, second)
    i = first.row
    j = first.col
    l = second.row
    p = second.col
    i != p || return nothing
    j != l || return nothing

    return (;
        rule_name = :disjoint_commutator_identity,
        replacement_records = (),
        metadata = (;
            indices = (; i, j, l, p),
            a = first.coefficient,
            b = second.coefficient,
        ),
    )
end

function _steinberg_commutator_window_candidate(window_records)
    all(record -> record.kind == :elementary, window_records) || return nothing
    first, second, third, fourth = window_records
    _steinberg_commutator_inverse_tail_matches(first, second, third, fourth) ||
        return nothing

    candidate = _steinberg_commutator_forward_candidate(first, second)
    candidate === nothing || return candidate
    candidate = _steinberg_commutator_reverse_candidate(first, second)
    candidate === nothing || return candidate
    return _steinberg_commutator_disjoint_candidate(first, second)
end
```

- [ ] **Step 2: Normalize candidate replacement records**

Immediately after the helper block from Step 1, add:

```julia
function _steinberg_commutator_replacement_factors(candidate)
    records = candidate.replacement_records isa Tuple &&
              hasproperty(candidate.replacement_records, :kind) ?
        (candidate.replacement_records,) :
        candidate.replacement_records
    return [_elementary_factor_record_matrix(record) for record in records]
end

function _steinberg_commutator_local_products_equal(
    original_window_factors,
    replacement_factors,
    R,
    n::Int,
)::Bool
    return _steinberg_factor_product(original_window_factors, R, n) ==
           _steinberg_factor_product(replacement_factors, R, n)
end
```

- [ ] **Step 3: Add the optimizer**

Insert this function after the Step 2 helpers and before `_steinberg_adjacent_rewrite_optimization_certificate`:

```julia
function _steinberg_commutator_rewrite_optimization_certificate(factors)
    original_context = _steinberg_sequence_context(factors, "original")
    optimized_factors = Any[]
    applied_rewrites = Any[]
    records = original_context.records

    index = 1
    while index <= length(records)
        if index + 3 <= length(records)
            window_records = records[index:(index + 3)]
            candidate = _steinberg_commutator_window_candidate(window_records)
            if candidate !== nothing
                original_window_factors = original_context.factors[index:(index + 3)]
                replacement_factors = _steinberg_commutator_replacement_factors(candidate)
                local_products_equal = _steinberg_commutator_local_products_equal(
                    original_window_factors,
                    replacement_factors,
                    original_context.ring,
                    original_context.n,
                )

                if local_products_equal
                    optimized_start = length(optimized_factors) + 1
                    append!(optimized_factors, replacement_factors)
                    optimized_stop = optimized_start + length(replacement_factors) - 1
                    push!(
                        applied_rewrites,
                        _steinberg_adjacent_rewrite_record(
                            candidate.rule_name,
                            index,
                            index + 3,
                            optimized_start,
                            optimized_stop,
                            merge(candidate.metadata, (; local_products_equal,)),
                        ),
                    )
                    index += 4
                    continue
                end
            end
        end

        push!(optimized_factors, original_context.factors[index])
        index += 1
    end

    return _steinberg_optimization_certificate(
        original_context.factors,
        optimized_factors,
        applied_rewrites,
    )
end
```

- [ ] **Step 4: Run the focused test to verify GREEN**

Run:

```bash
julia --project=. -e 'include("test/expert/steinberg_factor_count_optimization.jl")'
```

Expected: PASS.

- [ ] **Step 5: Commit the implementation**

Run:

```bash
git add src/algorithm/redundancy.jl
git commit -m "feat: add Steinberg commutator rewrites"
```

### Task 3: Verification And PR-Ready Cleanup

**Files:**
- Modify only if a verification failure identifies a bug in `src/algorithm/redundancy.jl` or `test/expert/steinberg_factor_count_optimization.jl`.

**Interfaces:**
- Consumes: the completed test and implementation commits.
- Produces: a verified branch ready for a pull request.

- [ ] **Step 1: Run issue-required focused verification**

Run:

```bash
julia --project=. -e 'include("test/expert/steinberg_factor_count_optimization.jl")'
```

Expected: command exits 0.

- [ ] **Step 2: Run full package verification**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: command exits 0.

- [ ] **Step 3: Inspect branch diff**

Run:

```bash
git status --short
git log --oneline --decorate -5
git diff origin/main...HEAD --stat
```

Expected: branch contains the design, test, and implementation commits; worktree is clean; changed files are limited to the spec, plan, `src/algorithm/redundancy.jl`, and `test/expert/steinberg_factor_count_optimization.jl`.

- [ ] **Step 4: Commit any verification fix if needed**

If Step 1 or Step 2 required a code or test fix, run:

```bash
git add src/algorithm/redundancy.jl test/expert/steinberg_factor_count_optimization.jl
git commit -m "fix: verify Steinberg commutator rewrites"
```

If no fix was needed, do not create an empty commit.

## Self-Review

The plan covers the spec's internal-only API boundary, exact four-factor
matching, local exact product checks, certificate verification, focused
positive cases, fixture-backed positives, negative controls, and both required
verification commands. It does not include a public optimizer, default
optimization, global search, or performance claims.
