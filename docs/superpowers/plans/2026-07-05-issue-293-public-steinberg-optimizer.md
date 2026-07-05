# Issue 293 Public Steinberg Optimizer Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Expose an opt-in public Steinberg optimizer for already-produced ordinary-polynomial elementary factor sequences.

**Architecture:** Reuse the existing internal Steinberg certificate, adjacent rewrite pass, commutator rewrite pass, and exact verifier from `src/algorithm/redundancy.jl`. Add only a public wrapper, public verifier alias, richer comparison-summary fields, and focused public/expert tests; do not touch `elementary_factorization(A)`.

**Tech Stack:** Julia, Oscar.jl exact ordinary polynomial matrix factors, Julia `Test`.

## Global Constraints

- `elementary_factorization(A)` must not optimize by default and must not gain an optimization keyword.
- Public optimizer input is an ordinary-polynomial elementary factor sequence.
- `rules = :safe` is the only supported public rule set in this issue.
- The safe rule set contains only #291 adjacent rules and #292 commutator rules.
- Non-elementary matrix sequences must throw a deliberate `ArgumentError`.
- The returned certificate must expose optimized factors, comparison summary, applied rule metadata, PR #179 metric summaries before and after optimization, and exact verification status.
- The public verifier must return `false` when returned optimized factors are tampered.
- Do not add Laurent, ToricBuilder, performance, or global-optimality claims.

---

## File Structure

- Modify `src/Suslin.jl`: export `optimize_elementary_factor_sequence` and `verify_steinberg_optimization_certificate`.
- Modify `src/algorithm/redundancy.jl`: add the public wrapper, public verifier, safe-rule validation, pass-summary metadata, and comparison-summary derived fields.
- Modify `test/public/api_surface.jl`: assert the optimizer and verifier are exported and bound to `Suslin`.
- Modify `test/expert/steinberg_factor_count_optimization.jl`: add public API behavior, negative controls, tamper verification, and unchanged factorization semantics tests.

### Task 1: Add Public Optimizer Tests

**Files:**
- Modify: `test/public/api_surface.jl`
- Modify: `test/expert/steinberg_factor_count_optimization.jl`

**Interfaces:**
- Consumes: desired public functions `optimize_elementary_factor_sequence(factors; rules = :safe)` and `verify_steinberg_optimization_certificate(certificate)::Bool`.
- Produces: failing tests for export surface, public certificate contents, product equality, negative controls, and unchanged `elementary_factorization(A)` semantics.

- [ ] **Step 1: Add API surface assertions**

In `test/public/api_surface.jl`, add these assertions near the other elementary exports:

```julia
@test isdefined(Suslin, :optimize_elementary_factor_sequence)
@test isdefined(Suslin, :verify_steinberg_optimization_certificate)
```

Add these identity assertions near the other `Suslin.foo === foo` checks:

```julia
@test Suslin.optimize_elementary_factor_sequence === optimize_elementary_factor_sequence
@test Suslin.verify_steinberg_optimization_certificate === verify_steinberg_optimization_certificate
```

- [ ] **Step 2: Add public optimizer behavior tests**

Append this testset to `test/expert/steinberg_factor_count_optimization.jl`:

```julia
@testset "public Steinberg elementary factor sequence optimizer" begin
    entries = SteinbergOptimizationFixtureCatalog.cases_by_id()
    entry = entries["steinberg-commutator-forward-qq"]
    original_factors = collect(entry.factors)

    certificate = optimize_elementary_factor_sequence(original_factors)
    summary = certificate.comparison_summary

    @test certificate isa Suslin.SteinbergOptimizationCertificate
    @test verify_steinberg_optimization_certificate(certificate)
    @test certificate.original_factors == original_factors
    @test certificate.optimized_factors == collect(entry.expected_rewrite_factors)
    @test length(certificate.optimized_factors) < length(original_factors)
    @test certificate.original_product == certificate.optimized_product
    @test summary.products_equal
    @test summary.verification_status
    @test summary.original_factor_count == length(original_factors)
    @test summary.optimized_factor_count == length(certificate.optimized_factors)
    @test summary.factor_count_delta ==
          length(certificate.optimized_factors) - length(original_factors)
    @test summary.factor_count == (;
        before = length(original_factors),
        after = length(certificate.optimized_factors),
        delta = length(certificate.optimized_factors) - length(original_factors),
    )
    @test summary.original_metrics.max_elementary_factor_monomial_degree ==
          max_elementary_factor_monomial_degree(original_factors)
    @test summary.optimized_metrics.max_elementary_factor_monomial_degree ==
          max_elementary_factor_monomial_degree(certificate.optimized_factors)
    @test summary.original_metrics.total_elementary_factor_offdiagonal_monomials ==
          total_elementary_factor_offdiagonal_monomials(original_factors)
    @test summary.optimized_metrics.total_elementary_factor_offdiagonal_monomials ==
          total_elementary_factor_offdiagonal_monomials(certificate.optimized_factors)
    @test summary.metric_deltas.max_elementary_factor_monomial_degree ==
          summary.optimized_metrics.max_elementary_factor_monomial_degree -
          summary.original_metrics.max_elementary_factor_monomial_degree
    @test summary.metric_deltas.total_elementary_factor_offdiagonal_monomials ==
          summary.optimized_metrics.total_elementary_factor_offdiagonal_monomials -
          summary.original_metrics.total_elementary_factor_offdiagonal_monomials
    @test summary.applied_rule_names == [:commutator_forward]
    @test length(summary.applied_rewrites) == 1
    @test summary.applied_rewrites[1].rule_name == :safe_steinberg_optimization
    @test summary.applied_rewrites[1].metadata.rules == :safe
    @test length(summary.applied_rewrites[1].metadata.passes) == 1
    @test summary.applied_rewrites[1].metadata.passes[1].pass_name == :commutator
    @test summary.applied_rewrites[1].metadata.passes[1].applied_rewrites[1].rule_name ==
          :commutator_forward

    explicit_safe_certificate = optimize_elementary_factor_sequence(original_factors; rules = :safe)
    @test verify_steinberg_optimization_certificate(explicit_safe_certificate)
    @test explicit_safe_certificate.optimized_factors == certificate.optimized_factors

    R = base_ring(first(original_factors))
    tampered_optimized = copy(certificate.optimized_factors)
    tampered_optimized[1] = elementary_matrix(3, 1, 3, one(R), R)
    tampered_certificate = Suslin.SteinbergOptimizationCertificate(
        certificate.original_factors,
        tampered_optimized,
        certificate.applied_rewrites,
        certificate.comparison_summary,
        certificate.original_product,
        certificate.optimized_product,
        certificate.verification,
    )
    @test !verify_steinberg_optimization_certificate(tampered_certificate)

    bad_factor = identity_matrix(R, 3)
    bad_factor[1, 1] = first(gens(R))
    @test_throws ArgumentError optimize_elementary_factor_sequence([bad_factor])
    @test_throws ArgumentError optimize_elementary_factor_sequence(original_factors; rules = :aggressive)

    A = matrix(R, [
        one(R) + first(gens(R))  one(R)                 zero(R);
        first(gens(R))           one(R)                 zero(R);
        zero(R)                  zero(R)                one(R)
    ])
    factors = elementary_factorization(A)
    @test verify_factorization(A, factors)
end
```

- [ ] **Step 3: Run tests to verify RED**

Run:

```bash
julia --project=. -e 'include("test/public/api_surface.jl")'
julia --project=. -e 'include("test/expert/steinberg_factor_count_optimization.jl")'
```

Expected: both commands fail because `optimize_elementary_factor_sequence` and `verify_steinberg_optimization_certificate` are not defined/exported yet.

- [ ] **Step 4: Commit failing tests**

Run:

```bash
git add test/public/api_surface.jl test/expert/steinberg_factor_count_optimization.jl
git commit -m "test: cover public Steinberg optimizer API"
```

Expected: commit succeeds with only test files changed.

### Task 2: Implement Public Optimizer API

**Files:**
- Modify: `src/Suslin.jl`
- Modify: `src/algorithm/redundancy.jl`

**Interfaces:**
- Consumes: failing tests from Task 1.
- Produces:
  - `optimize_elementary_factor_sequence(factors; rules = :safe)`
  - `verify_steinberg_optimization_certificate(certificate)::Bool`
  - comparison-summary fields `factor_count`, `metric_deltas`, and `applied_rule_names`.

- [ ] **Step 1: Export public functions**

In `src/Suslin.jl`, add these exports near the existing elementary and metric exports:

```julia
export optimize_elementary_factor_sequence
export verify_steinberg_optimization_certificate
```

- [ ] **Step 2: Add comparison-summary derived fields**

In `src/algorithm/redundancy.jl`, add these helpers before `_steinberg_comparison_summary`:

```julia
function _steinberg_metric_deltas(original_metrics, optimized_metrics)
    return (;
        max_elementary_factor_monomial_degree =
            optimized_metrics.max_elementary_factor_monomial_degree -
            original_metrics.max_elementary_factor_monomial_degree,
        total_elementary_factor_offdiagonal_monomials =
            optimized_metrics.total_elementary_factor_offdiagonal_monomials -
            original_metrics.total_elementary_factor_offdiagonal_monomials,
    )
end

function _steinberg_rewrite_rule_names(records)
    names = Symbol[]
    for record in records
        metadata = hasproperty(record, :metadata) ? record.metadata : (;)
        if hasproperty(metadata, :passes)
            for pass in metadata.passes
                for pass_record in pass.applied_rewrites
                    push!(names, pass_record.rule_name)
                end
            end
        else
            push!(names, record.rule_name)
        end
    end
    return names
end
```

Then update `_steinberg_comparison_summary` so it computes `original_metrics` and `optimized_metrics` once and returns these additional fields:

```julia
factor_count = (;
    before = original_factor_count,
    after = optimized_factor_count,
    delta = optimized_factor_count - original_factor_count,
),
metric_deltas = _steinberg_metric_deltas(original_metrics, optimized_metrics),
applied_rule_names = _steinberg_rewrite_rule_names(applied_rewrites),
```

Keep all existing summary fields unchanged.

- [ ] **Step 3: Add safe-rule public wrapper**

In `src/algorithm/redundancy.jl`, add these helpers after `_verify_steinberg_optimization_certificate`:

```julia
function _require_steinberg_optimization_rules(rules)
    rules isa Symbol ||
        throw(ArgumentError("Steinberg optimization rules must be a Symbol"))
    rules == :safe ||
        throw(ArgumentError("unsupported Steinberg optimization rule set; supported rule sets: :safe"))
    return rules
end

function _steinberg_pass_summary(pass_name::Symbol, certificate)
    return (;
        pass_name,
        original_factor_count = length(certificate.original_factors),
        optimized_factor_count = length(certificate.optimized_factors),
        factor_count_delta =
            length(certificate.optimized_factors) - length(certificate.original_factors),
        applied_rewrites = copy(certificate.applied_rewrites),
    )
end

function _empty_steinberg_pass_summary(pass_name::Symbol, factors)
    factor_count = length(factors)
    return (;
        pass_name,
        original_factor_count = factor_count,
        optimized_factor_count = factor_count,
        factor_count_delta = 0,
        applied_rewrites = Any[],
    )
end

function _steinberg_safe_public_rewrite_metadata(original_factors, optimized_factors, pass_summaries)
    active_passes = Tuple(pass for pass in pass_summaries if !isempty(pass.applied_rewrites))
    isempty(active_passes) && return Any[]

    return [(
        rule_name = :safe_steinberg_optimization,
        original_factor_count = length(original_factors),
        optimized_factor_count = length(optimized_factors),
        original_span = (start = 1, stop = length(original_factors)),
        optimized_span = (start = 1, stop = length(optimized_factors)),
        metadata = (;
            rules = :safe,
            passes = active_passes,
        ),
    )]
end

function optimize_elementary_factor_sequence(factors; rules = :safe)
    checked_rules = _require_steinberg_optimization_rules(rules)
    original_context = _steinberg_sequence_context(factors, "original")

    adjacent_certificate =
        _steinberg_adjacent_rewrite_optimization_certificate(original_context.factors)
    adjacent_summary = _steinberg_pass_summary(:adjacent, adjacent_certificate)

    if isempty(adjacent_certificate.optimized_factors)
        commutator_summary =
            _empty_steinberg_pass_summary(:commutator, adjacent_certificate.optimized_factors)
        optimized_factors = adjacent_certificate.optimized_factors
    else
        commutator_certificate =
            _steinberg_commutator_rewrite_optimization_certificate(
                adjacent_certificate.optimized_factors,
            )
        commutator_summary = _steinberg_pass_summary(:commutator, commutator_certificate)
        optimized_factors = commutator_certificate.optimized_factors
    end

    applied_rewrites = _steinberg_safe_public_rewrite_metadata(
        original_context.factors,
        optimized_factors,
        (adjacent_summary, commutator_summary),
    )
    return _steinberg_optimization_certificate(
        original_context.factors,
        optimized_factors,
        applied_rewrites,
    )
end

function verify_steinberg_optimization_certificate(certificate)::Bool
    return _verify_steinberg_optimization_certificate(certificate)
end
```

The `checked_rules` local is intentionally retained so unsupported rule sets are validated before any optimizer work; if Julia warns about an unused local in this codebase, replace it with `_require_steinberg_optimization_rules(rules)`.

- [ ] **Step 4: Run tests to verify GREEN**

Run:

```bash
julia --project=. -e 'include("test/public/api_surface.jl")'
julia --project=. -e 'include("test/expert/steinberg_factor_count_optimization.jl")'
```

Expected: both commands exit 0.

- [ ] **Step 5: Commit implementation**

Run:

```bash
git add src/Suslin.jl src/algorithm/redundancy.jl
git commit -m "feat: expose public Steinberg optimizer"
```

Expected: commit succeeds with only source files changed.

### Task 3: Verify Package and PR-Ready Diff

**Files:**
- Verify: `src/Suslin.jl`
- Verify: `src/algorithm/redundancy.jl`
- Verify: `test/public/api_surface.jl`
- Verify: `test/expert/steinberg_factor_count_optimization.jl`
- Verify: `docs/superpowers/specs/2026-07-05-issue-293-public-steinberg-optimizer-design.md`
- Verify: `docs/superpowers/plans/2026-07-05-issue-293-public-steinberg-optimizer.md`

**Interfaces:**
- Consumes: committed tests and implementation.
- Produces: verified branch ready to push and open as a pull request.

- [ ] **Step 1: Run issue verification commands**

Run:

```bash
julia --project=. -e 'include("test/public/api_surface.jl")'
julia --project=. -e 'include("test/expert/steinberg_factor_count_optimization.jl")'
```

Expected: both commands exit 0.

- [ ] **Step 2: Run package test command**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: command exits 0.

- [ ] **Step 3: Inspect diff and status**

Run:

```bash
git status --short
git log --oneline --decorate -6
git diff origin/main...HEAD --stat
```

Expected: branch contains only the design doc, plan doc, source changes, and tests for issue #293.

- [ ] **Step 4: Commit plan if needed**

If the plan file is still uncommitted, run:

```bash
git add docs/superpowers/plans/2026-07-05-issue-293-public-steinberg-optimizer.md
git commit -m "docs: plan public Steinberg optimizer"
```

Expected: commit succeeds. If the plan was already committed earlier, this step is a no-op with a clean status.
