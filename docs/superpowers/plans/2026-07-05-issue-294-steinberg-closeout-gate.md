# Issue 294 Steinberg Closeout Gate Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add parent-level #188 acceptance coverage and documentation for the optional Steinberg factor-sequence optimizer.

**Architecture:** Keep the optimizer implementation unchanged. Add acceptance tests around the public #293 API and documentation smoke tests around README, docs index, and a new #188 audit page. The docs state the optimizer is explicit opt-in, exact-product verified, and not part of the default `elementary_factorization(A)` contract.

**Tech Stack:** Julia, Oscar, Suslin public optimizer API, `Test`, Markdown documentation.

## Global Constraints

- Do not add new Park-Woodburn correctness modules.
- Do not broaden supported rings.
- Do not claim performance benchmarks or global optimality.
- Do not change `elementary_factorization(A)` semantics or add optimizer keywords to it.
- The safe rewrite set is exactly `:identity_removal`, `:same_position_merge`, `:inverse_cancellation`, `:commutator_forward`, `:commutator_reverse`, and `:disjoint_commutator_identity`.
- Documentation must cite Park-Woodburn Section 6.
- Documentation must say optimization is optional and opt-in through `optimize_elementary_factor_sequence(factors; rules = :safe)`.
- Documentation must say every optimized sequence is accepted only through exact product verification by `verify_steinberg_optimization_certificate`.
- Documentation must say #188 does not claim global minimum factor counts, Laurent `GL_n` support, or ToricBuilder support.
- Include a before/after table for the `steinberg-commutator-forward-qq` example with factor count `4 -> 1`, factor-count delta `-3`, max monomial degree `1 -> 2`, total off-diagonal monomial count `6 -> 2`, and applied rewrite `:commutator_forward`.

---

## File Structure

- Modify `test/expert/steinberg_factor_count_optimization.jl`: add public safe rewrite-set acceptance coverage over the #288 catalog.
- Modify `test/expert/documentation_smoke.jl`: add #188 documentation contract and overclaim negative-control assertions.
- Modify `README.md`: add the #188 optional optimizer scope bullet.
- Modify `docs/src/index.md`: mirror the README #188 optional optimizer scope bullet.
- Create `docs/audits/2026-07-04-issue-188-steinberg-optimization.md`: parent closeout audit with dependency map, safe rewrite set, before/after table, and non-claims.

### Task 1: Acceptance And Documentation Smoke Tests

**Files:**
- Modify: `test/expert/steinberg_factor_count_optimization.jl`
- Modify: `test/expert/documentation_smoke.jl`

**Interfaces:**
- Consumes: `SteinbergOptimizationFixtureCatalog.cases_by_id()`, `optimize_elementary_factor_sequence(factors; rules = :safe)`, `verify_steinberg_optimization_certificate(certificate)`.
- Produces: tests that fail before the #188 docs/audit content exists and prove the public safe rewrite set.

- [ ] **Step 1: Add the safe rewrite-set acceptance test**

Append this constant near the top of `test/expert/steinberg_factor_count_optimization.jl`, after the fixture include block:

```julia
const ACCEPTED_SAFE_STEINBERG_RULE_NAMES = Set([
    :identity_removal,
    :same_position_merge,
    :inverse_cancellation,
    :commutator_forward,
    :commutator_reverse,
    :disjoint_commutator_identity,
])
```

Append this testset after the existing `"public Steinberg elementary factor sequence optimizer"` testset:

```julia
@testset "public Steinberg safe rewrite set acceptance" begin
    entries = SteinbergOptimizationFixtureCatalog.cases_by_id()
    observed_rule_names = Set{Symbol}()

    for entry in values(entries)
        original_factors = collect(entry.factors)
        expected_factors = collect(entry.expected_rewrite_factors)
        certificate = optimize_elementary_factor_sequence(original_factors; rules = :safe)
        summary = certificate.comparison_summary

        @test verify_steinberg_optimization_certificate(certificate)
        @test certificate.original_factors == original_factors
        @test certificate.optimized_factors == expected_factors
        @test certificate.original_product == entry.original_product
        @test certificate.optimized_product == entry.rewritten_product
        @test certificate.original_product == certificate.optimized_product
        @test summary.products_equal
        @test summary.verification_status
        @test summary.factor_count == (;
            before = length(original_factors),
            after = length(expected_factors),
            delta = length(expected_factors) - length(original_factors),
        )
        @test summary.optimized_factor_count < summary.original_factor_count
        @test summary.applied_rule_names == [entry.rule_name]
        push!(observed_rule_names, entry.rule_name)
    end

    @test observed_rule_names == ACCEPTED_SAFE_STEINBERG_RULE_NAMES

    entry = entries["steinberg-commutator-forward-qq"]
    certificate = optimize_elementary_factor_sequence(collect(entry.factors))
    R = base_ring(first(certificate.optimized_factors))
    n = nrows(first(certificate.optimized_factors))
    corrupted_optimized_factors = copy(certificate.optimized_factors)
    corrupted_optimized_factors[1] = elementary_matrix(n, 1, 3, one(R), R)
    corrupted_certificate = Suslin.SteinbergOptimizationCertificate(
        certificate.original_factors,
        corrupted_optimized_factors,
        certificate.applied_rewrites,
        certificate.comparison_summary,
        certificate.original_product,
        certificate.optimized_product,
        certificate.verification,
    )

    @test !verify_steinberg_optimization_certificate(corrupted_certificate)
end
```

- [ ] **Step 2: Add #188 documentation smoke helpers**

In `test/expert/documentation_smoke.jl`, add this constant after `ISSUE187_ACCEPTANCE_AUDIT_PAGE`:

```julia
const ISSUE188_ACCEPTANCE_AUDIT_PAGE = joinpath(
    @__DIR__,
    "..",
    "..",
    "docs",
    "audits",
    "2026-07-04-issue-188-steinberg-optimization.md",
)
```

Add this reader after `_read_issue187_acceptance_audit()`:

```julia
function _read_issue188_acceptance_audit()
    @test isfile(ISSUE188_ACCEPTANCE_AUDIT_PAGE)
    return read(ISSUE188_ACCEPTANCE_AUDIT_PAGE, String)
end
```

Add these assertion helpers after `_assert_issue187_acceptance_audit(text)`:

```julia
function _assert_issue188_no_overclaims(text)
    for paragraph in _paragraphs(text)
        squashed = _squash_whitespace(paragraph)
        lower_squashed = lowercase(squashed)
        mentions_issue188 =
            occursin("#188", lower_squashed) || occursin("steinberg", lower_squashed)

        if mentions_issue188 && occursin("enabled by default", lower_squashed)
            @test occursin("not enabled by default", lower_squashed) ||
                  occursin("does not optimize by default", lower_squashed) ||
                  occursin("does not enable optimization by default", lower_squashed)
        end

        if mentions_issue188 && (
            occursin("global minimum", lower_squashed) ||
            occursin("globally minimal", lower_squashed) ||
            occursin("global optimal", lower_squashed)
        )
            @test occursin("does not claim", lower_squashed) ||
                  occursin("not claim", lower_squashed) ||
                  occursin("no claim", lower_squashed)
        end

        if mentions_issue188 &&
           (occursin("laurent", lower_squashed) || occursin("toricbuilder", lower_squashed)) &&
           (occursin("support", lower_squashed) || occursin("mainline", lower_squashed))
            @test occursin("does not add", lower_squashed) ||
                  occursin("not add", lower_squashed) ||
                  occursin("does not claim", lower_squashed) ||
                  occursin("not claim", lower_squashed) ||
                  occursin("separate", lower_squashed) ||
                  occursin("out of scope", lower_squashed)
        end
    end
end

function _assert_issue188_optimizer_contract(text)
    squashed = _squash_whitespace(text)
    @test occursin(
        "The optional Steinberg factor-count optimizer (#188) is available only through `optimize_elementary_factor_sequence(factors; rules = :safe)`",
        squashed,
    )
    @test occursin("It is not enabled by default", squashed)
    @test occursin(
        "every optimized sequence is accepted only through exact product verification by `verify_steinberg_optimization_certificate`",
        squashed,
    )
    for rule_name in (
        ":identity_removal",
        ":same_position_merge",
        ":inverse_cancellation",
        ":commutator_forward",
        ":commutator_reverse",
        ":disjoint_commutator_identity",
    )
        @test occursin(rule_name, squashed)
    end
    @test occursin(
        "#188 does not change the correctness contract of `elementary_factorization(A)`",
        squashed,
    )
    @test occursin("does not claim global minimum factor counts", squashed)
    @test occursin("does not add Laurent `GL_n` or ToricBuilder support", squashed)
    _assert_issue188_no_overclaims(text)
end

function _assert_issue188_acceptance_audit(text)
    squashed = _squash_whitespace(text)
    for issue_id in ("#288", "#289", "#290", "#291", "#292", "#293")
        @test occursin(issue_id, squashed)
    end
    @test occursin("Park-Woodburn Section 6", squashed)
    @test occursin("| Metric | Original | Optimized | Delta |", text)
    @test occursin("| Factor count | 4 | 1 | -3 |", text)
    @test occursin("| Max monomial degree | 1 | 2 | 1 |", text)
    @test occursin("| Total off-diagonal monomial count | 6 | 2 | -4 |", text)
    @test occursin("| Applied rewrites | `:commutator_forward` | `:commutator_forward` | accepted safe rewrite |", text)
    @test occursin("products_equal = true", squashed)
    @test occursin("verification_status = true", squashed)
    _assert_issue188_optimizer_contract(text)
end
```

- [ ] **Step 3: Wire the new smoke assertions into the existing testset**

Inside the `"support boundary evidence page"` nested testset, after the existing `audit = _read_issue187_acceptance_audit()` checks, add:

```julia
        issue188_audit = _read_issue188_acceptance_audit()
        _assert_issue188_acceptance_audit(issue188_audit)
```

Inside the `"ordinary-polynomial Park-Woodburn public contract"` nested testset, after the #187 README/docs index assertions, add:

```julia
        _assert_issue188_optimizer_contract(_read_repo_text(README_PATH))
        _assert_issue188_optimizer_contract(_read_repo_text(DOCS_INDEX_PATH))
```

- [ ] **Step 4: Run focused tests and capture red state**

Run:

```bash
julia --project=. -e 'include("test/expert/steinberg_factor_count_optimization.jl")'
julia --project=. -e 'include("test/expert/documentation_smoke.jl")'
```

Expected: optimizer acceptance coverage may already pass because #293 is merged; documentation smoke must fail because the #188 audit page and README/docs wording do not exist yet.

- [ ] **Step 5: Commit Task 1**

```bash
git add test/expert/steinberg_factor_count_optimization.jl test/expert/documentation_smoke.jl
git commit -m "test: cover Steinberg closeout acceptance"
```

### Task 2: Public Docs And Closeout Audit

**Files:**
- Modify: `README.md`
- Modify: `docs/src/index.md`
- Create: `docs/audits/2026-07-04-issue-188-steinberg-optimization.md`

**Interfaces:**
- Consumes: Task 1 documentation smoke assertions.
- Produces: public opt-in optimizer contract and parent #188 closeout audit.

- [ ] **Step 1: Add the #188 README scope bullet**

In `README.md`, insert this bullet in `## Current scope` after the `verify_factorization(A, factors)` bullet:

```markdown
- The optional Steinberg factor-count optimizer (#188) is available only
  through `optimize_elementary_factor_sequence(factors; rules = :safe)`. It is
  not enabled by default, and #188 does not change the correctness contract of
  `elementary_factorization(A)`: public factorization still returns the
  evidence-backed sequence for `A`, and callers may separately optimize that
  sequence only by making the opt-in call. The safe rewrite set is the
  Park-Woodburn Section 6 subset covered by the #288 catalog:
  `:identity_removal`, `:same_position_merge`, `:inverse_cancellation`,
  `:commutator_forward`, `:commutator_reverse`, and
  `:disjoint_commutator_identity`. Every optimized sequence is accepted only
  through exact product verification by
  `verify_steinberg_optimization_certificate`; #188 does not claim global
  minimum factor counts and does not add Laurent `GL_n` or ToricBuilder
  support. The closeout audit is
  `docs/audits/2026-07-04-issue-188-steinberg-optimization.md`.
```

- [ ] **Step 2: Mirror the #188 docs index scope bullet**

In `docs/src/index.md`, insert the same bullet in `## Scope` after the `verify_factorization(A, factors)` bullet.

- [ ] **Step 3: Create the #188 audit page**

Create `docs/audits/2026-07-04-issue-188-steinberg-optimization.md` with:

~~~markdown
# Issue 188 Steinberg Optimization Acceptance Audit

Date: 2026-07-04

Issue #294 is the parent-level closeout gate for #188. It records that the
optional Steinberg factor-count optimizer is accepted as an opt-in factor
sequence rewrite pass only. The optional Steinberg factor-count optimizer
(#188) is available only through
`optimize_elementary_factor_sequence(factors; rules = :safe)`. It is not
enabled by default, and #188 does not change the correctness contract of
`elementary_factorization(A)`.

The accepted rewrite set cites Park-Woodburn Section 6, "Eliminating
Redundancies", which lists the Steinberg relations for identity elementary
factors, same-position products, the two nontrivial commutators, and the
disjoint commutator identity.

## Upstream Layer Map

| Parent layer | Upstream closeout | Evidence boundary | #188 consumers |
| --- | --- | --- | --- |
| #288 fixture catalog | #288 | Ordinary-polynomial Section 6 Steinberg positive and negative catalog entries over exact field-backed rings | `steinberg-identity-removal-qq`, `steinberg-same-position-merge-qq`, `steinberg-inverse-cancellation-qq`, `steinberg-commutator-forward-qq`, `steinberg-commutator-reverse-qq`, and `steinberg-disjoint-commutator-identity-qq` |
| #289 canonical factors | #289 | Private canonical elementary-factor records for identity and single off-diagonal elementary matrices | #290 certificate validation and #291/#292 rewrite matching |
| #290 certificate replay | #290 | `SteinbergOptimizationCertificate` with exact original and optimized products, before/after counts, metrics, and verifier status | #291/#292 private rewrite certificates and #293 public optimizer certificates |
| #291 adjacent rewrites | #291 | Exact adjacent identity removal, same-position merge, and inverse cancellation | #293 `rules = :safe` adjacent pass |
| #292 commutator rewrites | #292 | Exact four-factor forward, reverse, and disjoint commutator windows with local product checks | #293 `rules = :safe` commutator pass |
| #293 public optimizer | #293 | Exported `optimize_elementary_factor_sequence` and `verify_steinberg_optimization_certificate` APIs | #294 closeout documentation and acceptance smoke |

## Accepted Safe Rewrite Set

The public safe rule set is exactly:

- `:identity_removal`
- `:same_position_merge`
- `:inverse_cancellation`
- `:commutator_forward`
- `:commutator_reverse`
- `:disjoint_commutator_identity`

Every optimized sequence is accepted only through exact product verification by
`verify_steinberg_optimization_certificate`. A certificate is accepted only
when the original and optimized products replay exactly, the summary agrees
with the replayed products and metrics, and `verification_status = true`.

## Documented Example

The compact example is the #288 `steinberg-commutator-forward-qq` sequence:

```julia
certificate = optimize_elementary_factor_sequence(factors; rules = :safe)
verify_steinberg_optimization_certificate(certificate)
```

| Metric | Original | Optimized | Delta |
| --- | --- | --- | --- |
| Factor count | 4 | 1 | -3 |
| Max monomial degree | 1 | 2 | 1 |
| Total off-diagonal monomial count | 6 | 2 | -4 |
| Applied rewrites | `:commutator_forward` | `:commutator_forward` | accepted safe rewrite |

For this documented example, `products_equal = true` and
`verification_status = true`.

## Closeout Note

#188 does not change the correctness contract of `elementary_factorization(A)`.
The public factorization call still returns the evidence-backed sequence for
`A` and remains accepted by `verify_factorization(A, factors)`. Users who want
factor-count cleanup must call `optimize_elementary_factor_sequence` explicitly
on the returned or supplied factor sequence and verify the returned certificate.

## Non-Claims

- #188 does not claim global minimum factor counts.
- #188 does not add Laurent `GL_n` or ToricBuilder support.
- Laurent/ToricBuilder mainline support remains separate from #188.
- This audit does not add algorithmic support beyond the landed #288-#293
  ordinary-polynomial Steinberg optimizer chain.

## Verification

The closeout gate is the focused optimizer/docs checks plus the full suite:

```bash
julia --project=. -e 'include("test/expert/steinberg_factor_count_optimization.jl")'
julia --project=. -e 'include("test/public/api_surface.jl")'
julia --project=. -e 'include("test/expert/documentation_smoke.jl")'
julia --project=. test/runtests.jl all
julia --project=. -e 'using Pkg; Pkg.test()'
```
~~~

- [ ] **Step 4: Run documentation smoke**

Run:

```bash
julia --project=. -e 'include("test/expert/documentation_smoke.jl")'
```

Expected: PASS.

- [ ] **Step 5: Commit Task 2**

```bash
git add README.md docs/src/index.md docs/audits/2026-07-04-issue-188-steinberg-optimization.md
git commit -m "docs: document Steinberg optimizer closeout"
```

### Task 3: Verification And Branch Review

**Files:**
- Modify only if verification exposes a defect in Task 1 or Task 2 files.

**Interfaces:**
- Consumes: Task 1 tests and Task 2 docs.
- Produces: clean verification evidence for the PR.

- [ ] **Step 1: Run issue-required focused commands**

```bash
julia --project=. -e 'include("test/expert/steinberg_factor_count_optimization.jl")'
julia --project=. -e 'include("test/public/api_surface.jl")'
julia --project=. -e 'include("test/expert/documentation_smoke.jl")'
```

Expected: all exit 0.

- [ ] **Step 2: Run full suite commands**

```bash
julia --project=. test/runtests.jl all
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: both exit 0.

- [ ] **Step 3: Review branch diff**

```bash
git diff --stat origin/main..HEAD
git diff --check origin/main..HEAD
git status --short
```

Expected: changed files are limited to the Superpowers spec/plan, README/docs/audit, and the two expert tests; `git diff --check` exits 0; working tree is clean after any final commit.

- [ ] **Step 4: Commit any verification fixes**

If Step 1 or Step 2 required edits, commit only those focused fixes:

```bash
git add README.md docs/src/index.md docs/audits/2026-07-04-issue-188-steinberg-optimization.md test/expert/documentation_smoke.jl test/expert/steinberg_factor_count_optimization.jl
git commit -m "fix: complete Steinberg closeout gate"
```
