# Issue 138 Laurent/ToricBuilder Evidence Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add one audit evidence page that reconciles #129 and #131 without widening the supported algorithm scope.

**Architecture:** Keep the evidence as a hand-written audit page under `docs/audits/` and pin its required claims with the existing expert documentation smoke test. The generated ToricBuilder Q-block report remains the raw route/timing report, while the new page summarizes support boundaries and verification commands.

**Tech Stack:** Julia, Test stdlib, Markdown documentation.

## Global Constraints

- Repository has no `AGENTS.md`, `CLAUDE.md`, `GEMINI.md`, or `CONVENTIONS.md`; follow `README.md` test commands and existing Suslin style.
- Add exactly one delivery evidence markdown page under `docs/audits/`.
- Do not add new algorithm support.
- Do not open or close #129 or #131.
- The page must distinguish `elementary_factorization` factor sequences from Laurent `GL_n` certificates.
- The page must not claim arbitrary Park-Woodburn support.
- The page must not claim arbitrary Laurent `GL_n` support.
- The page must name the supported `case_010` outcome.
- The page must name the bounded `case_008` outcome.
- The page must name the Laurent `GL_n` certificate route.
- The page must name the remaining staged boundary for original-input `elementary_factorization`.
- Required issue command: `julia --project=. test/runtests.jl`.
- Required issue command: `julia --project=. test/runtests.jl expert`.
- Required Agent Desk package command: `julia --project=. -e 'using Pkg; Pkg.test()'`.

---

## File Structure

- Create `docs/audits/2026-06-26-laurent-toricbuilder-support-boundary-evidence.md`: the audit evidence page with support matrix, route diagram, verification commands, and parent-issue reconciliation notes.
- Modify `test/expert/documentation_smoke.jl`: add a documentation smoke test that reads the evidence page and pins the required support-matrix rows and conservative scope text.
- Keep `docs/superpowers/specs/2026-06-26-issue-138-laurent-toricbuilder-evidence-design.md` and this plan in the PR.

---

### Task 1: Add Support-Boundary Evidence Page And Smoke Test

**Files:**
- Create: `docs/audits/2026-06-26-laurent-toricbuilder-support-boundary-evidence.md`
- Modify: `test/expert/documentation_smoke.jl`

**Interfaces:**
- Consumes: existing expert test runner entry `test/expert/documentation_smoke.jl`.
- Produces: `SUPPORT_BOUNDARY_EVIDENCE_PAGE`, `_read_support_boundary_evidence()`, and a smoke test that fails if required support-matrix rows or conservative boundary phrases are removed.

- [ ] **Step 1: Write the failing documentation smoke test**

In `test/expert/documentation_smoke.jl`, add this constant and helper after the
existing imports:

```julia
const SUPPORT_BOUNDARY_EVIDENCE_PAGE = joinpath(
    @__DIR__,
    "..",
    "..",
    "docs",
    "audits",
    "2026-06-26-laurent-toricbuilder-support-boundary-evidence.md",
)

function _read_support_boundary_evidence()
    @test isfile(SUPPORT_BOUNDARY_EVIDENCE_PAGE)
    return read(SUPPORT_BOUNDARY_EVIDENCE_PAGE, String)
end
```

Add this nested testset at the end of the existing `"documentation smoke"`
testset:

```julia
    @testset "support boundary evidence page" begin
        evidence = _read_support_boundary_evidence()
        @test occursin(
            "| `case_010` ToricBuilder Q-block | `laurent_gl_factorization_certificate` | `gl_certificate_pass`; verified `true`; decomposed base matrices `48` | public `elementary_factorization` remains `staged_boundary` for the original Laurent `GL_n` input |",
            evidence,
        )
        @test occursin(
            "| `case_008` bounded exercise | bounded Laurent `GL_n` certificate route | `certified_algorithm_boundary` at `certificate_construction` under explicit `--exercise=case_008 --timeout-seconds=120` | not a default report pass; remains a staged algorithm boundary |",
            evidence,
        )
        @test occursin(
            "ToricBuilder Q-block -> classify Laurent determinant -> normalize Laurent GL_n determinant -> factor determinant-one core -> verify Laurent GL_n certificate",
            evidence,
        )
        @test occursin(
            "elementary_factorization(A) -> exact elementary factor sequence -> verify_factorization(A, factors)",
            evidence,
        )
        @test occursin(
            "original-input `elementary_factorization` for Laurent `GL_n` remains a `staged boundary`",
            evidence,
        )
        @test occursin("not arbitrary Park-Woodburn", evidence)
        @test occursin("not arbitrary Laurent `GL_n`", evidence)
        @test occursin("julia --project=. test/runtests.jl", evidence)
        @test occursin("julia --project=. test/runtests.jl expert", evidence)
    end
```

- [ ] **Step 2: Run the focused test to verify RED**

Run:

```bash
julia --project=. -e 'include("test/expert/documentation_smoke.jl")'
```

Expected: FAIL because
`docs/audits/2026-06-26-laurent-toricbuilder-support-boundary-evidence.md` does
not exist yet. If the failure is a syntax error, fix the test and rerun until
the failure proves the missing evidence page.

- [ ] **Step 3: Create the evidence page**

Create `docs/audits/2026-06-26-laurent-toricbuilder-support-boundary-evidence.md`
with this content:

```markdown
# Laurent/ToricBuilder Support-Boundary Evidence

Date: 2026-06-26

This page reconciles the support-scope audit in #129 with the ToricBuilder
Q-block status in #131 after the `case_010` Laurent certificate route (#135)
and the bounded `case_008` triage (#137). It is evidence for review only; it
does not close the parent issues and does not add algorithm support.

## Support Matrix

| Evidence item | Route | Supported outcome | Boundary |
| --- | --- | --- | --- |
| `case_010` ToricBuilder Q-block | `laurent_gl_factorization_certificate` | `gl_certificate_pass`; verified `true`; decomposed base matrices `48` | public `elementary_factorization` remains `staged_boundary` for the original Laurent `GL_n` input |
| `case_008` bounded exercise | bounded Laurent `GL_n` certificate route | `certified_algorithm_boundary` at `certificate_construction` under explicit `--exercise=case_008 --timeout-seconds=120` | not a default report pass; remains a staged algorithm boundary |
| Default ToricBuilder Q-block report rows `case_001`-`case_006`, `case_010` | generated Q-block status report | `gl_certificate_pass` with public `staged_boundary` and Laurent monomial-unit determinants | evidence is the Laurent `GL_n` certificate route, not original-input elementary factor sequences |
| Laurent monomial-unit `GL_n` inputs in the staged certificate path | `normalize_laurent_gl_matrix` then `laurent_gl_factorization_certificate` | certificate verifies the normalized determinant-one core and exact reconstruction metadata | original-input `elementary_factorization` for Laurent `GL_n` remains a `staged boundary` |
| Ordinary-polynomial staged slices | `elementary_factorization` | exact elementary factor sequences that satisfy `verify_factorization(A, factors) == true` | not arbitrary Park-Woodburn support |
| Remaining Laurent scope | no broad public factor-sequence route | certificate-backed monomial-unit slices only where recorded tests exercise them | not arbitrary Laurent `GL_n` support |

## Route Diagram

`elementary_factorization(A) -> exact elementary factor sequence -> verify_factorization(A, factors)`

This is the public factor-sequence route for the supported ordinary-polynomial
and determinant-one Laurent `SL` staged slices. It returns factors only when the
current staged implementation can verify exact multiplication back to `A`.

`ToricBuilder Q-block -> classify Laurent determinant -> normalize Laurent GL_n determinant -> factor determinant-one core -> verify Laurent GL_n certificate`

This is the Laurent `GL_n` certificate route. It records determinant
normalization, factors the normalized determinant-one core, and verifies the
certificate metadata. It is not the same as returning an elementary factor
sequence for the original Laurent `GL_n` input.

## Verification Commands

Run the issue-required documentation and algorithm checks:

```bash
julia --project=. test/runtests.jl
julia --project=. test/runtests.jl expert
```

Run the package entry point used by Agent Desk workers:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Optional route evidence commands:

```bash
julia --project=. scripts/report_toricbuilder_cache_q_blocks.jl --exercise=case_010 --output=/tmp/case010-q-block-status.md
julia --project=. scripts/report_toricbuilder_cache_q_blocks.jl --exercise=case_008 --timeout-seconds=120 --output=/tmp/qblock-case008.md
```

Expected `case_010` evidence: the generated row contains
`gl_certificate_pass`, `verified true`, public `staged_boundary`, and a positive
decomposed base-matrix count.

Expected `case_008` evidence: the bounded generated row is structured, not
`not_exercised_in_default_report`, not `route_error`, and not an unstructured
timeout. The current observed bounded outcome is `certified_algorithm_boundary`
at `certificate_construction`.

## Parent-Issue Reconciliation

- #129 remains the support-scope audit. This page confirms the staged support
  boundary: exact factor sequences are available only on the supported
  `elementary_factorization` slices, while Laurent `GL_n` monomial-unit support
  is certificate evidence unless a later issue adds original-input factor
  sequences.
- #131 remains the ToricBuilder Q-block status thread. The default generated
  report now records `case_010` as a Laurent `GL_n` certificate pass, while
  `case_008` is documented through the explicit bounded exercise route.
- #135 is reflected as the supported `case_010` certificate outcome.
- #137 is reflected as the bounded `case_008` structured boundary outcome.

## Explicit Non-Claims

- This is not arbitrary Park-Woodburn `SL_n(k[x_1, ..., x_m])` support.
- This is not arbitrary Laurent `GL_n` support.
- This is not a claim that `elementary_factorization` returns factor sequences
  for original Laurent `GL_n` inputs with monomial-unit determinant.
- This is not a performance claim or benchmark table.
```

- [ ] **Step 4: Run the focused test to verify GREEN**

Run:

```bash
julia --project=. -e 'include("test/expert/documentation_smoke.jl")'
```

Expected: PASS, including the support-boundary evidence page checks.

- [ ] **Step 5: Commit Task 1**

Run:

```bash
git add docs/audits/2026-06-26-laurent-toricbuilder-support-boundary-evidence.md test/expert/documentation_smoke.jl docs/superpowers/plans/2026-06-26-issue-138-laurent-toricbuilder-evidence.md
git commit -m "docs: publish laurent toricbuilder evidence"
```

---

## Final Verification

After the task review is clean, run:

```bash
julia --project=. test/runtests.jl
julia --project=. test/runtests.jl expert
julia --project=. -e 'using Pkg; Pkg.test()'
git diff --check origin/main
```

The default and expert test commands must pass. The expert command must include
the documentation smoke test that pins `case_010`, `case_008`, the Laurent
`GL_n` certificate route, and the original-input `elementary_factorization`
staged boundary.

## Plan Self-Review

- Spec coverage: Task 1 creates the single evidence page, support matrix, route
  diagram, exact verification commands, parent-issue reconciliation notes, and
  documentation smoke coverage required by #138.
- Placeholder scan: no placeholder markers, deferred implementation, or copied
  step references remain.
- Type consistency: task names and test helpers consistently use
  `SUPPORT_BOUNDARY_EVIDENCE_PAGE`, `_read_support_boundary_evidence`, and the
  exact evidence-page path.
