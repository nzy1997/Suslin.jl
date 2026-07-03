# Issue 273 Ordinary-Polynomial Public Contract Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Publish the final ordinary-polynomial Park-Woodburn #187 public contract in README and Documenter without broadening implementation behavior.

**Architecture:** Add documentation-smoke assertions that read README and Documenter scope text directly, then update both public docs with the same narrow supported/out-of-scope language. Keep the existing runnable public example and do not touch factorization code or fixture behavior.

**Tech Stack:** Julia, Oscar, `Test`, Documenter.jl, Markdown documentation.

## Global Constraints

- Do not change factorization behavior, fixtures, catalog entries, or public APIs.
- Supported #187 scope is exact field-backed ordinary polynomial rings only.
- Supported #187 inputs are determinant-one `SL_3` and `SL_n`, `n > 3`, through the implemented evidence-backed route only.
- Unsupported coefficient rings, arbitrary Laurent `GL_n`, ToricBuilder mainline acceptance, and Steinberg factor-count optimization (#188) remain separate from #187.
- README and `docs/src/index.md` must not list full ordinary-polynomial Park-Woodburn acceptance (#187) as staged after #271/#272.
- Keep the README-style example runnable through `elementary_factorization(A)` and `verify_factorization(A, factors)`.

---

### Task 1: Documentation Contract Smoke And Public Scope Wording

**Files:**
- Modify: `test/expert/documentation_smoke.jl`
- Modify: `README.md`
- Modify: `docs/src/index.md`

**Interfaces:**
- Consumes: public docs text and the existing README-style ordinary-polynomial example.
- Produces: smoke assertions that fail on stale #187-staged wording or overclaimed Laurent/ToricBuilder/#188 wording.

- [ ] **Step 1: Write the failing documentation smoke assertions**

In `test/expert/documentation_smoke.jl`, add repository text paths after `SUPPORT_BOUNDARY_EVIDENCE_PAGE`:

```julia
const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const README_PATH = joinpath(REPO_ROOT, "README.md")
const DOCS_INDEX_PATH = joinpath(REPO_ROOT, "docs", "src", "index.md")
```

Add these helper functions before the top-level `@testset`:

```julia
function _read_repo_text(path)
    @test isfile(path)
    return read(path, String)
end

function _paragraphs(text)
    return split(replace(text, "\r\n" => "\n"), "\n\n")
end

function _assert_not_claimed_as_issue187(text, item)
    for paragraph in _paragraphs(text)
        if occursin("#187", paragraph) && occursin(item, paragraph)
            @test occursin("separate", paragraph) ||
                  occursin("out of scope", paragraph) ||
                  occursin("not part of", paragraph) ||
                  occursin("outside", lowercase(paragraph))
        end
    end
end

function _assert_issue187_public_contract(text)
    @test occursin(
        "The final ordinary-polynomial Park-Woodburn public contract (#187) is supported",
        text,
    )
    @test occursin(
        "exact field-backed ordinary-polynomial determinant-one `SL_3`",
        text,
    )
    @test occursin(
        "exact field-backed ordinary-polynomial determinant-one `SL_n`, `n > 3`",
        text,
    )
    @test occursin("implemented evidence-backed route", text)
    @test occursin(
        "Unsupported coefficient rings remain out of scope",
        text,
    )
    @test occursin(
        "arbitrary Laurent `GL_n`, ToricBuilder mainline acceptance, and Steinberg factor-count optimization (#188) remain separate from #187",
        text,
    )
    @test !occursin("full public Park-Woodburn acceptance (#187)", text)
    _assert_not_claimed_as_issue187(text, "unsupported coefficient rings")
    _assert_not_claimed_as_issue187(text, "arbitrary Laurent `GL_n`")
    _assert_not_claimed_as_issue187(text, "ToricBuilder")
    _assert_not_claimed_as_issue187(text, "factor-count optimization")
end
```

Inside the `"documentation smoke"` testset, after the existing README-style
example verification, add:

```julia
    @testset "ordinary-polynomial Park-Woodburn public contract" begin
        _assert_issue187_public_contract(_read_repo_text(README_PATH))
        _assert_issue187_public_contract(_read_repo_text(DOCS_INDEX_PATH))
    end
```

- [ ] **Step 2: Run the red smoke command**

Run:

```bash
julia --project=. -e 'include("test/expert/documentation_smoke.jl")'
```

Expected: FAIL because current README/docs still say `full public Park-Woodburn acceptance (#187)` remains staged and do not contain the new supported #187 contract phrase.

- [ ] **Step 3: Update README scope wording**

In `README.md`, keep the existing example unchanged. In `## Current scope`, insert this bullet immediately after the recursive #186 bullet:

```markdown
- The final ordinary-polynomial Park-Woodburn public contract (#187) is supported
  for exact field-backed ordinary-polynomial determinant-one `SL_3` and exact
  field-backed ordinary-polynomial determinant-one `SL_n`, `n > 3`, inputs
  through the implemented evidence-backed route. Public acceptance requires the
  factors returned by `elementary_factorization(A)` to satisfy
  `verify_factorization(A, factors)`, and the route certificate must replay the
  required #184 `SL_3` evidence plus #185/#186 ECP-backed recursive provenance
  for larger sizes. The example above is the README-style ordinary-polynomial
  public call shape covered by the #187 acceptance tests.
```

Then replace the final staged-boundary bullet with:

```markdown
- Staged ordinary-polynomial inputs include determinant-one matrices missing the
  required local-form, variable-change, normality/conjugation, Murthy, Quillen,
  ECP, or final `SL_3` evidence path. Those inputs fail before public factors
  are returned, with stable reason codes such as `:missing_ecp_evidence` and
  `:missing_final_sl3_route` for recursive staged failures. Unsupported
  coefficient rings remain out of scope. Arbitrary Laurent `GL_n`, ToricBuilder
  mainline acceptance, and Steinberg factor-count optimization (#188) remain
  separate from #187.
```

- [ ] **Step 4: Update Documenter scope wording**

Make the same two Markdown changes in `docs/src/index.md`, in its `## Scope`
section.

- [ ] **Step 5: Run the green focused smoke command**

Run:

```bash
julia --project=. -e 'include("test/expert/documentation_smoke.jl")'
```

Expected: PASS.

- [ ] **Step 6: Commit Task 1**

Run:

```bash
git add README.md docs/src/index.md test/expert/documentation_smoke.jl
git commit -m "docs: publish Park-Woodburn public contract"
```

### Task 2: Verification, Review, And Pull Request

**Files:**
- No additional source files expected.

**Interfaces:**
- Consumes: committed design, plan, docs, and smoke-test changes.
- Produces: verified worker branch pushed to a pull request against `main`.

- [ ] **Step 1: Run issue-required documentation verification**

Run:

```bash
julia --project=. -e 'include("test/expert/documentation_smoke.jl")'
julia --project=docs docs/make.jl
```

Expected: both commands exit 0.

- [ ] **Step 2: Run required package verification**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: command exits 0.

- [ ] **Step 3: Final review**

Dispatch a code reviewer subagent with the branch diff and these requirements:
README and Documenter must publish #187 only for exact field-backed
ordinary-polynomial determinant-one `SL_3` and `SL_n`, `n > 3`, through the
implemented evidence-backed route, while keeping unsupported coefficient rings,
arbitrary Laurent `GL_n`, ToricBuilder mainline acceptance, and #188
factor-count optimization separate.

- [ ] **Step 4: Finish branch**

Use `superpowers:finishing-a-development-branch`, choose "Push and create a Pull
Request", push `agent/issue-273-publish-the-full-ordinary-polynomial-park-woodbu-run-1`,
and create a PR against `main` with `Closes #273` in the body.
