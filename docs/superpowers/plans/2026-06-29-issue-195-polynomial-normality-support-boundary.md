# Issue 195 Polynomial Normality Support Boundary Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Document and gate the completed ordinary-polynomial normality/conjugation certificate support boundary from issue 181.

**Architecture:** Add a narrow public scope note to the README and Documenter index, then add an expert smoke gate that scans those docs and the test registry. Existing certificate suites remain the behavioral checks; this plan only makes their support boundary explicit and test-registered.

**Tech Stack:** Julia, Oscar, `Test`, existing `test/runtests.jl` group registry, Markdown docs.

## Global Constraints

- Ordinary-polynomial Cohn-type, rank-one, and conjugated-elementary normality/conjugation certificates are supported.
- The ECP induction/normality adapter may be described only as a staged adapter that replays a nested conjugated-elementary certificate.
- Arbitrary Park-Woodburn `SL_n(k[x_1, ..., x_m])` factorization remains staged until issues 182 through 187 are complete.
- Do not claim full issue 187 acceptance, Laurent `GL_n`, ToricBuilder `case_008`, or Steinberg factor-count optimization support as part of this issue.
- Do not implement Murthy local solving, Quillen patching, the general ECP reducer, recursive `SL_n` factorization, full Park-Woodburn public acceptance, Laurent/ToricBuilder support, or Steinberg factor-count optimization.

---

### Task 1: Add the Failing Support-Boundary Gate

**Files:**
- Create: `test/expert/polynomial_normality_support_boundary.jl`
- Modify: `test/runtests.jl`

**Interfaces:**
- Consumes: `README.md`, `docs/src/index.md`, and the literal expert file list in `test/runtests.jl`.
- Produces: `expert/polynomial_normality_support_boundary.jl` registered in the expert test group.

- [ ] **Step 1: Create the failing expert test file**

Create `test/expert/polynomial_normality_support_boundary.jl` with:

```julia
using Test

const POLYNOMIAL_NORMALITY_BOUNDARY_REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

function _normality_boundary_read(path_parts...)
    return read(joinpath(POLYNOMIAL_NORMALITY_BOUNDARY_REPO_ROOT, path_parts...), String)
end

function _normality_boundary_required_phrases()
    return [
        "ordinary-polynomial normality/conjugation certificates",
        "Cohn-type realization certificates",
        "rank-one normality certificates",
        "conjugated-elementary normality certificates",
        "staged ECP induction/normality adapter replays a nested conjugated-elementary certificate",
        "Park-Woodburn `SL_n(k[x_1, ..., x_m])` factorization remains staged",
        "Murthy local solving, general Quillen local realizability, the general ECP reducer",
        "Laurent/ToricBuilder mainline acceptance, and Steinberg factor-count optimization remain staged boundaries",
    ]
end

function _normality_boundary_contains(text::AbstractString, phrase::AbstractString)
    normalized_text = replace(text, r"\s+" => " ")
    normalized_phrase = replace(phrase, r"\s+" => " ")
    return occursin(normalized_phrase, normalized_text)
end

@testset "polynomial normality support boundary documentation" begin
    docs = Dict(
        "README.md" => _normality_boundary_read("README.md"),
        "docs/src/index.md" => _normality_boundary_read("docs", "src", "index.md"),
    )

    for (path, text) in docs
        @testset "$path" begin
            for phrase in _normality_boundary_required_phrases()
                @test _normality_boundary_contains(text, phrase)
            end
            @test !_normality_boundary_contains(text, "full #187 acceptance is supported")
            @test !_normality_boundary_contains(text, "ToricBuilder `case_008` is supported")
            @test !_normality_boundary_contains(text, "Steinberg optimization is supported")
            @test !_normality_boundary_contains(text, "Laurent `GL_n` is fully supported")
        end
    end
end

@testset "polynomial normality certificate expert gate registration" begin
    runtests = _normality_boundary_read("test", "runtests.jl")
    for expert_file in (
        "expert/cohn_type.jl",
        "expert/normality_rank_one.jl",
        "expert/normality.jl",
        "expert/ecp_induction_normality.jl",
    )
        @test occursin("\"$(expert_file)\"", runtests)
    end
    @test occursin("\"expert/polynomial_normality_support_boundary.jl\"", runtests)
end
```

- [ ] **Step 2: Register the new expert test**

In `test/runtests.jl`, add:

```julia
        "expert/polynomial_normality_support_boundary.jl",
```

immediately after:

```julia
        "expert/polynomial_normality_fixtures.jl",
```

- [ ] **Step 3: Verify RED**

Run:

```bash
julia --project=. -e 'include("test/expert/polynomial_normality_support_boundary.jl")'
```

Expected: FAIL because the new support-boundary phrases are not yet present in `README.md` and `docs/src/index.md`.

### Task 2: Document the Public Support Boundary

**Files:**
- Modify: `README.md`
- Modify: `docs/src/index.md`

**Interfaces:**
- Consumes: the support phrases enforced by Task 1.
- Produces: narrow public scope wording for issue 181's completed normality layer.

- [ ] **Step 1: Update README scope wording**

In `README.md`, add this bullet after the `laurent_gl_factorization_certificate(A)` bullet:

```markdown
- The supported ordinary-polynomial normality/conjugation certificates cover
  the completed #181 layer: Cohn-type realization certificates, rank-one
  normality certificates, and conjugated-elementary normality certificates. The
  staged ECP induction/normality adapter replays a nested
  conjugated-elementary certificate for its normality step.
```

Then replace the existing staged-boundary bullet with:

```markdown
- The implementation is not yet the full Park-Woodburn algorithm for arbitrary
  `SL_n(k[x_1, ..., x_m])`, `n >= 3`: Park-Woodburn `SL_n(k[x_1, ..., x_m])`
  factorization remains staged until #182 through #187 are complete. Murthy
  local solving, general Quillen local realizability, the general ECP reducer,
  coefficient-ring support beyond exact field-backed ordinary polynomial
  rings, arbitrary Laurent `GL_n` determinant correction, Laurent/ToricBuilder
  mainline acceptance, and Steinberg factor-count optimization remain staged
  boundaries.
```

- [ ] **Step 2: Update docs index scope wording**

Apply the same support bullet and staged-boundary bullet in `docs/src/index.md`
inside the `## Scope` list.

- [ ] **Step 3: Verify GREEN for the focused gate**

Run:

```bash
julia --project=. -e 'include("test/expert/polynomial_normality_support_boundary.jl")'
```

Expected: PASS.

### Task 3: Run Issue Verification

**Files:**
- No source edits.

**Interfaces:**
- Consumes: all changes from Tasks 1 and 2.
- Produces: verification evidence for the PR body and final Agent Desk result.

- [ ] **Step 1: Run the expert suite**

Run:

```bash
julia --project=. test/runtests.jl expert
```

Expected: exit 0, including the normality certificate suites and the new support-boundary gate.

- [ ] **Step 2: Run the package entry point**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: exit 0 for the default fast `public` and `internal` suites.

- [ ] **Step 3: Scan docs for unsupported claims**

Run:

```bash
rg -n "full #187 acceptance is supported|ToricBuilder `case_008` is supported|Steinberg optimization is supported|Laurent `GL_n` is fully supported" README.md docs/src/index.md
```

Expected: exit 1 with no matches.

## Plan Self-Review

- Every spec requirement maps to a task.
- The plan contains no placeholders.
- The only code change is a test gate; docs copy is exact and narrow.
- Verification includes the two issue-required Julia commands and a negative docs scan.
