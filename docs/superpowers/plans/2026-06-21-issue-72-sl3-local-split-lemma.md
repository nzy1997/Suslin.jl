# Issue 72 SL3 Local Split Lemma Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an internal Murthy split-lemma replay helper for local `SL_3` special-form targets and bridge it into replayable local certificates.

**Architecture:** Extend `src/algorithm/sl3_local.jl` with a structured split replay object, exact verifier, and certificate bridge. Add focused expert coverage that writes the split test first, verifies the RED failure, implements the helper, and then verifies focused and package tests.

**Tech Stack:** Julia, Oscar polynomial rings, Suslin elementary matrices, Test stdlib.

## Global Constraints

- Do not implement recursive Murthy solving.
- Do not implement q-degree normalization.
- Do not add Quillen local-to-global patching.
- Do not change `elementary_factorization` routing.
- Keep helper names internal and unexported.
- Preserve existing `realize_sl3_local(...)` behavior.
- Focused verification command is `julia --project=. -e 'include("test/expert/sl3_local_split_lemma.jl")'`.
- Required package verification command is `julia --project=. -e 'using Pkg; Pkg.test()'`.

---

## File Structure

- Modify `src/algorithm/sl3_local.jl`: add `SL3LocalSplitLemmaReplay`, split construction, split verification, certificate bridge, and `:murthy_split_lemma` certificate replay support.
- Create `test/expert/sl3_local_split_lemma.jl`: focused split-lemma examples, wrapper elementary checks, certificate bridge checks, and negative controls.
- Modify `test/runtests.jl`: register the new expert test.
- Keep `test/fixtures/sl3_murthy_gupta_cases.jl` unchanged and consume its `mg-split-lemma-x-square` fixture id.

---

### Task 1: Write the Split-Lemma Expert Test

**Files:**
- Create: `test/expert/sl3_local_split_lemma.jl`

**Interfaces:**
- Consumes: `Suslin.sl3_local_split_lemma_replay`, `Suslin.verify_sl3_local_split_lemma_replay`, `Suslin.sl3_local_split_lemma_certificate`, and existing `Suslin.realize_sl3_local_certificate`.
- Produces: failing expert coverage for the new helper.

- [x] **Step 1: Add the test file**

Create `test/expert/sl3_local_split_lemma.jl` with helper functions that build special-form targets, multiply factor lists, assert split replay fields, and check that wrapper factors have diagonal entries one, at most one off-diagonal nonzero entry, size `3 x 3`, and the same base ring.

- [x] **Step 2: Add hand-checkable examples**

Add one fixture-backed example using id `mg-split-lemma-x-square`, plus two certifiable examples:

```julia
a = X + 1; a_prime = X + 2; b = 1; d = 1
c = a * a_prime - 1
c1 = a - 1; d1 = 1
c2 = a_prime - 1; d2 = 1
```

and

```julia
a = X; a_prime = X + 1; b = 1; d = 1
c = a * a_prime - 1
c1 = 2 * X - 1; d1 = 2
c2 = 3 * (X + 1) - 1; d2 = 3
```

- [x] **Step 3: Add negative controls**

Assert that corrupting `c1` rejects construction and that swapping child
certificates rejects certificate construction.

- [x] **Step 4: Verify RED**

Run:

```bash
julia --project=. -e 'include("test/expert/sl3_local_split_lemma.jl")'
```

Expected: failure because the new split helper does not exist yet.

### Task 2: Implement Split Replay and Certificate Bridge

**Files:**
- Modify: `src/algorithm/sl3_local.jl`

**Interfaces:**
- Produces: `SL3LocalSplitLemmaReplay`.
- Produces: `sl3_local_split_lemma_replay(...)`.
- Produces: `verify_sl3_local_split_lemma_replay(replay)::Bool`.
- Produces: `sl3_local_split_lemma_certificate(replay, first_child_certificate, second_child_certificate, X)`.
- Extends: `verify_sl3_local_realization` for branch `:murthy_split_lemma`.

- [x] **Step 1: Add split replay construction**

Validate exact relations `a*a_prime*d - b*c == 1`, `a*d1 - b*c1 == 1`, and
`a_prime*d2 - b*c2 == 1`. Build original and child special-form targets.

- [x] **Step 2: Add wrapper factors**

Use the Park-Woodburn sequence:

```julia
E21(c*d1*d2 - d*(c2 + a_prime*c1*d2))
E23(d2 - 1) E32(1) E23(-1)
child1
E23(1) E32(-1) E23(1)
child2
E23(-1) E32(1) E23(a - 1) E31(-a_prime*c1) E32(-d1)
```

Record prefix, middle, suffix, flattened wrapper factors, and exact product.

- [x] **Step 3: Add replay verifier**

Recompute witness relations, child targets, wrapper factors, flattened wrapper
list, and reassembled product from the stored witness. Return `false` for
malformed data.

- [x] **Step 4: Add certificate bridge**

Require both child certificates to verify and match split child targets. Build
a local realization certificate with branch `:murthy_split_lemma` and factors
`prefix ++ child1.factors ++ middle ++ child2.factors ++ suffix`.

- [x] **Step 5: Extend certificate expected factors**

Teach `_sl3_local_certificate_expected_factors`, `_sl3_local_branch_witness_ok`,
and `_sl3_local_witness_keys_ok` to verify the `:murthy_split_lemma` witness and
recompute exact factors.

- [x] **Step 6: Verify GREEN**

Run:

```bash
julia --project=. -e 'include("test/expert/sl3_local_split_lemma.jl")'
```

Expected: pass.

### Task 3: Register and Verify

**Files:**
- Modify: `test/runtests.jl`

**Interfaces:**
- Produces: expert suite registration for `expert/sl3_local_split_lemma.jl`.

- [x] **Step 1: Register the expert test**

Add `"expert/sl3_local_split_lemma.jl"` after `"expert/sl3_local_certificate.jl"`.

- [x] **Step 2: Run required verification**

Run:

```bash
julia --project=. test/runtests.jl expert
julia --project=. -e 'using Pkg; Pkg.test()'
git diff --check
```

Expected: all commands exit 0.

- [x] **Step 3: Commit**

Stage only the files above plus this plan/spec and commit:

```bash
git add docs/superpowers/specs/2026-06-21-issue-72-sl3-local-split-lemma-design.md docs/superpowers/plans/2026-06-21-issue-72-sl3-local-split-lemma.md src/algorithm/sl3_local.jl test/expert/sl3_local_split_lemma.jl test/runtests.jl
git commit -m "Add SL3 local split lemma replay"
```

---

## Self-Review

- Spec coverage: the plan covers construction, exact replay, wrapper elementary checks, issue 69 fixture reuse, issue 70 certificate replay, negative controls, and registration.
- Placeholder scan: no deferred implementation markers remain.
- Type consistency: helper and verifier names match the design document.
