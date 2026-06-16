# Suslin Stability Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a Julia package that constructs explicit elementary-matrix factorizations for matrices in `SL_n(k[x_1, ..., x_m])` following Park-Woodburn's algorithmic proof of Suslin's stability theorem over polynomial rings.

**Architecture:** Use a layered design that separates exact algebra primitives, unimodular-column reduction, Quillen induction / patching, and the local `SL_3` realization algorithm. Rely on Oscar.jl, Singular.jl, and AbstractAlgebra.jl for polynomial-ring, ideal, Gröbner, module, normal-form, and resultant operations, while keeping Suslin-specific algorithm steps in small Julia modules with explicit invariants and focused tests.

**Tech Stack:** Julia 1.x, Oscar.jl, Singular.jl, AbstractAlgebra.jl, Test stdlib, Documenter.jl

---

## Scope

This plan is for a first serious implementation of the Park-Woodburn algorithm, not a polished end-user CAS replacement. The plan intentionally stages the work:

1. establish exact polynomial/algebra utilities,
2. implement constructive row/column and normality subroutines,
3. implement Quillen induction and local patching,
4. implement the local `SL_3` realization path,
5. compose these into a full elementary factorization API.

The plan avoids two traps:

- trying to implement the full theorem directly in one file,
- overcommitting to abstract APIs before verifying each algebraic step on small exact examples.

## Why Oscar/Singular

The package should use `Oscar.jl` and `Singular.jl` early instead of treating them as optional later additions.

Reasons:

- Park-Woodburn explicitly depends on effective ideal membership, Gröbner-style reductions, resultants, and constructive local-to-global glueing.
- The local reference repo `/Users/nzy/jcode/topological-code-decoupling/julia_code/ToricBuilder` already demonstrates that Julia + Oscar can support:
  - Laurent/polynomial ring construction,
  - Gröbner basis computation,
  - `normal_form`,
  - quotient/module coordinate extraction,
  - exact solver verification with layered tests.
- `ToricBuilder` also gives a proven project shape:
  - top-level module file,
  - focused `core/` submodules,
  - public/expert/internal test separation.

For SuslinStability, that suggests:

- use `src/core/` for algebra and matrix primitives,
- use `src/algorithm/` for theorem-specific steps,
- keep public APIs small,
- test mathematical kernels independently before end-to-end factorization.

## Proposed File Structure

### Package files

- Modify: `Project.toml`
- Modify: `src/SuslinStability.jl`
- Modify: `test/runtests.jl`
- Modify: `README.md`
- Modify: `docs/src/index.md`

### New source files

- Create: `src/core/rings.jl`
  - Polynomial-ring setup helpers.
  - Variable-ordering and monic-in-last-variable helpers.
- Create: `src/core/polynomials.jl`
  - Degree, leading coefficient, resultant, substitution, unit tests for ring conversions.
- Create: `src/core/elementary_matrices.jl`
  - `E_ij(a)` construction, multiplication helpers, factor sequence representation.
- Create: `src/core/matrix_normalization.jl`
  - Exact matrix checks: determinant `1`, ring compatibility, shape checks, block embeddings.
- Create: `src/core/unimodular.jl`
  - Unimodular row/column predicates and witness extraction.
- Create: `src/core/groebner_tools.jl`
  - Thin wrappers over Gröbner / ideal-membership / `normal_form` calls.

### Algorithm files

- Create: `src/algorithm/cohn_type.jl`
  - Realization of Cohn-type matrices for `n >= 3`.
- Create: `src/algorithm/normality.jl`
  - Constructive realization of conjugates `B E_ij(a) B^{-1}`.
- Create: `src/algorithm/column_reduction.jl`
  - Elementary Column Property path and reduction to `e_n`.
- Create: `src/algorithm/quillen_induction.jl`
  - Local realizability data, denominator tracking, patching steps.
- Create: `src/algorithm/sl3_local.jl`
  - Local-ring special-form realization path for
    `[[p,q,0],[r,s,0],[0,0,1]]`.
- Create: `src/algorithm/factorization.jl`
  - End-to-end driver from `SL_n` matrix to factor list.
- Create: `src/algorithm/redundancy.jl`
  - Optional Steinberg-style factor simplification pass.

### Test files

- Create: `test/public/api_surface.jl`
- Create: `test/expert/elementary_matrices.jl`
- Create: `test/expert/cohn_type.jl`
- Create: `test/expert/normality.jl`
- Create: `test/expert/unimodular_columns.jl`
- Create: `test/expert/quillen_induction.jl`
- Create: `test/expert/sl3_local.jl`
- Create: `test/expert/factorization_small_examples.jl`
- Create: `test/internal/rings.jl`
- Create: `test/internal/polynomials.jl`
- Create: `test/internal/groebner_tools.jl`

## Public API Target

The first stable public API should stay narrow:

```julia
using SuslinStability

R, vars = suslin_polynomial_ring(GF(2), ["x", "y"])
A = matrix(R, [...])
factors = elementary_factorization(A)
verify_factorization(A, factors)
```

Recommended initial exports:

- `suslin_polynomial_ring`
- `elementary_matrix`
- `elementary_factorization`
- `verify_factorization`
- `reduce_unimodular_column`

Expert-only APIs can be exported later if they prove useful:

- `realize_cohn_type`
- `realize_conjugate_elementary`
- `quillen_patch`
- `realize_sl3_local`

## Mathematical Implementation Strategy

### Layer 1: exact algebra kernel

Before theorem logic, the package needs exact wrappers for:

- polynomial ring creation over finite fields and small fields,
- checking whether a polynomial is monic in a chosen variable,
- substitution `X -> X + r^l g`,
- resultants,
- Gröbner basis construction,
- `normal_form`,
- ideal membership,
- extracting witnesses from unimodularity computations when available.

This layer should hide backend quirks from later code.

### Layer 2: elementary matrix calculus

Implement:

- elementary matrix constructors,
- composing factor sequences,
- applying factor sequences to matrices,
- embedding `SL_2` and `SL_3` blocks into larger `n x n` matrices,
- equality verification by exact multiplication.

This layer is the backbone for every later constructive proof.

### Layer 3: Cohn-type and normality

This is the first place where the Park-Woodburn proof becomes computational:

- realize Cohn-type matrices explicitly,
- use unimodular vector decomposition to realize matrices of the form `I + v w`,
- realize `B E_ij(a) B^{-1}` constructively,
- verify normality of `E_n` in `SL_n` by producing factors, not by proof only.

This layer should be implemented and tested before Quillen induction.

### Layer 4: Elementary Column Property

Implement reduction of a unimodular column to `e_n` using:

- change of variables to force monicity,
- reduction over smaller-variable rings,
- explicit `SL_2` bridge steps,
- normality to move `SL_2` actions through elementary factors.

This layer reduces general `SL_n` factorization to the `SL_3` special form.

### Layer 5: Quillen induction and patching

Implement constructive glueing:

- localize at maximal ideals,
- track denominators,
- choose finitely many `r_i`,
- use effective ideal membership / Nullstellensatz-style witness extraction,
- patch `A(X)` to `A(0)` via factorized correction matrices.

This is one of the most backend-sensitive parts; keep intermediate data types explicit.

### Layer 6: Local `SL_3` realization

Implement the local special-case algorithm for matrices:

```text
[ p q 0 ]
[ r s 0 ]
[ 0 0 1 ]
```

with `p` monic in the last variable over a local ring. This is the point where Murthy/Gupta becomes essential.

The first implementation can target:

- small finite fields,
- low variable count,
- correctness over speed.

### Layer 7: End-to-end factorization and simplification

Compose the previous layers into:

- `elementary_factorization(A)::Vector{ElementaryFactor}`,
- `verify_factorization(A, factors)::Bool`,
- optional redundancy elimination using Steinberg relations.

## Testing Strategy

Follow the `ToricBuilder` pattern and split tests into:

- `public`: exported API, docs snippets, end-to-end smoke tests,
- `expert`: algorithmic theorem steps on exact examples,
- `internal`: backend wrappers and helper invariants.

Key test classes:

1. **Sanity tests**
   - elementary matrices multiply as expected,
   - determinant remains `1`,
   - block embeddings preserve matrix size.

2. **Cohn-type examples**
   - reproduce the paper's explicit factorization pattern,
   - verify factor products exactly equal the target matrix.

3. **Normality examples**
   - for random small `B` and elementary `E`, construct factors for `BEB^{-1}` and verify.

4. **Unimodular-column examples**
   - reduce small unimodular columns over `GF(2)[x]`, `GF(2)[x,y]`.

5. **Quillen induction examples**
   - use tiny hand-constructed examples where local solutions and patching stay inspectable.

6. **`SL_3` local examples**
   - start with matrices close to triangular or known special-form examples.

7. **Full factorization examples**
   - tiny matrices over `GF(2)[x]` and `GF(2)[x,y]`,
   - compare `prod(factors)` with the original matrix.

8. **Negative tests**
   - reject `SL_2` cases where the theorem does not apply,
   - reject matrices with determinant not equal to `1`,
   - reject unsupported coefficient rings until explicitly implemented.

## Documentation Strategy

The package docs should initially contain:

- a short theorem-and-scope page,
- one worked small example over `GF(2)[x,y]`,
- a note on algorithmic layers and current limitations,
- a bibliography page with direct links.

Do not promise performance or broad coefficient-ring support in the first draft.

## Important References

### Primary target

1. Hyungju Park and Cynthia Woodburn, *An Algorithmic Proof of Suslin's Stability Theorem over Polynomial Rings*.
   - Local copy: [refs/arXiv-alg-geom9405003v1](/Users/nzy/.julia/dev/SuslinStability/refs/arXiv-alg-geom9405003v1:1)
   - arXiv ID: `alg-geom/9405003`
   - Use for the full algorithm skeleton and explicit Cohn-type factorization formulas.

### Core theorem lineage

2. A. A. Suslin, *On the structure of the special linear group over polynomial rings*.
   - Cited in the target paper bibliography as the 1977 source of the stability theorem.

3. A. A. Suslin, *Projective modules over a polynomial ring are free*.
   - Background for the Quillen-Suslin side of the reduction story.

### Local `SL_3` and constructive reductions

4. S. K. Gupta and M. P. Murthy, *Suslin's work on linear groups over polynomial rings and Serre problem*.
   - Explicitly cited by Park-Woodburn as the source of the local special-form algorithm.

5. Cynthia Woodburn, *An Algorithm for Suslin's Stability Theorem* (PhD dissertation, 1994).
   - Likely the best source for implementation-oriented details not fully spelled out in the shorter paper.

### Constructive Quillen-Suslin literature

6. A. Logar and B. Sturmfels, *Algorithms for the Quillen-Suslin theorem*, Journal of Algebra 145 (1992), 231-239.
   - Use for algorithm structure on unimodular rows/completions.

7. N. Fitchas, *Algorithmic aspects of Suslin's proof of Serre's conjecture*, Computational Complexity 3 (1993), 31-55.
   - Use for constructive/algebraic subroutines and proof strategy context.

8. N. Fitchas and A. Galligo, *Nullstellensatz effectif et conjecture de Serre (théorème de Quillen-Suslin) pour le calcul formel*, Math. Nachr. 149 (1990), 231-253.
   - Use for effective Nullstellensatz and witness extraction ideas.

### Background / edge cases

9. P. M. Cohn, *On the structure of the GL_2 of a ring*, Publ. Math. IHÉS 30 (1966), 365-413.
   - Important for `SL_2` non-realizability examples and boundary-case tests.

10. J. Milnor, *Introduction to Algebraic K-Theory*.
    - Background for Steinberg relations and redundancy elimination.

## Important Code Repositories and What to Borrow

### 1. Macaulay2 QuillenSuslin package

- Docs: <https://macaulay2.com/doc/Macaulay2/share/doc/Macaulay2/QuillenSuslin/html/index.html>
- Source: <https://github.com/Macaulay2/M2/blob/master/M2/Macaulay2/packages/QuillenSuslin.m2>
- Package paper: <https://msp.org/jsag/2013/5-1/p05.xhtml>
- DOI: `10.2140/jsag.2013.5.26`

Borrow from this:

- the staging `changeVar -> local solve -> patch -> global result`,
- documentation style for expert/public APIs,
- tests that assert exact equation correctness instead of one particular witness,
- examples that keep matrices small and exact.

Do not copy blindly:

- package-specific API naming,
- assumptions tied to projective module bases rather than elementary factorization.

### 2. `topological-code-decoupling` local reference

- Repo root: `/Users/nzy/jcode/topological-code-decoupling`
- Julia package: `/Users/nzy/jcode/topological-code-decoupling/julia_code/ToricBuilder`

Borrow from this:

- Julia package structure:
  - `src/core/`
  - layered exports in top-level module,
  - `test/public`, `test/expert`, `test/internal`.
- Oscar usage patterns:
  - `laurent_polynomial_ring`,
  - `ideal`,
  - `groebner_basis`,
  - `normal_form`,
  - module/quotient basis extraction.
- test style:
  - verify algebraic identities directly,
  - use exact finite-field examples,
  - keep public and research-level APIs distinct.

Be cautious with:

- Laurent-specific helpers, which are not automatically suitable for pure polynomial-ring Suslin steps,
- code that depends on toric-code semantics rather than generic commutative algebra.

### 3. Maple QuillenSuslin implementation

- Intro page: <https://who.rocq.inria.fr/Alban.Quadrat/QuillenSuslin/intro.html>

Notable signal:

- the page explicitly states that the package implements a constructive Quillen-Suslin theorem,
- cites Logar-Sturmfels inspiration,
- notes an extension with Park's algorithm for Laurent polynomial rings.

Use this as evidence that:

- Park-style constructive algebra can be packaged computationally,
- heuristic improvements may matter later.

### 4. Lean Quillen-Suslin formalization

- Repo: <https://github.com/mbkybky/QuillenSuslin>

Use only as:

- a proof-structure cross-check,
- a source of precise statements and invariants.

It is not a practical implementation reference for the first Julia version.

## Implementation Phases

### Phase 0: package skeleton and dependencies

- add `Oscar`, `Singular`, and `AbstractAlgebra` to `Project.toml`,
- create source/test directory structure,
- wire `test/runtests.jl` into `public/expert/internal` sets,
- add one smoke test that constructs a small polynomial ring.

### Phase 1: exact algebra and elementary matrix core

- implement ring creation and backend wrappers,
- implement elementary matrices and factor sequence type,
- verify exact multiplication and determinant checks.

### Phase 2: Cohn-type and normality

- encode the explicit Cohn-type factorization formulas from the target paper,
- implement `I + v w` realization for unimodular `v` and orthogonal `w`,
- expose expert verification helpers.

### Phase 3: unimodular column reduction

- implement monicity helper via variable changes,
- implement reduction to `e_n`,
- verify reduction for small examples.

### Phase 4: Quillen induction and patching

- represent local realizability certificates,
- implement patching denominators and substitution steps,
- prove correctness computationally on toy examples.

### Phase 5: local `SL_3` algorithm

- implement the special-form local solver,
- add exact examples from literature or derived toy cases,
- connect to global reduction path.

### Phase 6: end-to-end factorization

- expose `elementary_factorization`,
- verify on small matrices over `GF(2)[x]` and `GF(2)[x,y]`,
- document limitations.

### Phase 7: simplification and ergonomics

- add factor simplification,
- improve docs/examples,
- consider broader fields/rings only after the finite-field path is stable.

## Recommended First Milestone

The first milestone should **not** be the full theorem. It should be:

> given a small matrix in `SL_3(GF(2)[x,y])` that is already in a Cohn-type or near-Cohn-type form, return and verify an explicit elementary factorization.

Why:

- it exercises the exact factor-sequence machinery,
- it validates the backend stack,
- it gives a concrete green test early,
- it avoids getting trapped in Quillen induction before the low-level algebra is trustworthy.

## Risks and Mitigations

1. **Backend API mismatch**
   - Mitigation: isolate Oscar/Singular calls in `src/core/groebner_tools.jl`.

2. **Proof steps too non-constructive in practice**
   - Mitigation: start from the explicitly constructive portions: Cohn-type, normality, small unimodular reductions.

3. **Local `SL_3` step under-specified by the short paper**
   - Mitigation: consult Gupta-Murthy and Woodburn dissertation before locking the API.

4. **Test examples too large to debug**
   - Mitigation: require every layer to have tiny exact examples before any random or larger tests.

## Task Plan

### Task 1: Add algebra dependencies and test layout

**Files:**
- Modify: `Project.toml`
- Modify: `test/runtests.jl`
- Create: `test/public/api_surface.jl`
- Create: `test/internal/rings.jl`

- [ ] **Step 1: Write the failing test**

```julia
using SuslinStability
using Test

@testset "api surface" begin
    @test isdefined(SuslinStability, :suslin_polynomial_ring)
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project=. test/runtests.jl`  
Expected: `UndefVarError` or failed `isdefined` assertion for `suslin_polynomial_ring`.

- [ ] **Step 3: Add dependencies and test partition**

Update `Project.toml` to add:

```toml
[deps]
AbstractAlgebra = "c3fe647b-3220-5bb0-a1ea-a7954cac585d"
Oscar = "f1435218-dba5-11e9-1e4d-f1a5fab5fc13"
Singular = "bcd08a7b-43d2-5ff7-b6d4-c458787f915c"
```

Update `test/runtests.jl` to:

```julia
using SuslinStability
using Test

@testset "public" begin
    include("public/api_surface.jl")
end

@testset "expert" begin
end

@testset "internal" begin
    include("internal/rings.jl")
end
```

- [ ] **Step 4: Run test to verify it passes once the symbol exists**

Run: `julia --project=. test/runtests.jl`  
Expected: test tree runs; remaining failure should only be missing implementation.

- [ ] **Step 5: Commit**

```bash
git add Project.toml test/runtests.jl test/public/api_surface.jl test/internal/rings.jl
git commit -m "build: add algebra dependencies and layered test layout"
```

### Task 2: Implement ring and polynomial backend wrappers

**Files:**
- Create: `src/core/rings.jl`
- Create: `src/core/polynomials.jl`
- Modify: `src/SuslinStability.jl`
- Test: `test/internal/rings.jl`

- [ ] **Step 1: Write the failing test**

```julia
using Test
using SuslinStability
using Oscar

@testset "suslin polynomial ring" begin
    R, vars = suslin_polynomial_ring(GF(2), ["x", "y"])
    @test length(vars) == 2
    @test string(vars[1]) == "x"
    @test string(vars[2]) == "y"
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project=. -e 'include("test/internal/rings.jl")'`  
Expected: `UndefVarError: suslin_polynomial_ring not defined`.

- [ ] **Step 3: Write minimal implementation**

Add wrappers shaped like:

```julia
function suslin_polynomial_ring(F, names::Vector{String})
    R, vars = Oscar.polynomial_ring(F, names)
    return R, collect(vars)
end
```

- [ ] **Step 4: Run test to verify it passes**

Run: `julia --project=. -e 'include("test/internal/rings.jl")'`  
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/SuslinStability.jl src/core/rings.jl src/core/polynomials.jl test/internal/rings.jl
git commit -m "feat: add polynomial ring backend wrappers"
```

### Task 3: Implement elementary matrix primitives

**Files:**
- Create: `src/core/elementary_matrices.jl`
- Create: `test/expert/elementary_matrices.jl`
- Modify: `src/SuslinStability.jl`
- Modify: `test/runtests.jl`

- [ ] **Step 1: Write the failing test**

```julia
using Test
using SuslinStability
using Oscar

F = GF(2)
R, (x,) = Oscar.polynomial_ring(F, ["x"])
E = elementary_matrix(3, 1, 2, x + 1, R)

@test size(E) == (3, 3)
@test E[1, 2] == x + 1
@test det(E) == one(R)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project=. -e 'include("test/expert/elementary_matrices.jl")'`  
Expected: `UndefVarError: elementary_matrix not defined`.

- [ ] **Step 3: Write minimal implementation**

Implement:

```julia
function elementary_matrix(n::Int, i::Int, j::Int, a, R)
    i == j && throw(ArgumentError("i and j must differ"))
    E = identity_matrix(R, n)
    E[i, j] = a
    return E
end
```

- [ ] **Step 4: Run test to verify it passes**

Run: `julia --project=. -e 'include("test/expert/elementary_matrices.jl")'`  
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/core/elementary_matrices.jl src/SuslinStability.jl test/expert/elementary_matrices.jl test/runtests.jl
git commit -m "feat: add elementary matrix primitives"
```

### Task 4: Implement Cohn-type realization

**Files:**
- Create: `src/algorithm/cohn_type.jl`
- Create: `test/expert/cohn_type.jl`
- Modify: `src/SuslinStability.jl`

- [ ] **Step 1: Write the failing test**

Construct the `3 x 3` Cohn-type example from the paper and assert:

```julia
target == product_of_factors(factors)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project=. -e 'include("test/expert/cohn_type.jl")'`  
Expected: missing `realize_cohn_type` or equivalent helper.

- [ ] **Step 3: Write minimal implementation**

Encode the explicit eight-factor sequence from the paper for the base `i=1, j=2` case:

```julia
[
    E(1, 3, -v1),
    E(2, 3, -v2),
    E(3, 1, -a * v2),
    E(3, 2,  a * v1),
    E(1, 3,  v1),
    E(2, 3,  v2),
    E(3, 1,  a * v2),
    E(3, 2, -a * v1),
]
```

Then generalize by choosing an index `t != i, j` and embedding the same pattern into rows/columns `(i, j, t)`.

- [ ] **Step 4: Run test to verify it passes**

Run: `julia --project=. -e 'include("test/expert/cohn_type.jl")'`  
Expected: PASS with exact matrix equality.

- [ ] **Step 5: Commit**

```bash
git add src/algorithm/cohn_type.jl src/SuslinStability.jl test/expert/cohn_type.jl
git commit -m "feat: implement constructive Cohn-type realization"
```

### Task 5: Implement constructive normality

**Files:**
- Create: `src/algorithm/normality.jl`
- Create: `test/expert/normality.jl`
- Modify: `src/SuslinStability.jl`

- [ ] **Step 1: Write the failing test**

Choose a small `B` in `GL_3(F[x])` and one elementary `E`, then assert that the returned factor list multiplies to `B * E * inv(B)`.

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project=. -e 'include("test/expert/normality.jl")'`  
Expected: missing constructive conjugation helper.

- [ ] **Step 3: Write minimal implementation**

Implement the `I + v w` path using:

- extracted column `v`,
- extracted orthogonal row `w`,
- decomposition into Cohn-type pieces.

- [ ] **Step 4: Run test to verify it passes**

Run: `julia --project=. -e 'include("test/expert/normality.jl")'`  
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/algorithm/normality.jl src/SuslinStability.jl test/expert/normality.jl
git commit -m "feat: implement constructive normality step"
```

### Task 6: Implement unimodular-column reduction

**Files:**
- Create: `src/core/unimodular.jl`
- Create: `src/algorithm/column_reduction.jl`
- Create: `test/expert/unimodular_columns.jl`
- Modify: `src/SuslinStability.jl`

- [ ] **Step 1: Write the failing test**

Use a small unimodular column over `GF(2)[x,y]` and assert the returned elementary factors send it to `e_n`.

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project=. -e 'include("test/expert/unimodular_columns.jl")'`  
Expected: missing reduction helper.

- [ ] **Step 3: Write minimal implementation**

Start with the smallest supported case:

- detect unimodularity,
- apply monicity normalization in the last variable when needed,
- perform one-step reductions that are already covered by the previous normality primitives.

- [ ] **Step 4: Run test to verify it passes**

Run: `julia --project=. -e 'include("test/expert/unimodular_columns.jl")'`  
Expected: PASS on the small examples only.

- [ ] **Step 5: Commit**

```bash
git add src/core/unimodular.jl src/algorithm/column_reduction.jl src/SuslinStability.jl test/expert/unimodular_columns.jl
git commit -m "feat: add unimodular column reduction scaffolding"
```

### Task 7: Implement Quillen induction scaffolding

**Files:**
- Create: `src/core/groebner_tools.jl`
- Create: `src/algorithm/quillen_induction.jl`
- Create: `test/expert/quillen_induction.jl`
- Modify: `src/SuslinStability.jl`

- [ ] **Step 1: Write the failing test**

Write a toy case that exercises:

- denominator collection,
- substitution `X -> X + r^l g`,
- exact verification that each patch factor is elementary-realizable in the simplified setup.

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project=. -e 'include("test/expert/quillen_induction.jl")'`  
Expected: missing patching helpers.

- [ ] **Step 3: Write minimal implementation**

Implement three concrete helpers:

- `LocalCertificate(::Vector{Int}, ::Vector)` to store local row/column indices and denominators,
- `common_denominator_factor(entries)` to multiply or `lcm` denominators into one ring element,
- `patched_substitution(A, X, r, l, g)` to compute `A(X + r^l * g)` exactly.

- [ ] **Step 4: Run test to verify it passes**

Run: `julia --project=. -e 'include("test/expert/quillen_induction.jl")'`  
Expected: PASS on toy examples.

- [ ] **Step 5: Commit**

```bash
git add src/core/groebner_tools.jl src/algorithm/quillen_induction.jl src/SuslinStability.jl test/expert/quillen_induction.jl
git commit -m "feat: add quillen induction scaffolding"
```

### Task 8: Implement the local `SL_3` special-form solver

**Files:**
- Create: `src/algorithm/sl3_local.jl`
- Create: `test/expert/sl3_local.jl`
- Modify: `src/SuslinStability.jl`

- [ ] **Step 1: Write the failing test**

Prepare a small special-form matrix:

```text
[ p q 0 ]
[ r s 0 ]
[ 0 0 1 ]
```

and assert the solver returns factors whose product equals the matrix.

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project=. -e 'include("test/expert/sl3_local.jl")'`  
Expected: missing `realize_sl3_local`.

- [ ] **Step 3: Write minimal implementation**

Implement only a narrow API:

```julia
realize_sl3_local(p, q, r, s, X; check_monic=true)
```

with preconditions:

- inputs lie in one polynomial ring,
- `p` is monic in `X`,
- the constructed matrix has determinant `1`.

Reject everything outside this special form.

- [ ] **Step 4: Run test to verify it passes**

Run: `julia --project=. -e 'include("test/expert/sl3_local.jl")'`  
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/algorithm/sl3_local.jl src/SuslinStability.jl test/expert/sl3_local.jl
git commit -m "feat: implement local sl3 realization path"
```

### Task 9: Compose end-to-end elementary factorization

**Files:**
- Create: `src/algorithm/factorization.jl`
- Create: `test/expert/factorization_small_examples.jl`
- Modify: `src/SuslinStability.jl`
- Modify: `test/public/api_surface.jl`

- [ ] **Step 1: Write the failing test**

Write one end-to-end public test:

```julia
factors = elementary_factorization(A)
@test verify_factorization(A, factors)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project=. test/runtests.jl`  
Expected: missing `elementary_factorization` or verification helper.

- [ ] **Step 3: Write minimal implementation**

Compose:

- column reduction,
- reduction to `SL_3`,
- local `SL_3` realization,
- factor concatenation.

- [ ] **Step 4: Run test to verify it passes**

Run: `julia --project=. test/runtests.jl`  
Expected: PASS on the small supported examples.

- [ ] **Step 5: Commit**

```bash
git add src/algorithm/factorization.jl src/SuslinStability.jl test/expert/factorization_small_examples.jl test/public/api_surface.jl
git commit -m "feat: add end-to-end elementary factorization api"
```

### Task 10: Document supported scope and references

**Files:**
- Modify: `README.md`
- Modify: `docs/src/index.md`

- [ ] **Step 1: Write the failing docs hygiene expectation**

Add a doc test or smoke assertion that the README example names exported APIs that exist.

- [ ] **Step 2: Run docs-related test to verify it fails**

Run: `julia --project=. test/runtests.jl`  
Expected: docs example mismatch if docs mention missing APIs.

- [ ] **Step 3: Write minimal documentation**

Document:

- theorem scope,
- exact current limitations,
- one small example,
- key references from this plan.

- [ ] **Step 4: Run tests to verify they pass**

Run: `julia --project=. test/runtests.jl`  
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add README.md docs/src/index.md
git commit -m "docs: add suslin stability scope and references"
```

## Self-Review

### Spec coverage

Covered:

- implementation layers,
- file layout,
- dependency selection,
- test strategy,
- milestone ordering,
- reference literature,
- code repositories to borrow from.

Not covered in full detail:

- exact formulas for every Murthy local step,
- a fully specified citation list with verified DOI for every older source.

Those are acceptable gaps for the plan because the literature lookup belongs at implementation time for the exact local `SL_3` routines, and the plan already identifies the required sources.

### Placeholder scan

Checked for `TODO`, `TBD`, and vague "write tests later" phrasing. Remaining intentional abstractions are called out as staged implementation boundaries, not placeholders.

### Type consistency

The plan consistently assumes:

- matrices live over Oscar polynomial rings,
- factorization returns a factor sequence,
- public API centers on `elementary_factorization` and `verify_factorization`.
