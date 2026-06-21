# Issue 69 Murthy-Gupta Fixture Catalog Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a small shared Murthy-Gupta local `SL_3` fixture catalog and exact validator for later Issue 61 child work.

**Architecture:** Keep the catalog under `test/fixtures/` as test support, not public API. The fixture module builds exact `QQ[X]` special-form matrices and witness records, while the internal validator reconstructs targets, checks determinant and monicity, verifies every witness equality, checks current solver pass/staged-fail status, and carries negative controls.

**Tech Stack:** Julia, Oscar polynomial rings, Suslin local `SL_3` helpers, Test stdlib.

## Global Constraints

- Do not implement a Murthy-Gupta solver branch.
- Do not add Quillen induction, Elementary Column Property reduction, or a public Park-Woodburn driver.
- Do not add a public Suslin fixture API.
- Catalog file is `test/fixtures/sl3_murthy_gupta_cases.jl`.
- Validator file is `test/internal/sl3_murthy_gupta_fixtures.jl`.
- Catalog must validate at least five named entries.
- Required fixture ids are `mg-q-degree-normalization`, `mg-split-lemma-x-square`, `mg-q0-unit-recursion`, `mg-q0-nonunit-bezout-resultant`, and `mg-open-slice-control`.
- Required branches are `:q_degree_normalization`, `:split_lemma`, `:q0_unit_recursion`, `:q0_nonunit_bezout_resultant`, and `:open_slice_control`.
- All initial cases use `QQ[X]`.
- Every Murthy-path case must reconstruct the target matrix from `(p, q, r, s)`, satisfy `det(target) == 1`, and have monic `p` in `X`.
- At least two cases must have determinant one, monic nonconstant `p`, neither diagonal entry a unit, and current solver status `:staged_fail`.
- Every supplied witness relation must be replay-checked exactly.
- Negative controls must corrupt at least one split-lemma or Bezout equality and prove the validator rejects it.
- Focused validator command is `julia --project=. -e 'include("test/internal/sl3_murthy_gupta_fixtures.jl")'`.
- Package verification command is `julia --project=. -e 'using Pkg; Pkg.test()'`.
- Source basis is Park-Woodburn section 5, arXiv `alg-geom/9405003`, especially equations (54)-(58) and Lemma 4.

---

## File Structure

- Create `test/internal/sl3_murthy_gupta_fixtures.jl`: owns catalog validation, focused tests, branch coverage checks, current solver boundary checks, and negative controls.
- Create `test/fixtures/sl3_murthy_gupta_cases.jl`: owns exact fixture construction and witness metadata.
- Modify `test/runtests.jl`: add `internal/sl3_murthy_gupta_fixtures.jl` to the internal group after `internal/toricbuilder_problem_catalog.jl`.

---

### Task 1: Murthy-Gupta Fixture Catalog And Validator

**Files:**
- Create: `test/internal/sl3_murthy_gupta_fixtures.jl`
- Create: `test/fixtures/sl3_murthy_gupta_cases.jl`
- Modify: `test/runtests.jl`

**Interfaces:**
- Consumes: `Suslin.realize_sl3_local`, `Suslin.verify_factorization`, `Suslin.elementary_matrix`, Oscar matrix and polynomial APIs.
- Produces: `SL3MurthyGuptaFixtureCatalog.catalog()` returning `(; ring, cases)` where each case has `id`, `branch`, `ring_constructor`, `ring`, `variable`, `entries`, `target`, `murthy_path`, `expected_current_solver`, `witnesses`, `source_refs`, and `consumer_issue_ids`.
- Produces: internal validator functions `validate_sl3_murthy_gupta_fixture(entry)` and `validate_sl3_murthy_gupta_fixture_catalog(catalog)`.

- [ ] **Step 1: Write the failing validator**

Create `test/internal/sl3_murthy_gupta_fixtures.jl`.

The validator must:

```julia
using Test
using Suslin
using Oscar

const SL3_MURTHY_GUPTA_CATALOG_PATH =
    joinpath(@__DIR__, "..", "fixtures", "sl3_murthy_gupta_cases.jl")
const REQUIRED_SL3_MURTHY_GUPTA_FIELDS = (
    :id,
    :branch,
    :ring_constructor,
    :ring,
    :variable,
    :entries,
    :target,
    :murthy_path,
    :expected_current_solver,
    :witnesses,
    :source_refs,
    :consumer_issue_ids,
)
```

Implement helper functions with these exact names and behavior:

```julia
_sl3_mg_target(entry)
_sl3_mg_matrix(R, p, q, r, s)
_sl3_mg_product(factors, R)
_sl3_mg_field(entry, field::Symbol)
_sl3_mg_monic_in_variable(p, X)
_sl3_mg_constant_coefficient(value)
_sl3_mg_degree_in_variable(value, X)
_sl3_mg_assert_metadata(entry)
_sl3_mg_assert_target(entry)
_sl3_mg_assert_murthy_path(entry)
_sl3_mg_assert_current_solver_status(entry)
_sl3_mg_assert_q_degree_witness(entry, witness)
_sl3_mg_assert_split_lemma_witness(entry, witness)
_sl3_mg_assert_q0_unit_witness(entry, witness)
_sl3_mg_assert_q0_nonunit_bezout_witness(entry, witness)
_sl3_mg_assert_witnesses(entry)
validate_sl3_murthy_gupta_fixture(entry)
validate_sl3_murthy_gupta_fixture_catalog(catalog)
```

The split witness check must reconstruct the exact Lemma 4 elementary identity:

```julia
M = _sl3_mg_matrix(R, a * a_prime, b, c, d)
first = _sl3_mg_matrix(R, a, b, c1, d1)
second = _sl3_mg_matrix(R, a_prime, b, c2, d2)
rhs =
    elementary_matrix(3, 2, 1, c * d1 * d2 - d * (c2 + a_prime * c1 * d2), R) *
    elementary_matrix(3, 2, 3, d2 - one(R), R) *
    elementary_matrix(3, 3, 2, one(R), R) *
    elementary_matrix(3, 2, 3, -one(R), R) *
    first *
    elementary_matrix(3, 2, 3, one(R), R) *
    elementary_matrix(3, 3, 2, -one(R), R) *
    elementary_matrix(3, 2, 3, one(R), R) *
    second *
    elementary_matrix(3, 2, 3, -one(R), R) *
    elementary_matrix(3, 3, 2, one(R), R) *
    elementary_matrix(3, 2, 3, a - one(R), R) *
    elementary_matrix(3, 3, 1, -a_prime * c1, R) *
    elementary_matrix(3, 3, 2, -d1, R)
M == rhs || throw(ArgumentError("fixture $(entry.id) split-lemma elementary identity failed"))
```

The q(0)-nonunit witness check must reconstruct the exact Case 2 reduction:

```julia
bezout_matrix = _sl3_mg_matrix(R, p, q, q_prime, p_prime)
case1_matrix = _sl3_mg_matrix(R, p + q_prime, q + p_prime, q_prime, p_prime)
target ==
    elementary_matrix(3, 2, 1, r * p_prime - s * q_prime, R) * bezout_matrix ||
    throw(ArgumentError("fixture $(entry.id) Bezout reduction first equality failed"))
bezout_matrix ==
    elementary_matrix(3, 1, 2, -one(R), R) * case1_matrix ||
    throw(ArgumentError("fixture $(entry.id) Bezout reduction q0-unit equality failed"))
```

Add a testset named `"Murthy-Gupta local SL3 fixture catalog"` that includes the catalog file, validates the catalog, asserts all five required ids and branches are present, asserts at least two staged-fail nonunit-diagonal cases, and corrupts one split witness plus one Bezout witness with `@test_throws ArgumentError`.

- [ ] **Step 2: Run focused test to verify RED**

Run:

```bash
julia --project=. -e 'include("test/internal/sl3_murthy_gupta_fixtures.jl")'
```

Expected: FAIL because `test/fixtures/sl3_murthy_gupta_cases.jl` does not exist yet.

- [ ] **Step 3: Implement the catalog**

Create `test/fixtures/sl3_murthy_gupta_cases.jl` with `module SL3MurthyGuptaFixtureCatalog`.

Use these exact polynomial entries over `R, (X,) = Oscar.polynomial_ring(QQ, ["X"])`:

```julia
# mg-q-degree-normalization
p = X^2 + 1
q = X^3 + X + 1
r = -one(R)
s = -X
quotient = X
remainder = one(R)
normalized_s = zero(R)

# mg-split-lemma-x-square
p = X^2
q = one(R)
r = X^3 + X^2 - 1
s = X + 1
a = X
a_prime = X
c1 = r
c2 = r
d1 = X^2 + X
d2 = X^2 + X

# mg-q0-unit-recursion
p = X + 1
q = one(R)
r = X^2 + 2 * X
s = X + 1
p0 = one(R)
q0 = one(R)
q0_inverse = one(R)
right_e21_coefficient = -one(R)
normalized_p = X
normalized_r = X^2 + X - 1
normalized_s = X + 1
split = (;
    a = X,
    a_prime = one(R),
    b = one(R),
    c = X^2 + X - 1,
    c1 = X^2 + X - 1,
    c2 = X^2 + X - 1,
    d1 = X + 1,
    d2 = X^2 + X,
    d = X + 1,
)

# mg-q0-nonunit-bezout-resultant
p = X + 1
q = X
r = X + 2
s = X + 1
p0 = one(R)
q0 = zero(R)
p_prime = one(R)
q_prime = one(R)
resultant = one(R)
p_prime_degree = 0
q_prime_degree = 0
branch_unit = one(R)
case1_entries = (;
    p = X + 2,
    q = X + 1,
    r = one(R),
    s = one(R),
)

# mg-open-slice-control
p = X + 1
q = one(R)
r = X
s = one(R)
```

Each case must set `murthy_path = true`, `source_refs = ("Park-Woodburn arXiv:alg-geom/9405003 section 5",)`, and at least one of `consumer_issue_ids = ("#71",)`, `("#72",)`, `("#73",)`, `("#74",)`, or `("#75",)`.

The first four cases must set:

```julia
expected_current_solver = (; status = :staged_fail, message_substring = "staged local SL_3 solver failure")
```

The control must set:

```julia
expected_current_solver = (; status = :passes)
```

- [ ] **Step 4: Register the internal validator**

Modify `test/runtests.jl` and add:

```julia
"internal/sl3_murthy_gupta_fixtures.jl",
```

after `"internal/toricbuilder_problem_catalog.jl",`.

- [ ] **Step 5: Run focused test to verify GREEN**

Run:

```bash
julia --project=. -e 'include("test/internal/sl3_murthy_gupta_fixtures.jl")'
```

Expected: PASS with the Murthy-Gupta fixture testset, at least 18 passing checks, and no warnings.

- [ ] **Step 6: Run package verification**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: PASS.

- [ ] **Step 7: Commit**

Run:

```bash
git add docs/superpowers/plans/2026-06-21-issue-69-murthy-gupta-fixture-catalog.md \
    test/fixtures/sl3_murthy_gupta_cases.jl \
    test/internal/sl3_murthy_gupta_fixtures.jl \
    test/runtests.jl
git commit -m "Add Murthy-Gupta SL3 fixture catalog"
```

Expected: commit succeeds with only the plan, fixture, validator, and test registration staged.

---

## Self-Review

- Spec coverage: the task covers the fixture catalog, validator, branch ids, exact witness replay, negative controls, staged failure requirements, and verification commands.
- Placeholder scan: no TBD, TODO, or fill-in steps remain.
- Type consistency: `SL3MurthyGuptaFixtureCatalog.catalog()` and `validate_sl3_murthy_gupta_fixture_catalog(catalog)` are named consistently across the task.
