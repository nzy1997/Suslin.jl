# Murthy Resultant/Bezout Branch Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement the supported q(0)-nonunit Murthy resultant/Bezout branch for local `SL_3` special forms.

**Architecture:** Extend the existing local `SL_3` certificate path with an internal q(0)-nonunit reduction record. The branch verifies an extracted or supplied Bezout witness, rewrites the target with explicit elementary factors into a q(0)-unit child target, delegates that child to the existing q(0)-unit realization path, and replays every equation in certificate verification.

**Tech Stack:** Julia, Oscar/AbstractAlgebra univariate polynomial rings, existing Suslin elementary matrix, q-degree normalization, q(0)-unit, split-lemma, and certificate helpers.

## Global Constraints

- Preserve `realize_sl3_local(...)` returning a factor sequence.
- Preserve existing open-slice, unit-pivot, q-unit, q-degree normalization, split-lemma, and q(0)-unit branches.
- Supported extraction is ordinary univariate exact polynomial input where `gcdx(p, q)` verifies an exact unit Bezout relation.
- Supplied q(0)-nonunit witnesses must be checked before factors are returned.
- The q(0)-nonunit branch must transform into the existing q(0)-unit child path; do not duplicate q(0)-unit recursion.
- Record the Bezout relation, elementary transformation, q(0)-unit child certificate, and exact final product in replay metadata.
- Do not implement Quillen local-to-global patching, Elementary Column Property reduction, or the final public Park-Woodburn driver.
- Verification command required by issue: `julia --project=. -e 'include("test/expert/sl3_local_murthy_resultant.jl")'`.
- Repository verification command required by Agent Desk: `julia --project=. -e 'using Pkg; Pkg.test()'`.

---

### Task 1: Expert Coverage and Fixtures

**Files:**
- Modify: `test/fixtures/sl3_murthy_gupta_cases.jl`
- Modify: `test/internal/sl3_murthy_gupta_fixtures.jl`
- Create: `test/expert/sl3_local_murthy_resultant.jl`
- Modify: `test/runtests.jl`

**Interfaces:**
- Consumes: existing Murthy fixture catalog and `realize_sl3_local_certificate`.
- Produces: failing tests for supplied Bezout witness, extracted Bezout witness, child q(0)-unit routing, exact factorization replay, and corrupt witness rejection.

- [ ] **Step 1: Add and update fixture entries**

In `test/fixtures/sl3_murthy_gupta_cases.jl`, change `q0_nonunit_normalized_bezout_case` to use `now_supported` and add a second `:q0_nonunit_bezout_resultant` case:

```julia
    q0_nonunit_extracted_bezout_case = _case(
        "mg-q0-nonunit-extracted-bezout-resultant",
        :q0_nonunit_bezout_resultant,
        X,
        (;
            p = X^3 + X + 1,
            q = X^2,
            r = X^3 - X^2 + 2 * X,
            s = X^2 - X + 1,
        ),
        _target(R, X^3 + X + 1, X^2, X^3 - X^2 + 2 * X, X^2 - X + 1),
        ((
            p0 = one(R),
            q0 = zero(R),
            p_prime = one(R) - X,
            q_prime = -X^2 + X - 1,
            resultant = one(R),
            p_prime_degree = 1,
            q_prime_degree = 2,
            branch_unit = one(R),
            case1_entries = (;
                p = X^3 - X^2 + 2 * X,
                q = X^2 - X + 1,
                r = -X^2 + X - 1,
                s = one(R) - X,
            ),
        ),),
        now_supported,
        ("#74",),
    )
```

Add the new case to the returned `cases` array after the existing normalized Bezout case.

- [ ] **Step 2: Update internal catalog assertions**

In `test/internal/sl3_murthy_gupta_fixtures.jl`, assert the new fixture id is present and update the expected supported-status checks:

```julia
    @test haskey(by_id, "mg-q0-nonunit-extracted-bezout-resultant")
    @test by_id["mg-q0-nonunit-extracted-bezout-resultant"].branch == :q0_nonunit_bezout_resultant

    supported_nonunit_diagonal = [
        entry.id for entry in catalog.cases
        if entry.expected_current_solver.status == :passes && !is_unit(entry.target[1, 1]) &&
           !is_unit(entry.target[2, 2])
    ]
    @test "mg-q0-nonunit-normalized-bezout-resultant" in supported_nonunit_diagonal
    @test "mg-q0-nonunit-extracted-bezout-resultant" in supported_nonunit_diagonal
```

Keep the existing corrupted Bezout and branch-unit negative controls.

- [ ] **Step 3: Add the failing expert test**

Create `test/expert/sl3_local_murthy_resultant.jl` with assertions that:

```julia
using Test
using Suslin
using Oscar

const SL3_RESULTANT_FIXTURE_PATH =
    joinpath(@__DIR__, "..", "fixtures", "sl3_murthy_gupta_cases.jl")

function _resultant_product(factors, R)
    product = identity_matrix(R, 3)
    for factor in factors
        product *= factor
    end
    return product
end

function _degree_in_variable(value, X)
    var_idx = findfirst(isequal(X), collect(gens(parent(value))))
    return degree(value, var_idx)
end

function _assert_resultant_certificate(cert; expected_source::Symbol)
    @test cert.branch == :murthy_q0_nonunit_bezout_resultant
    @test Suslin.verify_sl3_local_realization(cert)
    @test Suslin.verify_factorization(cert.target, cert.factors)
    @test _resultant_product(cert.factors, base_ring(cert.target)) == cert.target

    reduction = cert.witness.reduction
    @test reduction.witness_source == expected_source
    @test Suslin.verify_sl3_local_murthy_q0_nonunit_reduction(reduction)
    @test reduction.child_certificate.branch == :murthy_q0_unit
    @test Suslin.verify_sl3_local_realization(reduction.child_certificate)
    @test reduction.bezout_target == reduction.first_elementary_factor * reduction.child_link_target
    @test reduction.target == reduction.left_factor * reduction.bezout_target
    @test reduction.branch_unit * reduction.branch_unit_inverse == one(base_ring(cert.target))
    @test _degree_in_variable(reduction.p_prime, reduction.selected_variable) < reduction.degree_q
    @test _degree_in_variable(reduction.q_prime, reduction.selected_variable) < reduction.degree_p
end
```

The testset should include the existing `mg-q0-nonunit-normalized-bezout-resultant` fixture with
`murthy_q0_nonunit_witness = first(fixture.witnesses)`, the new extracted fixture without a witness keyword, a corrupt supplied `p_prime`, and a tampered reduction whose Bezout relation or child target no longer verifies.

- [ ] **Step 4: Register the expert test**

Add `"expert/sl3_local_murthy_resultant.jl"` after `"expert/sl3_local_murthy_q_unit.jl"` in `test/runtests.jl`.

- [ ] **Step 5: Verify RED**

Run:

```bash
julia --project=. -e 'include("test/expert/sl3_local_murthy_resultant.jl")'
```

Expected: FAIL because `murthy_q0_nonunit_witness` and `verify_sl3_local_murthy_q0_nonunit_reduction` are not implemented.

- [ ] **Step 6: Commit Task 1**

```bash
git add test/fixtures/sl3_murthy_gupta_cases.jl test/internal/sl3_murthy_gupta_fixtures.jl test/expert/sl3_local_murthy_resultant.jl test/runtests.jl
git commit -m "test: cover Murthy resultant Bezout branch"
```

### Task 2: Solver and Replay Implementation

**Files:**
- Modify: `src/algorithm/sl3_local.jl`

**Interfaces:**
- Consumes: Task 1 tests, optional `murthy_q0_nonunit_witness`, existing q(0)-unit child realization.
- Produces: `SL3LocalMurthyQ0NonunitReduction`, `verify_sl3_local_murthy_q0_nonunit_reduction`, `:murthy_q0_nonunit_bezout_resultant` certificates, and exact witness extraction.

- [ ] **Step 1: Add reduction storage**

Add this struct near `SL3LocalMurthyQUnitReduction`:

```julia
struct SL3LocalMurthyQ0NonunitReduction
    target
    p0
    q0
    p_prime
    q_prime
    resultant
    bezout_target
    child_link_target
    left_factor
    first_elementary_factor
    child_certificate
    selected_variable
    degree_p::Int
    degree_q::Int
    degree_p_prime::Int
    degree_q_prime::Int
    branch_unit
    branch_unit_inverse
    witness_source::Symbol
end
```

- [ ] **Step 2: Thread optional witness data**

Add `murthy_q0_nonunit_witness = nothing` keywords to both `realize_sl3_local` methods, both `realize_sl3_local_certificate` methods, `_recognize_sl3_local_matrix`, `_recognize_sl3_local_parameters`, and recursive calls that may need to pass the data into the current recognition. Do not pass a parent witness to unrelated child certificates.

- [ ] **Step 3: Recognize the nonunit branch**

In `_recognize_sl3_local_parameters`, when ordinary univariate Murthy support applies, `deg(q) < deg(p)`, and `q(0)` is not a unit, return:

```julia
return (; family = :murthy_q0_nonunit_bezout_resultant, R, p, q, r, s, X, target, var_idx, murthy_q0_nonunit_witness)
```

If extraction is unsupported later inside the branch, throw a staged local `SL_3` failure. Keep the existing q(0)-unit and q-degree-normalization behavior unchanged.

- [ ] **Step 4: Build the reduction**

Implement `_sl3_local_murthy_q0_nonunit_reduction(form)` that:

```julia
p_prime, q_prime, witness_source = _sl3_local_murthy_q0_nonunit_bezout_pair(form)
resultant = p_prime * form.p - q_prime * form.q
resultant == one(form.R) || throw(ArgumentError("Murthy q(0)-nonunit Bezout equality p_prime*p - q_prime*q must equal 1"))
branch_unit = _sl3_local_constant_coefficient(form.q + p_prime, form.var_idx, form.R)
branch_unit_inverse = _unit_inverse_or_nothing(branch_unit)
branch_unit_inverse === nothing && _throw_staged_sl3_local_failure("Murthy q(0)-nonunit Bezout child q(0) is not a unit")
bezout_target = _sl3_local_special_form_target(form.R, form.p, form.q, q_prime, p_prime)
child_link_target = _sl3_local_special_form_target(form.R, form.p + q_prime, form.q + p_prime, q_prime, p_prime)
left_factor = elementary_matrix(3, 2, 1, form.r * p_prime - form.s * q_prime, form.R)
first_elementary_factor = elementary_matrix(3, 1, 2, -one(form.R), form.R)
child_certificate = realize_sl3_local_certificate(child_link_target, form.X)
```

Verify `form.target == left_factor * bezout_target`,
`bezout_target == first_elementary_factor * child_link_target`,
the child certificate verifies, and the degree guards
`degree(p_prime) < degree(q)` and `degree(q_prime) < degree(p)` hold.

- [ ] **Step 5: Implement supplied and extracted witnesses**

Implement `_sl3_local_murthy_q0_nonunit_bezout_pair(form)`:

```julia
if form.murthy_q0_nonunit_witness !== nothing
    witness = form.murthy_q0_nonunit_witness
    hasproperty(witness, :p_prime) || throw(ArgumentError("Murthy q(0)-nonunit witness must provide p_prime"))
    hasproperty(witness, :q_prime) || throw(ArgumentError("Murthy q(0)-nonunit witness must provide q_prime"))
    return (
        _coerce_into_ring(form.R, witness.p_prime, "Murthy q(0)-nonunit p_prime witness"),
        _coerce_into_ring(form.R, witness.q_prime, "Murthy q(0)-nonunit q_prime witness"),
        :supplied,
    )
end

g, a, b = gcdx(form.p, form.q)
g_inverse = _unit_inverse_or_nothing(g)
g_inverse === nothing && _throw_staged_sl3_local_failure("Murthy q(0)-nonunit Bezout extraction did not produce a unit gcd")
return (g_inverse * a, -g_inverse * b, :extracted)
```

- [ ] **Step 6: Add certificate and verifier integration**

Route `form.family == :murthy_q0_nonunit_bezout_resultant` in `_realize_sl3_local_certificate_form` to a new `_realize_sl3_local_murthy_q0_nonunit_certificate(form)`. Add:

```julia
function verify_sl3_local_murthy_q0_nonunit_reduction(reduction)::Bool
    try
        return _sl3_local_murthy_q0_nonunit_reduction_verification(reduction).overall_ok
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end
```

Update `_sl3_local_certificate_expected_factors`, `_sl3_local_branch_witness_ok`, and `_sl3_local_witness_keys_ok` for `:murthy_q0_nonunit_bezout_resultant`.

- [ ] **Step 7: Verify GREEN**

Run:

```bash
julia --project=. -e 'include("test/expert/sl3_local_murthy_resultant.jl")'
```

Expected: PASS.

- [ ] **Step 8: Commit Task 2**

```bash
git add src/algorithm/sl3_local.jl
git commit -m "feat: implement Murthy resultant Bezout branch"
```

### Task 3: Integration Verification

**Files:**
- Modify only if Task 2 reveals a real integration gap.

**Interfaces:**
- Consumes: Task 1 and Task 2 commits.
- Produces: full verification evidence and final cleanup.

- [ ] **Step 1: Run focused fixture validation**

Run:

```bash
julia --project=. -e 'include("test/internal/sl3_murthy_gupta_fixtures.jl")'
```

Expected: PASS.

- [ ] **Step 2: Run focused expert verification**

Run:

```bash
julia --project=. -e 'include("test/expert/sl3_local_murthy_resultant.jl")'
```

Expected: PASS.

- [ ] **Step 3: Run package verification**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: PASS for the default public and internal groups.

- [ ] **Step 4: Run whitespace check**

Run:

```bash
git diff --check
```

Expected: no output and exit code 0.

- [ ] **Step 5: Commit any integration corrections**

If files changed during Task 3:

```bash
git add <changed files>
git commit -m "test: verify Murthy resultant Bezout integration"
```

If no files changed, do not create an empty commit.
