# Issue 69 Murthy-Gupta Fixture Catalog Design

## Context

Issue 61 will add the Murthy-Gupta local `SL_3` solver for special-form
matrices

```text
[p q 0; r s 0; 0 0 1]
```

with `p` monic in the selected variable. The current `realize_sl3_local`
implementation only handles open-slice and unit-pivot cases, so the Murthy
work needs exact fixtures whose witnesses can be cited by later issues instead
of redefined in each solver test.

The source boundary is Park-Woodburn section 5, especially the q-degree
normalization, Lemma 4 split relation, q(0)-unit recursion, and q(0)-nonunit
Bezout/resultant branch. The issue comment narrows the scope: keep the catalog
small, accept only fixtures consumed by the planned Murthy child issues, and
replay every claimed equality exactly.

## Design Choice

Add a dedicated test fixture catalog at `test/fixtures/sl3_murthy_gupta_cases.jl`
and a focused internal validator at
`test/internal/sl3_murthy_gupta_fixtures.jl`.

The catalog returns structured named tuples. Each case records:

- `id`
- `branch`
- `ring_constructor`
- `ring`
- `variable`
- `entries`
- `target`
- `murthy_path`
- `expected_current_solver`
- `witnesses`
- `source_refs`
- `consumer_issue_ids`

This remains test support, not a public Suslin API. Later Murthy issues can
include the fixture file and select cases by `id`.

Alternatives considered:

- Add a public exported fixture API. That overstates a test-only contract and
  would create compatibility burden before the solver API is known.
- Put fixtures inside the future solver tests. That would repeat examples and
  lose the issue's single-source-of-truth goal.
- Build a broad catalog of many polynomial examples. The issue comment rejects
  coverage vanity; five small cases cover the named branches and current solver
  boundary.

## Fixture Scope

The initial catalog contains exactly the branch contracts needed by the next
Murthy issues:

- `mg-q-degree-normalization`: a degree-reduction example where `q = f*p + g`
  and right multiplication by `E12(-f)` is checked exactly.
- `mg-split-lemma-x-square`: a Lemma 4 split witness for `p = X * X`, including
  the two Bezout relations and the exact elementary-matrix identity from the
  lemma.
- `mg-q0-unit-recursion`: a q(0)-unit example that records the right `E21`
  normalization making `p(0) = 0`, plus the split data for the normalized
  matrix.
- `mg-q0-nonunit-bezout-resultant`: a q(0)-nonunit branch with the explicit
  Bezout/resultant certificate `p_prime*p - q_prime*q == 1` and the exact
  reduction to the q(0)-unit case.
- `mg-open-slice-control`: a current-pass control so the catalog records both
  staged failures and already-supported local `SL_3` behavior.

All cases live over `QQ[X]`. The first four are determinant-one, have monic
nonconstant `p`, have neither diagonal entry a unit, and should currently fail
through the staged local `SL_3` solver boundary. The control is intentionally
already supported.

## Validation Rules

Every case must have required metadata, unique ids, and at least one consumer
issue id. The validator reconstructs the `3 x 3` target from `(p, q, r, s)`,
checks the recorded target, and checks `det(target) == 1`.

For cases marked `murthy_path = true`, the validator checks that `p` is monic
in the selected variable. It also checks the current solver status: pass cases
must produce factors whose product verifies exactly, and staged-fail cases must
throw the existing staged local `SL_3` failure.

Witness validation is exact:

- q-degree witnesses must satisfy `q == quotient*p + remainder`, the degree
  bound for the remainder, and the recorded right elementary normalization.
- split-lemma witnesses must satisfy `p == a*a_prime`, both split Bezout
  equations, and the full elementary-matrix identity from Park-Woodburn Lemma 4.
- q(0)-unit witnesses must check the recorded constant terms, inverse, right
  `E21` normalization, and nested split witness.
- q(0)-nonunit witnesses must check `q(0)` is nonunit, the Bezout/resultant
  equality, degree bounds, the branch unit condition, and the exact Case 2
  reduction identity.

The test includes negative controls that deliberately corrupt witness data and
prove the validator rejects it.

## Files

- Create `test/fixtures/sl3_murthy_gupta_cases.jl`: exact catalog construction.
- Create `test/internal/sl3_murthy_gupta_fixtures.jl`: validator, focused tests,
  branch coverage checks, current solver status checks, and negative controls.
- Modify `test/runtests.jl`: include the focused validator in the internal group.

## Verification

Focused validator:

```bash
julia --project=. -e 'include("test/internal/sl3_murthy_gupta_fixtures.jl")'
```

Package entry point:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Full expert-inclusive suite remains available as:

```bash
julia --project=. test/runtests.jl all
```

## Spec Self-Review

- No incomplete markers remain.
- Scope is limited to test fixtures and internal validation.
- The source-derived relations map directly to Park-Woodburn section 5.
- The catalog has at least five named entries and at least two current staged
  Murthy targets with nonunit diagonal entries.
