# Issue 16 Extended Local SL3 Solver Design

## Goal

Extend `realize_sl3_local` from the two current open Gauss-cell slices to a
small, explicit family of embedded `SL_2` special forms inside `SL_3`, with
exact verification and staged failures for unsupported local inputs.

## Context

The current solver accepts parameter tuples `(p, q, r, s, X)` and builds the
embedded matrix

```text
[ p q 0 ]
[ r s 0 ]
[ 0 0 1 ]
```

over one exact parent ring. It only recognizes the unipotent products
`E12(q) * E21(r)` and `E21(r) * E12(q)`. Issue #16 asks for a broader local
solver but explicitly keeps general `SL_n` reduction and ToricBuilder
acceptance out of scope. Issues #12 and #21 are merged: Laurent elementary
routines and the documented full-suite command are available.

No `AGENTS.md`, `CLAUDE.md`, or `GEMINI.md` file is present in this checkout.
The README documents the test commands and states that the full suite is
`julia --project=. test/runtests.jl all`.

## Approaches Considered

1. Add recognition for embedded `SL_2` unit-pivot special forms. This is the
   chosen approach. It keeps the solver local, uses elementary formulas, and
   supports polynomial unit pivots and Laurent monomial-unit pivots without
   changing the public driver architecture.
2. Add a general `SL_2` Euclidean or Cohn-style reduction. This would be more
   powerful but is too broad for a local-solver issue and risks overlapping
   the later `SL_n` reduction work.
3. Route all unsupported local inputs through the existing Cohn/normality
   routines. Those routines are useful primitives, but using them as an
   unqualified local fallback would silently broaden the accepted domain
   without precise tests.

## Design

Keep the existing `realize_sl3_local(p, q, r, s, X; check_monic=true)` method
and add a matrix convenience method
`realize_sl3_local(A, X; check_monic=true)` for the exact embedded special
form. The matrix method checks that `A` is `3 x 3`, has zero third-row and
third-column off-block entries, and has `A[3, 3] == 1`; it then delegates to
the parameter method.

Split the implementation into three internal stages:

- Special-form recognition validates common parent rings, the selected
  generator `X`, determinant one, the optional polynomial monicity assumption,
  and the supported family tag.
- Constructive solving maps each family tag to an explicit factor sequence.
- Exact verification multiplies the returned factors and throws an internal
  verification error if the product differs from the target.

Preserve the two current paths exactly:

- `s == 1 && p == 1 + q*r` returns `[E12(q), E21(r)]`.
- `p == 1 && s == 1 + q*r` returns `[E21(r), E12(q)]`.

Add two broader unit-pivot families:

- If `s` is a unit, factor
  `A = E12(q*s^-1) * D(s^-1) * E21(r*s^-1)`.
- If `p` is a unit, factor
  `A = E21(r*p^-1) * D(p) * E12(q*p^-1)`.

Here `D(u) = diag(u, u^-1, 1)` is realized by elementary `SL_3` factors using

```text
D(u) = E12(u) E21(-u^-1) E12(u) E12(-1) E21(1) E12(-1)
```

The unit-pivot branches use `is_unit` and `inv` only after recognition has
confirmed that the pivot is a unit. If neither pivot is supported, the solver
throws an `ArgumentError` whose message starts with
`staged local SL_3 solver failure`.

## Tests

Add `test/expert/sl3_local_extended.jl` and register it in the expert group.
The focused test file covers:

- The two existing open Gauss-cell examples, including exact product checks.
- A new polynomial-ring form with `s = 2` and monic `p`, solved by the `s`
  unit-pivot branch.
- A new Laurent-ring form with `p = x` and nonunit `s`, solved by the `p`
  unit-pivot branch with `check_monic=false` to make the discharged assumption
  explicit.
- Matrix input for a supported form.
- An embedded determinant-one local `SL_3` input outside the implemented
  family, which must fail with a staged local-solver error.

All supported cases assert exact equality between the target matrix and the
factor product, and also use `verify_factorization` where the target and
factor parent rings match.

## Verification

Run these commands:

```bash
julia --project=. -e 'include("test/expert/sl3_local_extended.jl")'
julia --project=. -e 'using Pkg; Pkg.test()'
julia --project=. test/runtests.jl all
```

The first command is the issue-specific check. The second is required by Agent
Desk and runs the default public/internal groups. The third is the documented
full suite from issue #21.

## Scope Boundaries

This design does not implement general `SL_n` reduction, does not add a
ToricBuilder acceptance harness, and does not turn `realize_sl3_local` into a
general `SL_2` solver. It only accepts the explicit embedded local families
covered by the new expert tests.
