# Issue 75 Murthy-Gupta Acceptance Design

## Context

Issue #75 is the acceptance gate for parent issue #61. The preceding staged
issues already added the Murthy-Gupta fixture catalog, replayable local
certificates, q-degree normalization, split-lemma replay, q(0)-unit recursion,
and the q(0)-nonunit Bezout/resultant branch. The remaining work is to make the
complete local solver milestone checkable through `realize_sl3_local`.

No repository instruction file such as `AGENTS.md`, `CLAUDE.md`, or
`CONVENTIONS.md` is present in this checkout. The current worker branch is
already an Agent Desk linked worktree.

## Approaches Considered

1. Add a focused final acceptance test around the existing solver path.
   This is the chosen approach because the current main branch already routes
   the Murthy branches through `realize_sl3_local_certificate` and
   `realize_sl3_local`, while issue #75 asks for acceptance coverage rather
   than a new public driver.

2. Add a new expert wrapper API for the combined Murthy-Gupta solver.
   This would make the acceptance scenario explicit, but it would create an
   extra entry point without new algebraic responsibility and would not improve
   the existing public factor-returning path.

3. Route more cases through `elementary_factorization`.
   This is out of scope for #75. The issue asks to keep driver routing minimal
   and not expand into the later Park-Woodburn driver issues.

## Design

Add `test/expert/sl3_local_murthy_gupta.jl` as the final #61 local-solver
acceptance test. It will load the shared Murthy-Gupta fixture catalog and select
at least three hand-checkable ordinary `QQ[X]` special-form matrices whose
`p` entry is monic and whose diagonal entries are both non-units. For each
case, the test will call `realize_sl3_local_certificate` and `realize_sl3_local`
on the normal public local-solver path, assert a nonempty elementary sequence,
and verify exact product equality with `verify_factorization`.

The acceptance cases will include:

- a q-degree normalization case that previously could not be handled by only
  open-slice or unit-pivot branches;
- a q(0)-unit recursive case;
- two q(0)-nonunit Bezout/resultant cases, one with supplied witness metadata
  and one using extracted witness data.

The test will assert replayable certificate metadata instead of treating factor
success alone as sufficient. It will verify the top-level Murthy branch names,
normalization records, q(0)-unit reduction records, split-lemma replay records,
Bezout/resultant reduction records, witness-source metadata, degree guards, and
the exact elementary identities used by each branch.

Register the new file in the expert group in `test/runtests.jl`.

## Unsupported Boundary

The negative control will construct a determinant-one special-form matrix with
non-monic `p` and no supplied local witness. It must throw an `ArgumentError`
whose message is the staged local solver precondition failure, not return
factors and not surface an unstructured backend exception.

## Scope

This issue does not implement Quillen induction, the Elementary Column Property
pipeline, arbitrary local rings, factor-count optimization, or the final public
Park-Woodburn polynomial `SL_n` driver. Production solver changes are limited to
fixing any gap exposed by the acceptance test; if the existing Murthy branches
already satisfy the contract, the code change should remain test-only.

## Testing

Use test-first implementation:

- add the focused acceptance test and register it;
- run it before any solver code changes and confirm the expected red state is
  only the missing test registration or a real contract gap;
- if a solver gap appears, implement the minimal fix in
  `src/algorithm/sl3_local.jl`;
- run the focused acceptance command, the all-group test runner, and the package
  test command required by Agent Desk.
