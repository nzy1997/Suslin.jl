# Issue 116 Polynomial Public Acceptance Design

## Context

Issue 116 is the final public acceptance pass for parent issue 64. The driver
already has the staged route shell, recursive ordinary-polynomial column peel,
and deterministic Quillen fixture route from issues 113 and 115. This change
should prove the public `elementary_factorization(A)` contract through a single
focused acceptance file, register that file in the default public group, and
update user-facing docs to describe the conservative supported scope.

This checkout has no `AGENTS.md`. The branch is an isolated Agent Desk worktree
on `agent/issue-116-add-final-public-park-woodburn-polynomial-factor-run-1`.
Issues 113 and 115 are closed. Issue 64 remains open and lists issue 116 as the
M5 public acceptance/docs child.

## Approaches Considered

Recommended: add a public acceptance file that reuses the existing
Park-Woodburn polynomial fixture catalog. The test calls only public
factorization APIs for behavior, and uses route certificates as evidence that
each example exercises the intended conservative route. This avoids duplicating
large matrices by hand while keeping exact route coverage explicit.

Alternative: move all existing `factorization_driver_shell.jl` polynomial route
checks into the new file. That would centralize coverage, but it would make this
issue larger and risk weakening the existing shell regression file.

Alternative: create wholly new matrices in the public test. That would be more
independent from the catalog, but it would duplicate witness construction and
increase maintenance risk for no new driver behavior.

The chosen design is the recommended catalog-backed acceptance file because it
matches prior fixture-driven coverage, stays fail-closed, and gives the final
parent issue an exact public command.

## Design

Create `test/public/park_woodburn_polynomial_factorization.jl`. It loads
`test/fixtures/park_woodburn_polynomial_cases.jl` and tests three supported
ordinary-polynomial `SL_n` examples through public `elementary_factorization(A)`:

- `pw-poly-univariate-sl3-fast-local-qq`, expected automatic route
  `:fast_local_sl3`;
- `pw-poly-recursive-column-peel-sln-block-qq`, expected automatic route
  `:polynomial_column_peel` with a `PolynomialColumnPeelCertificate` whose final
  certificate route is `:disjoint_local_blocks`;
- `quillen-patched-substitution-witness-qq`, expected automatic route
  `:quillen_patch` with a verified `PolynomialQuillenPatchRouteAdapter`.

For every supported example, the returned factors must be a non-`nothing`
sequence and `verify_factorization(A, factors)` must return `true`. The route
certificate checks are private-test evidence, but the acceptance behavior stays
anchored on the public calls.

The negative controls use the catalog determinant-not-one matrix and the
determinant-one outside-witness matrix. Each is wrapped in a helper that records
whether a factor sequence was returned. The tests assert an `ArgumentError`,
clear staged error text, and `factors === nothing`.

Register the new file in the `public` group in `test/runtests.jl` so default
package tests exercise the public command.

Update `README.md` and `docs/src/index.md` so the documented scope now says:

- univariate local `SL_3` route is supported;
- selected `n > 3` ordinary-polynomial matrices are supported through
  block-local reduction and recursive polynomial column peel;
- deterministic multivariate Quillen/local-to-global fixture route is
  supported after the constructive Quillen path;
- arbitrary local realizability, broad coefficient-ring support, Laurent
  `GL_n` determinant corrections, and factor-count optimization remain staged
  boundaries.

## Verification

Required focused command:

```bash
julia --project=. -e 'include("test/public/park_woodburn_polynomial_factorization.jl")'
```

Required full command:

```bash
julia --project=. test/runtests.jl all
```

Required package command from the Agent Desk instructions:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Automatic Approval

This Agent Desk run is non-interactive. Under the Standing Answer Policy, the
design is approved automatically because it chooses the conservative
catalog-backed acceptance route, preserves public API shape, and documents
remaining staged boundaries rather than broadening support claims.

## Spec Self-Review

- No incomplete markers remain.
- The design includes the three issue-required supported public examples.
- Both required negative controls fail without returning factors.
- Documentation scope stays conservative and exact-field-backed.
