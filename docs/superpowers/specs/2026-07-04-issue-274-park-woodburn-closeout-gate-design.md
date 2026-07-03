# Issue 274 Park-Woodburn Closeout Gate Design

Issue #274 closes parent #187 by adding an audit-backed gate that ties the
already-landed final ordinary-polynomial Park-Woodburn catalog, public tests,
expert verification, and public documentation into one reviewable claim.

## Context

There is no repository `AGENTS.md`, `CLAUDE.md`, `GEMINI.md`, or
`CONVENTIONS.md`; the README, Documenter index, current tests, and existing
Superpowers specs/plans are the applicable local instructions. The worker branch
is already an isolated linked worktree on
`agent/issue-274-gate-parent-187-with-full-suite-park-woodburn-ac-run-1`.

GitHub issue #274 is open. Its required context is present on `main` through
merged PRs for #270 through #273:

- #270 added the final mainline acceptance catalog.
- #271 added public success coverage for multivariate `SL_3`, recursive
  `SL_n`, and README-style examples.
- #272 hardened final negative controls for determinant, unsupported
  coefficient ring, missing evidence, and Laurent boundaries.
- #273 published the narrow #187 public contract in README and Documenter.

Earlier parent closeout PRs show the upstream ordinary-polynomial layers are
complete: #195 for #181, #212 for #182, #220 for #183, #239 for #184, #249 for
#185, and #266 for #186.

## Approaches Considered

1. Add only a Markdown audit note. This is low risk, but it would not make the
   negative controls in #274 executable. Removing the unsupported coefficient
   negative case could still pass if no test names it directly.
2. Add a Markdown audit note plus documentation-smoke assertions. This proves
   the audit exists and that public docs keep #188, Laurent/ToricBuilder, and
   unsupported coefficient rings outside #187, but it still leaves final public
   case inventory implicit.
3. Add a Markdown audit note, documentation-smoke assertions, and public
   inventory assertions over the existing final acceptance catalog. This is the
   recommended approach because it is a closeout gate: it does not add new
   algorithmic support, but it makes the full-suite proof fail if the recursive
   `SL_n` acceptance case or unsupported coefficient-ring negative control is
   removed.

## Chosen Approach

Use approach 3. Keep implementation behavior unchanged. Add a short mechanical
audit note under `docs/audits/` mapping #181/#195, #182/#212, #183/#220,
#184/#239, #185/#249, and #186/#266 to the #270/#271/#272 public acceptance
and negative-control cases that consume each layer. Extend public tests to name
the required final acceptance case ids and negative-control ids directly. Extend
the documentation smoke test to read the audit note, verify every upstream row
is present, and verify the public docs still keep #188 and non-polynomial tracks
separate.

## Behavioral Requirements

- The full public suite must include final #187 success cases for:
  - multivariate exact field-backed ordinary-polynomial `SL_3`;
  - recursive exact field-backed ordinary-polynomial `SL_n`, `n > 3`;
  - README-style public `elementary_factorization(A)` and
    `verify_factorization(A, factors)` usage.
- The full public suite must include final #187 negative controls for:
  - determinant-not-one;
  - unsupported coefficient ring;
  - missing `SL_3` local-form evidence;
  - missing `SL_3` Quillen evidence;
  - missing ECP evidence;
  - missing final `SL_3` evidence;
  - Laurent boundary.
- The audit note must map every upstream parent layer #181-#186 to its upstream
  closeout issue and the final #187 case ids that consume that layer.
- README and `docs/src/index.md` must continue to claim only the narrow exact
  field-backed ordinary-polynomial #187 support boundary.
- README and `docs/src/index.md` must not claim Laurent/ToricBuilder mainline
  support, unsupported coefficient-ring support, or Steinberg factor-count
  optimization (#188) as part of #187.

## Scope Boundaries

Do not change factorization algorithms, fixture matrices, public APIs, or route
certificate semantics. Do not implement Laurent/ToricBuilder mainline support,
unsupported coefficient-ring support, or Steinberg factor-count optimization.
Do not close or claim #188.

## Tests

Focused red/green commands:

```bash
julia --project=. -e 'include("test/public/park_woodburn_polynomial_factorization.jl")'
julia --project=. -e 'include("test/expert/documentation_smoke.jl")'
```

Required final verification:

```bash
julia --project=. test/runtests.jl all
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Self-Review

This design is a closeout gate, not an implementation expansion. It names the
already-supported public cases and negative controls, adds the missing audit
artifact required by #274, and preserves all out-of-scope boundaries for
Laurent/ToricBuilder, unsupported coefficient rings, and #188.
