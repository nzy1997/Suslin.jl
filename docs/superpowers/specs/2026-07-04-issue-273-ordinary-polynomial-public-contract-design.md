# Issue 273 Ordinary-Polynomial Public Contract Design

Issue #273 publishes the final ordinary-polynomial Park-Woodburn #187 public
contract in README and Documenter after the final positive and negative
acceptance coverage landed in #271 and #272.

## Context

There is no repository `AGENTS.md`; the README, current Documenter index, and
existing Superpowers specs are the applicable local instructions. GitHub issue
#273 is open. Its dependencies #249, #266, #271, and #272 are closed as
completed, and current `main` includes merged PR #284 for #271 and PR #285 for
#272.

The current README and `docs/src/index.md` still say that "full public
Park-Woodburn acceptance (#187)" remains staged. That was correct before #271
and #272, but is now stale for the narrow public surface those issues proved:
exact field-backed ordinary polynomial rings, determinant-one `SL_3` and
`SL_n`, `n > 3`, inputs, and only the implemented evidence-backed route.

## Chosen Approach

Use a documentation-only change guarded by a focused documentation smoke test.

1. Update `README.md` and `docs/src/index.md` with the same support-boundary
   language: #187 is supported for exact field-backed ordinary-polynomial
   determinant-one `SL_3` and `SL_n`, `n > 3`, inputs when the implemented route
   can replay #184 final `SL_3` evidence and #185/#186 recursive ECP evidence.
2. Keep the runnable README-style ordinary-polynomial example and make the
   surrounding scope wording explain that it is part of the supported #187
   contract. The example should continue using public calls only:
   `elementary_factorization(A)` and `verify_factorization(A, factors)`.
3. Keep unsupported cases explicitly outside #187: missing evidence, unsupported
   coefficient rings, arbitrary Laurent `GL_n`, ToricBuilder mainline
   acceptance, and Steinberg factor-count optimization (#188).
4. Extend `test/expert/documentation_smoke.jl` so it reads README and Documenter
   text directly. The smoke test must fail if #187 is still listed as staged or
   if the docs imply arbitrary Laurent `GL_n`, ToricBuilder `case_008`,
   unsupported coefficient rings, or factor-count optimization are included in
   the #187 acceptance claim.

This is preferable to touching implementation or fixtures because #273 is a
public wording issue and the executable evidence already landed in #271/#272.
It is also preferable to README-only edits because the issue explicitly requires
Documenter scope language and negative-control protection.

## Behavioral Requirements

The public docs must state these supported inputs:

- exact field-backed ordinary polynomial rings;
- determinant-one `SL_3` inputs whose evidence-backed #184 route replays before
  factors are returned;
- determinant-one `SL_n`, `n > 3`, inputs whose recursive peel route verifies
  #185 ECP evidence, final #184 `SL_3` evidence, and #186 mainline provenance;
- returned public factors verified by `verify_factorization(A, factors)`.

The public docs must state these boundaries:

- determinant-one ordinary-polynomial inputs missing required evidence are
  staged failures before public factors are returned;
- unsupported coefficient rings remain out of scope;
- arbitrary Laurent `GL_n` and ToricBuilder mainline acceptance remain separate
  Laurent/ToricBuilder work;
- Steinberg factor-count optimization remains #188 and is not part of #187.

## Scope Boundaries

Do not change factorization behavior, fixtures, catalog entries, or public APIs.
Do not close or claim #188. Do not claim arbitrary Laurent determinant
correction, ToricBuilder `case_008` mainline acceptance, unsupported
coefficient-ring support, or optimized factor counts.

## Tests

Focused verification commands:

```bash
julia --project=. -e 'include("test/expert/documentation_smoke.jl")'
julia --project=docs docs/make.jl
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Self-Review

The design is scoped to public documentation and documentation smoke coverage.
It updates stale #187 staged wording only after #271/#272 are merged, preserves
the runnable public API example, and keeps Laurent/ToricBuilder, unsupported
coefficient rings, and #188 optimization outside the ordinary-polynomial public
contract.
