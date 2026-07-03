# Issue 272 Final Negative Controls Design

Issue #272 hardens the final Park-Woodburn #187 public boundary with negative
controls. The public route is credible only if invalid inputs fail before
factors are returned and if corrupted returned factors or nested evidence cannot
pass exact verification.

## Context

There is no repository `AGENTS.md`; the README, existing Superpowers specs, and
current test layout apply. Live GitHub issue fetching is unavailable in this
Agent Desk sandbox because `gh issue view` cannot reach the configured proxy, so
the supplied #272 issue body plus local #270/#271 merge history are the task
context.

The #270 catalog already records final #187 positive cases and three negative
controls. The #271 tests already prove public success for the final positive
cases and include earlier tamper coverage. This issue should make the final
negative surface explicit and required: determinant-not-one, unsupported
coefficient ring, missing `SL_3` local-form or Quillen evidence, missing ECP
peel evidence, missing final `SL_3` evidence in recursive `SL_n`, the
Laurent/ToricBuilder ordinary-public boundary, returned factor tampering, and
nested route evidence tampering.

## Chosen Approach

Extend tests and catalog metadata without broadening production support.

1. Extend `test/fixtures/park_woodburn_mainline_acceptance_cases.jl` with final
   negative-control entries for the missing evidence and Laurent boundary cases
   that #272 names explicitly.
2. Extend `test/internal/park_woodburn_mainline_acceptance_fixtures.jl` so the
   final negative-control ids are required and the validator still proves every
   negative catalog entry fails validation.
3. Extend `test/public/factorization_driver_shell.jl` and
   `test/public/park_woodburn_polynomial_factorization.jl` with helper assertions
   that call `elementary_factorization(A)`, prove no factors are returned, check
   staged route evidence when applicable, and check the intended staged-boundary
   phrases for each invalid public case.
4. Extend `test/expert/park_woodburn_route_certificate.jl` with one final
   nested tamper gate if the existing #271 corruption checks do not cover the
   exact #272 route shape.

This is preferable to changing route implementation first because the issue is a
negative-control coverage issue and the current code already has staged errors
and certificate verifiers. Production changes should only be made if a
test-first red check exposes a real leak.

## Behavioral Requirements

For unsupported public inputs, tests must call `elementary_factorization(A)` and
assert that it throws `ArgumentError` before returning factors. The assertions
must check both the staged boundary category and that the captured factor value
remains `nothing`.

For staged polynomial routes, tests must also call
`Suslin._polynomial_factorization_route_certificate(A)` when that internal route
is supported and prove it returns `:staged_failure` with the expected
`reason_code`. Determinant-not-one and Laurent cases can fail before a staged
certificate exists.

For tampering, tests must obtain factors from a final positive #187 case, corrupt
one returned elementary factor, and assert
`verify_factorization(A, corrupted_factors) == false`. A nested route
certificate must be rebuilt with corrupted #184/#185/#186 evidence while keeping
the public factor sequence/product intact, and the relevant route verifier must
return `false`.

## Scope Boundaries

Do not broaden accepted inputs. Do not add Laurent/ToricBuilder ordinary-public
support. Do not treat identity or legacy recursive examples as final #187
acceptance. Do not change public APIs. Keep helper code inside the named test
files unless a test-first failure exposes a production verifier bug.

## Tests

Focused verification commands:

```bash
julia --project=. -e 'include("test/public/factorization_driver_shell.jl")'
julia --project=. -e 'include("test/public/park_woodburn_polynomial_factorization.jl")'
julia --project=. -e 'include("test/expert/park_woodburn_route_certificate.jl")'
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Self-Review

The design is scoped to final negative-control fixtures and tests. It names each
negative case required by #272, requires real public API failure before factor
return, and keeps Laurent/ToricBuilder and unsupported coefficient rings outside
the accepted public surface.
