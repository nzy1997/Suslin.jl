# Issue 101 Quillen Denominator Cover Design

## Context

Issue #101 is the first cover-and-normalization child of the constructive
Quillen patching parent issue #63. Its dependency #99 is closed and merged, and
this checkout already contains the fixture catalog at
`test/fixtures/quillen_patch_cases.jl`. That catalog records exact ordinary
polynomial denominator data such as `r, 1-r` and nontrivial Bezout multipliers.

The existing public `construct_quillen_patch` path already checks a coverage
sum for supplied local contributions, but it does not expose a standalone cover
certificate that later normalization and final verification stages can replay.
This issue should add that reusable certificate without broadening into
Nullstellensatz solving or global Quillen factor assembly.

No `AGENTS.md`, `CLAUDE.md`, or `GEMINI.md` file is present in this checkout.

## Approaches Considered

1. Add a small expert/internal cover certificate constructor and replay helper
   in `src/algorithm/quillen_induction.jl`. Callers provide the ring,
   denominators, and coverage multipliers; the constructor normalizes them into
   the ring and verifies the exact identity `sum(c_i * d_i) == 1`.
2. Keep the certificate logic only in test validators. This would satisfy the
   fixture catalog but would not give later Quillen stages a reusable runtime
   object.
3. Add a solver that computes coverage multipliers from denominators. This is
   outside the issue and would create a general ideal-cover engine before the
   deterministic fixture-backed contract is stable.

The chosen design is approach 1. It is deterministic, exact, reusable by later
stages, and does not alter the exported public API.

## API Surface

The new names remain expert/internal and are intentionally not exported from
`src/Suslin.jl`:

- `QuillenDenominatorCoverVerification`
- `QuillenDenominatorCoverCertificate`
- `quillen_denominator_cover_certificate(R, denominators, coverage_multipliers)`
- `replay_quillen_denominator_cover(certificate)`
- `verify_quillen_denominator_cover(certificate)`

Expert tests use qualified `Suslin.<name>` access. `test/public/api_surface.jl`
does not change.

## Certificate Semantics

`quillen_denominator_cover_certificate` accepts an ordinary exact polynomial
ring, a finite nonempty list of denominators, and a same-length list of supplied
coverage multipliers. It coerces every denominator and multiplier into the
target ring, computes exact coverage terms `c_i * d_i`, and records:

- the normalized denominators,
- the normalized coverage multipliers,
- the exact coverage sum,
- verification metadata containing counts, parent/exact-ring flags, coverage
  terms, the replayed coverage sum, and a boolean coverage result.

Construction throws `ArgumentError` for unsupported or inexact rings,
mismatched lengths, empty covers, non-coercible inputs, or coverage sums not
equal to `one(R)`.

`replay_quillen_denominator_cover` recomputes the verification metadata from a
stored certificate. `verify_quillen_denominator_cover` returns `true` only when
the replayed metadata matches the stored metadata and the replay proves exact
coverage.

## Tests

Add `test/expert/quillen_denominator_cover.jl`. The test includes the #99
fixture catalog when needed and builds at least two positive certificates:

- `quillen-two-open-cover-qq`, whose denominators are exactly `r` and `1-r`,
  with multipliers `1, 1`.
- `quillen-nontrivial-multipliers-qq`, where at least one multiplier is not
  `1`.

For both certificates, the test asserts that replay recomputes the same exact
coverage sum equal to `one(R)`, that the certificate verifies, and that stored
denominators and multipliers equal the normalized ring elements with the
correct parent ring.

Negative controls cover removing a denominator, mutating a multiplier so the
sum is not one, building over an inexact polynomial ring, and replaying a
tampered certificate.

## Files

- Modify `src/algorithm/quillen_induction.jl`: add cover certificate structs,
  constructor, replay, and verifier.
- Create `test/expert/quillen_denominator_cover.jl`: focused positive and
  negative coverage.
- Modify `test/runtests.jl`: register the expert test.
- Do not modify `src/Suslin.jl` or `test/public/api_surface.jl`.

## Verification

Focused command:

```bash
julia --project=. -e 'include("test/expert/quillen_denominator_cover.jl")'
```

Expert group command:

```bash
julia --project=. test/runtests.jl expert
```

Package command required by Agent Desk:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Automatic Approval

This Agent Desk run is non-interactive. Under the Standing Answer Policy, the
design is approved automatically because it chooses the conservative internal
certificate API, reuses merged deterministic fixtures, avoids public API churn,
and avoids implementing a broad cover solver.

## Spec Self-Review

- No incomplete markers remain.
- Scope is limited to denominator cover certificates and replay.
- The API surface decision is explicit and avoids exported names.
- Negative controls cover uncovered, mutated, inexact-ring, and tampered replay
  cases.
