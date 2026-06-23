# Issue 103 Quillen Global Patch Assembly Design

## Context

Issue #103 is the assembly child of the constructive Quillen patching parent
#63. Its direct dependency #102 is merged in this checkout and provides
replayable `QuillenLocalContributionNormalization` records, constructed from
#100 local certificates and #101 denominator cover certificates.

The repository has no `AGENTS.md`, `CLAUDE.md`, or `GEMINI.md` file in this
checkout. The worker branch is already an isolated Agent Desk worktree on
`agent/issue-103-assemble-deterministic-quillen-patches-into-glob-run-1`.

## Approaches Considered

1. Add a new expert/internal global assembly record in
   `src/algorithm/quillen_induction.jl`. The record consumes verified local
   realization certificates, a verified denominator cover, and verified
   normalized local contributions. Verification replays cover, local
   certificate, normalization, global factor product, target equality, and
   replay metadata.
2. Extend the public `QuillenPatch` struct with normalization and replay
   fields. This matches the long-term shape, but it risks public field churn
   and direct-constructor compatibility before the final Park-Woodburn driver
   is ready.
3. Route assembly through `elementary_factorization` or solve new local
   certificates internally. This is out of scope because #103 is only supposed
   to assemble deterministic data supplied by earlier stages.

The chosen design is approach 1. It is the conservative path: keep the new
constructor expert/internal, reuse `verify_quillen_patch` by adding a method
for the new record, and leave the public API surface unchanged.

## API Surface

Add expert/internal names, intentionally not exported from `src/Suslin.jl`:

- `QuillenGlobalPatchAssembly`
- `QuillenGlobalPatchAssemblyVerification`
- `assemble_deterministic_quillen_patch(original_input, selected_variable, local_certificates, normalized_local_contributions, cover; target = original_input, ring = nothing, size = nothing)`
- `replay_deterministic_quillen_patch(patch)`

Extend the existing exported verifier with:

- `verify_quillen_patch(patch::QuillenGlobalPatchAssembly)::Bool`

Expert tests use qualified `Suslin.<name>` access for the new constructor and
record type. `test/public/api_surface.jl` does not change.

## Assembly Semantics

Construction requires:

- every local realization certificate verifies,
- the denominator cover verifies and sums to one,
- every normalized contribution verifies,
- the local certificates, normalized contributions, and cover are aligned by
  exact object content and cover pair,
- every normalized contribution records the requested original input and
  selected Quillen substitution variable,
- global elementary factors are exactly the concatenation of the normalized
  weighted global elementary factors in input order,
- the product of those global factors equals the requested target.

The returned record stores the ring, size, substitution variable, original
input, cover certificate, denominator coverage data, local certificates,
normalized local contributions, global elementary factors, patched product,
target, replay metadata, and verification data.

`replay_deterministic_quillen_patch` recomputes the same checks from stored
data. `verify_quillen_patch(patch::QuillenGlobalPatchAssembly)` returns `true`
only when the replayed verification matches the stored verification and all
component checks pass. Tampered records return `false`; construction throws
`ArgumentError` when supplied deterministic data is not already replayable or
does not assemble to the exact target.

## Tests

Add `test/expert/quillen_global_patch_assembly.jl`. It builds deterministic
ordinary-polynomial patches from two fixture catalog entries:

- `quillen-patched-substitution-witness-qq`
- `quillen-nontrivial-multipliers-qq`

For each fixture, the test builds #100 local certificates, a #101 denominator
cover, #102 normalized local contributions, and then the #103 global patch.
It checks that construction returns global elementary factors, exact factor
multiplication equals the fixture target matrix, and
`verify_quillen_patch(patch) == true`.

Negative controls cover a valid local certificate paired with the wrong
normalized denominator data and a tampered cover whose coverage sum no longer
equals one. Construction must throw `ArgumentError` or verification must return
`false`; it must not return a falsely verified patch.

Register the new expert test in `test/runtests.jl`.

## Files

- Modify `src/algorithm/quillen_induction.jl`: add the global assembly record,
  replay helper, constructor, verifier method, and exact alignment helpers.
- Create `test/expert/quillen_global_patch_assembly.jl`: focused positive and
  negative assembly coverage.
- Modify `test/runtests.jl`: register the expert test.
- Do not modify `src/Suslin.jl` or `test/public/api_surface.jl`.

## Verification

Focused command:

```bash
julia --project=. -e 'include("test/expert/quillen_global_patch_assembly.jl")'
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
design is approved automatically because it chooses the conservative
expert/internal assembly API, reuses merged deterministic replay schemas,
avoids public API churn, preserves exact product checks before returning
factors, and avoids the out-of-scope final public factorization driver.

## Spec Self-Review

- No incomplete markers remain.
- Scope is limited to deterministic global assembly and replay.
- The API surface decision is explicit and avoids exported names.
- Negative controls cover both wrong-denominator local pairing and broken
  cover coverage.
