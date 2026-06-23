# Issue 104 Quillen Patch Verification Hardening Design

## Context

Issue #104 is the verifier-hardening child of the constructive Quillen
patching parent #63. Its direct dependency #103 is merged in this checkout
and provides `QuillenGlobalPatchAssembly`, deterministic assembly from local
certificates, denominator covers, normalized local contributions, global
elementary factors, exact products, targets, replay metadata, and
`verify_quillen_patch(::QuillenGlobalPatchAssembly)`.

The repository has no `AGENTS.md`, `CLAUDE.md`, `GEMINI.md`, or `CODEX.md`
file in this checkout. The worker branch is already an isolated Agent Desk
worktree on
`agent/issue-104-harden-quillen-patch-verification-against-tamper-run-1`.

## Approaches Considered

1. Harden the existing #103 global patch replay in
   `src/algorithm/quillen_induction.jl`. Add exact replay helpers where the
   verifier needs to compare stored theorem-bearing fields against
   independently recomputed cover, local certificate, normalization,
   factor-product, target, and metadata values.
2. Move all Quillen verifier internals from `src/core/groebner_tools.jl` into
   `src/algorithm/quillen_induction.jl`. This could make long-term ownership
   clearer, but it risks broad churn across the legacy toy `QuillenPatch`
   tests and public exports while #104 only needs hardening.
3. Add new public verifier-related types or exported helpers. This is
   unnecessary because #103 already keeps the deterministic assembly names
   expert/internal and expert tests can use qualified `Suslin.<name>` access.

The chosen design is approach 1. It keeps the blast radius small, preserves
the existing toy exact patching tests, leaves the public API unchanged, and
hardens the #103 global patch verifier where the tamper surface exists.

## API Surface

No public API changes.

The existing exported `verify_quillen_patch` remains the verification entry
point. The #103 expert/internal names remain unexported:

- `Suslin.QuillenGlobalPatchAssembly`
- `Suslin.QuillenGlobalPatchAssemblyVerification`
- `Suslin.assemble_deterministic_quillen_patch`
- `Suslin.replay_deterministic_quillen_patch`

The new hardening tests use qualified `Suslin.<name>` access for those
expert/internal names. `src/Suslin.jl` and `test/public/api_surface.jl` do not
change.

## Verification Semantics

For `QuillenGlobalPatchAssembly`, `verify_quillen_patch` must return `true`
only when replay independently recomputes and matches all stored
theorem-bearing fields:

- the denominator cover certificate and its coverage sum,
- every local realization certificate and its stored verification summary,
- every patched-substitution witness summary,
- every normalized local contribution and its stored verification summary,
- denominator data derived from normalized local contributions,
- global elementary factors derived from normalized local contributions,
- the factor product derived from the stored global factors,
- target equality between the replayed product and the stored target,
- replay metadata derived from cover, local certificate, and normalization
  data,
- the stored global patch verification summary.

Malformed or tampered patch records return `false`. Constructors continue to
throw `ArgumentError` when supplied data is invalid.

## Implementation Notes

`replay_deterministic_quillen_patch` already recomputes most fields, but the
hardening pass should make the replay contract explicit and resistant to
tampered upstream data even when downstream fields are manually edited back to
look consistent. In particular:

- ensure product comparison uses a recomputed product from the stored global
  factor list, not only a stored product field,
- ensure target comparison is part of `overall_ok`,
- ensure stored global factors are exactly the concatenation of the replayed
  normalized weighted global elementary factors,
- ensure local certificates and normalized contributions are aligned by exact
  stored data and each child verifier matches its own replay,
- ensure replay metadata is exactly recomputed from child records,
- keep `try`/`catch` verifier behavior so malformed records return `false`
  except for `InterruptException`.

If tests reveal a specific field is already covered, keep the implementation
minimal and avoid rewriting stable verifier logic.

## Tests

Add `test/expert/quillen_patch_verification_hardening.jl`. It should build at
least one valid deterministic global patch from the existing fixture catalog
and assert `Suslin.verify_quillen_patch(patch) == true`.

The same file then rebuilds tampered copies of the patch and asserts
verification returns `false` for each of these changes:

- one cover multiplier,
- one local certificate factor,
- one patched-substitution witness field,
- one normalized contribution denominator,
- one global elementary factor,
- the stored patched product,
- the stored global verification summary.

Include a negative control where an upstream local certificate remains
corrupted but the final stored product field is manually set back to the
target. The verifier must still reject it.

Register the new expert test in `test/runtests.jl` immediately after
`expert/quillen_global_patch_assembly.jl`.

## Files

- Modify `src/algorithm/quillen_induction.jl`: harden replay/verifier logic
  only where tests expose a gap.
- Create `test/expert/quillen_patch_verification_hardening.jl`: focused
  tamper-resistance coverage for #104.
- Modify `test/runtests.jl`: register the expert test.
- Do not modify `src/Suslin.jl` or `test/public/api_surface.jl`.

## Verification

Focused command:

```bash
julia --project=. -e 'include("test/expert/quillen_patch_verification_hardening.jl")'
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
expert/internal verifier-hardening approach, avoids public API churn, preserves
existing toy patching tests, and stays within #104's requested tamper surfaces.

## Spec Self-Review

- No incomplete markers remain.
- Scope is limited to verifier hardening and focused expert tests.
- The API surface decision is explicit and avoids exported names.
- The negative control proves verification is not merely final matrix equality.
