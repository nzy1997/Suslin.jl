# Issue 105 Constructive Quillen Induction Acceptance Design

## Context

Issue #105 is the final acceptance issue for the constructive Quillen
patching parent #63. The checkout already contains the merged dependencies:

- #99 fixture catalog entries in `test/fixtures/quillen_patch_cases.jl`,
- #100 local realization certificates,
- #101 denominator cover certificates,
- #102 normalized local contributions,
- #103 deterministic global patch assembly, and
- #104 hardened global patch verification.

The repository has no `AGENTS.md`, `CLAUDE.md`, `GEMINI.md`, or `CODEX.md`
file in this checkout. The worker branch is already an isolated Agent Desk
worktree on
`agent/issue-105-add-the-final-constructive-quillen-induction-acc-run-1`.

## Approaches Considered

1. Add one focused expert acceptance test that builds the full constructive
   path from existing fixture entries: denominator cover, local certificates,
   normalized contributions, deterministic global patch assembly, replay, and
   final exact factor product. Register it in the expert group and keep all
   names expert/internal.
2. Add a new production wrapper named after constructive Quillen induction.
   This would make the acceptance test read like a single command, but it
   would introduce API surface and behavior beyond what #105 needs.
3. Route through the final public `elementary_factorization` driver. This is
   explicitly out of scope because the public driver belongs to #64.

The chosen design is approach 1. It proves the merged expert path works
together, avoids public API churn, does not solve new local realizability, and
does not route through the final public driver.

## API Surface

No public API changes.

The acceptance test uses qualified expert/internal names already present after
#103 and #104:

- `Suslin.quillen_denominator_cover_certificate`
- `Suslin.quillen_local_realization_certificate`
- `Suslin.normalize_quillen_local_contributions`
- `Suslin.assemble_deterministic_quillen_patch`
- `Suslin.replay_deterministic_quillen_patch`
- `Suslin.verify_quillen_patch`

`src/Suslin.jl` and `test/public/api_surface.jl` do not change.

## Acceptance Semantics

The new file `test/expert/quillen_induction_constructive.jl` loads the #99
fixture catalog and assembles at least two deterministic ordinary-polynomial
examples:

- `quillen-patched-substitution-witness-qq`, which exercises the
  patched-substitution witness replay chain over a two-open QQ cover, and
- `quillen-constructive-acceptance-gf2`, which exercises the constructive
  acceptance fixture over a verified finite GF(2) cover.

For each fixture, the test constructs:

- a verified denominator cover from the fixture denominator data,
- replayable local realization certificates from the fixture local factors,
- normalized local contributions from the cover and local certificates, and
- a deterministic global Quillen patch from the normalized data.

The test asserts exact replay at each theorem-bearing stage:

- denominator coverage equals one and matches the cover certificate,
- every local certificate verifies and its local product equals its stored
  correction,
- every patched-substitution witness summary replays,
- every normalized local contribution verifies and replays,
- normalized weighted factors concatenate to the patch global factor list,
- the global factor product equals the intended target matrix, and
- `Suslin.verify_quillen_patch(patch)` and replay `overall_ok` are true.

Negative controls must prove the acceptance harness rejects false positives:

- an invalid cover whose coverage sum no longer equals one must either fail
  construction or fail verification, and
- tampering one local factor or one patched-substitution witness field after
  construction must make `Suslin.verify_quillen_patch` return `false`.

## Test Runner Wiring

Register `expert/quillen_induction_constructive.jl` in the expert group in
`test/runtests.jl` immediately after the #104 hardening test. The current
expert group already includes the focused #99 through #104 files:

- `expert/quillen_denominator_cover.jl`,
- `expert/quillen_local_certificate.jl`,
- `expert/quillen_contribution_normalization.jl`,
- `expert/quillen_global_patch_assembly.jl`, and
- `expert/quillen_patch_verification_hardening.jl`.

## Files

- Create `test/expert/quillen_induction_constructive.jl`: final #63 expert
  acceptance test.
- Modify `test/runtests.jl`: register the new expert acceptance file.
- Do not modify `src/Suslin.jl` or `test/public/api_surface.jl`.
- Modify production code only if the acceptance test exposes a real gap in the
  already merged expert path.

## Verification

Focused command:

```bash
julia --project=. -e 'include("test/expert/quillen_induction_constructive.jl")'
```

Full suite command from #105:

```bash
julia --project=. test/runtests.jl all
```

Package command required by Agent Desk:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Automatic Approval

This Agent Desk run is non-interactive. Under the Standing Answer Policy, the
design is approved automatically because it chooses the conservative
expert/internal acceptance-test path, reuses the merged #99 through #104
fixtures and verifiers, avoids public API churn, and stays out of #64.

## Spec Self-Review

- No incomplete markers remain.
- Scope is limited to the final expert acceptance harness and runner wiring.
- The API surface decision is explicit and avoids exported names.
- The negative controls cover both broken cover data and tampered upstream
  replay data.
