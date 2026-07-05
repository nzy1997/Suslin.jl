# Issue 294 Steinberg Closeout Gate Design

Issue #294 is the parent-level closeout gate for #188. It must document the
optional Steinberg optimizer contract, prove the accepted safe rewrite set, and
make clear that `elementary_factorization(A)` keeps its existing correctness
contract and does not optimize by default.

## Context

There is no repository `AGENTS.md`; the README, docs index, existing audit
pages, and merged Steinberg optimizer issues are the working context. The
branch is already an isolated Agent Desk worktree. Live `gh issue view` failed
because the configured GitHub proxy is not reachable, and public web lookup did
not return accessible issue content, so the supplied issue body plus local
merge history and Superpowers specs for #288-#293 are the source of truth.

Dependencies are present on `main`:

- #288 added the Section 6 Steinberg fixture catalog.
- #289 added private canonical elementary-factor records.
- #290 added `SteinbergOptimizationCertificate` and exact replay verification.
- #291 added adjacent identity, merge, and cancellation rewrites.
- #292 added conservative commutator rewrites.
- #293 exported the opt-in public optimizer and public verifier.

Park-Woodburn Section 6 lists the Steinberg relations used here: identity
removal, same-position merge, two nontrivial commutator relations, and the
disjoint commutator identity. This issue documents that supported subset; it
does not add a new algorithmic route.

The run is non-interactive under Agent Desk, so the Standing Answer Policy is
used for design approval. No visual companion is needed because the design is
documentation and tests only.

## Approach Options

Recommended: update README and docs index with a short opt-in optimizer
contract, add one #188 audit page with an upstream map and before/after table,
and extend existing expert smoke tests to enforce wording and table content.
This matches the repository's #187 closeout pattern and keeps implementation
scoped to acceptance evidence.

Alternative: add a dedicated docs page under `docs/src/`. That would be useful
for a larger public guide, but #294 asks for a closeout gate and the existing
README/docs index scope section is enough for this contract.

Alternative: change `elementary_factorization(A)` to accept an optimization
keyword and document it. That conflicts with #293 and #294 because the optimizer
must remain explicitly opt-in through `optimize_elementary_factor_sequence`.

## Chosen Design

Update `README.md` and `docs/src/index.md` with a new scope bullet for #188.
The wording must state that:

- `optimize_elementary_factor_sequence(factors; rules = :safe)` is optional
  and opt-in;
- `elementary_factorization(A)` is unchanged and does not optimize by default;
- every optimized sequence is accepted only after exact product verification
  through `verify_steinberg_optimization_certificate`;
- the safe set is limited to `:identity_removal`, `:same_position_merge`,
  `:inverse_cancellation`, `:commutator_forward`, `:commutator_reverse`, and
  `:disjoint_commutator_identity`;
- #188 does not claim global minimum factor counts, Laurent `GL_n` support, or
  ToricBuilder mainline support.

Add `docs/audits/2026-07-04-issue-188-steinberg-optimization.md`. It should
cite Park-Woodburn Section 6, map #188 dependencies (#288-#293), list the safe
rewrite set, include a compact before/after comparison table for one documented
small example, and close with explicit non-claims. The documented example uses
the #288 `steinberg-commutator-forward-qq` sequence:

| Metric | Original | Optimized | Delta |
| --- | --- | --- | --- |
| Factor count | 4 | 1 | -3 |
| Max monomial degree | 1 | 2 | 1 |
| Total off-diagonal monomial count | 6 | 2 | -4 |

The applied rewrite is `:commutator_forward`; exact products must be equal and
the verification status must be true.

Extend `test/expert/steinberg_factor_count_optimization.jl` with public
acceptance coverage over all #288 positive catalog sequences. For each accepted
safe rule, call `optimize_elementary_factor_sequence`, require public verifier
success, exact product equality, expected optimized factors, expected rule
names, and a lower factor count. Add a corrupted optimized sequence negative
control that must fail verification.

Extend `test/expert/documentation_smoke.jl` to read README, docs index, and the
new audit page. It must require the opt-in wording, exact verification wording,
safe rewrite names, the before/after table values, and the closeout statement
that #188 does not change `elementary_factorization(A)`. It must fail if docs
claim optimization is enabled by default, claim global minimum factor counts, or
fold Laurent/ToricBuilder support into #188.

No source code change is needed unless the new acceptance tests expose a bug in
the already-public #293 optimizer.

## Tests

Focused red-green verification:

```bash
julia --project=. -e 'include("test/expert/steinberg_factor_count_optimization.jl")'
julia --project=. -e 'include("test/expert/documentation_smoke.jl")'
```

Issue-required verification:

```bash
julia --project=. -e 'include("test/expert/steinberg_factor_count_optimization.jl")'
julia --project=. -e 'include("test/public/api_surface.jl")'
julia --project=. -e 'include("test/expert/documentation_smoke.jl")'
julia --project=. test/runtests.jl all
julia --project=. -e 'using Pkg; Pkg.test()'
```

All commands must exit 0 before creating the pull request.

## Out Of Scope

Do not add new Park-Woodburn correctness modules. Do not broaden supported
rings. Do not claim performance benchmarks or global optimality. Do not change
`elementary_factorization(A)` semantics or add optimizer keywords to it.

## Self-Review

This design has no placeholders, follows the existing parent closeout audit
pattern, keeps #188 limited to optional factor-sequence optimization, and ties
every requirement to documentation, audit evidence, or focused expert tests.
