# Issue 293 Public Steinberg Optimizer Design

Issue #293 exposes an opt-in public API for the already-built Steinberg
factor-sequence optimizer. The API must optimize only when explicitly called,
return a certificate with clear before/after evidence, and leave
`elementary_factorization(A)` unchanged.

## Context

There is no repository `AGENTS.md`; the README, existing Superpowers docs, and
the merged Steinberg optimizer issues are the working context. The branch is an
isolated Agent Desk worktree. Baseline `julia --project=. -e 'using Pkg;
Pkg.test()'` passed before edits.

Dependencies are present on `main`:

- Issue #288 added the ordinary-polynomial Steinberg fixture catalog.
- Issue #290 added `SteinbergOptimizationCertificate`,
  `_steinberg_optimization_certificate`, and the private verifier.
- Issue #291 added the adjacent identity, merge, and cancellation pass.
- Issue #292 added the conservative commutator pass.
- PR #179 added the metric helpers used by the certificate summary.

The run is non-interactive under Agent Desk, so the Standing Answer Policy is
used for approvals. No visual companion is needed because the design has no UI
or diagram-dependent decision.

## Approach Options

Recommended: add a small public wrapper in `src/algorithm/redundancy.jl` that
validates `rules = :safe`, runs the existing adjacent pass and commutator pass,
and returns a `SteinbergOptimizationCertificate`. Export the wrapper and a
public verifier from `src/Suslin.jl`. This reuses the existing exact replay
certificate, avoids a second certificate shape, and keeps
`elementary_factorization(A)` semantics unchanged.

Alternative: export the private adjacent and commutator helpers directly. That
would expose implementation details and make it harder to adjust pass ordering
later.

Alternative: add an `optimize = false` keyword to `elementary_factorization(A)`.
That conflicts with the issue's requirement to leave factorization semantics
unchanged and not optimize by default.

## Chosen Design

Export:

```julia
optimize_elementary_factor_sequence(factors; rules = :safe)
verify_steinberg_optimization_certificate(certificate)::Bool
```

`optimize_elementary_factor_sequence` accepts an ordinary-polynomial elementary
factor sequence. It supports only `rules = :safe` for this issue. Any other
rule set throws `ArgumentError("unsupported Steinberg optimization rule set;
supported rule sets: :safe")`.

The safe rule set is conservative and contains only:

- the adjacent identity, same-position merge, and inverse-cancellation rules
  from #291;
- the exact four-factor commutator rules from #292.

The public wrapper first validates the input through the existing Steinberg
sequence context. This preserves the deliberate `ArgumentError` behavior for
non-elementary matrix sequences, Laurent rings, inexact coefficient rings, and
mixed sizes or rings. It then runs the adjacent pass, runs the commutator pass
on the adjacent result, and builds a final certificate from the original
factors, final optimized factors, and public rule metadata.

To keep the public certificate easy to read, the final applied rewrite log uses
one `:safe_steinberg_optimization` record with:

- `original_factor_count` and `optimized_factor_count` for the full sequence;
- `metadata.rules = :safe`;
- `metadata.passes`, a tuple of pass summaries that preserves each underlying
  adjacent and commutator rewrite record with its rule name, spans, counts, and
  metadata.

This avoids pretending that pass-local spans from the commutator pass refer to
the original sequence, while still exposing the detailed rules that fired. The
certificate's `comparison_summary.applied_rewrites` mirrors the same metadata,
and the summary keeps the #290 fields for before/after factor counts, factor
count delta, PR #179 metrics, exact products, product equality, and verification
status. The summary also gains lightweight derived fields for rule names and
metric deltas so users can answer "what changed?" without inspecting matrices.

`verify_steinberg_optimization_certificate` delegates to the existing private
verifier. It returns `false` for tampered optimized factors, stale summaries,
stale products, malformed logs, or invalid matrix sequences, and rethrows
interrupts.

## Tests

Extend `test/public/api_surface.jl` to prove the optimizer and verifier are
exported and bound to `Suslin`.

Extend `test/expert/steinberg_factor_count_optimization.jl` with public API
coverage:

- call `optimize_elementary_factor_sequence` on a #288 fixture catalog
  sequence and observe a lower factor count;
- verify the returned certificate with
  `verify_steinberg_optimization_certificate`;
- check exact product equality;
- assert the comparison summary reports before/after factor count, factor
  count delta, PR #179 before/after metric summaries, and safe-rule metadata;
- pass a non-elementary matrix sequence and require a deliberate
  `ArgumentError`;
- tamper returned optimized factors and require the public verifier to return
  `false`;
- call `elementary_factorization(A)` normally and verify the factors without
  requiring optimization.

## Verification

Focused verification:

```bash
julia --project=. -e 'include("test/public/api_surface.jl")'
julia --project=. -e 'include("test/expert/steinberg_factor_count_optimization.jl")'
```

Full verification:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

All commands must exit 0. The public factorization call must continue to verify
without invoking optimization.

## Out Of Scope

Do not optimize by default. Do not add a keyword to
`elementary_factorization(A)`. Do not add Laurent or ToricBuilder optimizer
claims. Do not make performance claims or search globally for an optimal factor
count.

## Self-Review

This design has no placeholders, keeps the API opt-in, reuses the existing
certificate and exact verifier, and ties every issue requirement to either the
public wrapper, exported verifier, summary fields, or focused tests. The only
public rule set is deliberately conservative and composed of the already-merged
#291 and #292 rules.
