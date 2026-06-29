# Issue 209 Local q-Degree and q(0)-Unit Replay Design

## Context

Issue #209 connects the existing Murthy q-degree normalization and q(0)-unit
recursive branch to the checked local context from #208 and the denominator-aware
local elementary factor replay from #207. Existing ordinary `QQ[X]` callers must
keep receiving the same materialized matrix factors. Nontrivial local-unit
inputs, such as the #206 `QQ[u, X]` q(0)-unit contract case with `q(0) = 1 + u`,
must be accepted by certificate replay rather than by ordinary
`verify_factorization` over the original polynomial ring.

The existing implementation already has:

- `SL3LocalMurthyInputContext`, with verified monicity, degrees, constants,
  global/local unit classification, local-unit witnesses, and optional split or
  Bezout witnesses;
- `SL3LocalElementaryFactor` and `SL3LocalElementaryFactorReplay`, which verify
  nontrivial local elementary factors by denominator-cleared replay;
- ordinary univariate q-degree and q(0)-unit certificates for `QQ[X]`.

## Approaches Considered

1. Add context-consuming internal helpers while preserving the existing public
   matrix-factor APIs. This is the chosen approach. It keeps ordinary
   denominator-one behavior unchanged, gives later #182 work a checked local
   replay object, and limits local acceptance to verifier-owned certificates.

2. Change `realize_sl3_local(A, X)` to return mixed ordinary matrices and local
   denominator records. This was rejected because existing callers expect
   ordinary matrix factors and the issue says nontrivial local-unit acceptance is
   localized replay, not ordinary factorization.

3. Extend `verify_factorization` to understand local denominators. This was
   rejected as a public API change outside this issue and would blur the
   boundary created by #207.

## Design

Add context methods in `src/algorithm/sl3_local.jl`:

```julia
sl3_local_q_degree_normalization(context::SL3LocalMurthyInputContext)
sl3_local_q_degree_normalization_certificate(context::SL3LocalMurthyInputContext)
realize_sl3_local_certificate(context::SL3LocalMurthyInputContext)
```

The q-degree context method first verifies the context and requires
`degree_q >= degree_p`. It reuses the existing q-degree record and certificate
types, so the local q-degree case remains an exact polynomial replay over the
original ring.

For q(0)-unit contexts with a global unit and ordinary univariate ring, the
context certificate delegates to the current ordinary branch. For a nontrivial
local unit, add an internal `SL3LocalMurthyQUnitLocalReduction` record. It
stores the checked context, `p0`, `q0`, a local representation of `q0^{-1}`, the
right `E21(-p0/q0)` elimination factor and inverse correction as
`SL3LocalElementaryFactor` records, the exact local factor replay, and the
source Murthy certificate used to align the q(0)-unit split children.

The local construction maps the checked `QQ[u, X]` style target to a univariate
polynomial ring over the fraction field of the coefficient variables. The
ordinary Murthy branch computes the recursive q(0)-unit and split-lemma replay
there. Each elementary factor coefficient is then translated back to a
numerator/denominator pair over the original ring, and every non-one denominator
receives a derived localization-at-maximal-ideal witness checked by the same
#208 local-unit verifier. The final acceptance check is the #207
denominator-cleared identity against the original target.

`verify_sl3_local_realization` remains unchanged for ordinary certificates. When
the q(0)-unit witness reduction is a local reduction, it verifies the local
reduction and uses `verify_sl3_local_elementary_factor_replay` instead of
ordinary `verify_factorization`.

Unsupported local contexts fail deliberately with `ArgumentError`: missing q0
local-unit evidence, q(0)-nonunit inputs, unsupported coefficient localization
shape, or denominators whose local-unit witness cannot be derived from the
checked context.

## Testing

Update:

- `test/expert/sl3_local_q_degree_normalization.jl`
- `test/expert/sl3_local_murthy_q_unit.jl`
- `test/expert/sl3_local_murthy_context.jl`

Coverage:

- existing ordinary `QQ[X]` q-degree and q(0)-unit cases still return ordinary
  matrix factors and pass `verify_factorization`;
- #206 local q-degree context returns a q-degree certificate whose replay
  verifies over `QQ[u, X]`;
- #206 local q(0)-unit context returns a `:murthy_q0_unit` certificate whose
  factors are local elementary records, whose local factor replay reconstructs
  the original target, and whose source split certificate verifies child
  alignment;
- corrupting `q0_inverse`, the right `E21` coefficient, the source split child
  target, or a local factor makes verification reject the certificate.

Focused verification:

```bash
julia --project=. -e 'include("test/expert/sl3_local_q_degree_normalization.jl")'
julia --project=. -e 'include("test/expert/sl3_local_murthy_q_unit.jl")'
```

Package verification:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Out Of Scope

- Do not implement the q(0)-nonunit local Bezout/resultant branch.
- Do not integrate with Quillen patching or public Park-Woodburn routing.
- Do not change `verify_factorization` or the ordinary matrix-factor contract.
- Do not expose new public APIs from `src/Suslin.jl`.

## Spec Self-Review

- No placeholders or unresolved questions remain.
- The scope is limited to context-consuming q-degree/q(0)-unit replay.
- Existing `QQ[X]` behavior is explicitly preserved.
- Nontrivial local-unit acceptance is tied to denominator-cleared replay, not
  ordinary base-ring factorization.
