# Issue 210 Local q(0)-Nonunit Resultant Replay Design

## Context

Issue #210 extends the existing Murthy q(0)-nonunit Bezout/resultant branch from
ordinary `QQ[X]` examples to the checked local context introduced by #208 and
the denominator-aware local factor replay from #207. Issue #209 already routes
q-degree and q(0)-unit local contexts through replayable certificates. The
existing q(0)-nonunit branch handles ordinary univariate inputs by accepting a
supplied Bezout witness or extracting one with `gcdx`, then reducing to a
q(0)-unit child.

The missing behavior is local-context routing for supplied #206-style
q(0)-nonunit witnesses. The verifier must not trust the branch tag or witness
shape. It must recompute the Bezout equality, degree guards, child q(0)-unit
condition, elementary identities, child replay, and final replay metadata.

## Approaches Considered

1. Extend the existing q(0)-nonunit reduction and certificate verifier. This is
   the chosen approach because it preserves ordinary `QQ[X]` matrix factors,
   reuses the checked context and local q(0)-unit child certificate, and keeps a
   single branch verifier responsible for replay.

2. Add a separate local-only q(0)-nonunit certificate type. This was rejected
   because it would duplicate most of the current reduction fields and make
   branch verification depend on which constructor produced the certificate.

3. Implement automatic Bezout/resultant extraction for all local coefficient
   contexts. This was rejected as out of scope. The issue requires broader local
   contexts to prefer supplied witness data and stage-fail if extraction has not
   been proven for that ring family.

## Design

Add context routing for `realize_sl3_local_certificate(context)` when
`deg(q) < deg(p)`, `q(0)` is not a global or local unit, and the checked context
contains supplied Bezout/resultant witness data. The constructor verifies the
context first, then delegates to a local q(0)-nonunit reduction helper.

The ordinary API path keeps the current behavior:

- supplied `murthy_q0_nonunit_witness` values remain accepted after exact
  verification;
- ordinary supported extraction still uses `gcdx`;
- returned `factors` remain ordinary matrices when the child q(0)-unit
  certificate is ordinary/materializable.

The local-context path does not try automatic extraction. It requires
`context.bezout_witness`, uses the normalized checked witness fields
`p_prime`, `q_prime`, `resultant`, and `branch_unit`, and constructs:

- `bezout_target = [p q; q_prime p_prime]`;
- `child_link_target = [p + q_prime, q + p_prime; q_prime, p_prime]`;
- `E21(r*p_prime - s*q_prime)`;
- `E12(-1)`;
- a q(0)-unit child context/certificate for `child_link_target`.

When `branch_unit` is only locally invertible, the child context receives the
branch-unit local witness as its q0 local-unit witness, so #209 owns all local
q(0)-unit replay. The final q(0)-nonunit certificate exposes prefix ordinary
factors followed by the child certificate factors. If the child factors are
local elementary records, final acceptance is through localized replay metadata
rather than ordinary base-ring `verify_factorization`.

## Verification

`verify_sl3_local_murthy_q0_nonunit_reduction` recomputes:

- target shape, determinant, selected variable, constants, and degree guards;
- `p_prime*p - q_prime*q == resultant == 1`;
- `branch_unit == q(0) + p_prime(0)` and either global inverse data or checked
  child q0 local-unit evidence;
- exact matrix identities
  `target == E21(r*p_prime - s*q_prime) * bezout_target` and
  `bezout_target == E12(-1) * child_link_target`;
- child q(0)-unit certificate replay and child target alignment;
- final factor replay, using ordinary multiplication when all factors are
  matrices and `SL3LocalElementaryFactorReplay` when local denominator factors
  appear.

The verifier continues to return `false` for malformed/tampered records and
rethrows `InterruptException`.

## Testing

Update `test/expert/sl3_local_murthy_resultant.jl` to keep the existing
ordinary supplied-witness and extracted `QQ[X]` cases, and add the #206 local
contract case `mg-local-q0-nonunit-bezout-at-u`.

Coverage must include:

- local context construction with the supplied Bezout/resultant witness;
- certificate branch `:murthy_q0_nonunit_bezout_resultant`;
- exact `A = E21(...)*B` and `B = E12(-1)*C` identities;
- child certificate branch `:murthy_q0_unit`, child q(0)-unit local status, and
  replay of the child local factor certificate;
- negative controls for corrupting `p_prime`, `q_prime`, resultant value,
  degree guard metadata, child local-unit evidence, and the child certificate.

Focused verification:

```bash
julia --project=. -e 'include("test/expert/sl3_local_murthy_resultant.jl")'
```

Package verification:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Out Of Scope

- No Quillen local-to-global patching or general `SL_3` driver.
- No automatic Bezout/resultant extraction for arbitrary local coefficient
  rings.
- No public API export changes.
- No change to ordinary `verify_factorization` semantics.

## Spec Self-Review

- No placeholder markers remain.
- Scope is limited to q(0)-nonunit local replay and focused verification.
- Existing ordinary univariate extraction and supplied-witness behavior are
  preserved.
- Unsupported automatic local extraction remains a staged failure.
