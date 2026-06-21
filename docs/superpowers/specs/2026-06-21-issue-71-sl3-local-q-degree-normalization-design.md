# Issue 71 SL3 Local Q-Degree Normalization Design

## Context

Issue 71 adds the first Murthy-Gupta local `SL_3` reduction step for
special-form matrices

```text
[p q 0; r s 0; 0 0 1]
```

with `p` monic in the selected variable `X`. The step writes
`q = f*p + g` with `deg(g, X) < deg(p, X)`, then right-multiplies the
normalized target by `E_12(f)` to reconstruct the original target:

```text
[p g       0;     [1 f 0;     [p q 0;
 r s-f*r   0;  *   0 1 0;  =   r s 0;
 0 0       1]      0 0 1]      0 0 1]
```

Issues 69 and 70 are merged. Issue 69 supplies a shared Murthy-Gupta fixture
catalog with a q-degree normalization case. Issue 70 supplies the local
`SL_3` certificate shape and replay verifier. This issue should extend that
certificate path, not solve the normalized target recursively.

## Design Choice

Add a narrow internal q-degree normalization helper in
`src/algorithm/sl3_local.jl` and teach the existing local certificate verifier
to replay a new branch tag, `:q_degree_normalization`.

The helper returns an `SL3LocalRealizationCertificate` whose `target` is the
original special-form matrix, whose `branch` is `:q_degree_normalization`, whose
`factors` field records the elementary correction `[E_12(f)]`, and whose
`witness` records:

- original entries `p`, `q`, `r`, `s`,
- quotient `f`,
- remainder `g`,
- normalized target,
- correction `E_12(f)`,
- selected variable degree data.

This keeps the implementation replayable without pretending to factor the
normalized target. The factor product is not expected to equal the original
target for this partial Murthy step; instead replay checks the exact
normalization identity `normalized_target * correction == target`.

Alternatives considered:

- Add the step directly to `realize_sl3_local`. That would change the public
  factor-returning path and imply recursive support that is explicitly out of
  scope.
- Return a new public normalization object. That would introduce API surface
  before later Murthy branches define the long-term public shape.
- Keep the proof only in the test fixture validator. That would miss the issue
  requirement to record the step in the local realization certificate path.

## Helper Contract

The implementation adds non-exported expert helpers:

- `normalize_sl3_local_q_degree(p, q, r, s, X; check_monic=true)`
- `normalize_sl3_local_q_degree(A, X; check_monic=true)`

The helpers reuse existing local special-form recognition for ring, shape,
monicity, and determinant checks. After recognition they require an ordinary
polynomial ring, perform polynomial division in the selected variable, and
construct the normalized target and correction.

Supported inputs are polynomial-only. Laurent inputs and wrong-ring selected
variables are rejected by the existing recognition or the polynomial-only
guard. Non-monic `p` is rejected when `check_monic=true`.

## Replay Rules

`verify_sl3_local_realization(cert)` gains a branch-specific path for
`:q_degree_normalization`. Replay checks:

- witness keys are exact, so extra decorative metadata is rejected;
- `target` has local special form and determinant one;
- `selected_variable` is a generator of the target base ring;
- witness entries match `target`;
- `q == f*p + g`;
- `deg(g, X) < deg(p, X)`;
- `normalized_target` equals `[p g 0; r s - f*r 0; 0 0 1]`;
- `correction == E_12(f)`;
- `normalized_target * correction == target`;
- `factors == [correction]`.

The existing complete factorization certificate branches continue to require
`verify_factorization(target, factors)`. The q-degree branch is intentionally a
partial Murthy certificate, so its replay success is the exact normalization
identity rather than full elementary factorization of the original target.

## Tests

Add `test/expert/sl3_local_q_degree_normalization.jl` and register it in the
expert group.

The focused test includes at least two exact `QQ[X]` examples with
`deg(q, X) >= deg(p, X)` and monic `p`. One example reuses the
`mg-q-degree-normalization` fixture from issue 69. A second hand-checkable
example covers a nonzero remainder and nontrivial `s - f*r`.

For each example the test asserts:

- quotient and remainder satisfy `q == f*p + g`;
- `deg(g, X) < deg(p, X)`;
- the normalized target and correction match the expected formulas;
- `verify_sl3_local_realization(cert)` returns true;
- `normalized_target * correction == target`.

Negative controls include a non-monic `p` rejection and a deliberately corrupted
q-degree certificate whose target still has determinant one but whose witness
fails replay.

## Verification

Focused expert command:

```bash
julia --project=. -e 'include("test/expert/sl3_local_q_degree_normalization.jl")'
```

Package command required by the Agent Desk run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Full expert-inclusive command remains available as:

```bash
julia --project=. test/runtests.jl expert
```

## Spec Self-Review

- No incomplete markers remain.
- Scope is limited to polynomial q-degree normalization and replay.
- The public `elementary_factorization` driver and recursive Murthy branches are
  untouched.
- The design distinguishes complete local factorizations from this partial
  normalization certificate.
