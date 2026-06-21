# Issue 71 SL3 Local q-Degree Normalization Design

## Context

Issue 71 adds only the first Murthy-Gupta local `SL_3` reduction step for
special-form polynomial inputs

```text
[p q 0; r s 0; 0 0 1]
```

with `p` monic in the selected variable `X`. The step writes
`q = f*p + g` with `deg(g, X) < deg(p, X)` and records the exact identity

```text
[p q 0; r s 0; 0 0 1] =
[p g 0; r s - f*r 0; 0 0 1] * E_12(f).
```

Issue 69 supplies the fixture catalog and exact witness style. Issue 70
supplies the local certificate replay path and requires every stored field to
participate in verification. The current branch also contains the split-lemma
replay work from issue 72, but this issue does not depend on that branch.

## Design Choice

Add a focused internal normalization record and replay path in
`src/algorithm/sl3_local.jl`:

- `SL3LocalQDegreeNormalization`
- `sl3_local_q_degree_normalization(p, q, r, s, X; check_monic=true)`
- `sl3_local_q_degree_normalization(A, X; check_monic=true)`
- `sl3_local_q_degree_normalization_certificate(...)`
- `verify_sl3_local_q_degree_normalization(record)::Bool`

The helper is intentionally not exported from `src/Suslin.jl`; expert tests can
call it as `Suslin.sl3_local_q_degree_normalization`. The certificate helper
returns an `SL3LocalRealizationCertificate` with branch
`:murthy_q_degree_normalization`. Its replay proves the stored normalization
step by exact matrix multiplication, but it does not recursively solve the
normalized target.

Alternatives considered:

- Change `realize_sl3_local` to run normalization before the staged failure.
  That would change the public driver and violate the issue's out-of-scope
  boundary.
- Store only a named tuple witness. That would miss the issue requirement to
  record the normalized target and elementary correction as replayed data.
- Add a broad Murthy result schema now. That conflicts with the issue 70
  guardrail to keep the certificate layer thin until later branches consume it.

## Normalization Record

The record stores only replayed fields:

- `target`: the original special-form target.
- `quotient`: `f`.
- `remainder`: `g`.
- `normalized_target`: `[p g 0; r s - f*r 0; 0 0 1]`.
- `elementary_correction`: `E_12(f)`.
- `selected_variable`: `X`.

Construction validates that the input is an ordinary polynomial-ring
special-form `SL_3` target, that `p` is monic in `X` when `check_monic=true`,
and that exact replay succeeds. Laurent inputs remain outside this helper.

Polynomial division uses an exponent-vector coefficient pass in the selected
variable. Because `p` is monic in `X`, the algorithm can repeatedly subtract
the leading-in-`X` coefficient times `X^(degree_gap) * p` without coefficient
division. This supports the same exact polynomial rings already handled by the
local monicity check and avoids string manipulation.

## Certificate Replay

The `:murthy_q_degree_normalization` certificate branch stores the
normalization record in `witness.normalization`. Verification checks:

- the witness key set is exact;
- the normalization record verifies independently;
- the record target equals the certificate target;
- the selected variable matches;
- the expected replay factors are exactly
  `[normalized_target, elementary_correction]`;
- the factor product equals the original target.

The certificate is a replay certificate for this reduction step, not a full
elementary factorization of the normalized target.

## Tests

Add `test/expert/sl3_local_q_degree_normalization.jl` and register it in the
expert group. The test covers:

- the issue 69 `mg-q-degree-normalization` fixture;
- a second exact polynomial example over `QQ[X]` with `deg(q) >= deg(p)`;
- exact quotient/remainder, degree bound, normalized target, and correction
  replay;
- certificate replay through `verify_sl3_local_realization`;
- matrix-input construction;
- rejection of non-monic `p`;
- replay failure for a deliberately corrupted quotient/remainder pair even
  when the original target still has determinant one.

## Verification

Focused command:

```bash
julia --project=. -e 'include("test/expert/sl3_local_q_degree_normalization.jl")'
```

Agent Desk package command:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Full expert-inclusive command:

```bash
julia --project=. test/runtests.jl all
```

## Spec Self-Review

- No incomplete placeholders remain.
- The helper is internal and polynomial-only.
- The design does not modify `elementary_factorization` or solve recursive
  Murthy branches.
- Every recorded field participates in replay verification.
