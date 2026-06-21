# Murthy q(0)-Unit Recursive Branch Design

## Context

Issue #73 adds the first recursive Murthy-Gupta path to the local `SL_3` solver. The current solver already has fast paths for open slices and unit diagonal pivots, plus replay helpers for q-degree normalization and the split lemma. The new branch should use those helpers without changing the public `realize_sl3_local(...)` factor-only API.

## Approaches Considered

1. Internal recursive branch in `src/algorithm/sl3_local.jl` using existing certificate machinery.
   This keeps the implementation near the current local solver, reuses q-degree normalization and split-lemma replay, and lets expert tests inspect nested certificates. This is the chosen approach.

2. Separate expert-only Murthy entry point.
   This would avoid touching the existing solver dispatch, but it would not cover the issue's requested `realize_sl3_local` path and would duplicate recognition logic.

3. Fixture-driven solver that consumes checked-in witness tuples.
   This would satisfy the catalog examples but would make the solver decorative outside fixtures. It is too narrow unless a future local-ring witness API is requested.

## Design

The branch remains internal. `realize_sl3_local_certificate` first preserves the existing open-slice and unit diagonal pivot fast paths. It then accepts ordinary univariate polynomial inputs where `p` is monic and the constant coefficient of normalized `q` is a unit. If `deg(q) >= deg(p)`, it records a q-degree normalization step, recursively realizes the normalized target, and appends the existing elementary correction.

For a normalized q(0)-unit target, the branch computes `lambda = -p(0) * inv(q(0))`, right-multiplies by `E_21(lambda)` to eliminate the constant term of `p`, verifies `p + lambda*q = X*p_prime`, and requires `degree(p_prime) < degree(p)`. It then builds split-lemma witnesses by exact univariate `gcdx` for `(X, q)` and `(p_prime, q)`, constructs the split replay, realizes both child certificates, and appends `E_21(-lambda)` to reconstruct the original target.

The first split child may have first entry `X`; if its `q` still has positive degree, the same q-degree normalization helper reduces it to a constant unit. A new internal q-unit base case handles matrices with unit top-right entry by the identity
`E_21((s - 1)/q) * E_12(q) * E_21((p - 1)/q)`.

## Certificate Replay

A new Murthy q(0)-unit reduction record stores the original target, `q(0)`, its inverse, `p(0)`, the elimination coefficient and inverse correction, the eliminated target, `p_prime`, the split certificate, the selected variable, and the degree guard. The `:murthy_q0_unit` certificate witness stores either a q-degree normalization child certificate or this reduction record. Replay verifies all stored equations and exact factor products.

## Scope

Supported inputs are ordinary univariate exact polynomial rings such as `QQ[X]`. Multivariate local coefficient witness extraction is intentionally staged because this issue does not define a public witness API. The q(0)-nonunit Bezout/resultant branch remains unsupported and must continue to throw a staged local `SL_3` failure.

## Testing

Add `test/expert/sl3_local_murthy_q_unit.jl` and register it in the expert group. Tests cover:

- the #69 q(0)-unit fixture;
- a degree-greater-than-one normalized q(0)-unit example;
- a degree-greater-than-one example that first uses q-degree normalization;
- nested certificate replay of normalization, constant-term elimination, split children, and recursive child certificates;
- exact final elementary factor products;
- a normalized nonunit `q(0)` negative control that remains staged.
