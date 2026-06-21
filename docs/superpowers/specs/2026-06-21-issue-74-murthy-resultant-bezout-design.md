# Murthy q(0)-Nonunit Resultant/Bezout Branch Design

## Context

Issue #74 completes the second normalized Murthy-Gupta local `SL_3` branch for
ordinary univariate polynomial special forms. The current solver already
supports open slices, unit pivots, replayable certificates, q-degree
normalization, split-lemma replay, and the q(0)-unit recursive branch. The
remaining gap is a normalized monic target whose `q(0)` is not a unit.

## Approaches Considered

1. Add an internal q(0)-nonunit reduction in `src/algorithm/sl3_local.jl`.
   This is the chosen approach because it reuses the existing certificate
   verifier, keeps branch metadata replayable, and routes the transformed
   target through the q(0)-unit solver instead of duplicating recursion.

2. Add a separate expert-only entry point that accepts Bezout data.
   This would keep the public solver dispatch unchanged, but it would not
   complete the requested local branch and would split certificate replay
   across two APIs.

3. Keep the branch entirely fixture-supplied.
   This would cover the checked examples but would make supported `QQ[X]`
   inputs depend on decorative metadata. Exact univariate `gcdx` can extract
   a Bezout relation for the fixture family, so extraction should be used when
   no explicit witness is supplied.

## Design

The existing `realize_sl3_local(...)` factor-only API remains intact. The
certificate-producing methods gain an optional `murthy_q0_nonunit_witness`
keyword that may supply `p_prime` and `q_prime` with
`p_prime*p - q_prime*q == 1`. When no witness is supplied for an ordinary
univariate polynomial ring, the branch tries exact `gcdx(p, q)` extraction and
normalizes a unit gcd to the exact Bezout equality. Unsupported extraction or
invalid witnesses raise staged `ArgumentError`s before factors are returned.

For a normalized q(0)-nonunit target `A = [p q; r s]`, the branch verifies
`p` is monic, `deg(q) < deg(p)`, `q(0)` is not a unit, the determinant is one,
and the Bezout equality holds. It then forms
`B = [p q; q_prime p_prime]` and
`C = [p + q_prime, q + p_prime; q_prime, p_prime]`. The exact identities are

```text
A = E_21(r*p_prime - s*q_prime) * B
B = E_12(-1) * C
```

The transformed child `C` has determinant one and top-right constant term
`q(0) + p_prime(0)`, which the branch requires to be a unit. The child is then
realized by `realize_sl3_local_certificate(C, X)`, so the existing q(0)-unit
path owns all recursive work.

## Certificate Replay

A new q(0)-nonunit reduction record stores the original target, `p(0)`,
`q(0)`, `p_prime`, `q_prime`, the exact resultant value, the Bezout matrix, the
q(0)-unit child target, both elementary transformation factors, the child
certificate, the selected variable, degree guards, and whether the witness was
supplied or extracted. Replay verifies the Bezout relation, determinant-one
child target, unit child `q(0)`, exact elementary transformation identities,
child certificate replay, and final factor product.

The `:murthy_q0_nonunit_bezout_resultant` certificate branch returns factors
`[E_21(...), E_12(-1), child factors...]`. The generic
`verify_sl3_local_realization` path reconstructs the expected sequence from
the reduction record and rejects tampered witnesses or factors.

## Scope

Supported extraction is limited to ordinary univariate exact polynomial rings
where `gcdx(p, q)` produces a unit gcd and the normalized Bezout relation
passes all branch checks. Supplied witnesses are accepted for the same
ordinary univariate local branch after exact verification. Multivariate
local-coefficient extraction, Quillen patching, ECP reduction, the final public
Park-Woodburn driver, and factor-count optimization stay out of scope.

## Testing

Add `test/expert/sl3_local_murthy_resultant.jl` and register it in the expert
group. Tests cover:

- the #69 q(0)-nonunit Bezout fixture with an explicitly supplied witness;
- a second normalized q(0)-nonunit example using extracted `gcdx` witness data;
- exact verification that each example transforms into a q(0)-unit child;
- replay metadata for the Bezout equality and elementary transformations;
- exact final factor products and `verify_factorization(target, factors)`;
- negative controls for a corrupt supplied witness and malformed replay data.
