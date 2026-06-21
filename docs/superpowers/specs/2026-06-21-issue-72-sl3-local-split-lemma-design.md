# Issue 72 SL3 Local Split Lemma Design

## Context

Issue 72 adds the Murthy-Gupta split-lemma replay helper for local `SL_3`
special-form matrices

```text
[a*a' b 0; c d 0; 0 0 1]
```

where `a*a'*d - b*c == 1`. Issue 69 supplies small fixture witnesses, and
issue 70 supplies the replayable local realization certificate machinery. This
issue should not add recursive Murthy solving or route `elementary_factorization`
through the Murthy branch.

## Design Choice

Add a focused internal replay object in `src/algorithm/sl3_local.jl`:

- `SL3LocalSplitLemmaReplay`
- `sl3_local_split_lemma_replay(a, a_prime, b, c, d, c1, d1, c2, d2; split_id=:murthy_split_lemma)`
- `verify_sl3_local_split_lemma_replay(replay)::Bool`
- `sl3_local_split_lemma_certificate(replay, first_child_certificate, second_child_certificate, X)`

The helper validates the parent determinant relation and the two child
relations before returning metadata. It records the original target, the two
child targets, prefix/middle/suffix elementary wrapper factors from the
Park-Woodburn split lemma, the flattened wrapper factor list, the exact
reassembled product, and the witness data used to recompute everything.

Alternatives considered:

- Return a raw tuple. This is terse but fragile for later recursive branches,
  which need named fields and a stable split identifier.
- Integrate directly into `realize_sl3_local`. That would imply recursive
  solver behavior, which is out of scope.
- Export a public API. The helper is expert/internal until Murthy recursion is
  implemented and compatibility requirements are clearer.

## Certificate Bridge

The split lemma itself rewrites the original target into wrapper factors around
two child targets. It does not produce an elementary factorization until the
children have certificates. The bridge function accepts a split replay plus two
already verified `SL3LocalRealizationCertificate` children, substitutes their
factor sequences into the prefix/child1/middle/child2/suffix expression, and
returns an `SL3LocalRealizationCertificate` with branch
`:murthy_split_lemma`.

`verify_sl3_local_realization` will learn this one new internal branch. The
branch verifier checks:

- the split replay verifies,
- the split original target equals the certificate target,
- both child certificates verify through the existing certificate path,
- child certificate targets match the split child targets,
- selected variables are in the same ring,
- the stored factor list exactly matches the recomputed wrapper/child factor
  sequence,
- the factor product equals the original target.

## Error Handling

Construction throws `ArgumentError` when any determinant relation fails, when
inputs cannot be coerced into one parent ring, or when exact reassembly fails.
The verifier catches malformed replay objects and returns `false`, rethrowing
only `InterruptException`. The certificate bridge throws if child certificates
do not verify or do not match the child targets.

## Tests

Add `test/expert/sl3_local_split_lemma.jl` and register it in the expert group.
The focused test will cover:

- the issue 69 fixture id `mg-split-lemma-x-square`,
- two hand-checkable examples whose children are currently certifiable,
- exact original and child determinant relations,
- exact reassembly from prefix/middle/suffix wrappers and child targets,
- every wrapper factor is an elementary `3 x 3` matrix over the same ring,
- certificate bridge replay through `verify_sl3_local_realization`,
- negative controls corrupting split child witness data and mismatching child
  certificates.

## Verification

Focused command:

```bash
julia --project=. -e 'include("test/expert/sl3_local_split_lemma.jl")'
```

Package command required by the Agent Desk run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```
