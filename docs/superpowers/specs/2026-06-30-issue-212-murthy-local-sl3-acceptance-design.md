# Issue 212 Murthy Local SL3 Acceptance Design

## Context

Issue #212 is the parent closeout gate for #182 after the Murthy local
children #206, #207, #208, #209, and #210. The current code already has the
checked Murthy fixture catalog, denominator-aware local factor replay, checked
input contexts, q-degree/q(0)-unit routing, and q(0)-nonunit
Bezout/resultant routing. The closeout needs acceptance tests and public
documentation that state exactly what is supported without claiming the later
Park-Woodburn mainline issues.

The important boundary is ordinary matrix factors versus localized replay.
Ordinary `realize_sl3_local`/`elementary_factorization` factor vectors are
supported when every factor materializes over the ordinary base ring. Local
Murthy witness cases are supported through `SL3LocalRealizationCertificate`
and denominator-cleared `SL3LocalElementaryFactorReplay`; they must not be
documented as ordinary factor vectors over the base ring.

## Approaches Considered

1. Extend the existing `test/expert/sl3_local_murthy_gupta.jl` acceptance file
   and docs. This is the chosen approach because the file is already the
   Murthy closeout anchor in `test/runtests.jl`, and it keeps the #212 gate
   focused on acceptance and wording.

2. Add a new separate #212 gate file. This was rejected because it would split
   the Murthy closeout assertions across two expert files while testing the
   same public/internal contracts.

3. Add new production APIs for Murthy-to-Quillen handoff. This was rejected as
   out of scope. The current tree has Quillen patch-route adapter code, but no
   Murthy-to-Quillen adapter callable. The optional #211 smoke check should be
   skipped unless that bridge exists.

## Design

Add an Issue #182 closeout testset to `test/expert/sl3_local_murthy_gupta.jl`.
It should construct non-fixture monic special-form `SL_3` inputs directly in
the test:

- an ordinary q-degree normalization case with `deg(q) >= deg(p)`;
- an ordinary q(0)-unit case whose factors materialize as ordinary matrices;
- an ordinary q(0)-nonunit case with a supplied Bezout/resultant witness;
- local q(0)-unit and local q(0)-nonunit contexts over `QQ[u, X]` using
  explicit local-unit/Bezout witnesses and denominator-cleared replay.

The ordinary cases must keep using `verify_factorization(A, factors)`. The
local cases must assert `SL3LocalElementaryFactorReplay.mode ==
:denominator_cleared`, `materialized_factors === nothing`, and verifier
acceptance through `verify_sl3_local_realization` and
`verify_sl3_local_elementary_factor_replay`.

Add negative controls to the same expert acceptance file for:

- determinant not one;
- non-monic `p`;
- missing local witness;
- unsupported automatic local Bezout extraction;
- corrupted certificate/factor sequence rejection.

Update `test/public/factorization_driver_shell.jl` with a non-fixture ordinary
Murthy case that proves the public ordinary factor path remains limited to
materializable factors.

Update `README.md` and `docs/src/index.md` so the scope says #182 Murthy local
certificates are supported for the proven ordinary/local-witness contract.
The wording must still say Quillen automatic patching (#183), general `SL_3`
(#184), ECP (#185), recursive `SL_n` (#186), full public Park-Woodburn
acceptance (#187), Laurent/ToricBuilder mainline acceptance, and
Steinberg factor-count optimization remain staged.

## Verification

Focused verification:

```bash
julia --project=. -e 'include("test/expert/sl3_local_murthy_gupta.jl")'
```

Full verification:

```bash
julia --project=. test/runtests.jl all
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Out Of Scope

- No implementation of #183 through #187.
- No Murthy-to-Quillen adapter unless a callable #211 bridge is already in
  this tree.
- No Steinberg factor-count optimization.
- No broader Laurent/ToricBuilder support.
- No public API export changes.

## Spec Self-Review

- No placeholder markers remain.
- Scope is limited to acceptance coverage and staged-boundary documentation.
- Ordinary and localized replay outputs are explicitly separated.
- The optional #211 smoke is intentionally skipped because no Murthy adapter
  callable is present in the current codebase.
