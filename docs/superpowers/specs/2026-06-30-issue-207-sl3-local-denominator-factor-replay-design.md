# Issue 207 SL3 Local Denominator Factor Replay Design

## Context

Issue #207 adds the replay format needed by later Park-Woodburn local SL3
branches. Existing `SL3LocalRealizationCertificate` values store ordinary
`3 x 3` elementary matrices, which is correct for the current `QQ[X]`
materializable solver. Issue #208 added the checked Murthy input context and
local-unit witness verifier, including the localization-at-maximal-ideal schema
from #206. This issue should consume that schema for factor replay only. It
must not route q(0)-unit or q(0)-nonunit local branches through the new record.

## Approaches Considered

1. Add an internal `SL3LocalElementaryFactor` record plus replay summary and
   keep public factorization APIs unchanged. This is the chosen approach. It
   gives later local branches a denominator-aware representation, reuses the
   #208 local-unit verifier, and keeps current solver outputs as ordinary
   matrices.

2. Replace certificate `factors` with a mixed vector of matrices and local
   records. This was rejected because verifiers would have to guess whether
   they are checking ordinary products or localized products, which the issue
   explicitly warns against.

3. Change public `elementary_matrix` or `verify_factorization` to accept local
   denominators. This was rejected as out of scope and would alter behavior for
   ordinary APIs that currently multiply exact matrices over the base ring.

## Design

Add two internal records in `src/algorithm/sl3_local.jl`:

- `SL3LocalElementaryFactor`: records `R`, `n`, `row`, `col`, `numerator`,
  `denominator`, `selected_variable`, and an optional `local_unit_witness`.
- `SL3LocalElementaryFactorReplay`: records a target matrix, a vector of local
  factors, replay mode, denominator product, denominator-cleared product, and
  the ordinary materialized factors when every denominator is one.

The constructor validates `n == 3`, `row != col`, index bounds, parent rings,
and denominator evidence. A denominator of one needs no local-unit witness and
can materialize as an ordinary elementary matrix. A non-one denominator must
carry a #208-compatible local-unit witness, verified with
`_sl3_local_murthy_verify_local_unit_witness(R, X, witness, denominator)`.

Replay uses an explicit denominator-cleared identity. For a factor
`I + numerator / denominator * E_ij`, define its cleared matrix as
`denominator * I + numerator * E_ij`. For a sequence with denominator product
`D`, verification checks:

```text
cleared_factor_1 * ... * cleared_factor_m == D * target
```

over the ordinary base ring. If all denominators are one, the replay mode is
`:ordinary`, the materialized factors are ordinary matrices, and the product is
checked with the existing exact `verify_factorization` semantics. Otherwise the
replay mode is `:denominator_cleared` and no ordinary materialized factor list
is exposed.

## Certificate Integration

Current `realize_sl3_local(...)` and `realize_sl3_local_certificate(...)`
continue to return ordinary matrix factors for supported materializable cases.
Add explicit adapter helpers so existing certificate factors can be represented
as denominator-one local records:

- `sl3_local_elementary_factor(row, col, numerator, denominator, X;
  local_unit_witness = nothing, n = 3)`
- `sl3_local_materialize_elementary_factor(record)`
- `sl3_local_elementary_factor_replay(target, records, X)`
- `sl3_local_denominator_one_records_from_matrices(factors, X)`
- `verify_sl3_local_elementary_factor(record)`
- `verify_sl3_local_elementary_factor_replay(replay)`

These helpers remain unexported. They avoid storing an ambiguous mixture of
raw matrices and local records while letting tests and future branches compare
ordinary certificate factors against denominator-one records through an
explicit adapter layer.

## Error Handling

Constructors throw `ArgumentError` for malformed row/column indexes, parent
ring mismatches, non-3 sizes, invalid denominator witnesses, and matrices that
cannot be unambiguously adapted as single elementary matrices. Verifier
functions catch ordinary exceptions and return `false`, matching the existing
SL3 local verifier style. `InterruptException` is rethrown.

## Tests

Add `test/expert/sl3_local_local_factors.jl` and register it in the expert
group. The focused test covers:

- denominator-one records materialize to the same matrices as existing
  `elementary_matrix` calls;
- adapting current certificate factors into denominator-one records produces an
  ordinary replay whose materialized product equals the certificate target;
- a nontrivial local denominator such as `1 + u` with the #206/#208 local-unit
  witness verifies by the denominator-cleared equality;
- corrupting numerator, denominator, row/column index, or unit witness makes
  construction throw `ArgumentError` or verification return `false`.

Update `test/expert/sl3_local_certificate.jl` with a small assertion that an
ordinary certificate can be replayed through denominator-one local records
without changing `realize_sl3_local` return values.

## Out Of Scope

- Do not implement q(0)-unit or q(0)-nonunit local Murthy factor production.
- Do not change exported public APIs.
- Do not make nontrivial local factors materialize as ordinary matrices.
- Do not alter ordinary `verify_factorization` behavior.

## Verification

Focused commands:

```bash
julia --project=. -e 'include("test/expert/sl3_local_local_factors.jl")'
julia --project=. -e 'include("test/expert/sl3_local_certificate.jl")'
```

Package command required by Agent Desk:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Spec Self-Review

- No placeholders or open questions remain.
- The design is limited to an internal replay format and verifier.
- Public ordinary factorization APIs keep their current return values.
- The #208 local-unit schema is reused instead of redefining local evidence.
