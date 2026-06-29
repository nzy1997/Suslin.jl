# Issue 208 Murthy Local Input Context Design

## Context

Issue #208 adds the checked input context needed before later Park-Woodburn
Section 5 local branches can consume denominator-aware local data. The current
`src/algorithm/sl3_local.jl` code recognizes special-form `SL_3` targets and
already records replayable q-degree, split-lemma, q(0)-unit, and q(0)-nonunit
reductions for ordinary `QQ[X]` inputs. Issue #206 added a local-contract
fixture catalog over `QQ[u, X]`, including local-unit witnesses for elements
such as `1 + u` that are not global polynomial-ring units.

This issue should add a replayable, internal context layer. It must validate
the shape and Murthy preconditions for both current ordinary cases and broader
local-contract cases, but it must not produce new elementary factors or claim
that multivariate local realization is implemented.

## Approaches Considered

1. Add a non-exported `SL3LocalMurthyInputContext` record in
   `src/algorithm/sl3_local.jl`, with constructors from either a matrix or
   `(p, q, r, s, X)` and a verifier that recomputes every stored field. This is
   the chosen approach because it lives beside the existing local solver records,
   reuses the #206 witness schema, and gives later branches one stable internal
   object to consume.

2. Extend the #206 fixture validator only. This would validate examples, but it
   would leave solver code without an internal context schema and would duplicate
   the same checks when later branches are added.

3. Change `realize_sl3_local` to accept local coefficients immediately. This
   would exceed the issue scope by starting factor production for local branches
   before q(0)-unit, q(0)-nonunit, denominator covers, and patching are ready.

## Design

Add `SL3LocalMurthyInputContext` in `src/algorithm/sl3_local.jl`. The record
stores:

- `R`, `X`, selected variable index, `(p, q, r, s)`, target matrix, and
  determinant;
- `degree_p`, `degree_q`, selected-variable constants `p0` and `q0`;
- `p_monic`, global-unit booleans for `p`, `q`, `r`, `s`, `p0`, and `q0`;
- local-unit booleans for `p0`, `q0`, `resultant`, and `branch_unit`;
- normalized local-unit witnesses, optional split witness, and optional Bezout
  witness data supplied by callers.

Constructors:

```julia
Suslin.sl3_local_murthy_input_context(A, X; witness = nothing,
    local_unit_witnesses = (;), split_witness = nothing, bezout_witness = nothing)

Suslin.sl3_local_murthy_input_context(p, q, r, s, X; witness = nothing,
    local_unit_witnesses = (;), split_witness = nothing, bezout_witness = nothing)
```

The functions are not exported. The broad `witness` keyword accepts the #206
fixture witness shape and extracts `local_unit_witness` as `q0`,
`branch_unit_witness` as `branch_unit`, `resultant_unit_witness` as
`resultant`, `split` as split data, and `p_prime`/`q_prime` as Bezout data.
Explicit keyword arguments remain available for later internal consumers.

Local-unit witnesses use the #206 schema exactly:

- `context.kind = :localization_at_maximal_ideal`;
- `context.selected_variable == X`;
- `context.maximal_ideal_generators == maximal_ideal_generators`;
- `unit`, `residue_unit`, `residue_inverse`,
  `maximal_ideal_generators`, and `residue_difference_coefficients`;
- exact equation `unit - residue_unit == sum(coeff_i * generator_i)`;
- exact equation `residue_unit * residue_inverse == 1`;
- optional `global_unit` must agree with `is_unit(unit)`.

The constructor rejects malformed target shape, parent-ring mismatches, a
non-generator selected variable, determinant not one, non-monic `p`, invalid
local-unit witnesses, invalid split witnesses, invalid Bezout witnesses, and
local q(0)-unit/q(0)-nonunit contract data that needs a local unit but does not
provide evidence. Ordinary univariate `QQ[X]` cases remain accepted without a
local witness.

## Verification

`verify_sl3_local_murthy_input_context(context)` recomputes every stored field
from the target and selected variable. It rejects tampered records by checking:

- target shape and entries;
- determinant one and stored determinant;
- selected variable parent and index;
- monicity and degree fields;
- `p0`, `q0`, global-unit classification, and local-unit classification;
- exact local-unit witness equations;
- exact split-lemma witness replay when supplied;
- exact Bezout equality, degree guards, resultant classification, and branch
  unit classification when supplied.

The verifier catches ordinary exceptions and returns `false`, following the
existing local solver verifier style.

## Testing

Add `test/expert/sl3_local_murthy_context.jl` and register it in the expert test
group. Tests cover:

- one existing ordinary `QQ[X]` Murthy case without a local witness;
- one #206 `QQ[u, X]` q(0)-unit local-contract case where `q0 = 1 + u` is not a
  global unit but is recorded as a local unit by the supplied witness;
- constructor rejection when the local-unit witness is removed;
- constructor rejection for wrong selected variable and non-monic `p`;
- verifier rejection for corrupted stored `q0`;
- constructor/verifier rejection for corrupted local-unit and Bezout witnesses.

Focused verification:

```bash
julia --project=. -e 'include("test/expert/sl3_local_murthy_context.jl")'
```

Package verification:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Out Of Scope

- No new elementary factor production.
- No q(0)-unit local denominator replay.
- No q(0)-nonunit local branch implementation.
- No Quillen denominator covers, global patching, ECP integration, or public
  Park-Woodburn driver.
- No exported public API.

## Spec Self-Review

- No placeholders or open questions remain.
- The design is scoped to one internal context record, constructor, verifier,
  and focused tests.
- The #206 local-unit witness schema is reused instead of duplicated.
- Existing solver behavior remains staged for unsupported local realization.
