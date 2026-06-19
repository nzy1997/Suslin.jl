# Issue 15 SL_n to Local SL_3 Reduction Design

## Context

`elementary_factorization(A)` currently accepts the narrow local `SL_3`
polynomial slice and explicitly rejects all larger `SL_n` matrices. The
dependency issues have landed the pieces needed for a first exact reduction
layer:

- block embeddings and factor-sequence composition,
- exact unimodular-column reduction for supported columns,
- a broader local `SL_3` solver with staged failures,
- Quillen patch metadata patterns,
- Laurent determinant normalization and exact reconstruction checks,
- the documented full-suite command `julia --project=. test/runtests.jl all`.

Issue #15 asks for the missing connection from supported larger `SL_n` inputs
to inspectable local `SL_3` obligations. It does not require a general Suslin
factorization algorithm in this step.

## Approaches Considered

1. Add a narrow block-local reduction layer. This recognizes matrices whose
   exact non-identity support decomposes into embedded 3x3 local blocks, solves
   each block through `realize_sl3_local`, embeds the returned factors, records
   obligation metadata, and verifies the assembled product. This is the chosen
   approach because it reuses the dependency helpers and keeps unsupported
   matrices staged.
2. Implement a broad constructive `SL_n` decomposition now. That would exceed
   the available local-solver support and blur issue #17's later factorization
   work.
3. Leave the behavior internal to `elementary_factorization` without an
   inspectable reduction object. That would make the requested local
   obligations and reassembly metadata hard to inspect and tamper-test.

## Selected Design

Add a small public reduction record layer:

- `SL3LocalObligation` records the block location, ring, target local matrix,
  required assumptions, embedded target, local factors, embedded factors, and
  exact reassembly data for one local block.
- `SLNToSL3Reduction` records the original input, any Laurent normalization
  metadata, the determinant-one core, obligations, assembled factors, product,
  and verification metadata.
- `reduce_sln_to_sl3(A)` validates and reduces a supported `SL_n` matrix.
- `verify_sln_to_sl3_reduction(reduction)` recomputes the exact metadata and
  verifies that the recorded factors multiply to the recorded target.

The supported polynomial family is intentionally exact and inspectable:

1. `A` must be square of size at least 3 over a supported ordinary polynomial
   ring or Laurent polynomial ring.
2. Ordinary polynomial inputs must already have determinant one. Laurent inputs
   pass through `normalize_laurent_gl_matrix(A)` first; only the normalized
   `SL_n` core is reduced by this issue.
3. The determinant-one core must decompose into disjoint 3x3 coordinate blocks
   plus identity coordinates. Every off-block entry must be exactly zero and
   every identity coordinate must have the expected identity row and column.
4. Each nontrivial block must be recognized by the local `SL_3` solver. For
   ordinary polynomial rings, the first generator is used as the local
   variable. For Laurent rings, the first generator is used with
   `check_monic=false`, because the caller has reached the normalized local
   layer and monicity is not currently meaningful there.
5. The reduction embeds each local factor sequence with `embed_factor_sequence`
   and concatenates the embedded sequences in block order.
6. Exact product verification is mandatory before returning. A valid but
   unsupported matrix fails with an `ArgumentError` containing
   `staged SL_n to local SL_3 reduction failure`, not with a misleading local
   solver message.

`elementary_factorization(A)` will call `reduce_sln_to_sl3(A)` when the current
direct local `SL_3` path does not apply. For normalized Laurent inputs with a
non-identity determinant correction, the reduction record remains available for
the normalized core, while `elementary_factorization` still returns only
factors for inputs whose exact assembled factors multiply to the original
matrix. If a Laurent correction would require non-elementary determinant
metadata in the returned factor sequence, the driver keeps a staged error
instead of claiming a complete elementary factorization.

## Reassembly Metadata

Each obligation records:

- `block_location`: a vector of global coordinate indices.
- `ring`: the base ring.
- `target_local_matrix`: the 3x3 local target sent to `realize_sl3_local`.
- `required_assumptions`: symbols such as `:disjoint_block_support`,
  `:determinant_one_core`, `:ordinary_polynomial_monicity` or
  `:laurent_normalized_check_monic_false`.
- `embedded_target`: the 3x3 target embedded into the full matrix.
- `local_factors` and `embedded_factors`.
- `reassembly_data`: a named tuple with local product, embedded product,
  expected embedded target, and exact product checks.

The top-level reduction verification records the assembled product and whether
the factor product equals the determinant-one core and, when applicable, the
original input.

## Tests

Add `test/expert/sln_to_sl3_reduction.jl` and register it in the expert group.
The focused tests cover:

- a 6x6 supported ordinary polynomial example with two local obligations,
- an 8x8 supported ordinary polynomial example with two local obligations and
  two identity coordinates,
- exact solving and reassembly through `verify_factorization`,
- `reduce_sln_to_sl3` metadata verification,
- dropping one recorded local obligation makes exact final verification fail,
- a determinant-one matrix outside the block-local support fails with the
  staged `SL_n` reduction error.

Also extend the API surface test for the new public records and functions.

## Verification

Run the issue-specific command:

```bash
julia --project=. -e 'include("test/expert/sln_to_sl3_reduction.jl")'
```

Run the Agent Desk package command:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Run the documented full suite from #21:

```bash
julia --project=. test/runtests.jl all
```

## Automatic Approval

This Agent Desk run is non-interactive. Under the standing answer policy, the
design and written spec are approved automatically because the selected
approach is the narrowest exact layer that satisfies issue #15, reuses the
landed dependency helpers, and avoids claiming a general `SL_n` algorithm.

## Spec Self-Review

- No incomplete markers remain.
- The supported input family is explicit and staged.
- Local obligations and reassembly data are inspectable.
- Negative controls are defined for omitted obligations and unsupported inputs.
- Laurent normalization is acknowledged without claiming unsupported determinant
  corrections as elementary factorizations of the original input.
