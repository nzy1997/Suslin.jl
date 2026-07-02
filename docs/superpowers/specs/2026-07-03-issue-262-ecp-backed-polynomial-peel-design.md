# Issue 262 ECP-Backed Polynomial Peel-Step Certificates Design

Issue #262 hardens ordinary-polynomial column-peel replay for parent #186. Each
real peel step must record the ECP evidence used to reduce the last column, the
left factors extracted from that evidence, deterministic right-clearing data,
the peeled block, determinant/descent metadata, and a verifier that recomputes
the step instead of trusting stored products.

## Context

There is no repository `AGENTS.md`; the README test instructions apply. Live
`gh issue view 262` is unavailable in this Agent Desk sandbox because the local
GitHub proxy connection is blocked, so the run-provided issue body and current
mainline code are the source of truth.

The relevant local predecessors are:

- `src/algorithm/column_reduction.jl`, which exposes
  `ecp_column_reduction_certificate`, `reduce_unimodular_column`, and the stable
  `verify_ecp_column_reduction` verifier from #185/#248.
- `src/algorithm/polynomial_column_peel.jl`, which already stores
  `left_certificate` on `PolynomialColumnPeelStep` but does not bind explicit
  peel-step ECP provenance, clearing coefficients, or determinant/descent
  metadata.
- `test/fixtures/park_woodburn_sln_driver_cases.jl`, which provides #260
  multivariate `SL_n` fixtures whose last columns reduce through replayed ECP
  certificates.
- `test/expert/park_woodburn_polynomial_column_peel.jl`, which covers the
  existing factor replay and should gain stronger negative controls.

## Chosen Approach

Extend `PolynomialColumnPeelStep` in place. This preserves the current
constructor that accepts only factor data while making the mainline constructor
record a richer certificate boundary for new #186 consumers.

The step will add these internal fields:

- `ecp_evidence`: the stored `ECPColumnReductionCertificate` or successor
  evidence object.
- `ecp_route_provenance`: a normalized `NamedTuple` derived from the evidence
  and stable verifier contract, not from fragile stage field assumptions.
- `right_clearing_coefficients`: the bottom-row coefficients of `B*A` that
  deterministically produce `E_{d,i}(-p_i)`.
- `block_embedding_indices`: the indices used to embed the next block back into
  the peeled matrix.
- `determinant_metadata`: input, peeled, next-block, and expected determinant
  values.
- `descent_metadata`: input and next dimensions, step index metadata, and the
  recorded dimension drop.
- `verification`: a per-step replay summary.

The existing `left_certificate` field remains as a compatibility alias for the
same ECP evidence. The legacy positional constructor will populate the new
fields from the supplied factor-only data when possible, but it will not invent
verified ECP evidence. Mainline `_polynomial_column_peel_step` will always
record a verified `ECPColumnReductionCertificate`.

## Verification Rules

The step verifier will recompute all evidence-sensitive data:

- The recorded last column must equal the input matrix last column.
- `ecp_evidence` must be an `ECPColumnReductionCertificate` accepted by
  `verify_ecp_column_reduction`.
- The ECP evidence must be over the same base ring, start from the recorded last
  column, end at `e_d`, and expose factors equal to the stored `left_factors`.
- `ecp_route_provenance` must match a normalized adapter result derived from the
  verified certificate. The adapter may read public/stable certificate fields
  and broad route metadata when present, but the verifier remains anchored on
  `verify_ecp_column_reduction`.
- `after_left_matrix` must be the product of the left factors with the input
  matrix, and its last column must be `e_d`.
- `right_clearing_coefficients` must equal the bottom row of `after_left_matrix`
  in columns `1:(d - 1)`.
- `right_factors` must equal the deterministic factors recomputed from those
  coefficients.
- `peeled_matrix` must equal `after_left_matrix * product(right_factors)`, must
  be the block embedding of `next_block`, and must isolate the final row/column.
- `next_block` must have determinant one.
- `determinant_metadata`, `descent_metadata`, `block_embedding_indices`, and
  stored per-step verification must match recomputation.

The certificate verifier will require all steps to pass this richer step
verifier and will keep the existing whole-certificate factor replay and final
route checks.

## Compatibility

Factor-only step construction remains available for tests and older internal
callers. Those steps may still satisfy exact factor replay, but they will fail
the mainline certificate verifier because no independent ECP evidence is
present. This keeps legacy behavior constructible without allowing #186 replay
to infer ECP support from multiplication alone.

No changes are needed in `src/algorithm/column_reduction.jl` unless a small
adapter proves necessary. The current stable verifier is sufficient.

## Tests

Add `test/expert/park_woodburn_sln_peel_step.jl` and register it in
`test/runtests.jl`. The new expert test will build the #260
`sln-driver-sl4-gf2-ecp-mainline` matrix, call
`_polynomial_column_peel_step`, and assert:

- the ECP certificate verifies independently through `verify_ecp_column_reduction`;
- route provenance identifies the stable ECP verifier and the public route when
  present;
- left factors send the last column to `e_d`;
- right-clearing factors are recomputed from the bottom row of `B*A`;
- the peeled matrix equals the upper-left block embedding;
- the next block has determinant one;
- determinant and descent metadata match recomputation.

Extend `test/expert/park_woodburn_polynomial_column_peel.jl` to use the new
constructor shape and add negative controls for corrupt ECP evidence, route
provenance, right-clearing coefficient, recorded last column, one left factor,
one right factor, and next block. Each negative control must fail even when the
stored final product is restored.

Verification commands:

```bash
julia --project=. -e 'include("test/expert/park_woodburn_sln_peel_step.jl")'
julia --project=. -e 'include("test/expert/park_woodburn_polynomial_column_peel.jl")'
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Non-Claims

This does not assemble multiple recursive steps into a public #186 driver,
implement `SL_3` route selection, broaden Laurent/ToricBuilder support, or
optimize factor counts. It only records and verifies the ECP-backed evidence for
ordinary-polynomial peel steps.

## Self-Review

The design is scoped to one internal step certificate boundary, one existing
certificate verifier, and two expert test files. The verifier rejects stored
products and temporary metadata unless they replay from the matrix and verified
ECP evidence.
