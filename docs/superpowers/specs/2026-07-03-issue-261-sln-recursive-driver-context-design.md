# Issue 261 SLn Recursive Driver Context Design

Issue #261 adds an internal checked input context for the staged
Park-Woodburn ordinary-polynomial `SL_n` recursive driver in parent #186. The
context records enough recomputable metadata to decide whether a matrix has the
ring, determinant, last-column, ECP, final `SL_3`, and route-provenance evidence
needed before any elementary factors are produced. It does not route public
`elementary_factorization` through this context.

## Context

There is no repository `AGENTS.md`; the README test instructions apply. The
branch starts from `main` after merged #260 fixture-catalog work. The relevant
local predecessors are:

- `src/algorithm/factorization.jl`, which already defines checked
  `SL3RealizationInputContext` construction and verification.
- `src/algorithm/polynomial_column_peel.jl`, which validates ordinary
  polynomial column-peel inputs and records replayed peel steps.
- `test/fixtures/park_woodburn_sln_driver_cases.jl` and
  `test/internal/park_woodburn_sln_driver_fixtures.jl`, which provide the #260
  supported, staged, legacy, and negative-control catalog data.

GitHub issue fetching is unavailable in this Agent Desk sandbox for GraphQL
issue calls because the configured proxy is blocked. The issue body supplied by
the run, the merged #260 PR metadata retrievable through `gh pr view 276`, local
merge history, and repository tests are the source of truth for this design.

## Chosen Approach

Add a private context type in `src/algorithm/polynomial_column_peel.jl`, near
the recursive column-peel code it will eventually front. This keeps the #186
staged boundary close to the recursive driver without changing public route
selection.

The new internal API will be:

- `SLnRecursiveDriverInputContext`
- `_sln_recursive_driver_input_context(A; variable_order, selected_variable,
  ecp_witness_metadata, final_route_metadata, route_provenance_metadata,
  catalog_id)`
- `_sln_recursive_driver_input_context_core_verification(context)`
- `_sln_recursive_driver_input_context_verification(context)`
- `_verify_sln_recursive_driver_input_context(context)`

The constructor recomputes fields from the matrix and hints, stores a replay
summary, and returns only a context whose verifier passes. The verifier
recomputes from stored inputs instead of trusting stored statuses or diagnostic
flags.

## Context Fields

The context records:

- The original matrix, base ring, coefficient ring, dimension, ring profile,
  determinant value, determinant status, and exact field-backed status.
- Generator metadata: generators, generator names, normalized variable order,
  variable-order status, selected variable, selected-variable index, and
  selected-variable status.
- Initial last-column data: entries, length, bottom entry, whether it is already
  the target column, and a canonical unimodularity witness when available.
- Evidence metadata supplied by the caller: ECP witness metadata for the first
  last-column peel, final `SL_3` route metadata, route-provenance metadata, and
  an optional #260 catalog id.
- Evidence statuses: `:replayed`, `:recorded`, or `:missing` for ECP and final
  route metadata; `:recorded`, `:catalog`, or `:missing` for route provenance.
- `support_classification`, `staged_reason_code`, and `staged_diagnostic`.

Reason-code priority is deterministic:

1. `:unsupported_coefficient_ring`
2. `:determinant_not_one`
3. `:missing_variable_metadata`
4. `:missing_ecp_evidence`
5. `:missing_final_sl3_route`
6. no reason code when all required evidence is replayed

The diagnostic is a `NamedTuple` with `status`, `reason_code`,
`missing_evidence`, `partial_evidence`, `determinant_status`,
`exact_field_status`, `variable_metadata_status`, route provenance fields, and a
human-readable `message`.

## Evidence Rules

The context stays conservative:

- Ordinary polynomial rings are required. Non-polynomial rings throw
  `ArgumentError`.
- Exact field-backed ordinary polynomial rings are required for support. A
  determinant-one `ZZ[x]` input may construct a staged context with
  `:unsupported_coefficient_ring`.
- Determinant is checked through `_require_polynomial_sl_determinant` and caught
  into `determinant_status = :not_one` for diagnostics.
- ECP evidence is `:replayed` only when metadata includes a verified
  `ECPColumnReductionCertificate` for the current last column over the base ring
  and the target column is `e_n`. Shell metadata without a verified certificate
  is `:recorded`.
- Final route evidence is `:replayed` only when metadata declares `status =
  :replayed`, has a target matrix over the same ring, and includes a route case
  id or replay payload. Shell metadata is `:recorded`; `status = :missing` is
  `:missing`.
- Route provenance may be recorded from #260 metadata or synthesized from
  `catalog_id`, but it never makes an unsupported or staged input supported by
  itself.

## Tests

Add `test/expert/park_woodburn_sln_driver_context.jl` and register it in the
expert test group. The test will build contexts for:

- A #260 supported multivariate `SL_4` case with replayed ECP and final route
  evidence.
- A #260 supported multistep `SL_5` case to verify the initial context records
  the first recursive boundary without peeling further.
- The legacy recursive column-peel regression case, which remains staged.
- A determinant-one staged case with missing final `SL_3` evidence.
- An unsupported coefficient-ring case over `ZZ[x]`.

Negative controls mutate stored determinant status, ring profile, generator
metadata, last column, route provenance, staged reason code, staged diagnostic,
and stored verification. Each mutation must make verification return `false` or
construction throw a deliberate `ArgumentError`.

## Non-Claims

This issue does not reduce the last column, produce elementary factors, call
#184 `SL_3` routes, change public route selection, broaden Laurent or
ToricBuilder support, or accept arbitrary recursive `SL_n` inputs. It only
builds and verifies the internal staged context that later #186 driver issues
can consume.

## Self-Review

The design is scoped to one internal context and one expert test. It reuses the
existing checked-context pattern from the `SL_3` driver, keeps route selection
unchanged, and uses #260 catalog metadata only as replayable input hints. The
diagnostic reason-code priority is explicit, stable, and covers every issue
reason code requested for this stage.
