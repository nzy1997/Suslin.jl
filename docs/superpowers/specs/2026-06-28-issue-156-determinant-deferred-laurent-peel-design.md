# Issue 156 Determinant-Deferred Laurent Peel Metadata Design

## Context

Issue #156 builds on the lazy Laurent determinant fixture catalog from #154 and
the lazy peel entry point from #155. The existing Laurent column-peel
certificate proves determinant-one input up front, then records enough
elementary left/right peel data to replay the full reduction to a final `2 x 2`
`SL` block. The lazy determinant route needs a sibling internal certificate
that starts with a Laurent `GL_n` matrix, records elementary peel operations
first, and exposes the smaller matrix where determinant classification is
deferred.

PR #170 for #155 added `_factor_laurent_gl_lazy_determinant_peel(A; ...)`.
That path completes one peel before probing the smaller block, then continues
only when the deferred block is already determinant-one. Issue #156 should not
change that public behavior. It should add metadata and verification for the
deferred point so later `GL_n` determinant correction work has exact replay
evidence.

## Approach Options

Recommended: extend `src/algorithm/laurent_column_peel.jl` with an internal
determinant-deferred certificate type and verifier. The new entry point records
one or more existing `LaurentColumnPeelStep` values, computes the accumulated
left and right products that transform the original input into
`blockdiag(deferred_submatrix, I)`, sets `determinant_source =
:deferred_submatrix`, and verifies the relation exactly. This reuses the
existing step data and keeps the public API unchanged.

Alternative: add fields to `LaurentColumnPeelFactorization`. This would avoid a
new type, but it would blur two different invariants: the `SL` certificate
proves a full factorization of the original matrix, while the deferred
certificate deliberately stops before determinant normalization.

Alternative: return only a named tuple from `_factor_laurent_gl_lazy_determinant_peel`.
This is small, but it would make tamper-resistant construction and replay less
consistent with the existing certificate style.

## Chosen Design

Add internal type `LaurentDeterminantDeferredPeelCertificate` in
`src/algorithm/laurent_column_peel.jl` with these fields:

- `original_matrix`;
- `peel_steps::Vector{LaurentColumnPeelStep}`;
- `deferred_submatrix`;
- `determinant_source::Symbol`;
- `left_factors::Vector`;
- `right_factors::Vector`;
- `left_product`;
- `right_product`;
- `target_matrix`;
- `verification`.

The constructor recomputes the replay artifacts from `original_matrix`,
`peel_steps`, and `deferred_submatrix`, ignoring caller-supplied replay values.
This mirrors `LaurentColumnPeelFactorization`, whose constructor also rebuilds
the factor sequence from step metadata.

For a single peel step, the stored step already satisfies:

```text
L_1 * A * R_1 = blockdiag(B_1, 1)
```

For multiple steps, the replay embeds each smaller step into the upper-left
corner of the original dimension:

```text
L_k' * ... * L_1 * A * R_1 * ... * R_k' =
blockdiag(B_k, I_k)
```

where `B_k` is the deferred submatrix. The helper should return both the
embedded left factors in multiplication order and the embedded right factors in
multiplication order. `target_matrix` is
`block_embedding(deferred_submatrix, n, 1:size(deferred_submatrix, 1))`.

Add internal entry point `_laurent_determinant_deferred_peel_certificate(A;
min_steps = 1, progress_callback = nothing)`. It validates the Laurent square
input shape and size, performs exactly `min_steps` peel steps without calling
`classify_laurent_determinant`, and returns the certificate. The default
`min_steps = 1` is the issue's lazy peel metadata surface. The keyword remains
internal and exists only to make the replay helper naturally support repeated
steps.

Add `_verify_laurent_determinant_deferred_peel_replay(certificate)::Bool`.
The verifier recomputes:

- input shape and determinant source metadata;
- step chain consistency from the original matrix to the deferred submatrix;
- each `LaurentColumnPeelStep` using existing `_is_valid_laurent_column_peel_step_data`;
- embedded left and right factor sequences;
- left and right products;
- exact target `blockdiag(deferred_submatrix, I)`;
- exact relation `left_product * original_matrix * right_product == target_matrix`.

Tampering with any recorded left or right elementary factor must make the
verifier return `false` or throw a verification error.

## Test Plan

Add `test/expert/laurent_lazy_peel_certificate.jl` and register it in the
expert test group.

The focused test should:

- load `test/fixtures/laurent_lazy_determinant_cases.jl`;
- use `determinant-one-triangular` to compare against the existing Laurent `SL`
  peel route;
- call `_laurent_determinant_deferred_peel_certificate(A)`;
- assert `determinant_source == :deferred_submatrix`;
- assert at least one peel step exists and the deferred submatrix is smaller
  than the original matrix;
- assert the deferred submatrix equals the first existing `SL` peel route
  block and that the first recorded left/right factors match exactly;
- assert the deferred replay verifier returns true;
- assert `left_product * A * right_product == target_matrix`;
- assert `target_matrix == block_embedding(deferred_submatrix, n, 1:(n - steps))`;
- tamper with one recorded left or right elementary factor and assert verifier
  rejection.

Focused verification:

```bash
julia --project=. -e 'include("test/expert/laurent_lazy_peel_certificate.jl")'
```

Package verification:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Out Of Scope

Do not add public exports. Do not normalize or correct a non-one deferred
determinant. Do not add row/column determinant-correction options. Do not change
`laurent_gl_factorization_certificate` or the public `elementary_factorization`
surface.

## Automatic Decisions

- Visual companion: skipped because this is an algebraic replay metadata design
  with no visual interface.
- Clarifying questions: resolved from the issue text because Agent Desk marked
  the run non-interactive.
- Approach: choose a sibling internal certificate type instead of extending the
  existing `SL` certificate because the stopped-at-deferred-block invariant is
  different from the full determinant-one factorization invariant.
- Execution path: choose the recommended Superpowers Subagent-Driven path when
  the plan offers it, falling back only if unavailable.

## Spec Self-Review

- No placeholders or incomplete sections remain.
- The design is scoped to an internal metadata path, verifier, expert test, and
  runner registration.
- The replay relation is explicit and covers the issue's required
  `blockdiag(deferred_submatrix, I)` evidence.
- The negative control checks behavior by mutating an elementary factor, not by
  checking field presence.
