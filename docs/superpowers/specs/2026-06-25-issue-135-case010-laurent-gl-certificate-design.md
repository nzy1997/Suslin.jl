# Issue 135 Case 010 Laurent GL Certificate Design

## Goal

Route the original ToricBuilder `case_010` Q-block through
`laurent_gl_factorization_certificate(A)` so the certificate verifies exactly
and the Q-block status report records `route_status == :gl_certificate_pass`.

## Context

Issue #134 added a certified `:laurent_unit_creation` column-reduction stage for
the preserved length-5 `case_010` Laurent boundary column. That change is
merged into `main` and the current branch starts after it.

Running the current full report probe still produces a route error, but the
remaining boundary is later in the same Laurent column-peel replay. The first
two peel dimensions are supported:

- dimension 6 uses the existing `:unit_entry` stage;
- dimension 5 uses the #134 `:laurent_unit_creation` stage;
- dimension 4 fails as an unsupported Laurent-normalized column.

The dimension-4 column is:

```julia
[
    u^-1 + u^-1*v^-1,
    v^3 + v^2 + v,
    v^4 + v^3 + v^2,
    v^2 + 1 + v^-1 + v^-2,
]
```

It has the same one-step exact unit-creation shape as the #134 boundary:

```julia
(v^3 + v^2 + v) + (u*v^3 + u*v) * (u^-1 + u^-1*v^-1) == 1
```

So this issue does not require a second independent Laurent algorithm boundary.
It can reuse the existing certificate/replay machinery with a narrower length
guard.

## Approach Options

Recommended: extend the certified Laurent unit-creation stage so it also
accepts the observed length-4 `case_010` boundary. Keep the algorithm otherwise
unchanged: exact division finds a coefficient, replay recomputes it from the
input column, and all factors still pass `_checked_reduction_factors`.

Alternative: add a second `case_010`-specific hard-coded factor sequence for the
length-4 column. This is smaller but worse, because it bypasses the reusable
stage and would be harder to replay safely.

Alternative: implement a more general Laurent unimodular-column algorithm. This
is out of scope for the issue and would require broader tests before claiming
support.

## Chosen Design

Update `_laurent_unit_creation_candidate(column, R)` in
`src/algorithm/column_reduction_case010.jl` so the existing stage can run for
length-4 and length-5 Laurent columns. The ring guard, exact-division guard,
nonzero coefficient guard, and exact unit equality check stay intact.

This keeps the public meaning of `elementary_factorization(case_010)`
unchanged. The public route remains `:staged_boundary`; only the explicit
Laurent GL certificate route is promoted.

## Certificate And Replay

No new certificate kind is needed. The existing `:laurent_unit_creation` stage
already records:

- the original input column;
- `pivot_index` and `source_index`;
- `target_unit == one(R)`;
- the exact `creation_coefficient`;
- the creation factor;
- the created column;
- the nested unit-entry stage;
- the complete factor sequence and output column.

Replay recomputes the coefficient and factor sequence from the stored indices
and the current input column. Extending the allowed column lengths does not
weaken replay: a tampered factor sequence or corrupted original matrix still
fails exact verification.

## Tests

Add `test/internal/toricbuilder_cache_case010_certificate.jl` for the issue
acceptance path. It should:

- materialize the original `case_010` from
  `test/fixtures/toricbuilder_cache_q_blocks.jl`;
- build `laurent_gl_factorization_certificate(A)`;
- assert `verify_laurent_gl_factorization_certificate(certificate)`;
- assert exact reconstruction of the original matrix;
- assert the normalized core determinant is `one(R)`;
- assert the decomposed base-matrix count is positive;
- assert the public `elementary_factorization(A)` route still reports the
  staged Laurent GL boundary;
- corrupt one sparse coordinate in a fresh materialized matrix and assert the
  original certificate does not verify for that corrupted matrix.

Update `test/internal/toricbuilder_cache_status_report.jl` so `case_010` expects
`gl_certificate_pass`, `verified == true`, a positive decomposed base-matrix
count, and no `Route Error Details` entry.

Refresh `docs/audits/2026-06-24-toricbuilder-cache-q-block-status.md` with the
report script after implementation.

## Verification

Issue verification:

```bash
julia --project=. -e 'include("test/internal/toricbuilder_cache_case010_certificate.jl")'
julia --project=. scripts/report_toricbuilder_cache_q_blocks.jl --exercise=case_010 --output=/tmp/case010-q-block-status.md
```

Package verification:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Out Of Scope

Do not make `elementary_factorization(case_010)` return factors for the original
Laurent `GL_n` input. Do not add a separate broad Laurent column algorithm. Do
not claim general Laurent unimodular-column support beyond the existing staged
boundaries.

## Automatic Decisions

- Clarifying questions were answered by the Standing Answer Policy because this
  is a non-interactive Agent Desk run.
- The visual companion was skipped because no visual decision would clarify the
  algebraic certificate path.
- The recommended approach was selected because the remaining length-4 boundary
  uses the same exact unit-creation relation as #134 and preserves certificate
  replay.
