# Issue 248 Public ECP Column Route Design

Issue #248 routes public ordinary-polynomial column reduction through the
checked #243-#247 ECP pipeline without fixture-id branching. The change keeps
the existing `reduce_unimodular_column(v, R)` factor-returning API and the
existing `ECPColumnReductionCertificate` public certificate boundary.

## Context

The repository now has replayable ECP input contexts, monicity normalization,
link witness extraction, link-step certificates, and induction/normality
certificates. Public dispatch still has one narrow fixture-preferred staged
path and otherwise returns the older exact reduction certificate. Dependency
#247 is merged on `main`.

There is no repository `AGENTS.md`; the README test instructions apply.

## Chosen Approach

Use the existing #243-#247 stages to build a new `:ecp_pipeline` stage inside
`ECPColumnReductionCertificate`. The public reducer returns that normal
certificate type, so existing consumers can keep using `certificate.factors`,
`certificate.final_column`, and `verify_ecp_column_reduction`.

The reducer keeps the low-risk legacy paths for:

- unit-entry reduction,
- witness-unit reduction,
- embedded three-block reduction,
- Laurent normalization and Laurent row preconditioning.

For ordinary-polynomial columns that are not handled by those preserved paths,
the reducer attempts the general ECP pipeline:

1. Build a checked ECP input context.
2. Produce a checked monicity normalization.
3. Extract or verify a link witness from the normalized column.
4. Build a link-step certificate.
5. Compose induction/normality.
6. Inverse-substitute the induction factors and append inverse-substituted
   coordinate-move factors so the final factors act on the original column.

Unsupported but unimodular ordinary-polynomial columns fail with an
`ArgumentError` whose diagnostics include the staged ECP pipeline failure.
Non-unimodular columns still fail during validation before route work.

## Link-Step Compatibility

The #246 link-step route can factor length-three obligations through the
general polynomial `SL_3` route. For length greater than three, the public ECP
pipeline will preserve exactness by recording a direct elementary endpoint
transport route. The verifier checks those endpoint factors directly and stores
route metadata, while length-three non-fixture obligations continue to use
polynomial `SL_3` route certificates.

This does not implement recursive matrix factorization for #186 or broaden
public Park-Woodburn matrix acceptance for #187.

## Polynomial Column Peel Boundary

`polynomial_column_peel` should stop treating public column reduction as only a
bare factor list. It should capture the ECP column certificate for each peel
step while keeping the older positional constructor shape available for tests
and older internal callers. This gives later #186 code a stable place to inspect
route metadata without knowing fixture families.

## Tests

Add expert coverage that:

- proves `ecp_column_reduction_certificate` and `reduce_unimodular_column`
  return factors reducing a non-fixture ordinary-polynomial length-four column
  to `e_n`;
- asserts that column is not handled by unit entry, witness unit, or embedded
  three-block reduction;
- checks route metadata on the public ECP certificate;
- confirms a non-unimodular column fails before route work;
- confirms an unsupported but unimodular ordinary column reports a staged ECP
  diagnostic;
- confirms tampering with an ECP certificate factor makes verification return
  `false`;
- preserves the existing Laurent and legacy exact-reduction coverage.

Verification commands:

```bash
julia --project=. -e 'include("test/expert/elementary_column_property.jl")'
julia --project=. -e 'include("test/expert/unimodular_reduction_exact.jl")'
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Self-Review

The scope is limited to public ordinary-polynomial column routing and metadata.
No new public API name is required. The design keeps explicit staged failures
instead of silently falling back after a failed general ECP attempt for columns
that reach that route. The out-of-scope #186 and #187 items remain out of scope.
