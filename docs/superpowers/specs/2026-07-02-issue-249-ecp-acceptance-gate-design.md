# Issue 249 ECP Acceptance Gate Design

Issue #249 closes parent #185 at the acceptance/documentation layer. The
implementation already routes public ordinary-polynomial unimodular column
reduction through the checked general ECP pipeline after #248. This gate proves
that route is usable by later matrix-peeling work and records the boundary
without broadening public Park-Woodburn claims.

## Context

There is no repository `AGENTS.md`; the README test instructions apply. The
current `main` includes #248, which added `:ecp_pipeline` stages to
`ECPColumnReductionCertificate` and taught polynomial column peel steps to store
the ECP left-column certificate. The existing expert tests already exercise a
representative length-four ordinary-polynomial ECP success case and staged
negative controls.

GitHub issue/PR fetching is unavailable in this Agent Desk sandbox because the
GitHub CLI cannot reach the configured proxy. The issue body supplied by the
run, local merge history for #248/#247, and repository docs/tests are the source
of truth for this design.

## Chosen Approach

Use a parent-level gate instead of algorithm expansion:

1. Strengthen existing ECP expert coverage with explicit determinant and route
   boundary checks.
2. Add a consumer smoke test in the `SL_n`/column-peel area that uses an
   `n > 3` ordinary-polynomial matrix and verifies its peel step contains a
   checked ECP left-column certificate.
3. Update README and Documenter scope text so #185 is no longer listed as a
   staged boundary, while #186 recursive `SL_n` matrix factorization and #187
   public Park-Woodburn acceptance remain staged.
4. Add a short parent coverage note mapping #185 child issues to the ECP
   algorithm stages and recording explicit non-claims.

This is preferable to changing reducer internals because the parent issue asks
for acceptance coverage and documentation after #248, not new algorithmic
support. If a new test exposes a real gap, fix that gap narrowly; otherwise keep
production code untouched.

## Acceptance Surface

The gate covers these public/expert calls:

- `reduce_unimodular_column(v, R)` returns factors reducing a representative
  ordinary-polynomial unimodular column to the last standard basis vector.
- `ecp_column_reduction_certificate(v, R)` returns a verified
  `ECPColumnReductionCertificate` whose final stage is `:ecp_pipeline`.
- ECP certificate verifiers reject tampered factor sequences.
- The existing polynomial column-peel consumer records and verifies the ECP
  certificate used for its last-column peel step.

Negative controls remain explicit:

- Non-unimodular ordinary-polynomial columns fail before route work.
- Unsupported but unimodular ordinary-polynomial columns still fail cleanly with
  staged ECP diagnostics.
- Determinant-not-one matrix inputs remain outside the polynomial
  column-peel/route certificate path.
- Laurent `GL_n`, ToricBuilder mainline support, #186 recursive matrix
  factorization, and #187 full public Park-Woodburn support are not claimed.

## Tests

Add focused tests in the suggested expert files:

- `test/expert/elementary_column_property.jl`: keep the representative ECP
  success case, assert determinant-one on its certificate factor product, and
  keep the tampered certificate verifier failure.
- `test/expert/unimodular_reduction_exact.jl`: keep staged non-unimodular and
  unsupported negative controls.
- `test/expert/sln_to_sl3_reduction.jl`: add a length `n > 3` consumer smoke
  case using the existing polynomial recursive column-peel certificate and
  check that the first peel step stores a verified ECP left certificate.

Verification commands:

```bash
julia --project=. -e 'include("test/expert/elementary_column_property.jl")'
julia --project=. -e 'include("test/expert/unimodular_reduction_exact.jl")'
julia --project=. -e 'include("test/expert/sln_to_sl3_reduction.jl")'
julia --project=. test/runtests.jl all
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Documentation

README and `docs/src/index.md` should say that #185 covers the
ordinary-polynomial unimodular-column reducer and certificate evidence needed by
later consumers. They must also say that recursive `SL_n` matrix factorization
(#186), final public Park-Woodburn acceptance (#187), arbitrary Laurent `GL_n`,
and ToricBuilder mainline support remain staged boundaries.

The parent coverage note belongs in `docs/audits/2026-07-02-issue-185-ecp-acceptance.md`.
It should map #242 through #248 to the implemented ECP stages and list #249 as
the parent acceptance/documentation gate.

## Self-Review

The design is scoped to tests and docs unless a test reveals a defect. It does
not add public APIs, change certificate structures, or reclassify unsupported
matrix routes. The negative controls keep the boundary wording precise and
prevent accidental claims about #186, #187, Laurent `GL_n`, or ToricBuilder.
