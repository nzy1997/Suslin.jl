# Issue 185 ECP Acceptance Coverage

Date: 2026-07-02

Issue #249 is the parent-level closeout gate for #185. It records that the
ordinary-polynomial ECP unimodular-column reducer is accepted for the staged
exact field-backed ordinary-polynomial route, and it records what remains out
of scope for later matrix-level work.

## Stage Map

| Issue | Stage | Evidence boundary |
| --- | --- | --- |
| #242 | ECP mainline catalog | Fixture/catalog entries for representative ordinary-polynomial columns |
| #243 | Input context | Checked ordinary-polynomial ring, column, variable order, and unimodularity context |
| #244 | Monicity normalization | Replayable coordinate changes that make the selected entry monic |
| #245 | Link witness extraction | Replayable Park-Woodburn link witness data and cover evidence |
| #246 | Link step | Exact link-step factors and route metadata, including direct endpoint transport for length greater than three |
| #247 | Induction/normality | Lower-column reduction plus conjugated-elementary normality replay |
| #248 | Public reducer route | `reduce_unimodular_column` and `ecp_column_reduction_certificate` dispatch through the general ECP pipeline |
| #249 | Parent acceptance gate | Tests and docs proving the #185 reducer boundary for later consumers |

## Acceptance Evidence

- `test/expert/elementary_column_property.jl` includes a representative
  ordinary-polynomial ECP success case, route metadata checks, determinant-one
  factor-product replay, and tampered certificate rejection.
- `test/expert/unimodular_reduction_exact.jl` keeps explicit non-unimodular and
  unsupported-but-unimodular negative controls.
- `test/expert/sln_to_sl3_reduction.jl` exercises a length `n > 3`
  polynomial column-peel consumer and verifies that the peel step stores a
  checked `ECPColumnReductionCertificate`.

## Non-Claims

- #186 recursive `SL_n` matrix factorization remains staged.
- #187 final public Park-Woodburn acceptance remains staged.
- Arbitrary Laurent `GL_n` determinant correction remains staged.
- ToricBuilder mainline support remains staged outside the documented
  certificate-backed Laurent slices.
- This gate does not optimize Steinberg factor counts.
