# Issue 206 Murthy Local Contract Catalog Design

## Context

Issue #206 extends the existing Murthy-Gupta fixture catalog from the completed
ordinary `QQ[X]` local `SL_3` solver into a shared Park-Woodburn Section 5
contract catalog for parent #182 and child issues #207 through #210. The work is
catalog-only. It must not claim that the full local-ring contract or the
Park-Woodburn driver is implemented.

The current catalog in `test/fixtures/sl3_murthy_gupta_cases.jl` already records
ordinary pass cases over `QQ[X]`. The validator in
`test/internal/sl3_murthy_gupta_fixtures.jl` reconstructs each special-form
target matrix, checks determinant one and monicity, validates q-degree,
split-lemma, q(0)-unit, and q(0)-nonunit witness equations, and checks current
solver behavior.

## Design Choice

Extend the same catalog and validator instead of creating a second fixture
module. Existing fixture ids stay stable. New ids are additive and marked with
`local_contract = true`, `expected_current_solver.status = :staged_fail`, and
`consumer_issue_ids = ("#182", "#208", "#207", "#209", "#210")`.

The catalog will add three Section 5 local-contract cases over `QQ[u, X]`:

- A multivariate q-degree case where `q = X*p + 1`, so the elementary
  q-degree normalization is a polynomial-ring identity over the current base
  ring even though the full solver remains staged for multivariate local input.
- A q(0)-unit local case where `q(0) = 1 + u` is not a global polynomial-ring
  unit but is a unit in the local context at the maximal ideal `(u)`.
- A q(0)-nonunit Bezout/resultant case with explicit `p_prime`, `q_prime`,
  degree guards, and a local branch-unit witness for `q(0) + p_prime(0) = 1 + u`.

The local unit witness schema is intentionally explicit and polynomial-only:

- `context.kind = :localization_at_maximal_ideal`
- `context.maximal_ideal_generators = (u,)`
- `unit`, `residue_unit`, `residue_inverse`
- `residue_difference_coefficients`
- exact equation `unit - residue_unit == sum(coeff_i * generator_i)`
- exact residue inverse equation `residue_unit * residue_inverse == 1`
- `global_unit = false` when the point is that the element is local but not a
  global polynomial-ring unit

This records the data later #208 can adopt without adding denominator-aware
replay in this issue.

## Validation Rules

The validator keeps all existing checks and adds:

- constants are evaluated in the selected Section 5 variable, so in `QQ[u, X]`
  the `X = 0` value may still depend on `u`;
- local-contract entries must declare staged current solver status and include
  #182 plus the child issue ids in their consumer metadata;
- local unit witnesses are checked by exact residue and maximal-ideal
  combination equations;
- q(0)-unit local witnesses can record that the required inverse is local,
  without pretending a polynomial inverse exists;
- q(0)-nonunit witnesses accept either a global unit or an explicit local-unit
  witness for resultant and branch-unit data;
- catalog-level negative controls must be rejected by the same positive-case
  validator.

Negative controls will cover non-monic `p`, determinant not one, corrupted
split witness, corrupted local-unit witness, and corrupted Bezout equality.

## Files

- Modify `test/fixtures/sl3_murthy_gupta_cases.jl`: add local-contract cases,
  local witness constructors, and named negative controls.
- Modify `test/internal/sl3_murthy_gupta_fixtures.jl`: extend validation to
  selected-variable constants, local witness records, local-contract metadata,
  and catalog negative controls.
- Add this design and a matching implementation plan under
  `docs/superpowers/`.

## Out Of Scope

- No new Murthy solving behavior.
- No denominator-aware replay.
- No Quillen patching, ECP, recursive `SL_n`, or public Park-Woodburn driver.
- No public Suslin fixture API.

## Verification

Focused validator:

```bash
julia --project=. -e 'include("test/internal/sl3_murthy_gupta_fixtures.jl")'
```

Package verification:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Spec Self-Review

- No open markers remain.
- Scope is limited to fixtures, internal validation, and Superpowers docs.
- Existing fixture ids and pass expectations stay unchanged.
- New local-contract cases are explicit staged failures for the current solver.
- Negative controls exercise malformed hypotheses rather than solver behavior.
