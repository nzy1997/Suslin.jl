# Issue 83 ECP Column Fixture Catalog Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a shared Elementary Column Property unimodular-column fixture catalog and exact validator for ordinary polynomial examples.

**Architecture:** Keep the catalog under `test/fixtures/` as test support, not public API. The fixture module builds exact `GF(2)[x,y]` and `QQ[x,y]` columns with named entries and metadata, while the internal validator reconstructs columns, checks unimodularity, monicity, witness equations, current reducer pass/staged-fail expectations, and negative controls.

**Tech Stack:** Julia, Oscar polynomial rings, Suslin column reduction helpers, Test stdlib.

## Global Constraints

- Do not implement a new reducer.
- Do not change `reduce_unimodular_column`.
- Do not add Quillen patching, a full public driver, or Laurent `GL_n` determinant correction.
- Do not add a public Suslin fixture API.
- Catalog file is `test/fixtures/ecp_column_cases.jl`.
- Validator file is `test/internal/ecp_column_fixtures.jl`.
- Register the validator in the `internal` group in `test/runtests.jl`.
- Catalog must validate at least eight named entries.
- Required valid fixture ids are `ecp-unit-entry-gf2`, `ecp-witness-unit-gf2`, `ecp-variable-change-monic-gf2`, `ecp-link-bezout-nonunit-witness-qq`, `ecp-longer-embedded-block-gf2`, `ecp-unsupported-unimodular-gf2`, `ecp-non-unimodular-gf2`, and `ecp-monic-first-entry-qq`.
- Required negative-control fixture ids are `ecp-corrupt-witness-control` and `ecp-corrupt-monicity-control`.
- Every valid case must have `id`, `kind`, `stage_coverage`, `ring_constructor`, `ring`, `variable_order`, `entries`, `column_order`, `monicity`, `witnesses`, `expected`, `source_refs`, and `consumer_issue_ids`.
- The validator must check `Suslin.is_unimodular_column(column, R) == true` for supported or staged-fail unimodular cases and `false` for the non-unimodular negative case.
- Cases with monicity metadata must check exact monicity in the selected variable; variable-change cases must reconstruct the substitution and compare the transformed entry.
- Cases with witness metadata must check exact equations such as `sum(w[i] * v[i]) == 1`.
- Negative controls must prove the validator rejects at least one corrupt witness coefficient and one corrupt monicity claim.
- The issue comment's link-theorem guardrail is binding: include a monic-first-entry case, a link-witness case with resultant/Bezout/coverage/path metadata, and a staged-failure case where those witnesses are missing.
- Focused validator command is `julia --project=. -e 'include("test/internal/ecp_column_fixtures.jl")'`.
- Package verification command is `julia --project=. -e 'using Pkg; Pkg.test()'`.

---

## File Structure

- Create `test/internal/ecp_column_fixtures.jl`: owns validation helpers, focused tests, current reducer checks, fixture id coverage checks, and negative controls.
- Create `test/fixtures/ecp_column_cases.jl`: owns exact ordinary polynomial fixture construction and witness/monicity metadata.
- Modify `test/runtests.jl`: add `internal/ecp_column_fixtures.jl` to the internal group after `internal/sl3_murthy_gupta_fixtures.jl`.

---

### Task 1: ECP Column Fixture Catalog And Validator

**Files:**
- Create: `test/internal/ecp_column_fixtures.jl`
- Create: `test/fixtures/ecp_column_cases.jl`
- Modify: `test/runtests.jl`

**Interfaces:**
- Consumes: `Suslin.is_unimodular_column`, `Suslin.reduce_unimodular_column`, `identity_matrix`, `matrix`, Oscar polynomial APIs, `Test`.
- Produces: `ECPColumnFixtureCatalog.catalog()` returning `(; cases, negative_controls)` where `cases` is the valid fixture list.
- Produces: `ECPColumnFixtureCatalog.cases_by_id()` returning `Dict(entry.id => entry for entry in catalog().cases)`.
- Produces: internal validator functions `validate_ecp_column_fixture(entry)` and `validate_ecp_column_fixture_catalog(catalog)`.

- [ ] **Step 1: Write the failing validator**

Create `test/internal/ecp_column_fixtures.jl` first. It must include the catalog path, define validation helpers, and include a testset named `"ECP column fixture catalog"`.

The file must start with:

```julia
using Test
using Suslin
using Oscar

const ECP_COLUMN_CATALOG_PATH = joinpath(@__DIR__, "..", "fixtures", "ecp_column_cases.jl")
const REQUIRED_ECP_COLUMN_FIELDS = (
    :id,
    :kind,
    :stage_coverage,
    :ring_constructor,
    :ring,
    :variable_order,
    :entries,
    :column_order,
    :monicity,
    :witnesses,
    :expected,
    :source_refs,
    :consumer_issue_ids,
)
```

Implement helper functions with these exact names and responsibilities:

```julia
_ecp_field(entry, field::Symbol)
_ecp_column(entry)
_ecp_column_matrix(column, R)
_ecp_target_column(R, n::Int)
_ecp_factor_product(factors, R, n::Int)
_ecp_apply_factors(factors, column, R)
_ecp_variable_index(entry, variable_name::Symbol)
_ecp_monic_in_variable(p, R, variable_name::Symbol)
_ecp_substitution_values(entry, substitution)
_ecp_assert_metadata(entry)
_ecp_assert_unimodularity(entry)
_ecp_assert_current_status(entry)
_ecp_assert_monicity(entry)
_ecp_assert_witness(entry, witness)
_ecp_assert_witnesses(entry)
validate_ecp_column_fixture(entry)
validate_ecp_column_fixture_catalog(catalog)
```

The witness validator must support these witness kinds:

```julia
:ideal_membership
:link_bezout
:missing_link_witness
```

For `:ideal_membership` and `:link_bezout`, compute the column in
`entry.column_order`, require `length(witness.coefficients) == length(column)`,
and require:

```julia
sum(witness.coefficients[idx] * column[idx] for idx in eachindex(column); init = zero(R)) == one(R)
```

If `witness.require_nonunit_coefficients == true`, require
`!any(is_unit, witness.coefficients)`. For `:link_bezout`, also require
`is_unit(witness.resultant)`, `witness.coverage.covers_unit_ideal == true`, and
nonempty `witness.path`. For `:missing_link_witness`, require
`entry.expected.current_status == :staged_fail` and a nonempty
`witness.missing` tuple.

The testset must include the catalog file, validate the catalog, assert all
required valid ids and negative-control ids are present, assert there are at
least eight valid cases, mutate one valid witness coefficient and one valid
monicity transformed entry at runtime, and check both mutations throw
`ArgumentError`.

- [ ] **Step 2: Run focused test to verify RED**

Run:

```bash
julia --project=. -e 'include("test/internal/ecp_column_fixtures.jl")'
```

Expected: FAIL because `test/fixtures/ecp_column_cases.jl` does not exist yet.

- [ ] **Step 3: Implement the catalog**

Create `test/fixtures/ecp_column_cases.jl` with `module ECPColumnFixtureCatalog`.

The catalog must define helper constructors:

```julia
_ring_metadata(description, R, generator_names, generators)
_ordinary_ring_constructor(coefficient, variables)
_case(; id, kind, stage_coverage, ring_constructor, ring, variable_order, entries, column_order, monicity, witnesses, expected, source_refs, consumer_issue_ids)
_negative_control(id, base_case_id, reason, entry)
catalog()
cases_by_id()
```

Use these exact seed columns:

```julia
R2, (x, y) = Oscar.polynomial_ring(GF(2), ["x", "y"])
RQ, (X, Y) = Oscar.polynomial_ring(QQ, ["x", "y"])

# ecp-unit-entry-gf2
(a = x, b = y, c = one(R2))

# ecp-witness-unit-gf2
(a = x, b = y, c = x + one(R2))
# witness coefficients (one(R2), zero(R2), one(R2))

# ecp-variable-change-monic-gf2
(a = x + y^2, b = x * y + x + one(R2), c = x^2 + x * y + y + one(R2))
# substitution x -> x + y, selected entry a, transformed entry x + y + y^2

# ecp-link-bezout-nonunit-witness-qq
(a = one(RQ) - X, b = X, c = X * Y)
# witness coefficients (one(RQ) + X, X, zero(RQ))
# equation (1 + X) * (1 - X) + X * X == 1

# ecp-longer-embedded-block-gf2
(a = x + y^2, b = x * y + x + one(R2), c = x^2 + x * y + y + one(R2),
 d = x^2, e = x * y, f = y^2 + x)

# ecp-unsupported-unimodular-gf2
(a = zero(R2), b = x^2, c = x * y + one(R2))

# ecp-non-unimodular-gf2
(a = x, b = y, c = x * y)

# ecp-monic-first-entry-qq
(a = X^2 + Y + one(RQ), b = X, c = Y)
```

Use `expected.current_status = :passes`, `:staged_fail`, or
`:rejects_non_unimodular` as appropriate. The staged-fail message substring is
`"unsupported exact unimodular column reduction"`. The non-unimodular rejection
message substring is `"v must be a unimodular column"`.

Set `consumer_issue_ids` to include `"#62"` for reducer-facing cases,
`"#87"` for monic-first-entry cases, and `"#88"` for link-witness or
missing-link-witness cases.

Create `negative_controls` by merging valid entries with corrupt metadata:

```julia
bad_witness = merge(
    witness_unit_case,
    (;
        id = "ecp-corrupt-witness-control",
        witnesses = (merge(only(witness_unit_case.witnesses), (;
            coefficients = (zero(R2), zero(R2), one(R2)),
        )),),
    ),
)

bad_monicity = merge(
    variable_change_case,
    (;
        id = "ecp-corrupt-monicity-control",
        monicity = merge(variable_change_case.monicity, (;
            transformed_entry = y,
        )),
    ),
)
```

- [ ] **Step 4: Register the internal validator**

Modify `test/runtests.jl` and add:

```julia
"internal/ecp_column_fixtures.jl",
```

after `"internal/sl3_murthy_gupta_fixtures.jl",`.

- [ ] **Step 5: Run focused test to verify GREEN**

Run:

```bash
julia --project=. -e 'include("test/internal/ecp_column_fixtures.jl")'
```

Expected: PASS and validate at least eight named catalog entries.

- [ ] **Step 6: Run package verification**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: PASS.

- [ ] **Step 7: Prepare commit or remote commit**

Stage these intended files only:

```bash
docs/superpowers/specs/2026-06-22-issue-83-ecp-column-fixture-catalog-design.md
docs/superpowers/plans/2026-06-22-issue-83-ecp-column-fixture-catalog.md
test/fixtures/ecp_column_cases.jl
test/internal/ecp_column_fixtures.jl
test/runtests.jl
```

If the sandbox permits local git writes, commit:

```bash
git commit -m "Add ECP column fixture catalog"
```

If local git writes are blocked by the managed worktree, create the equivalent
remote branch commit through the GitHub app with the same file set.

---

## Self-Review

- Spec coverage: the task covers fixture construction, validator functions,
  all required fixture ids, exact unimodularity checks, monicity checks,
  witness checks, current reducer pass/staged-fail checks, and negative
  controls.
- Placeholder scan: no TBD, TODO, fill-in, or placeholder steps remain.
- Type consistency: `ECPColumnFixtureCatalog.catalog()` and
  `validate_ecp_column_fixture_catalog(catalog)` are named consistently across
  the plan.
