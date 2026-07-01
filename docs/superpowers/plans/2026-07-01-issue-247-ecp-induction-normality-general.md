# Issue 247 ECP Induction Normality General Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Compose #246 ECP link steps with verified lower-variable induction and #181/#193 normality certificates.

**Architecture:** Extend the existing ECP induction/normality certificate in `src/algorithm/column_reduction.jl`. A verified link step supplies `v(X) -> v(0)`, a recorded descent measure bounds recursion, lower reduction is canonicalized to a verified `ECPColumnReductionCertificate`, and normality is constructed through `realize_conjugate_elementary_certificate` when no explicit witness is supplied.

**Tech Stack:** Julia, Oscar polynomial rings, existing ECP link-step certificates, existing ECP column-reduction certificates, existing #193/#181 conjugated-elementary normality certificates, Julia `Test`.

## Global Constraints

- Keep APIs internal and unexported; use `Suslin.<name>` from expert tests only.
- Do not change public reducer dispatch or implement recursive matrix factorization for #186.
- Input support remains ordinary polynomial rings only.
- A certificate must record and verify an explicit selected-variable descent measure.
- Passing from `v(X)` to `v(0)` must strictly reduce the selected-variable profile; a same-context recursive call must be rejected as a staged failure.
- The constructor must build a lower-variable reduction when no lower hint is supplied, and staged failures must distinguish missing lower-variable reduction from missing normality rewrite.
- The constructor must build a normality witness from the verified lower reduction and route it through `realize_conjugate_elementary_certificate` when no normality witness is supplied.
- The verifier must independently check the lower ECP certificate, lifted lower factors, #181/#193 normality certificate, final concatenated factors, and exact application to `e_n`.
- Focused verification commands are `julia --project=. -e 'include("test/expert/ecp_induction_normality_general.jl")'` and `julia --project=. -e 'include("test/expert/normality.jl")'`.
- Package verification command is `julia --project=. -e 'using Pkg; Pkg.test()'`.

---

## File Structure

- Modify `src/algorithm/column_reduction.jl`: add descent helpers, a descent field on `ECPInductionNormalityCertificate`, lower-reduction canonicalization, automatic normality witness construction, staged-failure message wrappers, and verifier checks.
- Create `test/expert/ecp_induction_normality_general.jl`: #242/#246 positive composition coverage and negative controls.
- Modify `test/expert/ecp_induction_normality.jl`: update legacy expectations for automatic normality construction and canonical lower certificates where needed.
- Modify `test/runtests.jl`: register the new expert test.
- Add this plan and `docs/superpowers/specs/2026-07-01-issue-247-ecp-induction-normality-general-design.md`.

### Task 1: Add Red General Composition Tests

**Files:**
- Create: `test/expert/ecp_induction_normality_general.jl`
- Modify: `test/expert/ecp_induction_normality.jl`
- Modify: `test/runtests.jl`

**Interfaces:**
- Consumes future `ECPInductionNormalityCertificate.descent_measure`.
- Consumes future automatic `normality_witness === nothing` behavior.
- Produces failing assertions for #242/#246 link-step composition, descent replay, automatic normality construction, and staged-failure diagnostics.

- [ ] **Step 1: Write the failing #242 general test**

Create `test/expert/ecp_induction_normality_general.jl` with helpers that load `test/fixtures/ecp_mainline_cases.jl`, build the `ecp-mainline-sl3-route-qq` column, call `ecp_link_step_certificate(column, R; link_witness = witness, route_mode = :polynomial_sl3)`, and call `ecp_induction_normality_certificate(column, R; link_step = link)` with no lower or normality hints.

The positive assertions must check:

```julia
@test cert isa Suslin.ECPInductionNormalityCertificate
@test cert.descent_measure.strict_descent
@test cert.descent_measure.parent_profile != cert.descent_measure.lower_profile
@test cert.lower_reduction_certificate isa Suslin.ECPColumnReductionCertificate
@test Suslin.verify_ecp_column_reduction(cert.lower_reduction_certificate)
@test cert.lower_reduction_certificate.original_column == collect(link.lower_variable_column)
@test cert.normality_witness.source == :constructed_normality_witness
@test cert.normality_certificate isa Suslin.ConjugatedElementaryNormalityCertificate
@test Suslin.verify_conjugate_elementary_certificate(cert.normality_certificate)
@test cert.normality_rewrite.normality_certificate == cert.normality_certificate
@test Suslin.verify_ecp_induction_normality_certificate(cert)
@test _apply_factors(cert.final_factors, column, R) == Suslin._target_reduced_column(R, length(column))
```

- [ ] **Step 2: Add negative controls**

In the new test, rebuild records by replacing struct fields and verify rejection for:

```julia
@test !Suslin.verify_ecp_induction_normality_certificate(tampered_descent)
@test !Suslin.verify_ecp_induction_normality_certificate(tampered_lower_certificate)
@test !Suslin.verify_ecp_induction_normality_certificate(tampered_lifted_factors)
@test !Suslin.verify_ecp_induction_normality_certificate(tampered_normality_witness)
@test !Suslin.verify_ecp_induction_normality_certificate(tampered_normality_certificate)
@test !Suslin.verify_ecp_induction_normality_certificate(tampered_final_factor)
```

Also assert staged-failure messages:

```julia
@test_throws Regex("same-context recursive lower-variable call") Suslin.ecp_induction_normality_certificate(
    constant_column,
    constant_R;
    link_step = constant_link,
)
@test_throws Regex("missing lower-variable reduction") Suslin.ecp_induction_normality_certificate(
    unsupported_column,
    unsupported_R;
    link_step = unsupported_link,
)
@test_throws Regex("missing normality rewrite") Suslin.ecp_induction_normality_certificate(
    column,
    R;
    link_step = link,
    lower_reduction = lower,
    normality_witness = bad_normality_witness,
)
```

- [ ] **Step 3: Update legacy test expectations**

In `test/expert/ecp_induction_normality.jl`, replace the no-normality-witness throw with an automatic construction assertion:

```julia
automatic_normality = Suslin.ecp_induction_normality_certificate(
    qq.column,
    qq.R;
    link_step = qq.link,
    lower_reduction = qq.lower,
)
@test automatic_normality.normality_witness.source == :constructed_normality_witness
@test Suslin.verify_ecp_induction_normality_certificate(automatic_normality)
```

Keep the malformed lower-reduction and identity-entry normality controls as throws.

- [ ] **Step 4: Register the new expert test**

Add `"expert/ecp_induction_normality_general.jl"` immediately before `"expert/ecp_induction_normality.jl"` in `test/runtests.jl`.

- [ ] **Step 5: Run RED verification**

Run:

```bash
julia --project=. -e 'include("test/expert/ecp_induction_normality_general.jl")'
```

Expected: FAIL because `descent_measure` is not yet a field and missing normality data is not yet constructed.

### Task 2: Implement Descent, Lower Canonicalization, And Automatic Normality

**Files:**
- Modify: `src/algorithm/column_reduction.jl`
- Test: `test/expert/ecp_induction_normality_general.jl`
- Test: `test/expert/ecp_induction_normality.jl`

**Interfaces:**
- Produces `_ecp_induction_descent_measure(link_step, R)`, `_ecp_selected_variable_profile(column, R, variable_index)`, `_ecp_descent_measure_strict(measure)`, `_ecp_construct_normality_witness(lower_factors, n, R, selected_variable)`, and extended `ECPInductionNormalityCertificate`.

- [ ] **Step 1: Add descent helpers and certificate field**

Add `descent_measure` after `lower_variable_column` in `ECPInductionNormalityCertificate`. Implement helpers that compute selected-variable degree profiles with `degree(entry, variable_index)`, store column length and variable count, and require componentwise non-increase plus at least one strict decrease.

- [ ] **Step 2: Reject non-descending link steps**

In `ecp_induction_normality_certificate`, compute descent immediately after extracting `lower_column`. Throw:

```julia
ArgumentError("ECP induction/normality staged failure: same-context recursive lower-variable call did not strictly reduce selected-variable profile")
```

when strict descent fails.

- [ ] **Step 3: Canonicalize lower reductions**

Update `_ecp_verified_lower_reduction` so it always returns a verified `ECPColumnReductionCertificate` and its factors. For raw factor-list hints, verify the factors are elementary and reduce `lower_column`, then construct the canonical lower certificate and require factor equality. Wrap automatic lower-construction failures in:

```julia
ArgumentError("ECP induction/normality staged failure: missing lower-variable reduction: " * sprint(showerror, err))
```

- [ ] **Step 4: Construct normality witness automatically**

Add `_ecp_construct_normality_witness` and use it when `normality_witness === nothing`. The witness must have `source = :constructed_normality_witness`, `conjugator = inv(_factor_sequence_product(lifted_lower_factors, R, n))`, `sl2_indices = (n, 1)`, and `sl2_entry = selected_variable + one(R)`.

Update `_ecp_induction_normality_rewrite` to accept `:supplied_normality_witness` and `:constructed_normality_witness`, and to return the stored source.

- [ ] **Step 5: Wrap normality failures distinctly**

In the constructor, wrap errors from `_ecp_induction_normality_rewrite` in:

```julia
ArgumentError("ECP induction/normality staged failure: missing normality rewrite: " * sprint(showerror, err))
```

while preserving `InterruptException`.

- [ ] **Step 6: Extend replay verification**

In `_ecp_induction_normality_replay_summary`, recompute descent, verify stored descent equality and strictness, require `lower_reduction_certificate isa ECPColumnReductionCertificate`, verify that certificate independently, verify lifted factors from the replayed lower certificate, and include `descent_measure_ok`, `descent_strict_ok`, and `lower_reduction_certificate_ok` in `overall_ok`.

- [ ] **Step 7: Run focused GREEN tests**

Run:

```bash
julia --project=. -e 'include("test/expert/ecp_induction_normality_general.jl")'
julia --project=. -e 'include("test/expert/ecp_induction_normality.jl")'
```

Expected: both PASS.

### Task 3: Verification, Review, And Commit

**Files:**
- Modify: `src/algorithm/column_reduction.jl`
- Modify: `test/expert/ecp_induction_normality.jl`
- Create: `test/expert/ecp_induction_normality_general.jl`
- Modify: `test/runtests.jl`
- Create: `docs/superpowers/specs/2026-07-01-issue-247-ecp-induction-normality-general-design.md`
- Create: `docs/superpowers/plans/2026-07-01-issue-247-ecp-induction-normality-general.md`

**Interfaces:**
- Consumes the finished implementation from Task 2.
- Produces a verified branch ready for PR.

- [ ] **Step 1: Run issue-required focused commands**

Run:

```bash
julia --project=. -e 'include("test/expert/ecp_induction_normality_general.jl")'
julia --project=. -e 'include("test/expert/normality.jl")'
```

Expected: both PASS.

- [ ] **Step 2: Run package verification**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: PASS.

- [ ] **Step 3: Review scope and commit**

Run:

```bash
git status --short
git diff --stat
git add docs/superpowers/specs/2026-07-01-issue-247-ecp-induction-normality-general-design.md docs/superpowers/plans/2026-07-01-issue-247-ecp-induction-normality-general.md src/algorithm/column_reduction.jl test/expert/ecp_induction_normality_general.jl test/expert/ecp_induction_normality.jl test/runtests.jl
git commit -m "Compose general ECP induction normality"
```

Expected: commit succeeds with only intended files staged.

## Plan Self-Review

- Every hard acceptance requirement maps to a focused assertion or replay check.
- The plan preserves TDD by adding the failing general test before production changes.
- The selected approach uses existing link-step, lower ECP, and #181/#193 normality certificate APIs.
- Staged-failure text is explicit for same-context descent, missing lower reduction, and missing normality rewrite.
- No incomplete markers remain.
