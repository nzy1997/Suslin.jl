# Issue 194 ECP Induction Normality Replayed Certificates Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Route the ECP induction/normality stage through #193-style replayed conjugated-elementary normality certificates.

**Architecture:** Preserve the existing staged ECP witness interface, convert legacy witness data into `ConjugatedElementaryNormalityCertificate` internally, and store that nested certificate on `ECPInductionNormalityCertificate`. Verification replays the lower reduction, reverifies or rebuilds the nested certificate, derives normality factors from the nested certificate, and then checks the final staged column reduction.

**Tech Stack:** Julia, Oscar, existing Suslin elementary matrix and normality certificate helpers.

## Global Constraints

- Keep backwards compatibility for callers passing `normality_witness` with `source`, `conjugator`, `sl2_indices`, and `sl2_entry`.
- Accept a supplied `normality_certificate` when it is a verified `ConjugatedElementaryNormalityCertificate` matching the lower-reduction conjugator and embedded `SL_2` contribution.
- Do not trust stored normality factor vectors without reverifying the nested #193 certificate.
- Final factors are `lifted_lower_variable_factors`, nested normality certificate factors, and `link_step.reduction_factors`, in that order.
- Verifiers must return `false` or throw a deliberate verification error for tampered nested normality data.
- Do not implement the general polynomial ECP reducer from issue 185 or any Murthy, Quillen, recursive `SL_n`, Laurent/ToricBuilder, or Steinberg optimization work.

---

### Task 1: Add Failing Nested-Certificate Tests

**Files:**
- Modify: `test/expert/ecp_induction_normality.jl`
- Modify: `test/expert/elementary_column_property.jl`

**Interfaces:**
- Consumes: existing `_fixture_certificate`, `_replace_record_field`, `_replace_rewrite_field`, and `_ecp_acceptance_tampered_staged_certificate` helpers.
- Produces: failing assertions that require `ECPInductionNormalityCertificate.normality_certificate`, `normality_rewrite.normality_certificate`, and the new `normality_certificate` constructor keyword.

- [ ] **Step 1: Add nested-certificate assertions to induction tests**

Add these checks after the first successful `qq` induction certificate assertions in `test/expert/ecp_induction_normality.jl`:

```julia
    @test qq.certificate.normality_certificate isa Suslin.ConjugatedElementaryNormalityCertificate
    @test qq.certificate.normality_rewrite.normality_certificate == qq.certificate.normality_certificate
    @test Suslin.verify_conjugate_elementary_certificate(qq.certificate.normality_certificate)
    @test qq.certificate.normality_rewrite.rewrite_factors == qq.certificate.normality_certificate.factors
    @test qq.certificate.normality_rewrite.rewrite_product == qq.certificate.normality_certificate.product
```

- [ ] **Step 2: Add supplied-certificate constructor coverage**

Add this block after `explicit_sequence` verifies:

```julia
    supplied_nested = Suslin.ecp_induction_normality_certificate(
        qq.column,
        qq.R;
        link_step = qq.link,
        lower_reduction = qq.lower,
        normality_witness = qq.witness,
        normality_certificate = qq.certificate.normality_certificate,
    )
    @test Suslin.verify_ecp_induction_normality_certificate(supplied_nested)
    @test supplied_nested.normality_certificate == qq.certificate.normality_certificate
    @test supplied_nested.final_factors == qq.certificate.final_factors
```

- [ ] **Step 3: Add nested tamper negative controls**

Add these helpers near the existing `_replace_rewrite_field` helper:

```julia
function _replace_conjugated_certificate_field(cert, field::Symbol, value)
    fields = fieldnames(typeof(cert))
    values = [getfield(cert, name) for name in fields]
    idx = findfirst(==(field), fields)
    idx === nothing && error("unknown conjugated certificate field: $(field)")
    values[idx] = value
    return typeof(cert)(values...)
end

function _tamper_nested_first_factor(cert, R, n::Int)
    factors = copy(cert.factors)
    factors[1] = identity_matrix(R, n)
    return _replace_conjugated_certificate_field(cert, :factors, factors)
end
```

Add these assertions after the existing tampered rewrite factor check:

```julia
    tampered_nested_factor = _replace_record_field(
        qq.certificate,
        :normality_certificate,
        _tamper_nested_first_factor(qq.certificate.normality_certificate, qq.R, length(qq.column)),
    )
    @test !Suslin.verify_ecp_induction_normality_certificate(tampered_nested_factor)

    bad_target = copy(qq.certificate.normality_certificate.conjugation_target)
    bad_target[1, 1] += one(qq.R)
    tampered_nested_target = _replace_record_field(
        qq.certificate,
        :normality_certificate,
        _replace_conjugated_certificate_field(
            qq.certificate.normality_certificate,
            :conjugation_target,
            bad_target,
        ),
    )
    @test !Suslin.verify_ecp_induction_normality_certificate(tampered_nested_target)
```

- [ ] **Step 4: Add public staged assertion**

In `test/expert/elementary_column_property.jl`, after
`@test staged.induction_normality.normality_rewrite.sl2_block != identity_matrix(R, 2)`, add:

```julia
    @test staged.induction_normality.normality_certificate isa Suslin.ConjugatedElementaryNormalityCertificate
    @test staged.induction_normality.normality_rewrite.normality_certificate ==
          staged.induction_normality.normality_certificate
    @test Suslin.verify_conjugate_elementary_certificate(staged.induction_normality.normality_certificate)
```

- [ ] **Step 5: Run tests to verify RED**

Run:

```bash
julia --project=. -e 'include("test/expert/ecp_induction_normality.jl")'
```

Expected: fail because `ECPInductionNormalityCertificate` has no `normality_certificate` field and `ecp_induction_normality_certificate` has no `normality_certificate` keyword.

### Task 2: Route ECP Normality Through Nested Certificates

**Files:**
- Modify: `src/algorithm/column_reduction.jl`

**Interfaces:**
- Consumes: `realize_conjugate_elementary_certificate`, `verify_conjugate_elementary_certificate`, and `ConjugatedElementaryNormalityCertificate` from `src/algorithm/normality.jl`.
- Produces: a new untyped `normality_certificate` field on `ECPInductionNormalityCertificate`, a `normality_certificate = nothing` keyword, and nested-certificate-aware replay summaries.

- [ ] **Step 1: Add the certificate field and keyword**

Change `ECPInductionNormalityCertificate` to store `normality_certificate` immediately after `normality_witness`. Add `normality_certificate = nothing` to `ecp_induction_normality_certificate`.

- [ ] **Step 2: Construct and store the nested certificate**

In `ecp_induction_normality_certificate`, pass `normality_certificate` into `_ecp_induction_normality_rewrite`. When constructing provisional and stored certificates, include `normality_rewrite.normality_certificate` in the new field.

- [ ] **Step 3: Validate legacy witnesses and supplied certificates in the rewrite helper**

Update `_ecp_induction_normality_rewrite` to accept `(normality_witness, normality_certificate, lower_column, lifted_lower_factors, R)`. Keep existing witness validation. Replace the raw call to `realize_conjugate_elementary` with:

```julia
    nested_certificate = _ecp_induction_normality_certificate_from_inputs(
        normality_certificate,
        conjugator,
        fixed_index,
        moving_index,
        entry,
        R,
        n,
    )
    rewrite_factors = nested_certificate.factors
    rewrite_product = nested_certificate.product
```

Add `_ecp_induction_normality_certificate_from_inputs` below `_ecp_normality_witness_keys_ok`. It should:

```julia
function _ecp_induction_normality_certificate_from_inputs(
    supplied_certificate,
    conjugator,
    fixed_index::Int,
    moving_index::Int,
    entry,
    R,
    n::Int,
)
    expected = realize_conjugate_elementary_certificate(conjugator, fixed_index, moving_index, entry)
    if supplied_certificate === nothing
        return expected
    end
    verify_conjugate_elementary_certificate(supplied_certificate) ||
        throw(ArgumentError("normality certificate does not verify"))
    _same_base_ring(supplied_certificate.ring, R) ||
        throw(ArgumentError("normality certificate ring must match the ECP ring"))
    supplied_certificate.n == n ||
        throw(ArgumentError("normality certificate dimension must match the ECP column"))
    supplied_certificate.A == expected.A &&
        supplied_certificate.i == expected.i &&
        supplied_certificate.j == expected.j &&
        supplied_certificate.a == expected.a &&
        supplied_certificate.conjugation_convention == expected.conjugation_convention &&
        supplied_certificate.conjugation_target == expected.conjugation_target &&
        supplied_certificate.factors == expected.factors &&
        supplied_certificate.product == expected.product &&
        supplied_certificate.verification == expected.verification ||
        throw(ArgumentError("normality certificate does not match the supplied witness data"))
    return supplied_certificate
end
```

The final implementation may compare additional #193 fields if useful, but it must compare at least the fields above.

- [ ] **Step 4: Store nested certificate in rewrite metadata**

The returned rewrite named tuple must include:

```julia
        normality_certificate = nested_certificate,
```

and preserve the existing keys `rewrite_factors`, `rewrite_product`,
`expected_rewrite_product`, `rewrite_product_ok`, and `fixed_lower_column_ok`.

- [ ] **Step 5: Replay nested certificates in the verifier**

In `_ecp_induction_normality_replay_summary`, pass `certificate.normality_certificate` into `_ecp_induction_normality_rewrite`. Add:

```julia
    normality_certificate_ok = normality_rewrite !== nothing &&
        certificate.normality_certificate == normality_rewrite.normality_certificate &&
        verify_conjugate_elementary_certificate(certificate.normality_certificate)
```

Include `normality_certificate_ok` in `overall_ok` and in the returned summary.

- [ ] **Step 6: Update staged verifier construction**

Update every `ECPInductionNormalityCertificate(...)` call to include the new field. Do not add a public export.

- [ ] **Step 7: Run focused GREEN tests**

Run:

```bash
julia --project=. -e 'include("test/expert/ecp_induction_normality.jl")'
```

Expected: pass.

### Task 3: Verify Staged ECP and Whole Package

**Files:**
- Modify only if Task 2 reveals a compile or verifier integration gap:
  `src/algorithm/column_reduction.jl`
  `test/expert/ecp_induction_normality.jl`
  `test/expert/elementary_column_property.jl`

**Interfaces:**
- Consumes: Task 1 tests and Task 2 implementation.
- Produces: verified ECP induction and public staged behavior.

- [ ] **Step 1: Run issue-required focused command**

Run:

```bash
julia --project=. -e 'include("test/expert/elementary_column_property.jl")'
```

Expected: pass.

- [ ] **Step 2: Run the required package test command**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: pass.

- [ ] **Step 3: Inspect the diff for scope**

Run:

```bash
git status --short
git diff --stat HEAD
git diff -- src/algorithm/column_reduction.jl test/expert/ecp_induction_normality.jl test/expert/elementary_column_property.jl
```

Expected: changes are limited to issue 194 docs, ECP induction normality code, and focused tests.

- [ ] **Step 4: Commit implementation**

Run:

```bash
git add src/algorithm/column_reduction.jl test/expert/ecp_induction_normality.jl test/expert/elementary_column_property.jl docs/superpowers/plans/2026-06-29-issue-194-ecp-induction-normality-replayed-certificates.md
git commit -m "Route ECP normality through replayed certificates"
```

Expected: commit succeeds.

## Plan Self-Review

- Spec coverage: Task 1 covers nested certificate acceptance and tamper rejection; Task 2 implements conversion, storage, and verifier replay; Task 3 runs the issue-required and package commands.
- Placeholder scan: no unresolved placeholder markers remain.
- Type consistency: `normality_certificate` is deliberately untyped in `ECPInductionNormalityCertificate` because `column_reduction.jl` is included before `normality.jl`; runtime calls already follow this pattern for normality helpers.
