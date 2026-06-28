# Issue 159 Lazy Laurent Row/Column Correction Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add row/left and column/right correction-side support to the internal lazy Laurent determinant hoist certificate.

**Architecture:** Extend `src/algorithm/laurent_gl_certificate.jl` with a side parser, right-side diagonal correction and rewrite helpers, and side-aware certificate fields and verification. Add a focused expert test that exercises both row and column certificates on the issue #38 lazy determinant fixture and covers the required negative controls.

**Tech Stack:** Julia, Oscar matrices over Laurent polynomial rings, existing Suslin Laurent determinant-deferred peel metadata, existing lazy Laurent hoist certificate internals, Test stdlib.

## Global Constraints

- Keep public exports unchanged.
- Accepted lazy hoist `correction_side` option names are exactly `:row` and `:column`.
- Internal metadata may record `:left` and `:right`.
- Row certificates verify `original_matrix == correction.factor * elementary_product`.
- Column certificates verify `original_matrix == elementary_product * correction.factor`.
- Row and column certificates on the same input report the same `overall_determinant`.
- Unsupported options such as `:diagonal` throw `ArgumentError` naming `:row` and `:column`.
- A certificate with row metadata but column reconstruction data must fail verification.
- Add focused test file `test/expert/laurent_lazy_row_column_correction.jl`.
- Register the focused expert test in `test/runtests.jl`.
- Focused verification command is `julia --project=. -e 'include("test/expert/laurent_lazy_row_column_correction.jl")'`.
- Package verification command is `julia --project=. -e 'using Pkg; Pkg.test()'`.
- Do not change `elementary_factorization(A)` to accept Laurent `GL_n` inputs.
- Do not commit `Manifest.toml`.

---

## File Structure

- Modify `src/algorithm/laurent_gl_certificate.jl`: add side parsing, side-aware lazy hoist certificate fields, column/right correction helpers, side-aware assembly, and stricter verification.
- Create `test/expert/laurent_lazy_row_column_correction.jl`: focused tests for row and column correction behavior, determinant equality, invalid option rejection, and row/column metadata mismatch rejection.
- Modify `test/runtests.jl`: register the new expert test near the existing lazy Laurent hoist tests.

---

### Task 1: Add Failing Row/Column Lazy Hoist Tests

**Files:**
- Create: `test/expert/laurent_lazy_row_column_correction.jl`
- Modify: `test/runtests.jl`

**Interfaces:**
- Consumes planned internal option `correction_side` on `Suslin._laurent_gl_lazy_deferred_correction_certificate`.
- Consumes planned certificate fields `correction_side` and `reconstruction_relation`.
- Produces a focused RED test that fails before implementation because the constructor does not accept `correction_side`.

- [ ] **Step 1: Create the focused failing test**

Create `test/expert/laurent_lazy_row_column_correction.jl` with tests that:

```julia
using Test
using Suslin
using Oscar

include("../fixtures/laurent_lazy_determinant_cases.jl")

function _issue159_fixture(id::AbstractString)
    catalog = LaurentLazyDeterminantCases.catalog()
    return only(filter(entry -> entry.id == id, catalog.cases))
end

function _issue159_product(factors, R, n::Int)
    product = identity_matrix(R, n)
    for factor in factors
        product *= factor
    end
    return product
end

function _issue159_rebuild(
    certificate;
    correction_side = certificate.correction_side,
    reconstruction_relation = certificate.reconstruction_relation,
    correction = certificate.correction,
    inverse_correction = certificate.inverse_correction,
    normalized_deferred_core = certificate.normalized_deferred_core,
    normalized_deferred_factorization = certificate.normalized_deferred_factorization,
    normalized_deferred_factors = certificate.normalized_deferred_factors,
    rewritten_left_factors = certificate.rewritten_left_factors,
    rewritten_right_factors = certificate.rewritten_right_factors,
    elementary_factors = certificate.elementary_factors,
    elementary_product = certificate.elementary_product,
    reconstructed_product = certificate.reconstructed_product,
    verification = certificate.verification,
)
    return Suslin.LaurentLazyGLHoistCertificate(
        certificate.original_matrix,
        certificate.deferred_metadata,
        certificate.overall_determinant,
        correction_side,
        reconstruction_relation,
        correction,
        inverse_correction,
        normalized_deferred_core,
        normalized_deferred_factorization,
        normalized_deferred_factors,
        rewritten_left_factors,
        rewritten_right_factors,
        elementary_factors,
        elementary_product,
        reconstructed_product,
        verification,
    )
end

@testset "lazy Laurent determinant correction side choices" begin
    entry = _issue159_fixture("issue-38-q-block-lazy-determinant")
    A = entry.inputs.matrix
    R = base_ring(A)
    n = nrows(A)

    row_certificate = Suslin._laurent_gl_lazy_deferred_correction_certificate(
        A;
        correction_side = :row,
    )
    column_certificate = Suslin._laurent_gl_lazy_deferred_correction_certificate(
        A;
        correction_side = :column,
    )

    @test row_certificate.correction_side == :row
    @test column_certificate.correction_side == :column
    @test row_certificate.correction.side == :left
    @test column_certificate.correction.side == :right
    @test row_certificate.reconstruction_relation == :left_correction_times_elementary_product
    @test column_certificate.reconstruction_relation == :elementary_product_times_right_correction
    @test row_certificate.overall_determinant == entry.determinant_profile.expected_determinant
    @test column_certificate.overall_determinant == entry.determinant_profile.expected_determinant
    @test row_certificate.overall_determinant == column_certificate.overall_determinant

    @test row_certificate.elementary_product ==
        _issue159_product(row_certificate.elementary_factors, R, n)
    @test column_certificate.elementary_product ==
        _issue159_product(column_certificate.elementary_factors, R, n)
    @test row_certificate.correction.factor * row_certificate.elementary_product == A
    @test column_certificate.elementary_product * column_certificate.correction.factor == A
    @test row_certificate.reconstructed_product == A
    @test column_certificate.reconstructed_product == A
    @test row_certificate.verification.overall_ok
    @test column_certificate.verification.overall_ok
    @test Suslin._verify_laurent_gl_lazy_deferred_correction_certificate(row_certificate)
    @test Suslin._verify_laurent_gl_lazy_deferred_correction_certificate(column_certificate)

    @test det(row_certificate.correction.factor) == row_certificate.overall_determinant
    @test det(column_certificate.correction.factor) == column_certificate.overall_determinant
    @test row_certificate.correction.factor * row_certificate.inverse_correction ==
        identity_matrix(R, n)
    @test column_certificate.correction.factor * column_certificate.inverse_correction ==
        identity_matrix(R, n)

    wrong_metadata = _issue159_rebuild(
        column_certificate;
        correction_side = :row,
        reconstruction_relation = :left_correction_times_elementary_product,
    )
    wrong_verification = Suslin._laurent_gl_lazy_deferred_correction_certificate_verification(
        wrong_metadata,
    )
    @test !wrong_verification.overall_ok
    @test !Suslin._verify_laurent_gl_lazy_deferred_correction_certificate(wrong_metadata)

    err = try
        Suslin._laurent_gl_lazy_deferred_correction_certificate(
            A;
            correction_side = :diagonal,
        )
        nothing
    catch caught
        caught
    end
    @test err isa ArgumentError
    @test occursin(":row", sprint(showerror, err))
    @test occursin(":column", sprint(showerror, err))
end
```

- [ ] **Step 2: Register the expert test**

In `test/runtests.jl`, add:

```julia
"expert/laurent_lazy_row_column_correction.jl",
```

immediately after `"expert/laurent_lazy_correction_hoist.jl",`.

- [ ] **Step 3: Run focused test to verify RED**

Run:

```bash
julia --project=. -e 'include("test/expert/laurent_lazy_row_column_correction.jl")'
```

Expected: FAIL with a method error or keyword error showing `correction_side` is not supported yet.

- [ ] **Step 4: Commit the failing tests**

Run:

```bash
git add test/expert/laurent_lazy_row_column_correction.jl test/runtests.jl
git commit -m "test: cover lazy laurent correction side choices"
```

---

### Task 2: Implement Side-Aware Lazy Laurent Hoist Certificate

**Files:**
- Modify: `src/algorithm/laurent_gl_certificate.jl`

**Interfaces:**
- Consumes existing #158 helpers and `metadata.deferred_correction`.
- Produces `correction_side`, `reconstruction_relation`, right-side correction helpers, side-aware factor assembly, and verification.

- [ ] **Step 1: Add option parsing and metadata helpers**

Add helpers that normalize `:row` to `:left` and `:column` to `:right`, expose the external side back as `:row` or `:column`, and throw `ArgumentError("correction_side must be :row or :column")` for unsupported values.

- [ ] **Step 2: Extend `LaurentLazyGLHoistCertificate`**

Insert fields after `overall_determinant`:

```julia
correction_side::Symbol
reconstruction_relation::Symbol
```

Insert a new `rewritten_right_factors::Vector` field after `rewritten_left_factors::Vector`.

Update every constructor call in `src/algorithm/laurent_gl_certificate.jl` and existing tests to pass the new fields.

- [ ] **Step 3: Add right-side correction and rewrite helpers**

Add a right diagonal correction that uses the same determinant and inverse as the left correction but records `side = :right` and `kind = :right_diagonal_determinant_correction`. Add:

```julia
_rewrite_right_elementary_factor_across_diagonal(factor, diagonal_factor)
```

with coefficient `diagonal_entries[row] * coefficient * inv(diagonal_entries[col])`, and a vector helper mirroring the left helper.

- [ ] **Step 4: Build normalized deferred cores by side**

For `:left`, use the existing `correction.inverse_factor * deferred` core. For `:right`, use `deferred * correction.inverse_factor`. For identity determinants, use a side-aware identity correction and keep the core equal to the deferred submatrix.

- [ ] **Step 5: Build side-aware original-dimension correction and elementary factors**

For `:left`, preserve the existing #158 factor sequence and reconstruction relation. For `:right`, build:

```text
left_product^-1 *
blockdiag(normalized_deferred_factors, I) *
rewrite(right_product^-1 across correction.factor)
```

and reconstruct with:

```text
elementary_product * correction.factor
```

- [ ] **Step 6: Update verification**

Recompute the expected side-aware correction, normalized deferred core, rewritten factors, factor sequence, elementary product, and reconstruction relation. Fail verification if any stored side metadata or relation does not match the recomputed side.

- [ ] **Step 7: Run focused test to verify GREEN**

Run:

```bash
julia --project=. -e 'include("test/expert/laurent_lazy_row_column_correction.jl")'
```

Expected: PASS.

- [ ] **Step 8: Commit implementation**

Run:

```bash
git add src/algorithm/laurent_gl_certificate.jl test/expert/laurent_lazy_correction_hoist.jl
git commit -m "feat: support lazy laurent correction side choices"
```

---

### Task 3: Regression Verification

**Files:**
- No planned source edits unless verification exposes a defect.

**Interfaces:**
- Consumes the new focused test and existing #158 test.
- Produces verified confidence that the new side-aware path preserves existing row/left behavior.

- [ ] **Step 1: Run the required issue verification**

Run:

```bash
julia --project=. -e 'include("test/expert/laurent_lazy_row_column_correction.jl")'
```

Expected: PASS.

- [ ] **Step 2: Run existing lazy hoist regression**

Run:

```bash
julia --project=. -e 'include("test/expert/laurent_lazy_correction_hoist.jl")'
```

Expected: PASS.

- [ ] **Step 3: Run package test gate**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: PASS.

- [ ] **Step 4: Confirm no forbidden artifacts are staged**

Run:

```bash
git status --short
git status --ignored --short -- Manifest.toml
```

Expected: no tracked `Manifest.toml`; only intended source, test, and docs changes are tracked or committed.
