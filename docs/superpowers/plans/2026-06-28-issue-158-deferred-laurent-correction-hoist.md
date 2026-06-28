# Issue 158 Deferred Laurent Correction Hoist Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Hoist supported deferred Laurent diagonal corrections into original matrix-level lazy `GL_n` certificate data with exact elementary coefficient rewriting.

**Architecture:** Add internal row/left hoist helpers and a focused `LaurentLazyGLHoistCertificate` in `src/algorithm/laurent_gl_certificate.jl`. The certificate consumes #157 deferred metadata, embeds the diagonal correction, rewrites inverse left peel factors using `E_ij(a) * D == D * E_ij(d_i^-1 * a * d_j)`, appends embedded normalized-core factors and inverse right factors, and verifies exact reconstruction.

**Tech Stack:** Julia, Oscar matrices over Laurent polynomial rings, existing Suslin Laurent column-peel and GL normalization helpers, Test stdlib.

## Global Constraints

- Keep public exports unchanged.
- Implement the internal row/left correction side only.
- Do not expose user-facing row/column options.
- Consume enriched #157 lazy determinant metadata.
- The certificate must store the original matrix, `overall_determinant`, an embedded diagonal correction and inverse correction, rewritten elementary factors, and exact reconstruction status.
- Moving a deferred diagonal correction across elementary factors must rewrite coefficients exactly with `E_ij(a) * D == D * E_ij(d_i^-1 * a * d_j)`.
- The focused expert test file is `test/expert/laurent_lazy_correction_hoist.jl`.
- Register the focused expert test in `test/runtests.jl`.
- Focused verification command is `julia --project=. -e 'include("test/expert/laurent_lazy_correction_hoist.jl")'`.
- Package verification command is `julia --project=. -e 'using Pkg; Pkg.test()'`.
- Do not commit `Manifest.toml`.

---

## File Structure

- Modify `src/algorithm/laurent_gl_certificate.jl`: add the internal hoist certificate struct, diagonal extraction and coefficient-rewrite helpers, embedded correction helper, constructor, and verifier.
- Create `test/expert/laurent_lazy_correction_hoist.jl`: focused tests for the elementary rewrite helper, supported monomial-unit hoist, exact reconstruction, determinant propagation, and wrong-unrewritten negative control.
- Modify `test/runtests.jl`: register the focused expert test near the existing lazy Laurent tests.

---

### Task 1: Add Failing Hoist Tests

**Files:**
- Create: `test/expert/laurent_lazy_correction_hoist.jl`
- Modify: `test/runtests.jl`

**Interfaces:**
- Consumes planned internal functions `Suslin._rewrite_left_elementary_factor_across_diagonal`, `Suslin._laurent_gl_lazy_deferred_correction_certificate`, and `Suslin._verify_laurent_gl_lazy_deferred_correction_certificate`.
- Produces a focused RED test that fails before implementation and proves that an unrewritten hoist is rejected.

- [ ] **Step 1: Create the focused failing test**

Create `test/expert/laurent_lazy_correction_hoist.jl` with:

```julia
using Test
using Suslin
using Oscar

include("../fixtures/laurent_lazy_determinant_cases.jl")

function _issue158_fixture(id::AbstractString)
    catalog = LaurentLazyDeterminantCases.catalog()
    return only(filter(entry -> entry.id == id, catalog.cases))
end

function _issue158_product(factors, R, n::Int)
    product = identity_matrix(R, n)
    for factor in factors
        product *= factor
    end
    return product
end

function _issue158_wrong_unrewritten_certificate(certificate)
    R = base_ring(certificate.original_matrix)
    n = nrows(certificate.original_matrix)
    deferred_n = nrows(certificate.normalized_deferred_core)
    peel_certificate = certificate.deferred_metadata.peel_certificate
    unrewritten_left = Suslin._inverse_elementary_sequence(peel_certificate.left_factors)
    embedded_core_factors = Suslin._embed_laurent_deferred_peel_factors(
        certificate.normalized_deferred_factors,
        R,
        n,
        deferred_n,
    )
    right_inverse = Suslin._inverse_elementary_sequence(peel_certificate.right_factors)
    wrong_factors = vcat(unrewritten_left, embedded_core_factors, right_inverse)
    wrong_product = _issue158_product(wrong_factors, R, n)

    return Suslin.LaurentLazyGLHoistCertificate(
        certificate.original_matrix,
        certificate.deferred_metadata,
        certificate.overall_determinant,
        certificate.correction,
        certificate.inverse_correction,
        certificate.normalized_deferred_core,
        certificate.normalized_deferred_factorization,
        certificate.normalized_deferred_factors,
        unrewritten_left,
        wrong_factors,
        wrong_product,
        certificate.correction.factor * wrong_product,
        certificate.verification,
    )
end

@testset "left elementary factors rewrite across Laurent diagonal correction" begin
    R, (x, y) = suslin_laurent_polynomial_ring(GF(2), ["x", "y"])
    diagonal = diagonal_matrix(R, [x, y, one(R)])
    coefficient = x + y
    factor = elementary_matrix(3, 2, 1, coefficient, R)

    rewritten = Suslin._rewrite_left_elementary_factor_across_diagonal(factor, diagonal)
    row, col, rewritten_coefficient = Suslin._elementary_factor_data(rewritten)

    @test (row, col) == (2, 1)
    @test rewritten_coefficient == y^-1 * coefficient * x
    @test factor * diagonal == diagonal * rewritten
end

@testset "lazy deferred Laurent correction hoists to original GL certificate" begin
    entry = _issue158_fixture("monomial-unit-row-column-cores")
    A = entry.inputs.matrix
    metadata = Suslin._factor_laurent_gl_lazy_determinant_peel(A)

    certificate = Suslin._laurent_gl_lazy_deferred_correction_certificate(metadata)

    @test certificate.original_matrix == A
    @test certificate.deferred_metadata == metadata
    @test certificate.overall_determinant == metadata.overall_determinant
    @test certificate.overall_determinant == entry.determinant_profile.expected_determinant
    @test certificate.correction.scope == :original_matrix
    @test certificate.correction.side == :left
    @test certificate.correction.factor * certificate.inverse_correction ==
        identity_matrix(base_ring(A), nrows(A))
    @test certificate.inverse_correction * certificate.correction.factor ==
        identity_matrix(base_ring(A), nrows(A))
    @test det(certificate.correction.factor) == certificate.overall_determinant
    @test certificate.normalized_deferred_core == metadata.normalized_deferred_core
    @test Suslin._verify_laurent_column_peel_replay(certificate.normalized_deferred_factorization)
    @test verify_factorization(
        certificate.normalized_deferred_core,
        certificate.normalized_deferred_factors,
    )

    expected_product = _issue158_product(
        certificate.elementary_factors,
        base_ring(A),
        nrows(A),
    )
    @test certificate.elementary_product == expected_product
    @test certificate.correction.factor * certificate.elementary_product == A
    @test certificate.reconstructed_product == A
    @test certificate.verification.overall_ok
    @test Suslin._verify_laurent_gl_lazy_deferred_correction_certificate(certificate)

    wrong = _issue158_wrong_unrewritten_certificate(certificate)
    @test !Suslin._verify_laurent_gl_lazy_deferred_correction_certificate(wrong)
    @test !Suslin._laurent_gl_lazy_deferred_correction_certificate_verification(wrong).rewritten_left_factors_ok
end
```

- [ ] **Step 2: Register the expert test**

In `test/runtests.jl`, add:

```julia
"expert/laurent_lazy_correction_hoist.jl",
```

immediately after `"expert/laurent_lazy_submatrix_normalization.jl",`.

- [ ] **Step 3: Run focused test to verify RED**

Run:

```bash
julia --project=. -e 'include("test/expert/laurent_lazy_correction_hoist.jl")'
```

Expected: FAIL with `UndefVarError` for `_rewrite_left_elementary_factor_across_diagonal` or `_laurent_gl_lazy_deferred_correction_certificate`.

- [ ] **Step 4: Commit the failing tests**

Run:

```bash
git add test/expert/laurent_lazy_correction_hoist.jl test/runtests.jl
git commit -m "test: cover deferred laurent correction hoist"
```

---

### Task 2: Implement Hoisted Lazy Laurent GL Certificate

**Files:**
- Modify: `src/algorithm/laurent_gl_certificate.jl`

**Interfaces:**
- Consumes: #157 metadata, `_verify_laurent_determinant_deferred_submatrix_normalization`, `_inverse_elementary_sequence`, `_embed_laurent_deferred_peel_factors`, `_factor_laurent_sl_column_peel`, `_verify_laurent_column_peel_replay`, `_factor_product`, and `_elementary_factor_data`.
- Produces: `LaurentLazyGLHoistCertificate`, `_rewrite_left_elementary_factor_across_diagonal`, `_laurent_gl_lazy_deferred_correction_certificate`, `_laurent_gl_lazy_deferred_correction_certificate_verification`, and `_verify_laurent_gl_lazy_deferred_correction_certificate`.

- [ ] **Step 1: Add the internal certificate struct**

At the top of `src/algorithm/laurent_gl_certificate.jl`, after `LaurentGLFactorizationCertificate`, add:

```julia
struct LaurentLazyGLHoistCertificate
    original_matrix
    deferred_metadata
    overall_determinant
    correction
    inverse_correction
    normalized_deferred_core
    normalized_deferred_factorization
    normalized_deferred_factors::Vector
    rewritten_left_factors::Vector
    elementary_factors::Vector
    elementary_product
    reconstructed_product
    verification
end
```

- [ ] **Step 2: Add diagonal and rewrite helpers**

After `_laurent_gl_factorization_certificate_from_normalization`, add:

```julia
function _laurent_diagonal_entries(diagonal_factor)
    n = _require_square_matrix(diagonal_factor, "diagonal factor")
    R = base_ring(diagonal_factor)
    entries = [diagonal_factor[i, i] for i in 1:n]
    for i in 1:n, j in 1:n
        i == j && continue
        diagonal_factor[i, j] == zero(R) ||
            throw(ArgumentError("diagonal factor must have zero off-diagonal entries"))
    end
    return entries
end

function _rewrite_left_elementary_factor_across_diagonal(factor, diagonal_factor)
    n = _require_square_matrix(diagonal_factor, "diagonal factor")
    nrows(factor) == n || throw(ArgumentError("elementary factor and diagonal factor sizes must match"))
    ncols(factor) == n || throw(ArgumentError("elementary factor and diagonal factor sizes must match"))
    _same_base_ring(base_ring(factor), base_ring(diagonal_factor)) ||
        throw(ArgumentError("elementary factor and diagonal factor must have the same base ring"))

    diagonal_entries = _laurent_diagonal_entries(diagonal_factor)
    row, col, coefficient = _elementary_factor_data(factor)
    rewritten_coefficient = inv(diagonal_entries[row]) * coefficient * diagonal_entries[col]
    return elementary_matrix(n, row, col, rewritten_coefficient, base_ring(factor))
end

function _rewrite_left_elementary_factors_across_diagonal(factors, diagonal_factor)
    R = base_ring(diagonal_factor)
    n = _require_square_matrix(diagonal_factor, "diagonal factor")
    rewritten = typeof(identity_matrix(R, n))[]
    for factor in factors
        push!(rewritten, _rewrite_left_elementary_factor_across_diagonal(factor, diagonal_factor))
    end
    return rewritten
end
```

- [ ] **Step 3: Add embedded correction and assembly helpers**

Add:

```julia
function _embed_laurent_deferred_correction(correction, original_dimension::Int)
    deferred_dimension = nrows(correction.factor)
    correction_side = correction.side
    correction_side == :left ||
        throw(ArgumentError("only left deferred Laurent corrections can be hoisted"))

    embedded_factor = block_embedding(
        correction.factor,
        original_dimension,
        collect(1:deferred_dimension),
    )
    embedded_inverse = block_embedding(
        correction.inverse_factor,
        original_dimension,
        collect(1:deferred_dimension),
    )
    return merge(correction, (;
        scope = :original_matrix,
        deferred_scope = correction.scope,
        factor = embedded_factor,
        inverse_factor = embedded_inverse,
        embedded_from_dimension = deferred_dimension,
    ))
end

function _laurent_lazy_hoist_elementary_factors(metadata, correction, normalized_deferred_factors)
    peel_certificate = metadata.peel_certificate
    R = base_ring(peel_certificate.original_matrix)
    original_dimension = nrows(peel_certificate.original_matrix)
    deferred_dimension = nrows(metadata.normalized_deferred_core)
    left_inverse = _inverse_elementary_sequence(peel_certificate.left_factors)
    rewritten_left = _rewrite_left_elementary_factors_across_diagonal(
        left_inverse,
        correction.factor,
    )
    embedded_core = _embed_laurent_deferred_peel_factors(
        normalized_deferred_factors,
        R,
        original_dimension,
        deferred_dimension,
    )
    right_inverse = _inverse_elementary_sequence(peel_certificate.right_factors)
    return (;
        rewritten_left,
        elementary_factors = vcat(rewritten_left, embedded_core, right_inverse),
    )
end
```

- [ ] **Step 4: Add the constructor**

Add:

```julia
function _laurent_gl_lazy_deferred_correction_certificate(
    metadata;
    progress_callback = nothing,
)
    _verify_laurent_determinant_deferred_submatrix_normalization(metadata) ||
        error("invalid deferred Laurent determinant normalization metadata")
    metadata.supported ||
        throw(ArgumentError("unsupported deferred Laurent determinant cannot be hoisted"))
    metadata.deferred_correction !== nothing ||
        throw(ArgumentError("supported deferred Laurent determinant metadata is missing a correction"))

    A = metadata.peel_certificate.original_matrix
    R = base_ring(A)
    n = nrows(A)
    correction = _embed_laurent_deferred_correction(metadata.deferred_correction, n)
    normalized_deferred_core = metadata.normalized_deferred_core
    normalized_deferred_factorization = _factor_laurent_sl_column_peel(
        normalized_deferred_core;
        progress_callback,
    )
    normalized_deferred_factors = normalized_deferred_factorization.factors
    assembly = _laurent_lazy_hoist_elementary_factors(
        metadata,
        correction,
        normalized_deferred_factors,
    )
    elementary_product = _factor_product(assembly.elementary_factors, R, n)
    reconstructed_product = correction.factor * elementary_product
    certificate = LaurentLazyGLHoistCertificate(
        A,
        metadata,
        metadata.overall_determinant,
        correction,
        correction.inverse_factor,
        normalized_deferred_core,
        normalized_deferred_factorization,
        normalized_deferred_factors,
        assembly.rewritten_left,
        assembly.elementary_factors,
        elementary_product,
        reconstructed_product,
        nothing,
    )
    verification = _laurent_gl_lazy_deferred_correction_certificate_verification(certificate)
    verification.overall_ok ||
        error("internal lazy Laurent GL correction hoist verification failed")
    return LaurentLazyGLHoistCertificate(
        certificate.original_matrix,
        certificate.deferred_metadata,
        certificate.overall_determinant,
        certificate.correction,
        certificate.inverse_correction,
        certificate.normalized_deferred_core,
        certificate.normalized_deferred_factorization,
        certificate.normalized_deferred_factors,
        certificate.rewritten_left_factors,
        certificate.elementary_factors,
        certificate.elementary_product,
        certificate.reconstructed_product,
        verification,
    )
end

function _laurent_gl_lazy_deferred_correction_certificate(
    A;
    determinant_probe = classify_laurent_determinant,
    progress_callback = nothing,
)
    deferred_certificate = _laurent_determinant_deferred_peel_certificate(
        A;
        progress_callback,
    )
    metadata = _normalize_laurent_determinant_deferred_submatrix(
        deferred_certificate;
        determinant_probe,
    )
    return _laurent_gl_lazy_deferred_correction_certificate(
        metadata;
        progress_callback,
    )
end
```

- [ ] **Step 5: Add the verifier**

Add:

```julia
function _verify_laurent_gl_lazy_deferred_correction_certificate(certificate)::Bool
    return _laurent_gl_lazy_deferred_correction_certificate_verification(certificate).overall_ok
end

function _laurent_gl_lazy_deferred_correction_certificate_verification(certificate)
    size_ok = try
        nrows(certificate.original_matrix) == ncols(certificate.original_matrix)
    catch err
        err isa InterruptException && rethrow()
        false
    end

    metadata_ok = false
    correction_ok = false
    normalized_core_ok = false
    normalized_factorization_ok = false
    rewritten_left_factors_ok = false
    elementary_factors_ok = false
    elementary_product_ok = false
    reconstructed_product_ok = false

    if size_ok
        try
            A = certificate.original_matrix
            R = base_ring(A)
            n = nrows(A)
            metadata = certificate.deferred_metadata
            metadata_ok =
                metadata.peel_certificate.original_matrix == A &&
                metadata.supported &&
                metadata.determinant_source == :deferred_submatrix &&
                _verify_laurent_determinant_deferred_submatrix_normalization(metadata)

            expected_correction = _embed_laurent_deferred_correction(
                metadata.deferred_correction,
                n,
            )
            identity = identity_matrix(R, n)
            correction_ok =
                certificate.overall_determinant == metadata.overall_determinant &&
                certificate.correction == expected_correction &&
                certificate.inverse_correction == expected_correction.inverse_factor &&
                certificate.correction.factor * certificate.inverse_correction == identity &&
                certificate.inverse_correction * certificate.correction.factor == identity &&
                det(certificate.correction.factor) == certificate.overall_determinant

            normalized_core_ok =
                certificate.normalized_deferred_core == metadata.normalized_deferred_core &&
                det(certificate.normalized_deferred_core) == one(base_ring(certificate.normalized_deferred_core))

            normalized_factorization_ok =
                certificate.normalized_deferred_factorization.original_matrix ==
                    certificate.normalized_deferred_core &&
                _verify_laurent_column_peel_replay(certificate.normalized_deferred_factorization) &&
                _factor_sequences_equal(
                    certificate.normalized_deferred_factors,
                    certificate.normalized_deferred_factorization.factors,
                ) &&
                verify_factorization(
                    certificate.normalized_deferred_core,
                    certificate.normalized_deferred_factors,
                )

            expected_left_inverse = _inverse_elementary_sequence(
                metadata.peel_certificate.left_factors,
            )
            expected_rewritten_left = _rewrite_left_elementary_factors_across_diagonal(
                expected_left_inverse,
                certificate.correction.factor,
            )
            rewritten_left_factors_ok = _factor_sequences_equal(
                certificate.rewritten_left_factors,
                expected_rewritten_left,
            )

            assembly = _laurent_lazy_hoist_elementary_factors(
                metadata,
                certificate.correction,
                certificate.normalized_deferred_factors,
            )
            elementary_factors_ok = _factor_sequences_equal(
                certificate.elementary_factors,
                assembly.elementary_factors,
            )
            rebuilt_product = _factor_product(certificate.elementary_factors, R, n)
            elementary_product_ok = certificate.elementary_product == rebuilt_product
            reconstructed_product_ok =
                certificate.reconstructed_product ==
                    certificate.correction.factor * certificate.elementary_product &&
                certificate.reconstructed_product == A
        catch err
            err isa InterruptException && rethrow()
        end
    end

    overall_ok = size_ok && metadata_ok && correction_ok && normalized_core_ok &&
        normalized_factorization_ok && rewritten_left_factors_ok &&
        elementary_factors_ok && elementary_product_ok && reconstructed_product_ok

    return (;
        overall_ok,
        size_ok,
        metadata_ok,
        correction_ok,
        normalized_core_ok,
        normalized_factorization_ok,
        rewritten_left_factors_ok,
        elementary_factors_ok,
        elementary_product_ok,
        reconstructed_product_ok,
    )
end
```

- [ ] **Step 6: Run focused test to verify GREEN**

Run:

```bash
julia --project=. -e 'include("test/expert/laurent_lazy_correction_hoist.jl")'
```

Expected: PASS.

- [ ] **Step 7: Commit implementation**

Run:

```bash
git add src/algorithm/laurent_gl_certificate.jl
git commit -m "feat: hoist deferred laurent correction certificates"
```

---

### Task 3: Verification And Regression Gate

**Files:**
- No source changes expected.

**Interfaces:**
- Consumes the test and implementation commits from Tasks 1 and 2.
- Produces verification evidence required before PR creation.

- [ ] **Step 1: Run the issue verification command**

Run:

```bash
julia --project=. -e 'include("test/expert/laurent_lazy_correction_hoist.jl")'
```

Expected: PASS.

- [ ] **Step 2: Run adjacent lazy metadata regression**

Run:

```bash
julia --project=. -e 'include("test/expert/laurent_lazy_submatrix_normalization.jl")'
```

Expected: PASS.

- [ ] **Step 3: Run package default tests**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: PASS.

- [ ] **Step 4: Check the final diff**

Run:

```bash
git status --short
git diff --stat origin/main...HEAD
```

Expected: only the issue #158 design, plan, test, registration, and implementation files are changed; `Manifest.toml` is not tracked.
