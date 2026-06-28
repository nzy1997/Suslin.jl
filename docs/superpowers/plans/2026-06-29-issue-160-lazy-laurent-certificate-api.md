# Issue 160 Lazy Laurent Certificate API Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Expose lazy determinant fields through the public Laurent GL certificate API while preserving the eager one-argument call and the Laurent GL boundary in `elementary_factorization`.

**Architecture:** Keep `laurent_gl_factorization_certificate(A)` on the eager normalization path. Add keyword routing for `determinant_strategy = :lazy`, expose the existing lazy hoist certificate as a public certificate type with a direct `determinant_source` field, and dispatch public verification to the lazy verifier for lazy certificates.

**Tech Stack:** Julia, Oscar matrices over Laurent polynomial rings, existing Suslin lazy Laurent determinant-deferred peel and hoist internals, Test stdlib.

## Global Constraints

- Keep `laurent_gl_factorization_certificate(A)` backward compatible and eager by default.
- Accepted `determinant_strategy` values are exactly `:eager` and `:lazy`.
- The public lazy route is `laurent_gl_factorization_certificate(A; determinant_strategy = :lazy, correction_side = :row)` with `correction_side = :row` or `:column`.
- `correction_side` without `determinant_strategy = :lazy` must throw `ArgumentError` instead of being ignored.
- Lazy certificates must expose direct public fields `overall_determinant`, `determinant_source`, and `correction_side`.
- Lazy route certificates must verify exact original-matrix reconstruction via `verify_laurent_gl_factorization_certificate(cert) == true`.
- On the issue #38 fixture, lazy row reports `overall_determinant == det(Q)`, `determinant_source == :deferred_submatrix`, and `correction_side == :row`; lazy column verifies and reports `correction_side == :column`.
- `elementary_factorization(Q)` for the original issue #38 Laurent `GL_n` input must still throw the staged Laurent `GL_n` boundary.
- Add focused public test file `test/public/laurent_gl_certificate_options.jl`.
- Register the focused public test in `test/runtests.jl`.
- Update `test/public/api_surface.jl` to cover the exported lazy certificate type, the existing one-argument call, and accepted/rejected keyword spellings.
- Focused verification command is `julia --project=. -e 'include("test/public/laurent_gl_certificate_options.jl")'`.
- Package verification command is `julia --project=. -e 'using Pkg; Pkg.test()'`.
- Do not update ToricBuilder Q-block reporting.
- Do not commit `Manifest.toml`.

---

## File Structure

- Modify `src/Suslin.jl`: export `LaurentLazyGLHoistCertificate`.
- Modify `src/algorithm/laurent_gl_certificate.jl`: add `determinant_source` to `LaurentLazyGLHoistCertificate`, add public keyword routing, and add lazy certificate verification dispatch.
- Create `test/public/laurent_gl_certificate_options.jl`: focused public lazy certificate contract and negative control tests.
- Modify `test/public/api_surface.jl`: export and keyword contract smoke tests.
- Modify `test/runtests.jl`: register the new public test.
- Modify `test/expert/laurent_lazy_correction_hoist.jl` and `test/expert/laurent_lazy_row_column_correction.jl`: update test helper constructors for the new lazy certificate field.

---

### Task 1: Add Failing Public API Tests

**Files:**
- Create: `test/public/laurent_gl_certificate_options.jl`
- Modify: `test/public/api_surface.jl`
- Modify: `test/runtests.jl`

**Interfaces:**
- Consumes planned public keyword route `laurent_gl_factorization_certificate(A; determinant_strategy = :lazy, correction_side = :row)`.
- Consumes planned exported type `LaurentLazyGLHoistCertificate`.
- Consumes planned direct lazy field `determinant_source`.
- Produces RED tests that fail before implementation because the public function currently has no keyword method and the lazy certificate type is not exported.

- [ ] **Step 1: Create the focused public test**

Create `test/public/laurent_gl_certificate_options.jl`:

```julia
using Test
using Suslin
using Oscar

include("../fixtures/toricbuilder_issue38_cases.jl")

function _issue160_caught_error(f)
    try
        f()
        return nothing
    catch err
        return err
    end
end

@testset "public Laurent GL certificate options" begin
    entry = only(ToricBuilderIssue38Cases.catalog().cases)
    Q = entry.inputs.matrix
    expected_determinant = det(Q)

    eager_certificate = laurent_gl_factorization_certificate(Q)
    @test eager_certificate isa LaurentGLFactorizationCertificate
    @test eager_certificate.reconstructed_product == Q
    @test verify_laurent_gl_factorization_certificate(eager_certificate)

    explicit_eager = laurent_gl_factorization_certificate(
        Q;
        determinant_strategy = :eager,
    )
    @test explicit_eager isa LaurentGLFactorizationCertificate
    @test explicit_eager.reconstructed_product == Q
    @test verify_laurent_gl_factorization_certificate(explicit_eager)

    row_certificate = laurent_gl_factorization_certificate(
        Q;
        determinant_strategy = :lazy,
        correction_side = :row,
    )
    @test row_certificate isa LaurentLazyGLHoistCertificate
    @test row_certificate.overall_determinant == expected_determinant
    @test row_certificate.determinant_source == :deferred_submatrix
    @test row_certificate.correction_side == :row
    @test row_certificate.reconstructed_product == Q
    @test row_certificate.verification.overall_ok
    @test verify_laurent_gl_factorization_certificate(row_certificate)

    default_row_certificate = laurent_gl_factorization_certificate(
        Q;
        determinant_strategy = :lazy,
    )
    @test default_row_certificate.correction_side == :row
    @test default_row_certificate.determinant_source == :deferred_submatrix
    @test verify_laurent_gl_factorization_certificate(default_row_certificate)

    column_certificate = laurent_gl_factorization_certificate(
        Q;
        determinant_strategy = :lazy,
        correction_side = :column,
    )
    @test column_certificate isa LaurentLazyGLHoistCertificate
    @test column_certificate.overall_determinant == expected_determinant
    @test column_certificate.overall_determinant == row_certificate.overall_determinant
    @test column_certificate.determinant_source == :deferred_submatrix
    @test column_certificate.correction_side == :column
    @test column_certificate.reconstructed_product == Q
    @test column_certificate.verification.overall_ok
    @test verify_laurent_gl_factorization_certificate(column_certificate)

    strategy_err = _issue160_caught_error(() ->
        laurent_gl_factorization_certificate(Q; determinant_strategy = :deferred))
    @test strategy_err isa ArgumentError
    @test occursin(":eager", sprint(showerror, strategy_err))
    @test occursin(":lazy", sprint(showerror, strategy_err))

    misplaced_side_err = _issue160_caught_error(() ->
        laurent_gl_factorization_certificate(Q; correction_side = :column))
    @test misplaced_side_err isa ArgumentError
    @test occursin("determinant_strategy = :lazy", sprint(showerror, misplaced_side_err))

    invalid_side_err = _issue160_caught_error(() ->
        laurent_gl_factorization_certificate(
            Q;
            determinant_strategy = :lazy,
            correction_side = :diagonal,
        ))
    @test invalid_side_err isa ArgumentError
    @test occursin(":row", sprint(showerror, invalid_side_err))
    @test occursin(":column", sprint(showerror, invalid_side_err))

    original_err = _issue160_caught_error(() -> elementary_factorization(Q))
    @test original_err isa ArgumentError
    @test occursin("Laurent GL_n normalization boundary", sprint(showerror, original_err))
end
```

- [ ] **Step 2: Update API surface tests**

In `test/public/api_surface.jl`, add `LaurentLazyGLHoistCertificate` next to the existing Laurent certificate type checks and add this keyword smoke test near the existing polynomial factorization smoke test:

```julia
    @test isdefined(Suslin, :LaurentLazyGLHoistCertificate)
    @test Suslin.LaurentLazyGLHoistCertificate === LaurentLazyGLHoistCertificate

    LR, (u,) = suslin_laurent_polynomial_ring(QQ, ["u"])
    L = matrix(LR, [
        one(LR) u       zero(LR);
        zero(LR) one(LR) zero(LR);
        zero(LR) zero(LR) one(LR)
    ])

    one_arg_certificate = laurent_gl_factorization_certificate(L)
    @test one_arg_certificate isa LaurentGLFactorizationCertificate
    @test verify_laurent_gl_factorization_certificate(one_arg_certificate)

    lazy_keyword_certificate = laurent_gl_factorization_certificate(
        L;
        determinant_strategy = :lazy,
        correction_side = :row,
    )
    @test lazy_keyword_certificate isa LaurentLazyGLHoistCertificate
    @test lazy_keyword_certificate.determinant_source == :deferred_submatrix
    @test lazy_keyword_certificate.correction_side == :row
    @test verify_laurent_gl_factorization_certificate(lazy_keyword_certificate)

    rejected_strategy = try
        laurent_gl_factorization_certificate(L; determinant_strategy = :unsupported)
        nothing
    catch err
        err
    end
    @test rejected_strategy isa ArgumentError
    @test occursin(":eager", sprint(showerror, rejected_strategy))
    @test occursin(":lazy", sprint(showerror, rejected_strategy))

    rejected_side_without_lazy = try
        laurent_gl_factorization_certificate(L; correction_side = :column)
        nothing
    catch err
        err
    end
    @test rejected_side_without_lazy isa ArgumentError
    @test occursin("determinant_strategy = :lazy", sprint(showerror, rejected_side_without_lazy))
```

- [ ] **Step 3: Register the focused public test**

In `test/runtests.jl`, add:

```julia
"public/laurent_gl_certificate_options.jl",
```

immediately after `"public/api_surface.jl",`.

- [ ] **Step 4: Run focused test to verify RED**

Run:

```bash
julia --project=. -e 'include("test/public/laurent_gl_certificate_options.jl")'
```

Expected: FAIL before implementation because the public constructor does not accept `determinant_strategy` and `LaurentLazyGLHoistCertificate` is not exported.

- [ ] **Step 5: Commit the failing tests**

Run:

```bash
git add test/public/laurent_gl_certificate_options.jl test/public/api_surface.jl test/runtests.jl
git commit -m "test: cover public lazy laurent certificate options"
```

---

### Task 2: Expose Lazy Certificate Fields Through the Public API

**Files:**
- Modify: `src/Suslin.jl`
- Modify: `src/algorithm/laurent_gl_certificate.jl`
- Modify: `test/expert/laurent_lazy_correction_hoist.jl`
- Modify: `test/expert/laurent_lazy_row_column_correction.jl`

**Interfaces:**
- Consumes existing internal `_laurent_gl_lazy_deferred_correction_certificate(A; correction_side)`.
- Produces exported `LaurentLazyGLHoistCertificate`, direct `determinant_source` field, public determinant strategy keyword routing, and public lazy verification dispatch.

- [ ] **Step 1: Export the lazy certificate type**

In `src/Suslin.jl`, add:

```julia
export LaurentLazyGLHoistCertificate
```

immediately after `export LaurentGLFactorizationCertificate`.

- [ ] **Step 2: Add the direct lazy determinant source field**

In `src/algorithm/laurent_gl_certificate.jl`, update `LaurentLazyGLHoistCertificate`:

```julia
struct LaurentLazyGLHoistCertificate
    original_matrix
    deferred_metadata
    overall_determinant
    determinant_source::Symbol
    correction_side::Symbol
    reconstruction_relation::Symbol
    correction
    inverse_correction
    normalized_deferred_core
    normalized_deferred_factorization
    normalized_deferred_factors::Vector
    rewritten_left_factors::Vector
    rewritten_right_factors::Vector
    elementary_factors::Vector
    elementary_product
    reconstructed_product
    verification
end
```

Update both `LaurentLazyGLHoistCertificate(...)` calls in `_laurent_gl_lazy_deferred_correction_certificate` to pass `metadata.determinant_source` and then `certificate.determinant_source` after `overall_determinant`.

- [ ] **Step 3: Verify the direct determinant source field**

In `_laurent_gl_lazy_deferred_correction_certificate_verification`, add a `determinant_source_ok` boolean. Set it with:

```julia
determinant_source_ok =
    certificate.determinant_source == :deferred_submatrix &&
    certificate.determinant_source == metadata.determinant_source
```

Include `determinant_source_ok` in `overall_ok` and in the returned named tuple.

- [ ] **Step 4: Add public determinant strategy routing**

Replace the one-argument-only public constructor with keyword routing:

```julia
function _laurent_gl_certificate_strategy(determinant_strategy::Symbol)::Symbol
    determinant_strategy in (:eager, :lazy) && return determinant_strategy
    throw(ArgumentError("determinant_strategy must be :eager or :lazy"))
end

function laurent_gl_factorization_certificate(
    A;
    determinant_strategy = :eager,
    correction_side = nothing,
    progress_callback = nothing,
)
    strategy = _laurent_gl_certificate_strategy(determinant_strategy)
    if strategy == :eager
        correction_side === nothing ||
            throw(ArgumentError("correction_side is supported only with determinant_strategy = :lazy"))
        return _laurent_gl_factorization_certificate(
            A;
            progress_callback,
        )
    end

    side = correction_side === nothing ? :row : correction_side
    return _laurent_gl_lazy_deferred_correction_certificate(
        A;
        correction_side = side,
        progress_callback,
    )
end
```

Leave `_laurent_gl_factorization_certificate(A; progress_callback = nothing)` as the eager internal route.

- [ ] **Step 5: Dispatch public verification for lazy certificates**

Add a specific verification method before or after the generic method:

```julia
function verify_laurent_gl_factorization_certificate(
    certificate::LaurentLazyGLHoistCertificate,
)::Bool
    return _laurent_gl_lazy_deferred_correction_certificate_verification(certificate).overall_ok
end
```

Keep the existing generic method for eager `LaurentGLFactorizationCertificate` and malformed tests.

- [ ] **Step 6: Update existing expert test helper constructors**

In `test/expert/laurent_lazy_correction_hoist.jl`, add `certificate.determinant_source` after `certificate.overall_determinant` in both `LaurentLazyGLHoistCertificate(...)` helper constructor calls.

In `test/expert/laurent_lazy_row_column_correction.jl`, add a keyword default:

```julia
determinant_source = certificate.determinant_source,
```

and pass `determinant_source` after `certificate.overall_determinant` in `_issue159_rebuild`.

- [ ] **Step 7: Run focused tests to verify GREEN**

Run:

```bash
julia --project=. -e 'include("test/public/laurent_gl_certificate_options.jl")'
julia --project=. -e 'include("test/public/api_surface.jl")'
julia --project=. -e 'include("test/expert/laurent_lazy_correction_hoist.jl"); include("test/expert/laurent_lazy_row_column_correction.jl")'
```

Expected: all commands exit 0.

- [ ] **Step 8: Commit implementation**

Run:

```bash
git add src/Suslin.jl src/algorithm/laurent_gl_certificate.jl test/expert/laurent_lazy_correction_hoist.jl test/expert/laurent_lazy_row_column_correction.jl
git commit -m "feat: expose lazy laurent certificate fields"
```

---

### Task 3: Full Verification and PR Readiness

**Files:**
- No planned source edits unless verification exposes a defect.

**Interfaces:**
- Consumes the public API from Task 2.
- Produces final verification evidence for PR creation.

- [ ] **Step 1: Run issue verification**

Run:

```bash
julia --project=. -e 'include("test/public/laurent_gl_certificate_options.jl")'
```

Expected: exits 0.

- [ ] **Step 2: Run package verification**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: exits 0.

- [ ] **Step 3: Check git status**

Run:

```bash
git status --short
```

Expected: no uncommitted files.
