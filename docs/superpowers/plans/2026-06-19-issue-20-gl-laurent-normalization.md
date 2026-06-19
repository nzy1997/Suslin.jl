# Issue 20 GL_n Laurent Normalization Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an explicit Laurent `GL_n` determinant normalization boundary before Suslin's staged `SL_n` factorization path.

**Architecture:** Put determinant classification, determinant correction construction, and exact reconstruction verification in a small core module. Export the named boundary API, and have `elementary_factorization` call it for Laurent matrices before continuing to the existing narrow `SL_3` checks.

**Tech Stack:** Julia, Oscar/AbstractAlgebra Laurent polynomial rings and matrices, Test stdlib, existing grouped test runner.

## Global Constraints

- Do not implement full Laurent elementary factorization in this issue.
- Do not hide determinant normalization inside algorithm internals; expose a clearly named validation/normalization layer.
- Supported determinant corrections are `:one`, `:permutation_sign_unit`, and `:laurent_monomial_unit`.
- Unsupported `:other_unit` and `:non_unit` determinants must fail with staged `ArgumentError` messages.
- Exact reconstruction must be checked as `correction.factor * normalized_matrix == original`.
- Register focused tests in the existing `internal` group.
- Run `julia --project=. -e 'using Pkg; Pkg.test()'` and `julia --project=. test/runtests.jl all` before finishing.
- Do not commit `Manifest.toml`.

---

## File Structure

- Create `src/core/gl_laurent_normalization.jl`: determinant classification, Laurent monomial metadata extraction, normalization construction, and exact reconstruction verification.
- Modify `src/Suslin.jl`: export `classify_laurent_determinant`, `normalize_laurent_gl_matrix`, and `verify_laurent_gl_normalization`; include the new core file before `algorithm/factorization.jl`.
- Modify `src/algorithm/factorization.jl`: call `normalize_laurent_gl_matrix(A)` for Laurent matrices before the current algorithm-specific checks.
- Create `test/internal/gl_laurent_normalization.jl`: focused tests for determinant-one, monomial-unit, permutation-sign, staged non-unit failure, tampered metadata, and driver boundary behavior.
- Modify `test/runtests.jl`: add the focused test file to the `internal` group.
- Modify `test/public/api_surface.jl`: assert the new public API functions are exported.
- Modify `docs/src/toricbuilder_contract.md`: document the new normalization boundary and staged behavior.

---

### Task 1: Laurent GL_n Normalization API

**Files:**
- Create: `src/core/gl_laurent_normalization.jl`
- Modify: `src/Suslin.jl`
- Create: `test/internal/gl_laurent_normalization.jl`
- Modify: `test/runtests.jl`

**Interfaces:**
- Consumes: `Suslin._is_laurent_polynomial_ring`, `Suslin._require_laurent_polynomial_ring`, Oscar `det`, `is_unit`, `exponents`, `coefficients`, `identity_matrix`, and exact matrix multiplication.
- Produces: `classify_laurent_determinant(A)`, `normalize_laurent_gl_matrix(A)`, and `verify_laurent_gl_normalization(A, normalization)::Bool`.

- [ ] **Step 1: Write the failing focused normalization tests**

Create `test/internal/gl_laurent_normalization.jl` with:

```julia
using Test
using Suslin
using Oscar

@testset "Laurent GL_n determinant classification and normalization" begin
    R, (x, y) = suslin_laurent_polynomial_ring(GF(2), ["x", "y"])

    determinant_one = matrix(R, [
        one(R) x;
        zero(R) one(R)
    ])
    one_normalization = normalize_laurent_gl_matrix(determinant_one)
    @test classify_laurent_determinant(determinant_one).classification == :one
    @test one_normalization.determinant_classification == :one
    @test one_normalization.normalized_matrix == determinant_one
    @test one_normalization.correction.kind == :identity
    @test verify_laurent_gl_normalization(determinant_one, one_normalization)

    monomial_unit = matrix(R, [
        x^-1 * y one(R);
        zero(R) one(R)
    ])
    monomial_profile = classify_laurent_determinant(monomial_unit)
    @test monomial_profile.classification == :laurent_monomial_unit
    @test monomial_profile.monomial_exponents == (-1, 1)
    monomial_normalization = normalize_laurent_gl_matrix(monomial_unit)
    @test monomial_normalization.determinant_classification == :laurent_monomial_unit
    @test det(monomial_normalization.normalized_matrix) == one(R)
    @test monomial_normalization.correction.kind == :left_diagonal_determinant_correction
    @test monomial_normalization.correction.side == :left
    @test verify_laurent_gl_normalization(monomial_unit, monomial_normalization)
    @test monomial_normalization.correction.factor * monomial_normalization.normalized_matrix == monomial_unit

    tampered_factor = copy(monomial_normalization.correction.factor)
    tampered_factor[1, 1] = one(R)
    tampered = merge(
        monomial_normalization,
        (;
            correction = merge(
                monomial_normalization.correction,
                (; factor = tampered_factor),
            ),
        ),
    )
    @test !verify_laurent_gl_normalization(monomial_unit, tampered)

    Q, (t,) = suslin_laurent_polynomial_ring(QQ, ["t"])
    transposition = matrix(Q, [
        zero(Q) one(Q);
        one(Q) zero(Q)
    ])
    sign_profile = classify_laurent_determinant(transposition)
    @test sign_profile.classification == :permutation_sign_unit
    sign_normalization = normalize_laurent_gl_matrix(transposition)
    @test sign_normalization.determinant_classification == :permutation_sign_unit
    @test det(sign_normalization.normalized_matrix) == one(Q)
    @test verify_laurent_gl_normalization(transposition, sign_normalization)

    non_unit = matrix(R, [
        x + one(R) zero(R);
        zero(R) one(R)
    ])
    non_unit_profile = classify_laurent_determinant(non_unit)
    @test non_unit_profile.classification == :non_unit
    err = try
        normalize_laurent_gl_matrix(non_unit)
        nothing
    catch caught
        caught
    end
    @test err isa ArgumentError
    @test occursin("unsupported Laurent GL_n determinant", sprint(showerror, err))
    @test occursin("outside the staged SL_n factorization path", sprint(showerror, err))
end
```

Add the new file to `test/runtests.jl`:

```julia
"internal/gl_laurent_normalization.jl",
```

- [ ] **Step 2: Run the focused test to verify it fails for the missing API**

Run:

```bash
julia --project=. test/runtests.jl internal
```

Expected: FAIL with `UndefVarError: normalize_laurent_gl_matrix not defined` or an earlier dependency-instantiation error in a fresh checkout.

- [ ] **Step 3: Implement the minimal normalization API**

Create `src/core/gl_laurent_normalization.jl` with:

```julia
function _laurent_monomial_metadata(value)
    try
        exponent_vectors = collect(exponents(value))
        coeffs = collect(coefficients(value))
        length(exponent_vectors) == 1 || return nothing
        length(coeffs) == 1 || return nothing
        return (;
            monomial_exponents = Tuple(Int.(collect(exponent_vectors[1]))),
            monomial_coefficient = coeffs[1],
        )
    catch err
        err isa MethodError || err isa ErrorException || rethrow()
        return nothing
    end
end

function classify_laurent_determinant(A)
    nrows(A) == ncols(A) || throw(ArgumentError("Laurent GL_n determinant classification requires a square matrix"))
    R = _require_laurent_polynomial_ring(base_ring(A); label="A base ring")
    determinant = det(A)
    determinant == one(R) && return (;
        determinant,
        classification = :one,
        monomial_exponents = ntuple(_ -> 0, ngens(R)),
        monomial_coefficient = one(R),
    )

    negative_one = -one(R)
    if negative_one != one(R) && determinant == negative_one
        return (;
            determinant,
            classification = :permutation_sign_unit,
            monomial_exponents = ntuple(_ -> 0, ngens(R)),
            monomial_coefficient = negative_one,
        )
    end

    is_unit(determinant) || return (;
        determinant,
        classification = :non_unit,
        monomial_exponents = nothing,
        monomial_coefficient = nothing,
    )

    monomial = _laurent_monomial_metadata(determinant)
    monomial !== nothing && return (;
        determinant,
        classification = :laurent_monomial_unit,
        monomial_exponents = monomial.monomial_exponents,
        monomial_coefficient = monomial.monomial_coefficient,
    )

    return (;
        determinant,
        classification = :other_unit,
        monomial_exponents = nothing,
        monomial_coefficient = nothing,
    )
end

function _identity_correction(R, n::Int, determinant)
    identity = identity_matrix(R, n)
    return (;
        kind = :identity,
        side = :left,
        factor = identity,
        inverse_factor = identity,
        determinant,
    )
end

function _left_diagonal_determinant_correction(R, n::Int, determinant)
    factor = identity_matrix(R, n)
    factor[1, 1] = determinant
    inverse_factor = identity_matrix(R, n)
    inverse_factor[1, 1] = inv(determinant)
    return (;
        kind = :left_diagonal_determinant_correction,
        side = :left,
        factor,
        inverse_factor,
        determinant,
    )
end

function normalize_laurent_gl_matrix(A)
    nrows(A) == ncols(A) || throw(ArgumentError("Laurent GL_n normalization requires a square matrix"))
    R = _require_laurent_polynomial_ring(base_ring(A); label="A base ring")
    n = nrows(A)
    determinant_profile = classify_laurent_determinant(A)
    classification = determinant_profile.classification

    if classification == :non_unit
        throw(ArgumentError("unsupported Laurent GL_n determinant: determinant is non-unit, so the input is outside the staged SL_n factorization path"))
    elseif classification == :other_unit
        throw(ArgumentError("unsupported Laurent GL_n determinant: non-monomial units are outside the staged SL_n factorization path"))
    end

    correction = classification == :one ?
        _identity_correction(R, n, determinant_profile.determinant) :
        _left_diagonal_determinant_correction(R, n, determinant_profile.determinant)
    normalized_matrix = correction.inverse_factor * A
    normalization = (;
        input_size = (n, n),
        ring = R,
        determinant = determinant_profile.determinant,
        determinant_classification = classification,
        determinant_profile,
        normalized_matrix,
        correction,
    )

    verify_laurent_gl_normalization(A, normalization) || throw(ArgumentError("Laurent GL_n normalization failed exact reconstruction verification"))
    return normalization
end

function verify_laurent_gl_normalization(A, normalization)::Bool
    try
        nrows(A) == ncols(A) || return false
        n = nrows(A)
        R = base_ring(A)
        _require_laurent_polynomial_ring(R; label="A base ring")
        normalization.input_size == (n, n) || return false
        normalization.ring === R || normalization.ring == R || return false
        nrows(normalization.normalized_matrix) == n || return false
        ncols(normalization.normalized_matrix) == n || return false
        base_ring(normalization.normalized_matrix) == R || return false
        correction = normalization.correction
        correction.side == :left || return false
        nrows(correction.factor) == n || return false
        ncols(correction.factor) == n || return false
        nrows(correction.inverse_factor) == n || return false
        ncols(correction.inverse_factor) == n || return false
        base_ring(correction.factor) == R || return false
        base_ring(correction.inverse_factor) == R || return false
        correction.factor * correction.inverse_factor == identity_matrix(R, n) || return false
        correction.inverse_factor * correction.factor == identity_matrix(R, n) || return false
        correction.factor * normalization.normalized_matrix == A || return false
        det(normalization.normalized_matrix) == one(R) || return false
        return true
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end
```

Modify `src/Suslin.jl`:

```julia
export classify_laurent_determinant
export normalize_laurent_gl_matrix
export verify_laurent_gl_normalization
```

and include the new file after `include("core/rings.jl")`:

```julia
include("core/gl_laurent_normalization.jl")
```

- [ ] **Step 4: Run the focused internal tests**

Run:

```bash
julia --project=. test/runtests.jl internal
```

Expected: PASS for `internal/gl_laurent_normalization.jl`, unless dependency instantiation is still missing.

- [ ] **Step 5: Commit**

Run:

```bash
git add src/Suslin.jl src/core/gl_laurent_normalization.jl test/internal/gl_laurent_normalization.jl test/runtests.jl
git commit -m "feat: add laurent gl determinant normalization"
```

Expected: commit succeeds and `Manifest.toml` is not staged.

---

### Task 2: Driver Boundary, Public API, and Docs

**Files:**
- Modify: `src/algorithm/factorization.jl`
- Modify: `test/internal/gl_laurent_normalization.jl`
- Modify: `test/public/api_surface.jl`
- Modify: `docs/src/toricbuilder_contract.md`

**Interfaces:**
- Consumes: `normalize_laurent_gl_matrix(A)` from Task 1.
- Produces: `elementary_factorization(A)` performs Laurent determinant normalization before the existing staged `SL_3` algorithm checks; public API tests cover the new exports.

- [ ] **Step 1: Add failing driver and public API tests**

Append to `test/internal/gl_laurent_normalization.jl`:

```julia
@testset "elementary factorization Laurent GL_n boundary" begin
    R, (x,) = suslin_laurent_polynomial_ring(GF(2), ["x"])

    normalized_then_rejected = matrix(R, [
        x zero(R);
        zero(R) one(R)
    ])
    err = try
        elementary_factorization(normalized_then_rejected)
        nothing
    catch caught
        caught
    end
    @test err isa ArgumentError
    @test occursin("currently supports only 3x3 matrices", sprint(showerror, err))

    non_unit = matrix(R, [
        x + one(R) zero(R) zero(R);
        zero(R) one(R) zero(R);
        zero(R) zero(R) one(R)
    ])
    non_unit_err = try
        elementary_factorization(non_unit)
        nothing
    catch caught
        caught
    end
    @test non_unit_err isa ArgumentError
    @test occursin("unsupported Laurent GL_n determinant", sprint(showerror, non_unit_err))
end
```

Add public API assertions to `test/public/api_surface.jl`:

```julia
@test isdefined(Suslin, :classify_laurent_determinant)
@test isdefined(Suslin, :normalize_laurent_gl_matrix)
@test isdefined(Suslin, :verify_laurent_gl_normalization)
@test Suslin.classify_laurent_determinant === classify_laurent_determinant
@test Suslin.normalize_laurent_gl_matrix === normalize_laurent_gl_matrix
@test Suslin.verify_laurent_gl_normalization === verify_laurent_gl_normalization
```

- [ ] **Step 2: Run focused tests to verify the driver test fails**

Run:

```bash
julia --project=. test/runtests.jl public internal
```

Expected: FAIL because `elementary_factorization` still rejects the Laurent `2 x 2` matrix before calling the normalization boundary, unless dependency instantiation is missing.

- [ ] **Step 3: Wire the public driver to the boundary**

Modify the top of `elementary_factorization(A)` in `src/algorithm/factorization.jl` to:

```julia
function elementary_factorization(A)
    if _is_laurent_polynomial_ring(base_ring(A))
        normalization = normalize_laurent_gl_matrix(A)
        A = normalization.normalized_matrix
    end

    nrows(A) == ncols(A) || throw(ArgumentError("A must be square"))
    nrows(A) == 3 || throw(ArgumentError("elementary_factorization currently supports only 3x3 matrices"))
```

Leave the existing univariate and local `SL_3` checks in place.

- [ ] **Step 4: Document the new boundary**

Add this section to `docs/src/toricbuilder_contract.md` before `## Suslin Output Contract`:

```markdown
## GL_n Laurent Normalization Boundary

Suslin exposes `normalize_laurent_gl_matrix(A)` for exact Laurent `GL_n`
inputs before any staged `SL_n` factorization attempt. The boundary computes
and classifies the determinant, then either returns a determinant-one core
with explicit correction metadata or throws a staged `ArgumentError`.

Supported corrections are determinant `1`, permutation/sign determinant `-1`
where the coefficient ring distinguishes it from `1`, and Laurent monomial
unit determinants such as `x^-1*y` over `GF(2)`. The correction metadata stores
a left diagonal factor `D` and verifies exact reconstruction as
`D * normalized_matrix == A`.

Non-unit determinants and non-monomial units remain outside the staged `SL_n`
path. `elementary_factorization` calls this boundary for Laurent matrices
before it continues to the current narrow `3 x 3` univariate algorithm checks.
```

- [ ] **Step 5: Run focused public and internal tests**

Run:

```bash
julia --project=. test/runtests.jl public internal
```

Expected: PASS for public and internal groups, unless dependency instantiation is missing.

- [ ] **Step 6: Commit**

Run:

```bash
git add src/algorithm/factorization.jl test/internal/gl_laurent_normalization.jl test/public/api_surface.jl docs/src/toricbuilder_contract.md
git commit -m "feat: wire laurent gl normalization boundary"
```

Expected: commit succeeds and `Manifest.toml` is not staged.

---

## Final Verification

- [ ] Run dependency setup if needed:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

- [ ] Run the required package test command:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

- [ ] Run the documented full suite command from Issue #21:

```bash
julia --project=. test/runtests.jl all
```

- [ ] Check git status:

```bash
git status --short
```

Expected: no uncommitted implementation files; no `Manifest.toml` staged.

## Plan Self-Review

- Spec coverage: both supported output modes and staged failure mode are covered by Task 1 and Task 2.
- Placeholder scan: no placeholder-only implementation steps remain.
- Type consistency: all task references use the same exported function names and named-tuple field names.
