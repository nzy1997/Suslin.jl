# Issue 8 Laurent Normalization Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add helpers that clear negative Laurent exponents into ordinary polynomial vectors/matrices and lift them back exactly.

**Architecture:** Keep the normalization layer in `src/core/polynomials.jl`, which is already the reserved split point for polynomial helpers. Normalize vectors with one monomial shift and matrices with one shift per column; store the original Laurent ring, ordinary polynomial ring, shift exponent vectors, monomial units, and determinant shift metadata. Expose a small public helper API plus focused fixture-based tests.

**Tech Stack:** Julia, Oscar/AbstractAlgebra Laurent and ordinary polynomial rings, Oscar matrices, Test stdlib, existing grouped test runner.

## Global Constraints

- Reuse the shared Laurent fixtures from #6 instead of inventing ad hoc examples.
- Inputs whose entries do not all belong to the same Laurent ring must throw `ArgumentError`.
- Matrix normalization uses one monomial shift per column.
- Column/vector normalization uses one monomial shift for the whole column.
- Metadata must record monomial shift exponents and monomial units for later determinant/unit accounting.
- For square matrices, metadata must record `determinant_shift_exponents = sum(column_shifts)`.
- Lift-back must reconstruct the original Laurent object exactly.
- Tampering with shift metadata must make exact verification fail.
- Register focused tests in the existing `internal` group.
- Run `julia --project=. -e 'include("test/internal/laurent_normalization.jl")'`, `julia --project=. -e 'using Pkg; Pkg.test()'`, and `julia --project=. test/runtests.jl all` before finishing.
- Do not commit `Manifest.toml`.

---

## File Structure

- Modify `src/core/polynomials.jl`: Laurent normalization, metadata construction, lift-back, and verification helpers.
- Modify `src/Suslin.jl`: export `normalize_laurent_object`, `lift_laurent_normalization`, and `verify_laurent_normalization`.
- Create `test/internal/laurent_normalization.jl`: focused tests using `test/fixtures/laurent_cases.jl`.
- Modify `test/runtests.jl`: include `internal/laurent_normalization.jl` in the `internal` group.
- Modify `test/public/api_surface.jl`: assert the new public helper names are exported.

---

### Task 1: Laurent Normalization and Lift-Back Helpers

**Files:**
- Modify: `src/core/polynomials.jl`
- Modify: `src/Suslin.jl`
- Create: `test/internal/laurent_normalization.jl`
- Modify: `test/runtests.jl`
- Modify: `test/public/api_surface.jl`

**Interfaces:**
- Consumes: `Suslin._require_laurent_polynomial_ring`, `Suslin._require_laurent_element`, `Suslin._require_same_laurent_parent`, Oscar `coefficient_ring`, `symbols`, `gens`, `exponents`, `coefficients`, `zero_matrix`, `base_ring`, `nrows`, and `ncols`.
- Produces: `normalize_laurent_object(obj)`, `lift_laurent_normalization(normalization)`, `lift_laurent_normalization(polynomial_object, metadata)`, and `verify_laurent_normalization(original, normalization)::Bool`.

- [ ] **Step 1: Write the failing focused tests**

Create `test/internal/laurent_normalization.jl` with:

```julia
using Test
using Suslin
using Oscar

const LAURENT_NORMALIZATION_FIXTURES = joinpath(@__DIR__, "..", "fixtures", "laurent_cases.jl")

function _fixture_by_id(catalog, id)
    return only(filter(entry -> entry.id == id, catalog.cases))
end

function _all_polynomial_matrix_entries_have_nonnegative_exponents(A)
    for value in A
        for exponent_vector in collect(exponents(value))
            all(exponent -> exponent >= 0, exponent_vector) || return false
        end
    end
    return true
end

function _all_polynomial_vector_entries_have_nonnegative_exponents(values)
    for value in values
        for exponent_vector in collect(exponents(value))
            all(exponent -> exponent >= 0, exponent_vector) || return false
        end
    end
    return true
end

@testset "Laurent normalization and lift-back helpers" begin
    include(LAURENT_NORMALIZATION_FIXTURES)
    catalog = LaurentFixtureCatalog.catalog()

    column_fixture = _fixture_by_id(catalog, "laurent-negative-exponent-normalization")
    column = column_fixture.inputs.vector
    column_normalization = normalize_laurent_object(column)

    @test column_normalization.metadata.kind == :matrix
    @test column_normalization.metadata.shape == (2, 1)
    @test column_normalization.metadata.column_shifts == ((2, 0),)
    @test column_normalization.metadata.determinant_shift_exponents === nothing
    @test base_ring(column_normalization.normalized_object) == column_normalization.metadata.polynomial_ring
    @test _all_polynomial_matrix_entries_have_nonnegative_exponents(column_normalization.normalized_object)
    @test lift_laurent_normalization(column_normalization) == column
    @test lift_laurent_normalization(column_normalization.normalized_object, column_normalization.metadata) == column
    @test verify_laurent_normalization(column, column_normalization)

    matrix_fixture = _fixture_by_id(catalog, "toricbuilder-factor-toric-block-3-pinv")
    A = matrix_fixture.inputs.matrix
    matrix_normalization = normalize_laurent_object(A)

    @test matrix_normalization.metadata.kind == :matrix
    @test matrix_normalization.metadata.shape == (8, 8)
    @test length(matrix_normalization.metadata.column_shifts) == 8
    @test matrix_normalization.metadata.determinant_shift_exponents ==
        ntuple(i -> sum(shift[i] for shift in matrix_normalization.metadata.column_shifts), 2)
    @test base_ring(matrix_normalization.normalized_object) == matrix_normalization.metadata.polynomial_ring
    @test _all_polynomial_matrix_entries_have_nonnegative_exponents(matrix_normalization.normalized_object)
    @test lift_laurent_normalization(matrix_normalization) == A
    @test verify_laurent_normalization(A, matrix_normalization)

    R = column_fixture.ring.object
    x, y = column_fixture.ring.generators
    vector_normalization = normalize_laurent_object([x^-1, x^-2 * y])
    @test vector_normalization.metadata.kind == :vector
    @test vector_normalization.metadata.shape == (2,)
    @test vector_normalization.metadata.column_shifts == ((2, 0),)
    @test _all_polynomial_vector_entries_have_nonnegative_exponents(vector_normalization.normalized_object)
    @test lift_laurent_normalization(vector_normalization) == [x^-1, x^-2 * y]

    S, (u, v) = suslin_laurent_polynomial_ring(GF(2), ["u", "v"])
    @test_throws ArgumentError normalize_laurent_object([x, u])

    tampered_shift = ((1, 0),)
    tampered_metadata = merge(
        column_normalization.metadata,
        (;
            column_shifts = tampered_shift,
            shift_monomials = (x,),
            inverse_shift_monomials = (x^-1,),
        ),
    )
    tampered = (;
        normalized_object = column_normalization.normalized_object,
        metadata = tampered_metadata,
    )
    @test lift_laurent_normalization(tampered) != column
    @test !verify_laurent_normalization(column, tampered)
end
```

Modify `test/public/api_surface.jl` inside the existing testset:

```julia
@test isdefined(Suslin, :normalize_laurent_object)
@test isdefined(Suslin, :lift_laurent_normalization)
@test isdefined(Suslin, :verify_laurent_normalization)
@test Suslin.normalize_laurent_object === normalize_laurent_object
@test Suslin.lift_laurent_normalization === lift_laurent_normalization
@test Suslin.verify_laurent_normalization === verify_laurent_normalization
```

Modify `test/runtests.jl` by adding the focused file to the `internal` list after `internal/laurent_fixtures.jl`:

```julia
"internal/laurent_normalization.jl",
```

- [ ] **Step 2: Run the focused test to verify it fails for the missing API**

Run:

```bash
julia --project=. -e 'include("test/internal/laurent_normalization.jl")'
```

Expected: FAIL with `UndefVarError: normalize_laurent_object not defined`.

- [ ] **Step 3: Implement the minimal helper API**

Modify `src/Suslin.jl` by adding these exports near the existing Laurent exports:

```julia
export normalize_laurent_object
export lift_laurent_normalization
export verify_laurent_normalization
```

Replace `src/core/polynomials.jl` with the existing reserved comment plus the helper implementation. Keep helper names underscore-prefixed except for the public API functions:

```julia
# Reserved split point for polynomial helpers that the later algorithm layers build on.

function _laurent_normalization_variable_names(R)
    return Tuple(string.(symbols(R)))
end

function _ordinary_polynomial_ring_for_laurent(R)
    P, _ = suslin_polynomial_ring(coefficient_ring(R), collect(_laurent_normalization_variable_names(R)))
    return P
end

function _laurent_monomial_from_exponents(R, exponent_vector)
    length(exponent_vector) == ngens(R) || throw(ArgumentError("shift exponent vector length must match the Laurent ring generators"))
    term = one(R)
    for (j, exponent) in enumerate(exponent_vector)
        exponent == 0 && continue
        term *= gen(R, j)^exponent
    end
    return term
end

function _column_shift_for_laurent_entries(entries, R)
    min_exponents = zeros(Int, ngens(R))
    for value in entries
        _require_laurent_element(value, R; label="normalization entry")
        for raw_exponents in collect(exponents(value))
            exponent_vector = Int.(collect(raw_exponents))
            length(exponent_vector) == ngens(R) || throw(ArgumentError("Laurent exponent vector length must match the parent ring"))
            for j in 1:ngens(R)
                min_exponents[j] = min(min_exponents[j], exponent_vector[j])
            end
        end
    end
    return ntuple(j -> max(0, -min_exponents[j]), ngens(R))
end

function _polynomial_term_from_exponents(P, coeff, exponent_vector)
    term = P(coeff)
    for (j, exponent) in enumerate(exponent_vector)
        exponent < 0 && throw(ArgumentError("normalized Laurent term still has a negative exponent"))
        exponent == 0 && continue
        term *= gen(P, j)^exponent
    end
    return term
end

function _laurent_entry_to_polynomial(value, R, P, shift)
    _require_laurent_element(value, R; label="normalization entry")
    length(shift) == ngens(R) || throw(ArgumentError("shift exponent vector length must match the Laurent ring generators"))

    result = zero(P)
    for (coeff, raw_exponents) in zip(collect(coefficients(value)), collect(exponents(value)))
        exponent_vector = Int.(collect(raw_exponents))
        shifted_exponents = ntuple(j -> exponent_vector[j] + shift[j], ngens(R))
        result += _polynomial_term_from_exponents(P, coeff, shifted_exponents)
    end
    return result
end

function _sum_laurent_shift_exponents(column_shifts, R)
    return ntuple(j -> sum(shift[j] for shift in column_shifts), ngens(R))
end

function _laurent_normalization_metadata(kind::Symbol, shape, R, P, column_shifts)
    shift_monomials = map(shift -> _laurent_monomial_from_exponents(R, shift), column_shifts)
    inverse_shift_monomials = map(
        shift -> _laurent_monomial_from_exponents(R, ntuple(j -> -shift[j], ngens(R))),
        column_shifts,
    )
    determinant_shift_exponents = kind == :matrix && shape[1] == shape[2] ?
        _sum_laurent_shift_exponents(column_shifts, R) :
        nothing
    return (;
        kind,
        shape,
        laurent_ring = R,
        polynomial_ring = P,
        variable_names = _laurent_normalization_variable_names(R),
        column_shifts,
        shift_monomials,
        inverse_shift_monomials,
        determinant_shift_exponents,
    )
end

function _has_oscar_matrix_interface(obj)::Bool
    try
        nrows(obj)
        ncols(obj)
        base_ring(obj)
        return true
    catch err
        err isa MethodError || err isa ErrorException || rethrow()
        return false
    end
end

function _normalize_laurent_matrix(A)
    R = _require_laurent_polynomial_ring(base_ring(A); label="input base ring")
    rows, cols = nrows(A), ncols(A)
    P = _ordinary_polynomial_ring_for_laurent(R)
    column_shifts = ntuple(
        j -> _column_shift_for_laurent_entries((A[i, j] for i in 1:rows), R),
        cols,
    )

    normalized = zero_matrix(P, rows, cols)
    for j in 1:cols
        for i in 1:rows
            normalized[i, j] = _laurent_entry_to_polynomial(A[i, j], R, P, column_shifts[j])
        end
    end

    metadata = _laurent_normalization_metadata(:matrix, (rows, cols), R, P, column_shifts)
    return (; normalized_object = normalized, metadata)
end

function _normalize_laurent_vector(values::AbstractVector)
    R = _require_same_laurent_parent(values; label="vector entries")
    P = _ordinary_polynomial_ring_for_laurent(R)
    shift = _column_shift_for_laurent_entries(values, R)
    normalized = [_laurent_entry_to_polynomial(value, R, P, shift) for value in values]
    metadata = _laurent_normalization_metadata(:vector, (length(values),), R, P, (shift,))
    return (; normalized_object = normalized, metadata)
end

function normalize_laurent_object(obj)
    _has_oscar_matrix_interface(obj) && return _normalize_laurent_matrix(obj)
    obj isa AbstractVector && return _normalize_laurent_vector(obj)
    throw(ArgumentError("input must be an Oscar matrix over a Laurent polynomial ring or a vector of Laurent elements"))
end

function _require_laurent_normalization_metadata(metadata)
    required = (
        :kind,
        :shape,
        :laurent_ring,
        :polynomial_ring,
        :variable_names,
        :column_shifts,
        :shift_monomials,
        :inverse_shift_monomials,
        :determinant_shift_exponents,
    )
    for field in required
        hasproperty(metadata, field) || throw(ArgumentError("normalization metadata missing field $(field)"))
    end

    R = _require_laurent_polynomial_ring(metadata.laurent_ring; label="metadata Laurent ring")
    metadata.variable_names == _laurent_normalization_variable_names(R) || throw(ArgumentError("normalization metadata variable names do not match the Laurent ring"))
    metadata.column_shifts isa Tuple || throw(ArgumentError("normalization metadata column shifts must be a tuple"))

    for shift in metadata.column_shifts
        length(shift) == ngens(R) || throw(ArgumentError("shift exponent vector length must match the Laurent ring generators"))
        all(exponent -> exponent >= 0, shift) || throw(ArgumentError("shift exponents must be nonnegative"))
    end

    expected_shift_monomials = map(shift -> _laurent_monomial_from_exponents(R, shift), metadata.column_shifts)
    expected_inverse_shift_monomials = map(
        shift -> _laurent_monomial_from_exponents(R, ntuple(j -> -shift[j], ngens(R))),
        metadata.column_shifts,
    )
    metadata.shift_monomials == expected_shift_monomials || throw(ArgumentError("normalization shift monomials do not match shift exponents"))
    metadata.inverse_shift_monomials == expected_inverse_shift_monomials || throw(ArgumentError("normalization inverse shift monomials do not match shift exponents"))

    if metadata.kind == :matrix && metadata.shape[1] == metadata.shape[2]
        metadata.determinant_shift_exponents == _sum_laurent_shift_exponents(metadata.column_shifts, R) ||
            throw(ArgumentError("normalization determinant shift metadata does not match column shifts"))
    elseif metadata.determinant_shift_exponents !== nothing
        throw(ArgumentError("normalization determinant shift metadata is only defined for square matrices"))
    end

    return metadata
end

function _lift_laurent_matrix(polynomial_object, metadata)
    R = metadata.laurent_ring
    P = metadata.polynomial_ring
    rows, cols = metadata.shape
    nrows(polynomial_object) == rows || throw(ArgumentError("normalized matrix row count does not match metadata"))
    ncols(polynomial_object) == cols || throw(ArgumentError("normalized matrix column count does not match metadata"))
    base_ring(polynomial_object) == P || throw(ArgumentError("normalized matrix base ring does not match metadata"))
    length(metadata.column_shifts) == cols || throw(ArgumentError("matrix metadata must record one shift per column"))

    lifted = zero_matrix(R, rows, cols)
    for j in 1:cols
        inverse_shift = metadata.inverse_shift_monomials[j]
        for i in 1:rows
            lifted[i, j] = inverse_shift * R(polynomial_object[i, j])
        end
    end
    return lifted
end

function _lift_laurent_vector(polynomial_object, metadata)
    R = metadata.laurent_ring
    P = metadata.polynomial_ring
    length(polynomial_object) == metadata.shape[1] || throw(ArgumentError("normalized vector length does not match metadata"))
    length(metadata.column_shifts) == 1 || throw(ArgumentError("vector metadata must record exactly one shift"))
    inverse_shift = only(metadata.inverse_shift_monomials)

    lifted = Vector{Any}(undef, length(polynomial_object))
    for (i, value) in enumerate(polynomial_object)
        parent(value) == P || throw(ArgumentError("normalized vector entry parent does not match metadata"))
        lifted[i] = inverse_shift * R(value)
    end
    return lifted
end

function lift_laurent_normalization(polynomial_object, metadata)
    metadata = _require_laurent_normalization_metadata(metadata)
    if metadata.kind == :matrix
        return _lift_laurent_matrix(polynomial_object, metadata)
    elseif metadata.kind == :vector
        return _lift_laurent_vector(polynomial_object, metadata)
    end
    throw(ArgumentError("unsupported Laurent normalization kind $(metadata.kind)"))
end

function lift_laurent_normalization(normalization)
    hasproperty(normalization, :normalized_object) || throw(ArgumentError("normalization missing normalized_object"))
    hasproperty(normalization, :metadata) || throw(ArgumentError("normalization missing metadata"))
    return lift_laurent_normalization(normalization.normalized_object, normalization.metadata)
end

function verify_laurent_normalization(original, normalization)::Bool
    try
        return lift_laurent_normalization(normalization) == original
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end
```

- [ ] **Step 4: Run the focused test to verify it passes**

Run:

```bash
julia --project=. -e 'include("test/internal/laurent_normalization.jl")'
```

Expected: PASS, including exact lift-back for the negative column and ToricBuilder matrix fixtures.

- [ ] **Step 5: Run the default package tests**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: PASS for the default `public` and `internal` groups.

- [ ] **Step 6: Commit**

Run:

```bash
git add src/core/polynomials.jl src/Suslin.jl test/internal/laurent_normalization.jl test/runtests.jl test/public/api_surface.jl
git commit -m "feat: add laurent normalization lift-back helpers"
```

Expected: commit succeeds and `Manifest.toml` is not staged.

---

## Plan Self-Review

- Spec coverage: the task covers helper API, metadata, lift-back, mixed-ring rejection, tampered metadata, fixture reuse, test registration, exports, and required focused/default/full verification commands.
- Placeholder scan: no incomplete-marker patterns remain.
- Type consistency: the produced API names and metadata fields match between tests, implementation, and public API checks.
