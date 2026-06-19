# Issue 8 Laurent Normalization Design

## Context

Suslin now has exact Laurent ring constructors, shared Laurent fixtures, and a
documented test command contract. Later Groebner, normal-form, and module
workflows need to move Laurent vectors and matrices into ordinary polynomial
rings by clearing negative exponents, while retaining enough metadata to lift
the result back to the original Laurent parent exactly.

Issue #8 depends on #6, #7, and #21. The shared fixture catalog already
contains a small negative-exponent column and a ToricBuilder-derived Laurent
matrix with negative exponents, so the normalization tests should reuse those
fixtures instead of inventing unrelated examples.

## Design Choice

Add a small public helper layer in `src/core/polynomials.jl`:

- `normalize_laurent_object(obj)`: accepts an Oscar matrix or a Julia vector of
  Laurent elements, returns an ordinary polynomial representation and metadata.
- `lift_laurent_normalization(normalization)`: lifts the normalized object back
  using its metadata.
- `lift_laurent_normalization(polynomial_object, metadata)`: lifts an explicit
  ordinary object and metadata pair.
- `verify_laurent_normalization(original, normalization)`: checks exact
  reconstruction and returns `false` for invalid or tampered metadata.

The helper names are exported because the issue defines an input/output
interface that later issues should be able to call directly. The names remain
Laurent-specific to avoid colliding with the existing determinant-only
`normalize_laurent_gl_matrix` boundary.

Alternatives considered:

- Keep the helpers private. This would make future issues rely on
  underscore-prefixed internals despite the issue asking for a reusable
  interface.
- Apply one global monomial shift to every matrix entry. That clears negative
  exponents, but it records less useful determinant metadata for square
  matrices.
- Shift each matrix entry independently. That minimizes each polynomial entry,
  but it does not correspond to multiplying by a manageable monomial matrix and
  is less useful for later determinant/unit accounting.

## Normalization Semantics

For a column vector, compute one componentwise monomial shift from all terms in
the column. The ordinary polynomial column is `m * v`, and lift-back returns
`m^-1 * normalized`.

For a matrix, compute shifts column by column. Each normalized column is
`m_j * A[:, j]`, where `m_j` is the smallest nonnegative monomial that clears
negative exponents in that column. Lift-back divides each column by the same
recorded monomial. For a square matrix, metadata records
`determinant_shift_exponents = sum(column_shifts)`, since
`det(normalized) = prod(m_j) * det(original)`.

Zero entries do not force a shift. Columns with no negative exponents record
the zero exponent vector and are unchanged.

## Metadata Contract

Normalization returns:

- `normalized_object`: an Oscar matrix over an ordinary polynomial ring for
  matrix inputs, or a Julia vector of ordinary polynomial elements for vector
  inputs.
- `metadata`: a named tuple containing `kind`, `shape`, `laurent_ring`,
  `polynomial_ring`, `variable_names`, `column_shifts`, `shift_monomials`,
  `inverse_shift_monomials`, and `determinant_shift_exponents`.

Lift-back validates that the ordinary object shape and parent ring match the
metadata before reconstructing Laurent entries. It also checks that recorded
monomials agree with the recorded exponent shifts. Tampering with shift
metadata must therefore either make lift-back throw or reconstruct a different
object; `verify_laurent_normalization` catches those failures and returns
`false`.

Inputs whose entries do not all belong to the same Laurent parent are rejected
with `ArgumentError`. Matrix inputs use their `base_ring` and vector inputs use
the existing `_require_same_laurent_parent` validator.

## Files

- Modify `src/core/polynomials.jl`: add Laurent-to-polynomial normalization,
  metadata construction, lift-back, and verification helpers.
- Modify `src/Suslin.jl`: export the three helper functions.
- Create `test/internal/laurent_normalization.jl`: focused fixture-based tests
  for a negative column, a ToricBuilder-derived matrix, mixed-ring rejection,
  and tampered shift metadata.
- Modify `test/runtests.jl`: include the focused test in the default internal
  group.
- Modify `test/public/api_surface.jl`: cover the new exported helper names.

## Verification

Issue-specific focused command:

```bash
julia --project=. -e 'include("test/internal/laurent_normalization.jl")'
```

Package entry point required by Agent Desk:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Documented full suite from #21:

```bash
julia --project=. test/runtests.jl all
```

## Automatic Approval

This Agent Desk run is non-interactive. Under the standing answer policy, the
design is approved automatically because it selects the narrowest helper layer
that satisfies issue #8, reuses #6 fixtures, preserves determinant shift
metadata for later accounting, and avoids unrelated factorization changes.

## Spec Self-Review

- No incomplete markers remain.
- The scope is limited to Laurent normalization, metadata, lift-back,
  verification, exports, and focused tests.
- Matrix determinant metadata follows directly from the chosen column-shift
  normalization.
- Negative controls cover mixed Laurent parents and shift-metadata tampering.
