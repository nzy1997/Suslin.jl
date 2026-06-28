# Elementary Factor Analysis Design

## Goal

Add two public analysis helpers for elementary factor sequences:

- the largest monomial exponent weight appearing in any elementary factor;
- the total number of monomial terms appearing in non-diagonal entries across
  all elementary factors.

For a Laurent monomial such as `x^5*y^-7`, the exponent weight is
`abs(5) + abs(-7) == 12`.

## Context

`elementary_factorization(A)` returns a `Vector` of elementary matrices. Several
certificate paths also expose already-computed factor vectors, including
polynomial route certificates, Laurent column-peel certificates, and Laurent GL
certificates through `core_factors`.

The analysis should not recompute a factorization from `A`; it should inspect
the factor sequence the caller chooses to analyze. This keeps the helpers cheap,
unambiguous, and reusable across all current factorization routes.

## Approach Options

Recommended: add two public pure functions that accept a factor sequence:

```julia
max_elementary_factor_monomial_degree(factors)
total_elementary_factor_offdiagonal_monomials(factors)
```

This is the chosen approach because it works uniformly on
`elementary_factorization(A)`, `.factors`, and `.core_factors` without binding
the API to a specific certificate type.

Alternative: add one aggregate stats function returning a named tuple. That is
convenient, but it adds an extra public API before there is a clear need.

Alternative: add helpers that accept `A` and call `elementary_factorization(A)`.
That is short for users, but it hides potentially expensive factorization work
and makes Laurent GL/core-factor boundaries less explicit.

## Chosen Design

Add the two public functions to the elementary-matrix core layer and export
them from `Suslin.jl`.

Both functions accept any iterable of Oscar matrices. They validate each factor
as a square matrix and scan only non-diagonal entries. Diagonal entries are
ignored even if a caller passes a non-elementary matrix; this matches the
requested "non-diagonal elements" analysis.

`max_elementary_factor_monomial_degree(factors)` returns the maximum
`sum(abs, exponent_vector)` among every monomial term in every nonzero
non-diagonal entry. Empty factor sequences, zero non-diagonal entries, and
constant non-diagonal terms all contribute `0`, so the function returns `0` if
no positive exponent weight appears.

`total_elementary_factor_offdiagonal_monomials(factors)` returns the sum of the
number of monomial terms in all nonzero non-diagonal entries. A coefficient-only
nonzero entry counts as one monomial. A zero entry counts as zero.

## Data Flow

Each helper iterates over factors, then over `(row, column)` pairs with
`row != column`.

For each nonzero entry, it uses Oscar/AbstractAlgebra term metadata:

- `exponents(entry)` supplies exponent vectors for Laurent polynomial entries;
- `AbstractAlgebra.exponent_vectors(entry)` is the ordinary polynomial fallback;
- `coefficients(entry)` supplies term counts when needed.

The implementation should keep the term-extraction logic private so both public
helpers share the same interpretation of monomials.

## Error Handling

If a factor is not square, throw `ArgumentError("factor must be square")`.

If a nonzero entry does not expose polynomial/Laurent term metadata, throw an
`ArgumentError` explaining that elementary factor analysis requires polynomial
or Laurent polynomial entries. Interrupts must still be rethrown.

The helpers do not verify that each matrix is elementary. They are analysis
tools for factor sequences and should avoid the extra cost and ambiguity of
proving elementary shape.

## Tests

Add focused tests under `test/expert/elementary_matrices.jl` because the feature
belongs to the elementary-matrix API.

Test cases should cover:

- Laurent terms with negative exponents, including `x^5*y^-7` contributing
  maximum degree `12`;
- a non-diagonal polynomial with multiple monomial terms contributing the
  correct term count;
- diagonal entries ignored for both metrics;
- empty factor sequences returning `0`;
- public API surface exports in `test/public/api_surface.jl`.

## Verification

Run:

```bash
julia --project=. -e 'include("test/expert/elementary_matrices.jl")'
julia --project=. -e 'include("test/public/api_surface.jl")'
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected result: all commands exit 0.

## Out Of Scope

Do not add an aggregate stats function yet. Do not add helpers that accept and
factorize an input matrix. Do not change existing factorization or certificate
logic.
