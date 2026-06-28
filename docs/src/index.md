```@meta
CurrentModule = Suslin
```

# Suslin

Constructive elementary-matrix factorizations for staged supported slices over polynomial and Laurent polynomial rings.

## Example

```julia
using Suslin, Oscar

R, (X,) = Oscar.polynomial_ring(QQ, ["X"])
A = matrix(R, [
    one(R)      one(R) + X      zero(R);
    X           one(R) + X + X^2 zero(R);
    zero(R)     zero(R)         one(R)
])

factors = elementary_factorization(A)
verify_factorization(A, factors)
```

## Scope

- `elementary_factorization(A)` is staged. It supports univariate local `SL_3`
  ordinary-polynomial matrices, selected `n > 3` ordinary-polynomial matrices
  through block-local reduction and recursive polynomial column peel,
  deterministic multivariate Quillen/local-to-global fixture-backed matrices,
  and determinant-one Laurent inputs through the existing Laurent `SL` path.
- `verify_factorization(A, factors)` checks exact multiplication against `A`.
- `laurent_gl_factorization_certificate(A)` defaults to the eager Laurent
  normalization/core certificate. With `determinant_strategy = :lazy`, it
  records the supported monomial-unit deferred determinant correction path for
  original Laurent `GL_n` inputs such as the issue #38 fixture.
- The implementation is not yet the full Park-Woodburn algorithm for arbitrary
  `SL_n(k[x_1, ..., x_m])`, `n >= 3`. General Quillen local realizability,
  coefficient-ring support beyond exact field-backed ordinary polynomial
  rings, arbitrary Laurent `GL_n` determinant correction, and factor-count
  optimization remain staged boundaries.

See [ToricBuilder Integration Contract](@ref) for the first recorded
consumer-boundary fixture contract.

## References

- Park and Woodburn, *An algorithmic proof of Suslin's stability theorem for polynomial rings*.
- Logar and Sturmfels, *Algorithms for the Quillen-Suslin theorem*.
- Fitchas and Galligo, *Nullstellensatz effectif et conjecture de Serre*.

```@index
```

```@autodocs
Modules = [Suslin]
```
