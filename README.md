# SuslinStability

Constructive elementary-matrix factorizations for small supported `SL_3` slices over polynomial rings.

## Example

```julia
using SuslinStability, Oscar

R, (X,) = Oscar.polynomial_ring(QQ, ["X"])
A = matrix(R, [
    one(R)      one(R) + X      zero(R);
    X           one(R) + X + X^2 zero(R);
    zero(R)     zero(R)         one(R)
])

factors = elementary_factorization(A)
verify_factorization(A, factors)
```

## Current scope

- `elementary_factorization(A)` currently supports only a narrow 3x3 univariate `SL_3` slice.
- `verify_factorization(A, factors)` checks exact multiplication against `A`.
- The implementation is staged; it is not yet the full Park-Woodburn algorithm.

## References

- Park and Woodburn, *An algorithmic proof of Suslin's stability theorem for polynomial rings*.
- Logar and Sturmfels, *Algorithms for the Quillen-Suslin theorem*.
- Fitchas and Galligo, *Nullstellensatz effectif et conjecture de Serre*.
