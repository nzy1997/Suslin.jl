# Suslin

[![CI](https://github.com/nzy1997/Suslin.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/nzy1997/Suslin.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/nzy1997/Suslin.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/nzy1997/Suslin.jl)

Constructive elementary-matrix factorizations for small supported `SL_3` slices over polynomial rings.

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

## Current scope

- `elementary_factorization(A)` currently supports only a narrow 3x3 univariate `SL_3` slice.
- `verify_factorization(A, factors)` checks exact multiplication against `A`.
- The implementation is staged; it is not yet the full Park-Woodburn algorithm.

## Testing

This repository does not commit a `Manifest.toml`. In a fresh checkout,
instantiate dependencies before running tests:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

The test runner separates routine checks from expert algorithm checks:

| Command | Coverage |
| --- | --- |
| `julia --project=. test/runtests.jl` | Default fast tests: `public` and `internal` groups |
| `julia --project=. test/runtests.jl expert` | Expert-only algorithm and documentation checks |
| `julia --project=. test/runtests.jl all` | Full suite: `public`, `internal`, and `expert` groups |
| `julia --project=. -e 'using Pkg; Pkg.test()'` | Package test entry point; runs the default fast tests |

CI uploads full-suite coverage to Codecov for the `main` branch and pull
requests. Private repositories need a `CODECOV_TOKEN` GitHub Actions secret for
coverage uploads.

## References

- Park and Woodburn, *An algorithmic proof of Suslin's stability theorem for polynomial rings*.
- Logar and Sturmfels, *Algorithms for the Quillen-Suslin theorem*.
- Fitchas and Galligo, *Nullstellensatz effectif et conjecture de Serre*.
