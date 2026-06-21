# Suslin

[![CI](https://github.com/nzy1997/Suslin.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/nzy1997/Suslin.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/nzy1997/Suslin.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/nzy1997/Suslin.jl)

Constructive elementary-matrix factorizations for staged Suslin stability
experiments over polynomial and Laurent polynomial rings.

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

- `elementary_factorization(A)` is staged. It supports the local embedded
  `SL_3` families used by the current tests, block-local reductions for
  selected `SL_n` inputs, and determinant-one Laurent inputs handled by the
  recursive Laurent column-peel path.
- `verify_factorization(A, factors)` checks exact multiplication against `A`.
- `laurent_gl_factorization_certificate(A)` records Laurent normalization and
  determinant-one core factorization data. General Laurent `GL_n` determinant
  correction is still a staged boundary.
- The implementation is not yet the full Park-Woodburn algorithm for arbitrary
  `SL_n(k[x_1, ..., x_m])`, `n >= 3`.

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

## Related Implementations

No maintained public implementation of the full Park-Woodburn elementary
factorization algorithm was located in this search pass (2026-06-21). The
closest public implementations found are for the neighboring
Quillen-Suslin/unimodular-completion problem:

- [M2-QuillenSuslin](https://github.com/bbarwick/M2-QuillenSuslin) is the
  Macaulay2 package source by Brett Barwick and Branden Stone. It computes free
  bases for projective modules and completions of unimodular matrices over
  polynomial and Laurent polynomial rings. See also the package paper:
  [arXiv:1107.4383](https://arxiv.org/abs/1107.4383),
  [local PDF](refs/barwick-stone-2011-quillensuslin-macaulay2.pdf), and
  [JSAG DOI](https://doi.org/10.2140/jsag.2013.5.26).
- [QuillenSuslin-M2](https://github.com/bstone/QuillenSuslin-M2) is an older
  Macaulay2 repository/mirror. Its README says the package follows
  Logar-Sturmfels and Fabianska-Quadrat for free bases and unimodular
  completion; it is related background rather than a direct Park-Woodburn
  `SL_n` elementary factorization implementation.

## References

Local open PDFs are kept in [`refs/`](refs/README.md).

Primary algorithmic references:

- H. Park and C. Woodburn, "An algorithmic proof of Suslin's stability theorem
  for polynomial rings," *Journal of Algebra* 178 (1995), 277-298.
  [arXiv](https://arxiv.org/abs/alg-geom/9405003),
  [local PDF](refs/park-woodburn-1994-suslin-stability.pdf),
  [DOI](https://doi.org/10.1006/jabr.1995.1349).
- A. Logar and B. Sturmfels, "Algorithms for the Quillen-Suslin theorem,"
  *Journal of Algebra* 145 (1992), 231-239.
  [DOI](https://doi.org/10.1016/0021-8693(92)90189-s).
- N. Fitchas and A. Galligo, "Nullstellensatz effectif et Conjecture de Serre
  (Theoreme de Quillen-Suslin) pour le Calcul Formel,"
  *Mathematische Nachrichten* 149 (1990), 231-253.
  [DOI](https://doi.org/10.1002/mana.19901490118).
- L. Caniglia, G. Cortinas, S. Danon, J. Heintz, T. Krick, and P. Solerno,
  "Algorithmic aspects of Suslin's proof of Serre's conjecture,"
  *Computational Complexity* 3 (1993), 31-55.
  [DOI](https://doi.org/10.1007/bf01200406).

Related constructive and Laurent-polynomial background:

- H. Park, "Symbolic computations and signal processing,"
  *Journal of Symbolic Computation* 37 (2004), 209-226.
  [DOI](https://doi.org/10.1016/j.jsc.2002.06.003).
- B. Barwick and B. Stone, "The QuillenSuslin Package for Macaulay2."
  [arXiv](https://arxiv.org/abs/1107.4383),
  [local PDF](refs/barwick-stone-2011-quillensuslin-macaulay2.pdf).
  See also the related JSAG package article, "Computing free bases for
  projective modules," *Journal of Software for Algebra and Geometry* 5
  (2013), 26-32. [DOI](https://doi.org/10.2140/jsag.2013.5.26).
- A. Fabianska and A. Quadrat, "Applications of the Quillen-Suslin theorem to
  multidimensional systems theory," in *Groebner Bases in Control Theory and
  Signal Processing* (2007), 23-106.
  [DOI](https://doi.org/10.1515/9783110909746.23).
- P. M. Cohn, "On the structure of the `GL_2` of a ring,"
  *Publications Mathematiques de l'IHES* 30 (1966), 5-54.
  [DOI](https://doi.org/10.1007/BF02684355).

Classical theorem references:

- D. Quillen, "Projective modules over polynomial rings,"
  *Inventiones Mathematicae* 36 (1976), 167-171.
  [DOI](https://doi.org/10.1007/bf01390008).
- A. A. Suslin, "Projective modules over a polynomial ring are free,"
  *Soviet Mathematics Doklady* 17 (1976), 1160-1164.
- A. A. Suslin, "On the structure of the special linear group over polynomial
  rings," *Mathematics of the USSR-Izvestiya* 11 (1977), 221-238.
  [DOI](https://doi.org/10.1070/IM1977v011n02ABEH001709).
- T. Y. Lam, *Serre's Conjecture*, Lecture Notes in Mathematics 635,
  Springer, 1978. [DOI](https://doi.org/10.1007/BFb0068340).
