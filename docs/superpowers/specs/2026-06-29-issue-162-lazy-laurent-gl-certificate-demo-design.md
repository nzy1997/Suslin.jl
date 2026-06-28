# Issue 162 Lazy Laurent GL Certificate Demo Design

## Context

Issue #38 is the canonical small ToricBuilder Laurent `GL_n` boundary. The
original `Q` matrix has determinant `u*v`; row- and column-normalized cores
factor, while `elementary_factorization(Q)` remains a staged boundary for
returning elementary factors of the original `GL_n` input. Issue #160 added a
public lazy certificate route:

```julia
laurent_gl_factorization_certificate(
    Q;
    determinant_strategy = :lazy,
    correction_side = :row,
)
```

That route returns a `LaurentLazyGLHoistCertificate` with public
`overall_determinant`, `determinant_source`, and `correction_side` fields and
public verification dispatch.

## Goal

Demonstrate the supported lazy determinant Laurent `GL_n` certificate on the
original issue #38 `Q` input, without expanding `elementary_factorization(Q)`
beyond its staged boundary.

## Chosen Approach

Update the existing issue #38 example and expert certificate test in place.
This keeps all reviewer context next to the fixture that already owns the
boundary and certificate checks. Add a concise `lazy_gl_certificate` PASS line
to the example, and extend the expert test to assert exact original-input
reconstruction plus negative controls for the lazy certificate fields.

Alternatives considered:

- Add a separate expert test file for the lazy original-input certificate. This
  would avoid touching existing eager-certificate tests, but it would duplicate
  the issue #38 fixture setup and split a single acceptance story across files.
- Change `elementary_factorization(Q)` to return the lazy certificate factors.
  This is out of scope because issue #162 explicitly preserves the staged
  boundary for original Laurent `GL_n` inputs.

## Example Behavior

`example/toric_decoupling/issue38_unit_correction.jl` will keep the existing
original-boundary, row-normalized, and column-normalized output. It will add one
lazy original-input certificate line:

```text
lazy_gl_certificate status=PASS det=u*v correction_side=row determinant_source=deferred_submatrix factors=<N> verified=true
```

The example will construct the certificate through the public API with
`determinant_strategy = :lazy` and `correction_side = :row`, verify it with
`verify_laurent_gl_factorization_certificate`, and require exact
`reconstructed_product == Q`.

## Tests

`test/expert/issue38_laurent_gl_certificate.jl` will extend its issue #38 test
with a lazy original-input certificate block that asserts:

- the certificate is a `LaurentLazyGLHoistCertificate`;
- `overall_determinant == u*v`;
- `determinant_source == :deferred_submatrix`;
- `correction_side == :row`;
- `elementary_product` equals the product of the recorded hoisted elementary
  factors;
- `correction.factor * elementary_product == Q` for the row correction side;
- `reconstructed_product == Q`;
- public verification returns `true`.

Negative controls will rebuild the immutable lazy certificate with one field
tampered at a time and assert public verification returns `false` for:

- determinant correction tampering;
- reported correction side tampering;
- one hoisted elementary factor tampering.

The same test will keep asserting that `elementary_factorization(Q)` throws the
existing Laurent `GL_n` normalization boundary error.

## Documentation

The README and generated docs index currently describe
`laurent_gl_factorization_certificate(A)` as core-only and say Laurent `GL_n`
determinant correction is still wholly staged. Update that wording narrowly:
the eager default records normalization and determinant-one core factorization;
the lazy strategy records the supported monomial-unit deferred determinant
correction path. Keep arbitrary Laurent `GL_n` determinant correction out of
scope.

## Verification

Run the issue commands:

```bash
julia --project=. example/toric_decoupling/issue38_unit_correction.jl
julia --project=. -e 'include("test/expert/issue38_laurent_gl_certificate.jl")'
```

Run the package gate when applicable:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

The example must exit 0 and print the `lazy_gl_certificate status=PASS ...`
line. The expert test must include exact reconstruction of the original issue
#38 `Q` matrix and the required negative controls.
