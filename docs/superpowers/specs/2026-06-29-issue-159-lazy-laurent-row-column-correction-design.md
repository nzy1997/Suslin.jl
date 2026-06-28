# Issue 159 Lazy Laurent Row/Column Correction Design

## Context

Issue #158 added an internal lazy Laurent GL hoist certificate for the row/left
correction path. It consumes determinant-deferred metadata from #157, embeds
the deferred diagonal correction into the original matrix dimension, rewrites
the inverse left peel factors across that diagonal correction, and verifies:

```text
original_matrix == correction.factor * elementary_product
```

Issue #159 extends that internal certificate route so callers can choose
whether the Laurent monomial-unit determinant correction is represented on the
left as a row correction or on the right as a column correction. The same input
must report the same `overall_determinant` for both choices, while the
certificate metadata and reconstruction relation must identify the chosen side.

## Approach Options

Recommended: add a narrow side-aware layer inside the lazy hoist certificate
constructor. It accepts `correction_side = :row` or `:column`, normalizes those
to internal `:left` and `:right` sides, builds the matching deferred
determinant-one core, hoists the chosen correction to original dimension, and
verifies the side-specific reconstruction relation. This keeps the feature at
the certificate internals requested by the issue and reuses the existing #158
left-side algebra.

Alternative: extend `normalize_laurent_gl_matrix` and the public Laurent GL
certificate API with row/column options. This is too broad for #159 because the
issue targets lazy determinant certificate internals and explicitly keeps
general Laurent `GL_n` public factorization out of scope.

Alternative: store a `correction_side` field but keep always using the
left/row reconstruction algebra. This is rejected because the negative control
requires row metadata paired with column reconstruction data to fail
verification.

## Chosen Design

Add an internal option parser:

```text
:row    -> :left
:column -> :right
```

Unsupported values throw `ArgumentError` with a message naming `:row` and
`:column`. The certificate stores both `correction_side` and
`reconstruction_relation` so verification can distinguish the requested side
from the algebra it claims to satisfy.

For row/left, preserve the #158 relation:

```text
D_left * normalized_deferred_core == deferred_submatrix
original_matrix == correction.factor * elementary_product
```

The left-side elementary product remains:

```text
rewrite(left_product^-1 across D_left) *
blockdiag(core factors, I) *
right_product^-1
```

For column/right, build the determinant-one deferred core as:

```text
normalized_deferred_core * D_right == deferred_submatrix
```

where `D_right` is a diagonal correction with determinant equal to the deferred
submatrix determinant. The original reconstruction is:

```text
original_matrix == elementary_product * correction.factor
```

The column-side elementary product moves the diagonal correction past the
inverse right peel factors:

```text
left_product^-1 *
blockdiag(core factors, I) *
rewrite(D_right across right_product^-1)
```

For a diagonal `D = diag(d_i)` and an elementary factor `E_ij(a)`, the new
right-side rewrite uses:

```text
D * E_ij(a) == E_ij(d_i * a * d_j^-1) * D
```

The existing left-side helper still uses:

```text
E_ij(a) * D == D * E_ij(d_i^-1 * a * d_j)
```

## Verification

The lazy hoist verifier recomputes the expected side from
`certificate.correction_side`, rebuilds the expected correction and normalized
core from the verified deferred metadata, and checks:

- `correction_side` is exactly `:row` or `:column`;
- `reconstruction_relation` is `:left_correction_times_elementary_product` for
  row certificates and `:elementary_product_times_right_correction` for column
  certificates;
- the correction and inverse correction are original-dimension matrices over
  the correct ring and multiply to identity in both orders;
- the correction determinant equals `overall_determinant`;
- the normalized deferred core has determinant one and satisfies the
  side-specific deferred relation;
- the normalized deferred core factorization verifies exactly;
- the recorded rewritten left or right factors match the side-specific rewrite;
- `elementary_product == product(elementary_factors)`;
- the side-specific reconstructed product equals the original matrix.

A certificate with row metadata and column reconstruction data must fail
verification because the expected relation and expected factor assembly are
derived from `correction_side`, not trusted from stored products.

## Testing

Add `test/expert/laurent_lazy_row_column_correction.jl` and register it in the
expert group. The focused test covers:

- row and column certificates on the #154 issue #38 lazy determinant fixture;
- both certificates verify exact reconstruction of the original matrix;
- both certificates report `overall_determinant == u*v`;
- `correction_side` differs as requested;
- row reconstruction verifies `correction.factor * elementary_product == A`;
- column reconstruction verifies `elementary_product * correction.factor == A`;
- unsupported `correction_side = :diagonal` throws an `ArgumentError` naming
  accepted `:row` and `:column` options;
- a tampered certificate with row metadata and column reconstruction data fails
  verification.

Focused verification command:

```bash
julia --project=. -e 'include("test/expert/laurent_lazy_row_column_correction.jl")'
```

Package verification command:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Scope

Do not export new APIs. Do not change `elementary_factorization(A)` to accept
Laurent `GL_n` inputs. Do not broaden unsupported determinant classes. Do not
commit `Manifest.toml`.

## Automatic Decisions

- Clarifying questions and design approval were resolved by the Standing Answer
  Policy because this is a non-interactive Agent Desk run.
- The visual companion was skipped because no visual decision would clarify the
  algebraic certificate design.
- The recommended narrow internal side-aware certificate option was selected
  because it directly satisfies #159 and avoids public API expansion.
