# Issue 160 Lazy Laurent Certificate API Design

## Context

Issue #159 added the internal lazy Laurent determinant hoist certificate with
row/left and column/right correction choices. The public
`laurent_gl_factorization_certificate(A)` entry point still only exposes the
eager Laurent normalization certificate, and public callers cannot request the
lazy route or read the lazy determinant metadata through a stable certificate
contract.

Issue #160 asks for that public contract while preserving the existing
`elementary_factorization(A)` boundary: original Laurent `GL_n` inputs must
still throw the staged Laurent boundary instead of returning elementary factors
for only a normalized core.

## Approach Options

Recommended: keep the existing one-argument public call on the eager
normalization path and add explicit keyword routing:
`laurent_gl_factorization_certificate(A; determinant_strategy = :lazy,
correction_side = :row)`. The eager default preserves compatibility. The lazy
route returns the existing lazy hoist certificate type, promoted to public API
and extended with a direct `determinant_source` field. Public verification
dispatches to the lazy verifier for lazy certificates.

Alternative: make lazy determinant correction the default and add
`determinant_strategy = :eager` for regression comparison. This is a larger
behavior change for existing callers and is unnecessary for the issue's stated
lazy opt-in interface.

Alternative: keep the lazy route internal and expose only accessor functions
for selected metadata. This avoids exporting another certificate type, but it
does not give callers a concrete public certificate contract with direct fields
and would leave the public constructor unable to choose row or column
correction.

## Chosen Design

The public constructor keeps this compatibility rule:

```julia
laurent_gl_factorization_certificate(A)
```

continues to return `LaurentGLFactorizationCertificate` from the eager
normalization route.

The new lazy public route is:

```julia
laurent_gl_factorization_certificate(
    A;
    determinant_strategy = :lazy,
    correction_side = :row,
)
```

with `correction_side = :column` supported by the #159 column/right algebra.
`determinant_strategy = :eager` explicitly selects the old route. Unsupported
strategy values throw `ArgumentError` naming `:eager` and `:lazy`.
`correction_side` is accepted only when `determinant_strategy = :lazy`; using it
without the lazy strategy throws a specific `ArgumentError` instead of silently
ignoring the caller's option.

`LaurentLazyGLHoistCertificate` becomes part of the exported certificate API.
It gains a direct field:

```julia
determinant_source::Symbol
```

set from the verified deferred metadata. Lazy certificates already expose
`overall_determinant` and `correction_side`; after this change public callers
can assert:

```julia
cert.overall_determinant == det(A)
cert.determinant_source == :deferred_submatrix
cert.correction_side == :row # or :column
verify_laurent_gl_factorization_certificate(cert)
cert.reconstructed_product == A
```

The lazy verifier checks the direct `determinant_source` field against the
verified metadata so a certificate with mismatched public determinant source
does not verify.

## Testing

Add `test/public/laurent_gl_certificate_options.jl` and register it in the
public test group. The focused public test uses the issue #38 fixture and
checks:

- the existing one-argument call remains valid and verifies;
- the explicit lazy row route reports `overall_determinant == det(Q)`,
  `determinant_source == :deferred_submatrix`, `correction_side == :row`, exact
  reconstruction, and public verification success;
- the explicit lazy column route verifies and reports
  `correction_side == :column`;
- invalid strategy and invalid correction-side combinations are rejected with
  `ArgumentError`;
- `elementary_factorization(Q)` for the original Laurent `GL_n` input still
  throws the staged Laurent `GL_n` boundary.

Update `test/public/api_surface.jl` to assert the exported lazy certificate
type and to exercise the existing one-argument call plus accepted/rejected
keyword spellings on a small Laurent determinant-one matrix.

Focused verification command:

```bash
julia --project=. -e 'include("test/public/laurent_gl_certificate_options.jl")'
```

Package verification command:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Scope

Do not change `elementary_factorization(A)` to accept original Laurent `GL_n`
inputs. Do not update ToricBuilder Q-block reporting. Do not commit
`Manifest.toml`.

## Automatic Decisions

- Clarifying questions, design approval, and written spec review were resolved
  by the Standing Answer Policy because this is a non-interactive Agent Desk
  run.
- The visual companion was skipped because no visual question would clarify
  this algebraic API change.
- The explicit lazy strategy with eager default was selected because it is the
  conservative compatibility-preserving option.
- `LaurentLazyGLHoistCertificate` is exported because the issue asks for a
  public certificate carrying lazy determinant fields.
