# Issue 70 SL3 Local Certificates Design

## Context

`realize_sl3_local` recognizes local `SL_3` special-form matrices

```text
[p q 0; r s 0; 0 0 1]
```

and returns elementary factors for the currently supported open-slice and
unit-pivot families. Issue 70 needs an expert path that can replay the same
realization and inspect the algebraic witness used by the branch, without
changing the legacy public behavior where `realize_sl3_local(...)` returns only
the factor sequence.

Issue 69 is already merged and supplies the Murthy-Gupta fixture style:
metadata is useful only when a verifier checks the claimed equality exactly.
The issue 70 guardrail repeats that constraint for local `SL_3` certificates:
keep the layer thin and internal, and reject branch tags or witness fields that
do not participate in replay.

## Design Choice

Add a non-exported expert certificate path in `src/algorithm/sl3_local.jl`:

- `SL3LocalRealizationCertificate`
- `realize_sl3_local_certificate(A, X; check_monic=true)`
- `realize_sl3_local_certificate(p, q, r, s, X; check_monic=true)`
- `verify_sl3_local_realization(cert)::Bool`

The names are available as `Suslin.<name>` for expert tests, but they are not
exported through `src/Suslin.jl`. This keeps the API surface conservative while
giving later Murthy-Gupta issues a stable internal result shape to extend.

Alternatives considered:

- Export a public certificate schema. This creates compatibility pressure before
  the Murthy branches exist, and the issue asks for the thinnest useful layer.
- Change `realize_sl3_local` to return a rich result. This breaks existing
  callers and violates the objective.
- Add a parallel certificate solver. This would duplicate recognition and risk
  drift from the factor-returning implementation.

## Certificate Shape

The certificate records only replayed fields:

- `target`: the exact `3 x 3` matrix being realized.
- `branch`: the recognized branch tag, one of `:open_s_one`, `:open_p_one`,
  `:s_unit`, or `:p_unit`.
- `factors`: the exact factor sequence returned by the branch.
- `selected_variable`: the generator passed as `X`.
- `witness`: a branch-specific named tuple whose fields are all checked by
  verification.

Current witness contracts:

- `:open_s_one`: `q` and `r`, checked against `s == 1`, `p == 1 + q*r`, and
  the expected `[E12(q), E21(r)]` factor sequence.
- `:open_p_one`: `q` and `r`, checked against `p == 1`, `s == 1 + q*r`, and
  the expected `[E21(r), E12(q)]` factor sequence.
- `:s_unit`: `pivot = s` and `pivot_inverse`, checked against
  `s*pivot_inverse == 1` and the exact unit-pivot factor sequence.
- `:p_unit`: `pivot = p` and `pivot_inverse`, checked against
  `p*pivot_inverse == 1` and the exact unit-pivot factor sequence.

The branch tag is not trusted by itself. Replay first reconstructs the target
entries from the stored target matrix, checks the selected variable and branch
witness relation, recomputes the expected factors for that branch, and then
checks exact multiplication.

## Data Flow

Recognition remains the single entry point:

1. `realize_sl3_local_certificate` calls the existing recognition helper for
   matrix or `(p, q, r, s, X)` inputs.
2. `_realize_sl3_local_certificate_form(form)` computes factors with the same
   branch logic as the legacy path and builds a witness from the recognized
   form.
3. The constructor path immediately verifies the new certificate. Internal
   verification failure is an error, matching the existing exact-factorization
   guard.
4. `realize_sl3_local` delegates to the certificate path and returns
   `certificate.factors`, preserving legacy behavior.

## Error Handling

Recognition errors remain unchanged. A verifier returns `false` for malformed
or tampered certificates, except `InterruptException` is rethrown. The internal
constructor path throws if it builds a certificate that does not replay.

Unknown branch tags fail verification. Witness fields are accessed explicitly
per branch; missing or corrupted witness fields fail replay instead of being
ignored.

## Tests

Add `test/expert/sl3_local_certificate.jl` and register it in the expert group.
The focused test builds certificates for:

- an open `s == 1` slice,
- an open `p == 1` slice,
- an `s`-unit pivot,
- a `p`-unit pivot.

For each certificate the test checks:

- the factor product equals `cert.target`,
- `verify_factorization(cert.target, cert.factors) == true`,
- `verify_sl3_local_realization(cert) == true`.

The test also checks that legacy `realize_sl3_local(p, q, r, s, X)` still
returns factors with the same multiplication behavior as before.

Negative controls mutate one factor and one witness field and prove replay
rejects each tampered certificate.

## Verification

Focused expert command:

```bash
julia --project=. -e 'include("test/expert/sl3_local_certificate.jl")'
```

Package command required by the Agent Desk run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Full expert-inclusive command:

```bash
julia --project=. test/runtests.jl all
```

## Spec Self-Review

- No incomplete placeholders remain.
- The design preserves the legacy factor-returning API.
- Every stored field is consumed by `verify_sl3_local_realization`.
- The implementation is limited to `src/algorithm/sl3_local.jl`, expert tests,
  and expert test registration.
