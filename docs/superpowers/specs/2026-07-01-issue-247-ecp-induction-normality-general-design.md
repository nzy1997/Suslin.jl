# Issue 247 ECP Induction Normality General Design

## Context

Issue #247 builds on the merged #246 link-step route, #245 link-witness
extraction, and #193/#181 conjugated-elementary normality certificates. The
current `ecp_induction_normality_certificate` can replay an explicit
lower-variable reduction and an explicit normality witness, but #185 needs the
same stage to assemble those pieces from a verified link step when possible.

There is no `AGENTS.md` in this checkout. Live GitHub issue/PR metadata could
not be fetched because the sandbox blocks the configured GitHub proxy; the
issue body supplied by Agent Desk plus the checked-in #246/#193/#194 workflow
docs and merged code are the binding context.

## Approaches Considered

1. Extend the existing `ecp_induction_normality_certificate` constructor. This
   is the chosen approach. It preserves the current internal API, keeps replay
   verification in one place, and lets the #246 link-step certificate drive
   lower-variable induction and normality construction.
2. Add a separate `ecp_induction_normality_general_certificate` wrapper. This
   would avoid touching the replayed constructor, but it would split verifier
   behavior and make #248 diagnostics choose between two certificate shapes.
3. Change the public reducer dispatch to call the recursive route directly.
   This is out of scope for #247 and would couple this proof-composition work
   to #186/#185 reducer policy.

## Design

Add an explicit `descent_measure` field to
`ECPInductionNormalityCertificate`. The measure records the selected variable
from the link witness, the parent selected-variable degree profile, the
lower-variable degree profile, the column length, the variable count, and a
`strict_descent` flag. Construction rejects a same-context recursive call with
an `ArgumentError` staged-failure message when the lower profile is not
componentwise no larger and strictly smaller in at least one entry. The verifier
recomputes the measure from the stored link step and requires exact equality
and strict descent.

Lower reductions remain accepted as either absent, an
`ECPColumnReductionCertificate`, or a factor-list hint. The constructor always
stores a verified lower ECP certificate. For a factor-list hint, it verifies the
sequence is elementary, reduces the lower column, and matches the canonical
`ecp_column_reduction_certificate(lower_column, R)` factors before storing that
canonical certificate. If no lower certificate can be constructed, the staged
failure message starts with `missing lower-variable reduction`.

Normality data becomes optional. When no witness is supplied, the constructor
derives one from the verified lower reduction: the conjugator is the inverse of
the lifted lower factor product, the fixed coordinate is `e_n`, the moving
coordinate is `1`, and the nonzero elementary entry is `selected_variable + 1`.
The rewrite then routes through `realize_conjugate_elementary_certificate`; a
supplied normality certificate is still accepted only if
`verify_conjugate_elementary_certificate` succeeds and it matches the derived
certificate. If normality replay cannot be produced, the staged failure message
starts with `missing normality rewrite`.

The final factor sequence remains:

```julia
lifted_lower_variable_factors
normality_rewrite.rewrite_factors
link_step.reduction_factors
```

The verifier independently checks the link step, descent measure, lower ECP
certificate, lifted lower factors, nested #193/#181 normality certificate,
stored normality rewrite, final factor sequence, elementary-factor shape, and
exact reduction of the original column to `e_n`.

## Tests

Add `test/expert/ecp_induction_normality_general.jl` and register it in the
expert group. The focused test builds the #242 `ecp-mainline-sl3-route-qq`
case through `route_mode = :polynomial_sl3`, then calls
`ecp_induction_normality_certificate` with no lower reduction and no normality
witness. Assertions cover:

- strict selected-variable descent from `v(X)` to `v(0)`;
- an actual stored `ECPColumnReductionCertificate` for the lower-variable
  column;
- automatically constructed normality witness data;
- independently verified `ConjugatedElementaryNormalityCertificate` data;
- exact final reduction of the original column to `e_n`;
- negative controls for corrupted descent metadata, lower reduction, lifted
  factor list, normality witness, nested normality certificate, and one final
  factor; and
- staged-failure messages that distinguish same-context descent failure,
  missing lower-variable reduction, and missing normality rewrite.

Keep `test/expert/normality.jl` as the independent #181/#193 normality
verification suite. Required verification commands:

```bash
julia --project=. -e 'include("test/expert/ecp_induction_normality_general.jl")'
julia --project=. -e 'include("test/expert/normality.jl")'
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Out Of Scope

This change does not alter public reducer dispatch, implement recursive matrix
factorization for #186, broaden Laurent/ToricBuilder support, or teach the
link-step route to factor new classes of `SL_n` matrices.

## Automatic Approval

This Agent Desk run is non-interactive. Under the Standing Answer Policy, the
design is approved automatically because it follows the existing internal API,
uses verified certificates as authorities, records descent to bound recursion,
and fails with explicit staged messages instead of guessing.

## Self-Review

- No incomplete markers remain.
- The design covers every hard acceptance requirement from #247.
- Same-context recursion is rejected before attempting lower reduction.
- The normality witness is constructed from the verified lower reduction and
  checked through the #193/#181 certificate API.
- The change stays scoped to `column_reduction.jl`, focused expert tests, and
  test registration.
