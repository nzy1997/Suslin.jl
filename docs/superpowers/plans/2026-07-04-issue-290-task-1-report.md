# Task 1 Report: Certificate Replay And Verification

## Scope

Implemented Task 1 from `docs/superpowers/plans/2026-07-04-issue-290-steinberg-optimization-certificate.md`:

- added internal `SteinbergOptimizationCertificate`
- added `_steinberg_optimization_certificate`
- added `_verify_steinberg_optimization_certificate`
- included `src/algorithm/redundancy.jl` from `src/Suslin.jl`
- extended `test/expert/steinberg_factor_count_optimization.jl` with no-op and negative-control certificate tests

No public exports were added. Factor shortening was not implemented.

## TDD Evidence

### RED

After adding the new expert testset and before implementing the certificate support, I ran:

```bash
julia --project=. -e 'include("test/expert/steinberg_factor_count_optimization.jl")'
```

Observed failure:

```text
Test Summary:                                 | Pass  Total  Time
Steinberg canonical elementary factor records |   15     15  0.4s
Steinberg optimization certificate replay: Error During Test at .../test/expert/steinberg_factor_count_optimization.jl:48
  Got exception outside of a @test
  UndefVarError: `_steinberg_optimization_certificate` not defined in `Suslin`
...
ERROR: LoadError: Some tests did not pass: 0 passed, 0 failed, 1 errored, 0 broken.
```

This matched the plan’s expected missing-symbol failure and proved the new test was covering unimplemented behavior.

### GREEN

After implementing `src/algorithm/redundancy.jl` and wiring the include, I reran:

```bash
julia --project=. -e 'include("test/expert/steinberg_factor_count_optimization.jl")'
```

Observed success:

```text
Test Summary:                                 | Pass  Total  Time
Steinberg canonical elementary factor records |   15     15  0.5s
Test Summary:                             | Pass  Total  Time
Steinberg optimization certificate replay |   22     22  0.7s
```

## Implementation Notes

### New internal module file

Created `src/algorithm/redundancy.jl` with:

- `SteinbergOptimizationCertificate`
- ordinary-polynomial ring validation
- factor sequence normalization and validation
- exact factor-product replay
- rewrite-record normalization and span/count checks
- comparison summary generation using existing factor metrics
- exact certificate verification with stored-summary consistency checks

### Integration

Added:

```julia
include("algorithm/redundancy.jl")
```

to `src/Suslin.jl` immediately after `include("algorithm/quillen_induction.jl")`.

### Tests added

Extended `test/expert/steinberg_factor_count_optimization.jl` to cover:

- no-op certificate construction on catalog case `steinberg-same-position-merge-qq`
- exact replay verification success
- summary metric consistency
- tampered optimized factor rejection
- stale rewrite-log delta rejection
- mixed-ring rejection with `ArgumentError`

## Full Verification

Ran:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Observed success:

```text
Test Summary: | Pass  Total     Time
public        |  725    725  2m45.5s
Test Summary: | Pass  Total     Time
internal      |  720    720  4m56.7s
Testing Suslin tests passed
```

The full run also emitted pre-existing Julia 1.12 world-age warnings in Quillen fixture tests, but the suite completed with exit code `0`.

## Self-Review

### Correctness

- certificate verification replays exact matrix products instead of relying on counts alone
- rewrite-log validation checks factor-count deltas against normalized metadata
- mixed base-ring inputs are rejected deliberately
- certificate data remains internal to the module and accessible only through `Suslin.<name>`

### Scope adherence

- no public API export changes
- no factor-shortening logic added
- only the required source and test surfaces were touched, plus this report

### Changed files

- `src/algorithm/redundancy.jl`
- `src/Suslin.jl`
- `test/expert/steinberg_factor_count_optimization.jl`
- `docs/superpowers/plans/2026-07-04-issue-290-task-1-report.md`

## Task 1 Review Fix: PolyRing Acceptance

Addressed the Task 1 review finding that `_require_steinberg_ordinary_polynomial_ring` accepted only `MPolyRing` and rejected supported exact ordinary univariate `PolyRing`.

### TDD Evidence

#### RED

After adding a univariate no-op certificate regression over `Oscar.polynomial_ring(QQ, "x")`, I ran:

```bash
julia --project=. -e 'include("test/expert/steinberg_factor_count_optimization.jl")'
```

Observed expected failure:

```text
Steinberg optimization certificate accepts univariate ordinary polynomial rings: Error During Test
  ArgumentError: Steinberg optimization certificates require an ordinary polynomial ring
  Stacktrace:
    [1] _require_steinberg_ordinary_polynomial_ring(R::QQPolyRing)
```

This confirmed the test was exercising the missing `PolyRing` support rather than an unrelated issue.

#### GREEN

After updating `src/algorithm/redundancy.jl` to accept `PolyRing` alongside `MPolyRing`, I reran:

```bash
julia --project=. -e 'include("test/expert/steinberg_factor_count_optimization.jl")'
```

Observed success:

```text
Test Summary:                                 | Pass  Total  Time
Steinberg canonical elementary factor records |   15     15  0.5s
Test Summary:                             | Pass  Total  Time
Steinberg optimization certificate replay |   22     22  0.7s
Test Summary:                                                                   | Pass  Total  Time
Steinberg optimization certificate accepts univariate ordinary polynomial rings |    7      7  0.1s
```

### Scope Notes

- kept the construction internal; no public exports changed
- added only no-op certificate coverage for the supported univariate ordinary polynomial ring case
- did not implement any shortening rules
