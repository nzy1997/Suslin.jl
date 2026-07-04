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

## Task 1 Review Fix: Impossible Same-Delta Rewrite Log Counts

Addressed the final whole-branch review finding for issue #290: certificate verification previously accepted a coherently forged rewrite log when its aggregate factor-count delta matched the sequence delta, even if the logged counts were impossible for the actual factor sequences.

### TDD Evidence

#### RED

After adding a focused expert regression for a same-delta but impossible rewrite log, I ran:

```bash
julia --project=. -e 'include("test/expert/steinberg_factor_count_optimization.jl")'
```

Observed expected failure:

```text
Steinberg optimization certificate replay: Test Failed at .../test/expert/steinberg_factor_count_optimization.jl:114
  Expression: !(Suslin._verify_steinberg_optimization_certificate(impossible_same_delta_certificate))
...
Test Summary:                             | Pass  Fail  Total  Time
Steinberg optimization certificate replay |   22     1     23  1.5s
ERROR: LoadError: Some tests did not pass: 22 passed, 1 failed, 0 errored, 0 broken.
```

This showed the regression was real: the verifier still accepted a forged certificate whose rewrite log claimed counts far larger than either sequence length while preserving net delta zero.

#### GREEN

After tightening rewrite-log validation in `src/algorithm/redundancy.jl`, I reran:

```bash
julia --project=. -e 'include("test/expert/steinberg_factor_count_optimization.jl")'
```

Observed success:

```text
Test Summary:                                 | Pass  Total  Time
Steinberg canonical elementary factor records |   15     15  0.5s
Test Summary:                             | Pass  Total  Time
Steinberg optimization certificate replay |   23     23  0.9s
Test Summary:                                                                   | Pass  Total  Time
Steinberg optimization certificate accepts univariate ordinary polynomial rings |    7      7  0.2s
```

### Scope Notes

- verification now rejects rewrite logs whose summed `original_factor_count` exceeds the original factor sequence length
- verification now rejects rewrite logs whose summed `optimized_factor_count` exceeds the optimized factor sequence length
- malformed rewrite records still raise deliberate `ArgumentError`s during normalization, while tampered certificates return `false` from `_verify_steinberg_optimization_certificate`
- verification remained focused to `src/algorithm/redundancy.jl`, the expert regression file, and this report
