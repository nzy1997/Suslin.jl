# ToricBuilder Case 007 And Case 009 Bounded Peel Evidence

Date: 2026-06-27

This artifact records issue #147 bounded evidence for `case_007` and
`case_009`. These cases remain outside `DEFAULT_EXERCISED_CASE_IDS`; the
evidence comes only from explicit slow commands.

## Commands

```text
julia --project=. scripts/report_toricbuilder_cache_q_blocks.jl --exercise=case_007 --timeout-seconds=180 --output=/tmp/qblock-case007-peel-progress.md
julia --project=. scripts/report_toricbuilder_cache_q_blocks.jl --exercise=case_009 --timeout-seconds=650 --output=/tmp/qblock-case009-progress.md
julia --project=. -e 'include("test/internal/toricbuilder_cache_status_report.jl")'
```

## Generated Reports

- `/tmp/qblock-case007-peel-progress.md`
- `/tmp/qblock-case009-progress.md`

## Observed Evidence

### `case_007`

Case table row:

```text
| case_007 | 42x42 | 546 | default_contract | timed_out | not_run | timed_out | 0 | 180.020 |
```

Stage timing row:

```text
| case_007 | pass (11.953s) | pass (45.232s) | timed_out (83.990s) | not_run |
```

Route Error Details:

```text
- `case_007`: timed out after 180.000 seconds; peel progress: current d=30, completed steps=12, last completed d=31 (elapsed 2.680s, left factors=24, right factors=27), last-column nnz=25, max entry terms=8 while running certificate_construction
```

### `case_009`

Case table row:

```text
| case_009 | 62x62 | 739 | default_contract | timed_out | not_run | timed_out | 0 | 650.025 |
```

Stage timing row:

```text
| case_009 | pass (101.627s) | pass (312.326s) | timed_out (178.979s) | not_run |
```

Route Error Details:

```text
- `case_009`: timed out after 650.000 seconds; peel progress: current d=61, completed steps=1, last completed d=62 (elapsed 10.423s, left factors=9, right factors=9), last-column nnz=23, max entry terms=3 while running certificate_construction
```

## Interpretation

The generated rows are structured bounded evidence rows, not default report
passes. A `certificate_construction` timeout is counted as issue #147
per-peel evidence only when its Route Error Details include `peel progress:`,
`current d=`, `completed steps=`, and at least one last-column statistic.
