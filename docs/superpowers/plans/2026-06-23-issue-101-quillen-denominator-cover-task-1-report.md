# Issue 101 Task 1 Report: Quillen Denominator Cover RED Test

## Files Changed

- `test/expert/quillen_denominator_cover.jl`
- `test/runtests.jl`
- `docs/superpowers/plans/2026-06-23-issue-101-quillen-denominator-cover-task-1-report.md`

## RED Verification

Command:

```bash
julia --project=. -e 'include("test/expert/quillen_denominator_cover.jl")'
```

Result: exited 1 as expected after the RED test was added.

Output summary:

- Testset: `Quillen denominator cover certificates`
- Failure mode: `UndefVarError: quillen_denominator_cover_certificate not defined in Suslin`
- Location: `test/expert/quillen_denominator_cover.jl:18`
- Test summary: `0 passed, 0 failed, 1 errored`

## Concerns

- The first focused run could not load `Suslin` because the local Julia environment had not installed `Oscar`. Running `julia --project=. -e 'using Pkg; Pkg.instantiate()'` installed the recorded dependencies; it left no `Project.toml` or `Manifest.toml` diff.
- The RED failure is an error outside `@test` rather than an assertion failure because the planned non-exported API does not exist yet. This matches Task 1's expected undefined-helper failure.
