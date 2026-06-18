# Issue 21 CI and Test Contract Design

## Context

`test/runtests.jl` already separates tests into `public`, `internal`, and
`expert` groups. With no arguments it runs `public` and `internal`; with
`all` it runs all three groups. The current CI workflow uses Julia package
helpers and only runs expert tests from a manual dispatch input, so it does not
show the exact commands named in issue 21.

## Design Choice

Use the existing test runner semantics and make the command contract explicit
in CI and README documentation.

Alternative considered: change the default test runner to run every group.
That would make default testing slower and contradict the existing separation
between routine and expert coverage.

Alternative considered: keep CI on `Pkg.test` only and describe the mapping in
documentation. That is less clear because CI logs would not show the requested
fresh-checkout commands.

## Test Contract

- Dependency instantiation: `julia --project=. -e 'using Pkg; Pkg.instantiate()'`.
- Default fast tests: `julia --project=. test/runtests.jl`, covering `public`
  and `internal`.
- Expert-only tests: `julia --project=. test/runtests.jl expert`.
- Full suite: `julia --project=. test/runtests.jl all`, covering `public`,
  `internal`, and `expert`.
- Julia package test entry point: `julia --project=. -e 'using Pkg; Pkg.test()'`
  remains valid and exercises the default fast groups because it invokes
  `test/runtests.jl` without extra arguments.

## CI Shape

Keep one CI workflow and split the named checks by contract:

- `Instantiate dependencies` runs the explicit `Pkg.instantiate()` command.
- `Default Fast Tests` runs `julia --project=. test/runtests.jl`.
- `Full Suite Tests` runs `julia --project=. test/runtests.jl all`.

The full-suite check must run on pull requests and pushes so an expert-test
failure fails CI. CI labels must not describe the default fast check as full
suite coverage.

## Documentation

Add a README test section that lists the commands above and states exactly
which groups each command covers. No Manifest is committed; contributors should
instantiate dependencies from `Project.toml`.

## Issue Metadata

Repository documentation is the durable source for this PR. Editing existing
GitHub issue bodies is a public metadata change, so this implementation will
avoid issue-body rewrites unless a specific stale reference is identified and
the final PR can point maintainers at the corrected README contract instead.
