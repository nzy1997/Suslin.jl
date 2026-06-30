# Task 1 Report

Status: red validator contract added and verified.

Test summary: `julia --project=. -e 'include("test/internal/ecp_mainline_fixtures.jl")'` now fails for the intended reason after instantiating dependencies: `test/fixtures/ecp_mainline_cases.jl` does not exist yet.

Relevant output:

```text
ECP mainline fixture catalog: Error During Test at .../test/internal/ecp_mainline_fixtures.jl:290
  Got exception outside of a @test
  SystemError: opening file ".../test/fixtures/ecp_mainline_cases.jl": No such file or directory
```

Note: the first attempt failed earlier because `Oscar` was not yet instantiated in this checkout; `Pkg.instantiate()` resolved that environment issue before the required red run.

Update for review findings:

- The ring contract now requires an ordinary Oscar polynomial ring (`MPolyRing` or `PolyRing`), an exact coefficient type, and a field coefficient ring.
- The stage-evidence contract now requires all status symbols, replays supported `link_step` and `lower_variable` evidence when replayable metadata is present, and rejects supported entries that still carry `missing_evidence`.
- Catalog validation now checks ID uniqueness across `cases` and `negative_controls` together.
- Supported entries now have to present real replayable support evidence; `missing_evidence = ()` no longer softens the supported/staged boundary.

Covering command:

```text
julia --project=. -e 'include("test/internal/ecp_mainline_fixtures.jl")'
```

Output summary:

```text
Testset reached `include("test/fixtures/ecp_mainline_cases.jl")` and failed only because that file is absent.
```

Files changed:

- `test/internal/ecp_mainline_fixtures.jl`
- `.agent-desk/sdd/task-1-report.md`

Commit SHA:

- `e0354078a5cebd31ef457767749ac49363da01de`
