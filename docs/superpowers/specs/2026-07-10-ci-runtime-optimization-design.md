# CI Runtime Optimization Design

## Goal

Reduce routine pull-request feedback from more than one hour to 5–15 minutes
for ordinary single-domain changes while preserving strong coverage enforcement
and leaving GitHub Actions concurrency available for other pull requests and
repositories owned by the same account.

The design replaces per-PR full-suite coverage with conservative affected-test
selection and fresh patch coverage. A complete coverage run refreshes the
project baseline once per day at 06:00 Asia/Shanghai.

## Current Evidence

The current workflow runs three independent jobs on every pull request:

- `Instantiate Dependencies`;
- `Default Fast Tests`, which runs the `public` and `internal` groups without
  coverage; and
- `Full Suite Tests`, which runs `public`, `internal`, and `expert` again with
  coverage.

In GitHub Actions run `29081337398`, dependency instantiation and Codecov
processing took seconds, while the covered test command took about 111 minutes:

- `public`: 7m20s;
- `internal`: 23m33s; and
- `expert`: 79m58s.

The bottleneck is therefore the sequential full test suite, especially the
expert group, rather than checkout, dependency caching, coverage conversion, or
Codecov upload.

## Non-goals

- Do not reduce the set of tests included in the daily complete suite.
- Do not treat carried-forward coverage as fresh coverage for changed lines.
- Do not depend on developers remembering to run a local command before push.
- Do not introduce larger paid runners or self-hosted runners.
- Do not automatically create GitHub issues or send external notifications for
  nightly failures.

## Workflow Architecture

### Pull requests

The pull-request workflow in `.github/workflows/CI.yml` computes the files
changed from the PR base to the head, selects the affected test shards, and runs
those shards as a matrix with `max-parallel: 4`. It then combines their LCOV
reports and exposes one stable `PR Gate` check.

The PR workflow does not run the complete suite. At most four selected shards
run concurrently, so two simultaneous pull requests consume at most eight
standard runner slots. On the GitHub Free concurrency limit of 20 standard
jobs, this leaves at least twelve slots for other work when the nightly workflow
is not running.

Superseded commits on the same pull request retain the existing
`cancel-in-progress` behavior.

### Complete coverage

A separate `.github/workflows/Nightly.yml` complete-coverage workflow runs:

- daily at `0 22 * * *` UTC, which is 06:00 Asia/Shanghai;
- on manual dispatch; and
- on tag pushes, preserving a full test signal for releases.

It runs every registered test shard with `max-parallel: 4` and refreshes the
complete Codecov baseline. The scheduled workflow may start a few minutes late
under GitHub load.

The complete suite does not run on every push to `main`. Consequently, project
coverage between nightly runs is partly carried forward and can be up to 24
hours stale. Fresh patch coverage remains mandatory on every source-changing
pull request.

The existing manually dispatched documentation job retains its behavior.

### Existing jobs

After the replacement gate is verified:

- remove the standalone `Instantiate Dependencies` job because its runner does
  not prepare the independent test runners;
- remove the duplicated `Default Fast Tests` job; and
- remove the per-PR `Full Suite Tests` job.

Each selected shard remains responsible for checking out the repository,
setting up Julia, restoring the Julia cache, and instantiating dependencies.

## Test Shards

### Single source of truth

Add `test/ci/shards.toml` as the single source of truth for local selection,
pull-request selection, and the complete workflow. It records:

- ten stable coverage shard identifiers;
- the test files contained in each shard;
- source-path and fixture impact rules; and
- the special non-coverage documentation-smoke target.

The initial coverage shards are:

1. `public`;
2. `internal-core`;
3. `internal-fixtures`;
4. `expert-core`;
5. `expert-laurent-a`;
6. `expert-laurent-b`;
7. `expert-sl3`;
8. `expert-quillen`;
9. `expert-ecp`; and
10. `expert-integration`.

The implementation assigns every file currently registered in
`TEST_GROUP_FILES` to exactly one shard. Assignment follows the logical domain
first and uses measured per-file timings to balance the two Laurent shards and
the two internal shards. The target is 10–15 minutes per shard on a warm
GitHub-hosted runner. `test/expert/documentation_smoke.jl` is also exposed as a
standalone `documentation-smoke` target for documentation-only pull requests;
it remains registered in exactly one coverage shard for the complete run.

### Test runner interface

Extend `test/runtests.jl` without removing its current public interface:

- no arguments still run `public` and `internal`;
- `public`, `internal`, `expert`, and `all` retain their current meanings;
- `shard:<id>` runs one configured coverage shard; and
- `documentation-smoke` runs only the documentation smoke test.

Every included test file reports its elapsed time. A manifest validation mode
checks that all configured files exist, every file previously registered in
`TEST_GROUP_FILES` belongs to exactly one coverage shard, no shard is empty,
and the shard union equals the existing full suite. Test files that were not
part of the existing full suite are not silently added by this migration.

## Affected-test Selection

Add `test/ci/select_shards.jl`. It accepts explicit base and head revisions and
emits the selected matrix plus a human-readable reason for every selection.
The local coverage command and GitHub workflow call the same selector.

Selection is deliberately fail-closed:

- a specific algorithm change selects `public` plus its mapped internal and
  expert shards;
- a changed test file selects the shard that owns it;
- a changed fixture selects every shard that consumes it;
- documentation-only changes select `documentation-smoke` and do not request
  source coverage;
- changes under `src/core/`, changes to `src/Suslin.jl`, dependency metadata,
  `test/runtests.jl`, or the CI selection configuration select all ten coverage
  shards;
- an unknown changed source path selects all ten coverage shards; and
- a Git diff or selector failure falls back to all ten coverage shards.

A source-changing pull request is never allowed to produce an empty selection.
The workflow contains an independent fallback step so that a syntax or startup
failure in the selector itself still expands the matrix to all shards.

## Coverage Data Flow

### Pull-request reports

Each selected shard runs Julia coverage scoped to `src`, generates a uniquely
named LCOV report, uploads the small report as a short-retention workflow
artifact, and uploads it to Codecov under the `pr-selected` flag. The aggregate
`PR Gate` downloads all selected reports and calculates diff coverage against
the PR base.

The gate fails when changed executable source lines have less than 99% line
coverage. Documentation-only changes have no changed executable lines and pass
the coverage portion of the gate.

### Complete reports

Every complete-workflow shard uploads its report under the `full-suite` flag
with a unique upload name. Codecov merges the reports for the commit into one
complete result.

Configure Codecov flags as follows:

- `full-suite`: `carryforward: true`;
- `pr-selected`: `carryforward: false`.

The most recent successful complete report supplies unchanged coverage between
nightly runs. Coverage for current changed lines must come from the current
`pr-selected` reports and is independently enforced by `PR Gate`; it cannot be
satisfied by carryforward data.

The existing project status retains `target: auto` and a 1% threshold. The
nightly workflow treats failed tests or a failed Codecov upload as a failure so
that an unrefreshed baseline is visible. PR Codecov upload errors remain
non-blocking because the repository-owned diff-coverage gate is authoritative.

### Local changed coverage

Add a documented local entry point that accepts `--base`, defaults to
`origin/main`, invokes the shared selector, runs the selected targets, writes
coverage to a unique temporary LCOV trace file, filters the report to `src`, and
runs `diff-cover` through `uvx`. The command checks for `uvx` and reports a
clear installation prerequisite when it is unavailable. Temporary files
prevent stale ignored `.cov` files from contaminating the result.

The local result is explicitly called patch coverage. It must not be presented
as a fresh complete-project percentage.

## Failure Handling

The fixed-name `PR Gate` runs with `if: always()` semantics and fails if:

- shard selection fails without producing the all-shards fallback;
- any selected test shard fails or times out;
- any required LCOV artifact is missing; or
- changed-line coverage is below 99%.

Individual shard jobs have a 30-minute timeout. The expected runtime is 10–15
minutes; reaching the timeout indicates a hung test or a shard that needs
rebalancing.

A failed complete workflow does not retroactively block merged pull requests.
It remains red in GitHub Actions, and the previous successful `full-suite`
report remains the carryforward baseline. No issue or external notification is
created automatically.

`main` currently has no branch protection. Repository changes will provide the
stable `PR Gate` check and document that it should be the required check.
Changing repository protection settings is an explicit operational follow-up,
not an automatic part of the code change. Until protection is enabled, a direct
push to `main` can bypass the PR gate and will not receive a complete test
signal until the next scheduled run.

## Verification

### Selector and manifest tests

Automated tests cover at least these cases:

- a Laurent algorithm change selects the public and Laurent-related shards;
- an SL3, Quillen, or ECP change selects its mapped shards;
- a changed test selects its owning shard;
- a changed fixture selects all consumers;
- a documentation-only change selects only `documentation-smoke`;
- a `src/core` change selects all ten coverage shards;
- an unknown source file selects all ten coverage shards;
- a diff failure selects all ten coverage shards;
- an empty source selection is rejected; and
- a missing or duplicated existing full-suite test file fails manifest
  validation.

### Coverage-gate tests

Small synthetic LCOV and diff fixtures prove that:

- at least 99% changed-line coverage passes;
- less than 99% fails;
- combining multiple selected reports counts a line as covered when any report
  covers it;
- a missing selected report fails; and
- a documentation-only diff passes without inventing source coverage.

### Workflow verification

Before removing the old per-PR full jobs:

1. run all ten shards once and verify their union passes;
2. confirm the complete reports merge into the expected Codecov baseline;
3. exercise a source-changing selection and its negative coverage control;
4. exercise a documentation-only selection; and
5. confirm that the aggregate check is named exactly `PR Gate`.

After these checks pass, remove the duplicated jobs and update the README test
contract.

## Acceptance Criteria

- An ordinary single-domain pull request completes in 5–15 minutes under normal
  runner availability and no later than 20 minutes after a warm cache.
- A pull request starts no more than four test jobs concurrently.
- Two simultaneous pull requests consume no more than eight concurrent test
  jobs.
- The complete workflow starts no more than four test jobs concurrently.
- Two simultaneous pull requests plus the complete workflow consume no more
  than twelve concurrent test jobs, leaving at least eight Free-plan slots for
  other repositories.
- Shared-core and unknown changes safely fall back to all ten shards even when
  that exceeds the ordinary 20-minute target.
- Changed executable lines have at least 99% fresh patch coverage.
- The daily complete suite contains every test that the current `all` command
  contains, with no omissions or duplicates.
- The complete Codecov project status retains its current automatic target and
  1% tolerance.
- The README documents routine commands, changed-coverage usage, complete-suite
  behavior, the 06:00 Asia/Shanghai schedule, and the carryforward freshness
  limitation.
