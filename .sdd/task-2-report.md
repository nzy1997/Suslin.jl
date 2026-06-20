# Task 2 Report

## Implementation Summary

Created [`test/fixtures/toricbuilder_issue38_cases.jl`](/Users/nzy/pycode/agent-desk/config/.agent-desk/worktrees/nzy1997-suslin.jl/issue-39-run-1-agent-issue-39-add-an-offline-fixture-for-the-issue-38-toricbui-run-1/test/fixtures/toricbuilder_issue38_cases.jl) with module `ToricBuilderIssue38Cases`.

The fixture:
- constructs the exact Issue #38 `6x6` `Q` block over `GF(2)[u^+/-1, v^+/-1]`
- records stable fixture metadata with id `toricbuilder-issue-38-q-block`
- captures determinant metadata from `classify_laurent_determinant(Q)` with determinant exactly `u*v`
- stores both the row normalization from `normalize_laurent_gl_matrix(Q)` and the explicit column normalization `Q * Dcol`
- records the current staged factorization failure substrings for both normalized cores

I did not modify [`test/internal/toricbuilder_issue38_fixture.jl`](/Users/nzy/pycode/agent-desk/config/.agent-desk/worktrees/nzy1997-suslin.jl/issue-39-run-1-agent-issue-39-add-an-offline-fixture-for-the-issue-38-toricbui-run-1/test/internal/toricbuilder_issue38_fixture.jl). The existing validator already matched the concrete fixture shape.

## TDD Evidence

### RED

Command:

```bash
julia --project=. -e 'include("test/internal/toricbuilder_issue38_fixture.jl")'
```

Observed failure before implementation:

```text
Expression: isfile(TORICBUILDER_ISSUE38_FIXTURE_PATH)
SystemError: opening file ".../test/fixtures/toricbuilder_issue38_cases.jl": No such file or directory
```

### GREEN

Same command after adding the fixture:

```text
Test Summary:                         | Pass  Total  Time
ToricBuilder Issue 38 Q block fixture |    2      2  2.7s
```

## Tests

Ran:

```bash
julia --project=. -e 'include("test/internal/toricbuilder_issue38_fixture.jl")'
```

Result: PASS.

## Files Changed

- Added [`test/fixtures/toricbuilder_issue38_cases.jl`](/Users/nzy/pycode/agent-desk/config/.agent-desk/worktrees/nzy1997-suslin.jl/issue-39-run-1-agent-issue-39-add-an-offline-fixture-for-the-issue-38-toricbui-run-1/test/fixtures/toricbuilder_issue38_cases.jl)
- Added [`.sdd/task-2-report.md`](/Users/nzy/pycode/agent-desk/config/.agent-desk/worktrees/nzy1997-suslin.jl/issue-39-run-1-agent-issue-39-add-an-offline-fixture-for-the-issue-38-toricbui-run-1/.sdd/task-2-report.md)

## Self-Review

- The fixture is offline and does not depend on a ToricBuilder checkout.
- The determinant metadata is derived from runtime computation, not duplicated by hand.
- `normalizations.row.core == normalizations.row.normalization.normalized_matrix`.
- Column normalization uses the exact right diagonal correction implied by determinant `u*v`.
- Both normalized cores reproduce the required staged failure strings under `elementary_factorization`.

## Concerns

- Only the focused validator command was required and run in this task. `test/runtests.jl` registration is intentionally left for Task 3.
