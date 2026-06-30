Status: DONE
Commits created: `f480251` (`test: add park woodburn ecp mainline catalog`)
One-line test summary: `julia --project=. -e 'include("test/internal/ecp_mainline_fixtures.jl")'` passed and `julia --project=. -e 'using Pkg; Pkg.test()'` passed.
Concerns, if any: None.
Report file path: `/Users/nzy/pycode/agent-desk/config/.agent-desk/worktrees/nzy1997-suslin.jl/issue-242-run-1-agent-issue-242-add-a-park-woodburn-ecp-mainline-problem-catalog-run-1/.agent-desk/sdd/task-2-report.md`

Fix note: rebuilt `negative_supported_without_evidence` without `:missing_evidence` so the control now depends only on the missing replayable link-step and lower-variable evidence expected by the review finding.
Verification: `julia --project=. -e 'include("test/internal/ecp_mainline_fixtures.jl")'` passed.
