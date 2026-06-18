# Issue 21 CI Test Contract Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make Suslin.jl's dependency instantiation, default fast tests, expert tests, and full-suite CI contract explicit.

**Architecture:** Keep the existing `test/runtests.jl` group semantics. Update CI so named jobs run the exact fresh-checkout commands, and update README so contributors can see which command maps to which test groups.

**Tech Stack:** Julia package project, GitHub Actions, `julia-actions/setup-julia`, `julia-actions/cache`, Markdown.

## Global Constraints

- Do not commit a `Manifest.toml`; dependency instantiation uses `Project.toml`.
- Default fast tests are `public` plus `internal`.
- Expert tests are selected with `expert`.
- Full suite means `julia --project=. test/runtests.jl all`, covering `public`, `internal`, and `expert`.
- CI must not describe default fast coverage as full-suite coverage.
- Full-suite CI must run on pull requests and pushes so expert-test failures fail CI.

---

## File Structure

- Modify `.github/workflows/CI.yml`: define clear jobs for dependency instantiation, default fast tests, full-suite tests, and keep manual docs deployment.
- Modify `README.md`: add a Testing section with the exact commands and coverage groups.
- No changes to `test/runtests.jl`: its current default and `all` behavior already matches the contract.

---

### Task 1: Document the Local Test Contract

**Files:**
- Modify: `README.md`

**Interfaces:**
- Consumes: existing `test/runtests.jl` command-line groups: no args, `expert`, `all`.
- Produces: contributor-facing command table for CI and local verification.

- [ ] **Step 1: Add the Testing section after Current scope**

Update `README.md` so the content after the Current scope list includes:

~~~markdown
## Testing

This repository does not commit a `Manifest.toml`. In a fresh checkout,
instantiate dependencies before running tests:

```julia
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

The test runner separates routine checks from expert algorithm checks:

| Command | Coverage |
| --- | --- |
| `julia --project=. test/runtests.jl` | Default fast tests: `public` and `internal` groups |
| `julia --project=. test/runtests.jl expert` | Expert-only algorithm and documentation checks |
| `julia --project=. test/runtests.jl all` | Full suite: `public`, `internal`, and `expert` groups |
| `julia --project=. -e 'using Pkg; Pkg.test()'` | Package test entry point; runs the default fast tests |
~~~

- [ ] **Step 2: Review the rendered Markdown structure**

Run:

```bash
sed -n '1,220p' README.md
```

Expected: README shows exactly one `## Testing` section, the code fence is closed, and the command table appears before `## References`.

- [ ] **Step 3: Commit**

Run:

```bash
git add README.md
git commit -m "docs: document test suite commands"
```

---

### Task 2: Make CI Jobs Match the Command Contract

**Files:**
- Modify: `.github/workflows/CI.yml`

**Interfaces:**
- Consumes: README command contract from Task 1.
- Produces: CI jobs named `Instantiate Dependencies`, `Default Fast Tests`, and `Full Suite Tests`.

- [ ] **Step 1: Replace the workflow test jobs**

Update `.github/workflows/CI.yml` so the workflow keeps the existing triggers and `docs` job, but the test section is:

```yaml
jobs:
  instantiate:
    name: Instantiate Dependencies
    runs-on: ubuntu-latest
    timeout-minutes: 45
    permissions:
      actions: write
      contents: read
    steps:
      - uses: actions/checkout@v6
      - uses: julia-actions/setup-julia@v3
        with:
          version: '1.12'
          arch: x64
      - uses: julia-actions/cache@v3
      - name: julia --project=. -e 'using Pkg; Pkg.instantiate()'
        run: julia --project=. -e 'using Pkg; Pkg.instantiate()'

  default-tests:
    name: Default Fast Tests
    runs-on: ubuntu-latest
    timeout-minutes: 45
    permissions:
      actions: write
      contents: read
    steps:
      - uses: actions/checkout@v6
      - uses: julia-actions/setup-julia@v3
        with:
          version: '1.12'
          arch: x64
      - uses: julia-actions/cache@v3
      - name: Instantiate dependencies
        run: julia --project=. -e 'using Pkg; Pkg.instantiate()'
      - name: julia --project=. test/runtests.jl
        run: julia --project=. test/runtests.jl

  full-suite-tests:
    name: Full Suite Tests
    runs-on: ubuntu-latest
    timeout-minutes: 60
    permissions:
      actions: write
      contents: read
    steps:
      - uses: actions/checkout@v6
      - uses: julia-actions/setup-julia@v3
        with:
          version: '1.12'
          arch: x64
      - uses: julia-actions/cache@v3
      - name: Instantiate dependencies
        run: julia --project=. -e 'using Pkg; Pkg.instantiate()'
      - name: julia --project=. test/runtests.jl all
        run: julia --project=. test/runtests.jl all
```

Keep the existing manual `docs` job after these jobs. Remove the old manual-only `expert` job because the full-suite job now covers expert tests on every pull request and push.

- [ ] **Step 2: Check workflow text for stale labels**

Run:

```bash
rg -n "lightweight|expert test suite|run_expert|Pkg.test|coverage|Full Suite|Default Fast|Instantiate" .github/workflows/CI.yml
```

Expected: no stale `lightweight`, manual `run_expert`, `Pkg.test`, or coverage labels remain in CI. `Instantiate`, `Default Fast`, and `Full Suite` labels are present.

- [ ] **Step 3: Commit**

Run:

```bash
git add .github/workflows/CI.yml
git commit -m "ci: run explicit default and full test suites"
```

---

### Task 3: Verify the Contract End to End

**Files:**
- Read: `README.md`
- Read: `.github/workflows/CI.yml`
- Read: `test/runtests.jl`

**Interfaces:**
- Consumes: documented and CI commands from Tasks 1 and 2.
- Produces: verification evidence for the PR.

- [ ] **Step 1: Verify dependency instantiation**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

Expected: command exits 0.

- [ ] **Step 2: Verify default fast tests**

Run:

```bash
julia --project=. test/runtests.jl
```

Expected: command exits 0 and testsets include `public` and `internal`, not `expert`.

- [ ] **Step 3: Verify full suite tests**

Run:

```bash
julia --project=. test/runtests.jl all
```

Expected: command exits 0 and testsets include `public`, `internal`, and `expert`.

- [ ] **Step 4: Verify the package test entry point**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: command exits 0 and runs the default fast tests.

- [ ] **Step 5: Check GitHub Actions syntax enough for review**

Run:

```bash
sed -n '1,260p' .github/workflows/CI.yml
```

Expected: YAML indentation is consistent and each named test job runs the documented command.

- [ ] **Step 6: Commit any verification-only doc fixes**

If verification finds wording or label drift, update the relevant file and run:

```bash
git add README.md .github/workflows/CI.yml
git commit -m "docs: align CI labels with test contract"
```
