# Issue 249 ECP Acceptance Gate Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add parent-level #185 ECP acceptance coverage and scope documentation without expanding #186 or #187 support.

**Architecture:** Keep production reducer code unchanged unless tests expose a defect. Add focused expert tests around the existing ECP public certificate and polynomial column-peel consumer, then update README/Documenter scope and add a parent coverage audit note.

**Tech Stack:** Julia, Oscar, Suslin expert tests, Documenter Markdown, Test stdlib.

## Global Constraints

- Repository has no `AGENTS.md`; follow README test commands.
- Preserve `reduce_unimodular_column(v, R)` returning a factor sequence.
- Preserve `ecp_column_reduction_certificate(v, R)` returning `ECPColumnReductionCertificate`.
- Do not implement recursive matrix factorization for #186.
- Do not implement final public Park-Woodburn acceptance for #187.
- Do not claim arbitrary Laurent `GL_n`, ToricBuilder mainline, Laurent determinant correction, or Steinberg factor-count optimization support.
- Non-unimodular columns must still fail before route work.
- Unsupported but unimodular ordinary-polynomial columns must still fail cleanly with staged ECP diagnostics.

---

## File Structure

- Modify `test/expert/elementary_column_property.jl`: add determinant-one replay assertions for the representative general ECP success case.
- Modify `test/expert/sln_to_sl3_reduction.jl`: add a length `n > 3` polynomial column-peel consumer smoke test and a determinant-not-one column-peel negative control.
- Modify `README.md`: update current scope language for #185 and keep #186/#187/Laurent/ToricBuilder boundaries staged.
- Modify `docs/src/index.md`: mirror README scope language for generated docs.
- Create `docs/audits/2026-07-02-issue-185-ecp-acceptance.md`: map #185 child issues to ECP algorithm stages and list non-claims.

### Task 1: Parent Acceptance Tests

**Files:**
- Modify: `test/expert/elementary_column_property.jl`
- Modify: `test/expert/sln_to_sl3_reduction.jl`

**Interfaces:**
- Consumes: `Suslin.reduce_unimodular_column`, `Suslin.ecp_column_reduction_certificate`, `Suslin.verify_ecp_column_reduction`, `Suslin._polynomial_column_peel_certificate`, `Suslin._verify_polynomial_column_peel_certificate`.
- Produces: acceptance checks proving representative ordinary-polynomial ECP success, tamper rejection remains present, and the polynomial column-peel consumer stores a verified ECP left certificate for a length `n > 3` input.

- [ ] **Step 1: Add the failing determinant replay assertion**

In `test/expert/elementary_column_property.jl`, inside the existing `"public ECP unimodular-column pipeline"` testset immediately after:

```julia
    @test _ecp_acceptance_apply(general_cert.factors, general_column, general_R) ==
          _ecp_acceptance_target(general_R, length(general_column))
```

add:

```julia
    @test det(_ecp_acceptance_product(general_cert.factors, general_R, length(general_column))) ==
          one(general_R)
```

- [ ] **Step 2: Add consumer smoke helpers**

In `test/expert/sln_to_sl3_reduction.jl`, after `_issue15_supported_matrix`, add:

```julia
function _issue15_wrap_column_peel_matrix(final_block, tail_entries)
    R = base_ring(final_block)
    n = nrows(final_block) + 1
    length(tail_entries) == n - 1 || throw(ArgumentError("tail_entries must match final block size"))
    wrapped = block_embedding(final_block, n, collect(1:(n - 1)))
    for row in 1:(n - 1)
        wrapped[row, n] = tail_entries[row]
    end
    return wrapped
end
```

- [ ] **Step 3: Add the length `n > 3` column-peel consumer smoke test**

In `test/expert/sln_to_sl3_reduction.jl`, inside `"SL_n to local SL3 reduction supported examples"` after the `custom_reduction` assertions, add:

```julia
    peel_matrix = _issue15_wrap_column_peel_matrix(
        block_a,
        [X, X + one(R), X^2 + X],
    )
    @test nrows(peel_matrix) > 3
    peel_cert = Suslin._polynomial_column_peel_certificate(peel_matrix)
    @test Suslin._verify_polynomial_column_peel_certificate(peel_cert)
    @test length(peel_cert.peel_steps) == 1
    first_step = only(peel_cert.peel_steps)
    @test first_step.dimension == nrows(peel_matrix)
    @test first_step.left_certificate isa Suslin.ECPColumnReductionCertificate
    @test Suslin.verify_ecp_column_reduction(first_step.left_certificate)
    @test first_step.left_certificate.original_column == first_step.last_column
    @test first_step.left_certificate.factors == first_step.left_factors
    @test first_step.left_certificate.final_column ==
          Suslin._column_peel_target_column(R, first_step.dimension)
```

- [ ] **Step 4: Add the determinant-not-one negative control**

In `test/expert/sln_to_sl3_reduction.jl`, inside `"SL_n to local SL3 reduction staged failures"` after `unsupported_err` assertions, add:

```julia
    det_not_one = identity_matrix(R, 4)
    det_not_one[1, 1] = one(R) + X
    det_err = _issue15_captured_error(() -> Suslin._polynomial_column_peel_certificate(det_not_one))
    @test det_err isa ArgumentError
    @test occursin("determinant-one input", sprint(showerror, det_err))
```

- [ ] **Step 5: Run RED tests**

Run:

```bash
julia --project=. -e 'include("test/expert/elementary_column_property.jl")'
julia --project=. -e 'include("test/expert/sln_to_sl3_reduction.jl")'
```

Expected: `sln_to_sl3_reduction.jl` fails before implementation because the new helper/test is not present yet; if an assertion already passes after insertion because #248 implemented the consumer metadata, record that this task is a characterization addition rather than a production-code change.

- [ ] **Step 6: Implement minimal test additions**

Apply the exact changes from Steps 1 through 4. Do not edit production code unless these tests reveal a genuine product defect.

- [ ] **Step 7: Run GREEN tests**

Run:

```bash
julia --project=. -e 'include("test/expert/elementary_column_property.jl")'
julia --project=. -e 'include("test/expert/sln_to_sl3_reduction.jl")'
```

Expected: both commands exit 0.

- [ ] **Step 8: Commit**

```bash
git add test/expert/elementary_column_property.jl test/expert/sln_to_sl3_reduction.jl
git commit -m "test: add issue 185 ECP acceptance gate"
```

### Task 2: Scope Documentation and Parent Coverage Note

**Files:**
- Modify: `README.md`
- Modify: `docs/src/index.md`
- Create: `docs/audits/2026-07-02-issue-185-ecp-acceptance.md`

**Interfaces:**
- Consumes: README current-scope bullets, Documenter index scope bullets, issue #249 scope language.
- Produces: documentation stating #185 ordinary-polynomial column reduction is accepted while #186/#187 and Laurent/ToricBuilder broad support remain staged.

- [ ] **Step 1: Update README scope**

In `README.md`, after the Murthy local `SL_3` solver bullet, add:

```markdown
- The ordinary-polynomial ECP unimodular-column reducer (#185) is accepted for
  exact field-backed ordinary polynomial rings through
  `reduce_unimodular_column`, `ecp_column_reduction_certificate`, and replaying
  ECP certificate verifiers. The route covers the checked input context (#243),
  monicity normalization (#244), link witness extraction (#245), link-step
  replay (#246), induction/normality composition (#247), and public reducer
  dispatch (#248). Polynomial column peel records the verified ECP certificate
  used for each last-column peel step, which gives later #186 work a stable
  consumer boundary without claiming recursive `SL_n` matrix factorization.
```

Then remove `the general ECP reducer (#185),` from the staged-boundaries bullet while leaving `recursive SL_n (#186)` and `full public Park-Woodburn acceptance (#187)` in that boundary list.

- [ ] **Step 2: Mirror the README scope in Documenter index**

Apply the same scope addition and staged-boundary removal in `docs/src/index.md`.

- [ ] **Step 3: Add the audit note**

Create `docs/audits/2026-07-02-issue-185-ecp-acceptance.md` with:

```markdown
# Issue 185 ECP Acceptance Coverage

Date: 2026-07-02

Issue #249 is the parent-level closeout gate for #185. It records that the
ordinary-polynomial ECP unimodular-column reducer is accepted for the staged
exact field-backed ordinary-polynomial route, and it records what remains out
of scope for later matrix-level work.

## Stage Map

| Issue | Stage | Evidence boundary |
| --- | --- | --- |
| #242 | ECP mainline catalog | Fixture/catalog entries for representative ordinary-polynomial columns |
| #243 | Input context | Checked ordinary-polynomial ring, column, variable order, and unimodularity context |
| #244 | Monicity normalization | Replayable coordinate changes that make the selected entry monic |
| #245 | Link witness extraction | Replayable Park-Woodburn link witness data and cover evidence |
| #246 | Link step | Exact link-step factors and route metadata, including direct endpoint transport for length greater than three |
| #247 | Induction/normality | Lower-column reduction plus conjugated-elementary normality replay |
| #248 | Public reducer route | `reduce_unimodular_column` and `ecp_column_reduction_certificate` dispatch through the general ECP pipeline |
| #249 | Parent acceptance gate | Tests and docs proving the #185 reducer boundary for later consumers |

## Acceptance Evidence

- `test/expert/elementary_column_property.jl` includes a representative
  ordinary-polynomial ECP success case, route metadata checks, determinant-one
  factor-product replay, and tampered certificate rejection.
- `test/expert/unimodular_reduction_exact.jl` keeps explicit non-unimodular and
  unsupported-but-unimodular negative controls.
- `test/expert/sln_to_sl3_reduction.jl` exercises a length `n > 3`
  polynomial column-peel consumer and verifies that the peel step stores a
  checked `ECPColumnReductionCertificate`.

## Non-Claims

- #186 recursive `SL_n` matrix factorization remains staged.
- #187 final public Park-Woodburn acceptance remains staged.
- Arbitrary Laurent `GL_n` determinant correction remains staged.
- ToricBuilder mainline support remains staged outside the documented
  certificate-backed Laurent slices.
- This gate does not optimize Steinberg factor counts.
```

- [ ] **Step 4: Check documentation wording**

Run:

```bash
rg -n "#185|#186|#187|Laurent `GL_n`|ToricBuilder|general ECP reducer|recursive `SL_n`" README.md docs/src/index.md docs/audits/2026-07-02-issue-185-ecp-acceptance.md
```

Expected: README and docs say #185 is accepted, and #186/#187/Laurent/ToricBuilder boundaries remain staged.

- [ ] **Step 5: Commit**

```bash
git add README.md docs/src/index.md docs/audits/2026-07-02-issue-185-ecp-acceptance.md
git commit -m "docs: record issue 185 ECP acceptance boundary"
```

### Task 3: Full Verification and Review

**Files:**
- No planned file edits.

**Interfaces:**
- Consumes: all changes from Tasks 1 and 2.
- Produces: verified branch ready for PR creation.

- [ ] **Step 1: Run required issue verification**

Run:

```bash
julia --project=. -e 'include("test/expert/elementary_column_property.jl")'
julia --project=. -e 'include("test/expert/unimodular_reduction_exact.jl")'
julia --project=. -e 'include("test/expert/sln_to_sl3_reduction.jl")'
julia --project=. test/runtests.jl all
```

Expected: all commands exit 0.

- [ ] **Step 2: Run Agent Desk extra verification**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: command exits 0.

- [ ] **Step 3: Request code review**

Use `superpowers:requesting-code-review` with the branch diff against `origin/main`. Fix Critical or Important findings before proceeding.

- [ ] **Step 4: Finish branch**

Use `superpowers:verification-before-completion` and `superpowers:finishing-a-development-branch`. When the finishing skill asks what to do, choose option 2, "Push and create a Pull Request", per the Agent Desk standing instructions.
