### Task 3: Expert Test Registration And Verification

**Files:**
- Modify: `test/runtests.jl`
- Modify: `docs/superpowers/plans/2026-06-21-issue-71-sl3-local-q-degree-normalization.md`

**Interfaces:**
- Consumes: focused test from Task 1 and helper/verifier implementation from Task 2.
- Produces: registered expert test and committed implementation.

- [ ] **Step 1: Register the expert test**

In `test/runtests.jl`, add this entry immediately after
`"expert/sl3_local_certificate.jl",`:

```julia
"expert/sl3_local_q_degree_normalization.jl",
```

- [ ] **Step 2: Run expert and package verification**

Run:

```bash
julia --project=. -e 'include("test/expert/sl3_local_q_degree_normalization.jl")'
julia --project=. test/runtests.jl expert
julia --project=. -e 'using Pkg; Pkg.test()'
git diff --check
```

Expected:

- Focused q-degree test exits 0.
- Expert suite exits 0.
- Package test exits 0.
- `git diff --check` exits 0.

- [ ] **Step 3: Commit implementation**

Stage and commit:

```bash
git add \
    docs/superpowers/plans/2026-06-21-issue-71-sl3-local-q-degree-normalization.md \
    src/algorithm/sl3_local.jl \
    test/expert/sl3_local_q_degree_normalization.jl \
    test/runtests.jl
git commit -m "Add SL3 local q-degree normalization"
```

Expected: commit succeeds with only the plan, source, focused expert test, and test registration staged.

---

## Self-Review

- Spec coverage: the tasks cover a polynomial-only helper, fixture reuse, two exact q-degree examples, quotient/remainder degree checks, normalized target and correction recording, exact replay through `verify_sl3_local_realization`, non-monic rejection, corrupted witness rejection, focused verification, package verification, and no public driver change.
- Placeholder scan: no deferred implementation markers remain.
- Type consistency: `normalize_sl3_local_q_degree`, `SL3LocalRealizationCertificate`, `verify_sl3_local_realization`, `:q_degree_normalization`, and witness field names are used consistently across test and implementation steps.
