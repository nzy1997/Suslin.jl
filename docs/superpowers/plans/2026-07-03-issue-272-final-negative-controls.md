# Issue 272 Final Negative Controls Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Prove final Park-Woodburn unsupported inputs and corrupted evidence fail before unverifiable factors can pass.

**Architecture:** Extend the final #270 fixture catalog with explicit negative-control records and public failure metadata, then add public tests that loop over those records through `elementary_factorization(A)`. Add one expert route-certificate tamper check for final recursive `SL_n` nested ECP evidence while preserving the returned public factor sequence.

**Tech Stack:** Julia, Oscar, `Test`, existing Suslin Park-Woodburn route certificate internals.

## Global Constraints

- Do not broaden the accepted input class.
- Do not treat Laurent/ToricBuilder or Steinberg optimization as part of #187.
- Unsupported public inputs must throw `ArgumentError` before factors are returned.
- Verifier tamper checks must return `false`, not throw as the expected success path.
- Required negative cases: determinant-not-one, unsupported coefficient ring, missing `SL_3` local-form evidence, missing `SL_3` ordinary Quillen evidence, missing ECP evidence, missing final `SL_3` evidence for recursive `SL_n`, and Laurent/ToricBuilder ordinary-public boundary.
- Required tamper cases: corrupted returned factor sequence and corrupted nested route evidence.

---

### Task 1: Final Negative Catalog Entries

**Files:**
- Modify: `test/fixtures/park_woodburn_mainline_acceptance_cases.jl`
- Modify: `test/internal/park_woodburn_mainline_acceptance_fixtures.jl`

**Interfaces:**
- Consumes: `ParkWoodburnMainlineAcceptanceFixtureCatalog.catalog()`.
- Produces: final negative controls with a `public_failure` NamedTuple containing `terms`, `staged_route`, and optional `reason_code`.

- [ ] **Step 1: Write the failing validator requirements**

In `test/internal/park_woodburn_mainline_acceptance_fixtures.jl`, expand `REQUIRED_PW_MAINLINE_NEGATIVE_IDS` to include:

```julia
"pw-mainline-negative-missing-sl3-local-form-evidence"
"pw-mainline-negative-missing-sl3-quillen-evidence"
"pw-mainline-negative-missing-ecp-evidence"
"pw-mainline-negative-missing-final-sl3-evidence"
"pw-mainline-negative-laurent-boundary"
```

Add `REQUIRED_PW_MAINLINE_NEGATIVE_KINDS` with one exact kind per required id:

```julia
const REQUIRED_PW_MAINLINE_NEGATIVE_KINDS = Dict(
    "pw-mainline-negative-det-not-one" => :determinant_not_one,
    "pw-mainline-negative-unsupported-coefficient-ring" => :unsupported_coefficient_ring,
    "pw-mainline-negative-missing-evidence" => :missing_final_sl3_evidence,
    "pw-mainline-negative-missing-sl3-local-form-evidence" => :missing_sl3_local_form_evidence,
    "pw-mainline-negative-missing-sl3-quillen-evidence" => :missing_sl3_quillen_evidence,
    "pw-mainline-negative-missing-ecp-evidence" => :missing_ecp_evidence,
    "pw-mainline-negative-missing-final-sl3-evidence" => :missing_final_sl3_evidence,
    "pw-mainline-negative-laurent-boundary" => :laurent_boundary,
)
```

Inside the negative-control loop, require `negative_kind` and `public_failure.terms`, and check the kind matches `REQUIRED_PW_MAINLINE_NEGATIVE_KINDS[entry.id]`.

- [ ] **Step 2: Run the red validator command**

Run:

```bash
julia --project=. -e 'include("test/internal/park_woodburn_mainline_acceptance_fixtures.jl")'
```

Expected: FAIL because the newly required negative ids are absent from the catalog.

- [ ] **Step 3: Add catalog negative controls**

In `test/fixtures/park_woodburn_mainline_acceptance_cases.jl`, add `_negative_failure`:

```julia
function _negative_failure(terms; staged_route = true, reason_code = nothing)
    data = (; terms = Tuple(terms), staged_route = staged_route)
    reason_code === nothing || (data = merge(data, (; reason_code = reason_code)))
    return data
end
```

Update `_negative_control` to attach `negative_kind` and `public_failure`.

Add matrices for:

```julia
missing_local_form_matrix =
    elementary_matrix(3, 1, 3, issue238_X, issue238_R) *
    elementary_matrix(3, 2, 1, issue238_r, issue238_R)

missing_quillen_matrix =
    elementary_matrix(3, 1, 2, one(issue238_R), issue238_R) *
    elementary_matrix(3, 2, 1, issue238_X, issue238_R) *
    elementary_matrix(3, 1, 2, one(issue238_R), issue238_R)

missing_ecp_matrix = matrix(missing_ecp_R, [
    one(missing_ecp_R)  zero(missing_ecp_R) zero(missing_ecp_R) zero(missing_ecp_R);
    zero(missing_ecp_R) zero(missing_ecp_R) missing_ecp_p     missing_ecp_a;
    zero(missing_ecp_R) zero(missing_ecp_R) missing_ecp_q     missing_ecp_b;
    zero(missing_ecp_R) one(missing_ecp_R)  zero(missing_ecp_R) zero(missing_ecp_R)
])

normalizable_laurent = matrix(L, [
    lx      zero(L) zero(L);
    zero(L) one(L)  zero(L);
    zero(L) zero(L) one(L)
])
```

Change `"pw-mainline-negative-missing-evidence"` so its matrix is the staged
missing-final-`SL_3` matrix, not the supported recursive positive matrix.

- [ ] **Step 4: Run the green validator command**

Run:

```bash
julia --project=. -e 'include("test/internal/park_woodburn_mainline_acceptance_fixtures.jl")'
```

Expected: PASS.

- [ ] **Step 5: Commit Task 1**

Run:

```bash
git add test/fixtures/park_woodburn_mainline_acceptance_cases.jl test/internal/park_woodburn_mainline_acceptance_fixtures.jl docs/superpowers/plans/2026-07-03-issue-272-final-negative-controls.md
git commit -m "test: catalog final Park-Woodburn negative controls"
```

### Task 2: Public Unsupported-Input Failure Checks

**Files:**
- Modify: `test/public/factorization_driver_shell.jl`
- Modify: `test/public/park_woodburn_polynomial_factorization.jl`

**Interfaces:**
- Consumes: final catalog `negative_controls` and their `public_failure` metadata.
- Produces: helper assertions that prove each negative public call throws before assigning factors.

- [ ] **Step 1: Write the failing public loop**

Add helper functions in both public test files:

```julia
function _pw_failure_field(failure, field::Symbol, default)
    return hasproperty(failure, field) ? getproperty(failure, field) : default
end

function _pw_assert_mainline_negative_public_failure(entry)
    factors = nothing
    err = _captured_error(() -> begin
        factors = elementary_factorization(entry.matrix)
        nothing
    end)
    @test factors === nothing
    @test err isa ArgumentError
    msg = sprint(showerror, err)
    for term in entry.public_failure.terms
        @test occursin(term, msg)
    end

    if _pw_failure_field(entry.public_failure, :staged_route, false)
        cert = Suslin._polynomial_factorization_route_certificate(entry.matrix)
        @test cert.route == :staged_failure
        @test cert.status == :staged
        @test isempty(cert.factors)
        @test Suslin._verify_polynomial_factorization_route_certificate(cert)
        expected_reason = _pw_failure_field(entry.public_failure, :reason_code, nothing)
        if expected_reason !== nothing
            @test cert.evidence.reason_code == expected_reason
        end
    end
end
```

Use `_pw_acceptance_result_or_error` instead of `_captured_error` in
`park_woodburn_polynomial_factorization.jl`, where that helper already returns
both factors and errors.

Loop over every `ParkWoodburnMainlineAcceptanceFixtureCatalog.catalog().negative_controls`.

- [ ] **Step 2: Run the red public commands**

Run:

```bash
julia --project=. -e 'include("test/public/factorization_driver_shell.jl")'
julia --project=. -e 'include("test/public/park_woodburn_polynomial_factorization.jl")'
```

Expected before Task 1 is complete: FAIL because new catalog metadata or helper
expectations are absent. Expected after Task 1 and before helper adaptation is
complete: FAIL if any negative catalog entry returns factors or lacks expected
failure metadata.

- [ ] **Step 3: Complete public helpers and loop placement**

In `factorization_driver_shell.jl`, place the loop after `mainline_entries` are
loaded:

```julia
for entry in ParkWoodburnMainlineAcceptanceFixtureCatalog.catalog().negative_controls
    _pw_assert_mainline_negative_public_failure(entry)
end
```

In `park_woodburn_polynomial_factorization.jl`, add the same loop after final
positive catalog cases are checked and keep the existing older polynomial
negative controls.

- [ ] **Step 4: Run green public commands**

Run:

```bash
julia --project=. -e 'include("test/public/factorization_driver_shell.jl")'
julia --project=. -e 'include("test/public/park_woodburn_polynomial_factorization.jl")'
```

Expected: both commands exit 0.

- [ ] **Step 5: Commit Task 2**

Run:

```bash
git add test/public/factorization_driver_shell.jl test/public/park_woodburn_polynomial_factorization.jl
git commit -m "test: prove final Park-Woodburn public negatives fail"
```

### Task 3: Final Nested Evidence Tamper Gate

**Files:**
- Modify: `test/expert/park_woodburn_route_certificate.jl`

**Interfaces:**
- Consumes: `PolynomialColumnPeelCertificate` from the final recursive #187 catalog case.
- Produces: explicit verifier-false coverage for corrupted nested ECP evidence in the final recursive route.

- [ ] **Step 1: Write the failing expert tamper assertion**

In the final mainline acceptance certificate block, rebuild the first peel
step's ECP evidence with `verification = nothing`, rebuild the peel certificate
with that step, and assert:

```julia
@test !Suslin.verify_ecp_column_reduction(bad_ecp_evidence)
@test verify_factorization(bad_ecp_mainline_cert.matrix, bad_ecp_mainline_cert.factors)
@test !Suslin._verify_polynomial_column_peel_certificate(bad_ecp_mainline_cert.evidence)
@test !Suslin._verify_polynomial_factorization_route_certificate(bad_ecp_mainline_cert)
```

- [ ] **Step 2: Run the red expert command**

Run:

```bash
julia --project=. -e 'include("test/expert/park_woodburn_route_certificate.jl")'
```

Expected: FAIL until the exact rebuild uses the current `PolynomialColumnPeelStep`
fields correctly.

- [ ] **Step 3: Implement the rebuild using existing `_pw_rebuild` helpers**

Use:

```julia
bad_ecp_evidence = _pw_rebuild(first_step.ecp_evidence; verification = nothing)
bad_ecp_step = _pw_rebuild(first_step; ecp_evidence = bad_ecp_evidence)
bad_ecp_steps = copy(mainline_acceptance_cert.evidence.peel_steps)
bad_ecp_steps[1] = bad_ecp_step
```

Then rebuild the column-peel certificate and route certificate through
`_pw_replace_column_peel_certificate` and `_pw_replace_certificate`.

- [ ] **Step 4: Run the green expert command**

Run:

```bash
julia --project=. -e 'include("test/expert/park_woodburn_route_certificate.jl")'
```

Expected: command exits 0.

- [ ] **Step 5: Commit Task 3**

Run:

```bash
git add test/expert/park_woodburn_route_certificate.jl
git commit -m "test: reject final nested Park-Woodburn evidence tampering"
```

### Task 4: Full Verification And PR

**Files:**
- No source files expected beyond Tasks 1-3.

**Interfaces:**
- Consumes: all committed task changes.
- Produces: verified branch pushed to a pull request.

- [ ] **Step 1: Run required verification**

Run:

```bash
julia --project=. -e 'include("test/public/factorization_driver_shell.jl")'
julia --project=. -e 'include("test/public/park_woodburn_polynomial_factorization.jl")'
julia --project=. -e 'include("test/expert/park_woodburn_route_certificate.jl")'
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: all commands exit 0.

- [ ] **Step 2: Final review**

Dispatch a code reviewer subagent with the branch diff and fix Critical or
Important findings before continuing.

- [ ] **Step 3: Finish branch**

Use `superpowers:finishing-a-development-branch`, choose "Push and create a Pull
Request", push `agent/issue-272-harden-final-park-woodburn-unsupported-input-and-run-1`,
and create a PR against `main`.
