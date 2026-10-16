# Issue 355 Park LaurentToPoly Conversion Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Land the generic Park LaurentToPoly preprocessing certificate boundary for non-unit Noether columns.

**Architecture:** Extend `LaurentToPolynomialColumnCertificate` with an explicit polynomialization power, polynomial Bezout witness, witness metadata, and replay fields. Keep the ordinary ECP bridge unchanged except for using the stricter conversion validator.

**Tech Stack:** Julia, Oscar, existing Suslin internal certificate helpers.

## Global Constraints

- Close only #355.
- Do not implement #366 monicity planning, #367 sparse planned transformation, #368 ordinary child construction, or #369 full d14 bridge acceptance.
- Do not change production reducer dispatch.
- Do not include the run-193 ECP verifier weakening based on `stage.is_unimodular === true`.
- Keep `case_008` fixture choices outside generic algorithm code.

---

### Task 1: Focus the Branch Scope

**Files:**
- Delete: `test/expert/case008_d14_laurent_to_polynomial_bridge.jl`
- Modify: `test/ci/runtests.jl`
- Modify: `test/ci/shards.toml`

**Interfaces:**
- Consumes: updated #355 issue body.
- Produces: branch diff without full-bridge test coverage.

- [ ] Remove old full-bridge test and stale full-bridge Superpowers artifacts.
- [ ] Replace CI entries with the focused conversion expert test.

### Task 2: Certificate Witness Replay

**Files:**
- Modify: `src/core/laurent_to_polynomial.jl`
- Modify: `test/internal/laurent_to_polynomial_certificate.jl`

**Interfaces:**
- Produces: `LaurentToPolynomialColumnCertificate.polynomialization_power`.
- Produces: `LaurentToPolynomialColumnCertificate.polynomial_witness`.
- Produces: `LaurentToPolynomialColumnCertificate.polynomial_witness_metadata`.

- [ ] Add helper functions to count terms, construct witness metadata, and verify the polynomial Bezout identity.
- [ ] Use a unit-entry witness for unit fast paths, a Park-shaped witness when the selected entry has constant term one and other entries are selected-variable multiples, and the existing coordinates-backed witness for small ordinary polynomial columns.
- [ ] Update replay so `polynomial_unimodular` derives from witness identity replay, not from a trusted flag.
- [ ] Restore ECP bridge validation to call `_validate_laurent_to_polynomial_certificate` directly.
- [ ] Add internal tests for witness replay and negative witness/power controls.

### Task 3: Bounded d14 Conversion Report

**Files:**
- Create: `scripts/report_case008_d14_laurent_to_polynomial_conversion.jl`
- Create: `test/expert/case008_d14_laurent_to_polynomial_conversion.jl`

**Interfaces:**
- Produces: `case008_d14_laurent_to_polynomial_conversion_report(; timeout_seconds)`.
- Produces: `validate_case008_d14_laurent_to_polynomial_conversion_report(report)::Symbol`.

- [ ] Load the existing validated d14 boundary fixture.
- [ ] Run the generic conversion with selected entry 9 and selected generator `v`.
- [ ] Record conversion status, replay status, factor count, source terms, polynomial terms, polynomialization power, polynomial witness replay, witness terms, elapsed time, and final status.
- [ ] Add negative controls for non-unimodular source mutation and tampered conversion/witness payloads.

### Task 4: Verification and PR

**Files:**
- No new implementation files.

**Interfaces:**
- Consumes: focused verification commands from #355.
- Produces: pushed branch and PR containing `Closes #355`.

- [ ] Run `julia --project=. -e 'include("test/internal/laurent_to_polynomial_certificate.jl")'`.
- [ ] Run `julia --project=. scripts/report_case008_d14_laurent_to_polynomial_conversion.jl --timeout-seconds=300 --output=/tmp/case008-d14-ltp-conversion.md`.
- [ ] Run touched CI manifest verification if CI manifest files changed.
- [ ] Commit, push, and create the PR.
