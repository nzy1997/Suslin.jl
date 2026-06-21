# Issue 60 ToricBuilder Cache Smoke Gate Design

## Goal

Add an optional local smoke gate for selected ToricBuilder cache cases that
reports deterministic, reviewer-friendly route results for the recursive
Laurent column-peel path and the Laurent `GL_n` certificate path.

## Context

The Suslin package and CI must not depend on a local ToricBuilder checkout.
The local cache environment lives outside this repository at the ToricBuilder
checkout and contains serialized cache cases such as `case_001.jls`,
`case_004.jls`, and failed/unavailable cases such as `case_013.jls`.

The relevant existing Suslin routes are:

- `Suslin._factor_laurent_sl_column_peel(A)` for determinant-one Laurent
  `SL_n` blocks.
- `Suslin.laurent_gl_factorization_certificate(A)` and
  `Suslin.verify_laurent_gl_factorization_certificate(certificate)` for
  Laurent monomial-unit `GL_n` blocks after determinant normalization.

Issue #60 should exercise those routes locally and summarize stable statuses.
It should not add ToricBuilder to the main package dependencies, main tests, or
GitHub Actions.

## Approaches Considered

1. Add a focused optional example project under `example/toric_decoupling`.

   This is the selected approach. It keeps ToricBuilder loading dynamic and
   local, exposes the requested command path, and provides a local test command
   without registering it in Suslin's package tests.

2. Add the smoke gate to the main test suite behind environment variables.

   This is rejected because it risks turning the optional cache checkout into a
   CI concern and widens the routine test surface.

3. Expand the checked-in fixture catalog with serialized local cache cases.

   This is out of scope. The issue asks for an optional local smoke gate, not
   for copying cache artifacts into Suslin.

## Design

Create `example/toric_decoupling/try_column_block_decoupling.jl` as a guarded
script with reusable functions and a `main(args=ARGS; io=stdout)` entry point.
The script will:

- Accept `--case=case_001,case_004`, defaulting to those two cases.
- Accept optional `--toricbuilder-dir=...` and `--cache-dir=...` overrides.
- Default `toricbuilder-dir` from `ENV["TORICBUILDER_DIR"]` or the local
  checkout path already recorded in the ToricBuilder fixture provenance.
- Load Suslin and ToricBuilder dynamically by pushing the Suslin repository
  root and ToricBuilder checkout onto `LOAD_PATH`.
- Never throw raw cache or factorization stack traces during normal smoke
  reporting; cache/package/load failures become stable `CACHE_ERROR` rows.

For each available `:ok` cache case, the initial selected blocks are:

- `column_Q`: the upper-left decoder basis-change block from
  `transfer_result.column_transformation`. Its determinant is commonly a
  Laurent monomial unit, so it exercises the `GL_n` certificate path.
- `pair_mix_2_1`: the deterministic toric pair mix used by the existing local
  cache inspection script. Its determinant is `1`, so it exercises the direct
  `SL_n` column-peel path.

Each output row is a single stable line beginning with `TORIC_SMOKE` and using
`key=value` tokens:

```text
TORIC_SMOKE case=case_001 block=pair_mix_2_1 size=6x6 det=one normalization=NORMALIZATION_SKIP sl_core=SL_CORE_PASS gl_cert=GL_CERT_SKIP factors=41 verified=true status=PASS failure=NONE
```

The row includes the requested case id, block role, size, determinant
classification, normalization status, `SL` core status, `GL` certificate
status, factor count, exact verification result, and stable failure code.

Status tokens:

- `SL_CORE_PASS`: determinant-one block factored and exactly verified.
- `GL_CERT_PASS`: monomial-unit block certified and exactly verified.
- `UNSUPPORTED_STAGED`: Suslin reached a staged unsupported algorithm boundary.
- `CACHE_ERROR`: cache case, ToricBuilder checkout, package loading, or cache
  payload is unavailable or unusable.
- `NORMALIZATION_PASS`, `NORMALIZATION_SKIP`, `NORMALIZATION_FAIL`.
- `GL_CERT_SKIP`, `SL_CORE_SKIP`, and `NONE` for non-applicable fields.

`status=PASS` means exact verification succeeded. `status=WARN` means a stable
unsupported staged boundary was reached without claiming a factorization.
`status=FAIL` means cache or unexpected factorization failure.

Create `example/toric_decoupling/runtests.jl` as an optional local test command.
It will always verify that missing cache inputs produce `CACHE_ERROR` without
throwing. When the local ToricBuilder checkout and cache are available, it will
also run `case_001` and `case_004`, checking that deterministic rows are
printed and that supported determinant-one blocks report `SL_CORE_PASS` with
`verified=true`.

## Verification

Focused optional commands:

```bash
julia --project=example/toric_decoupling example/toric_decoupling/runtests.jl
julia --project=example/toric_decoupling example/toric_decoupling/try_column_block_decoupling.jl --case=case_001,case_004
julia --project=example/toric_decoupling example/toric_decoupling/try_column_block_decoupling.jl --case=case_013
```

Package verification remains:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Out Of Scope

Do not add ToricBuilder to Suslin's main `Project.toml`, CI, or main test
suite. Do not require every local cache case or every pair mix to pass. Do not
copy local cache artifacts into this repository.

## Spec Self-Review

- No placeholders or incomplete sections remain.
- The design keeps ToricBuilder optional and local.
- The command surface is limited to the two requested example files plus the
  optional example project file needed by `--project=example/toric_decoupling`.
- The output row has stable tokens for supported, unsupported, and unavailable
  paths.
