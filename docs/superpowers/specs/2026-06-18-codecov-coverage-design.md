# Codecov Coverage Design

## Goal

Show test coverage for Suslin.jl in README and on pull requests through Codecov.

## Approach

The existing full-suite GitHub Actions job remains the source of truth for
coverage. It runs `test/runtests.jl all` with Julia coverage enabled, converts
Julia `.cov` files under `src` to `lcov.info` with
`julia-actions/julia-processcoverage@v1`, and uploads that report with
`codecov/codecov-action@v6`.

README displays the existing CI status and the Codecov badge for the `main`
branch. `codecov.yml` enables a PR comment using Codecov's standard diff,
flags, and files layout.

## Operational Notes

Private repositories must define the GitHub Actions secret `CODECOV_TOKEN`.
Public repositories may be able to upload without a token depending on Codecov
repository settings, but this workflow always passes the secret when present.

## Verification

Local verification covers YAML parseability and confirms that the workflow
contains the coverage generation, Julia coverage processing, and Codecov upload
steps. The actual upload is verified by GitHub Actions after the secret is set.
