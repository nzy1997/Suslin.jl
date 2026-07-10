#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$repo_root"

base="origin/main"
if [[ "${1:-}" == --base=* ]]; then
  base="${1#--base=}"
  shift
fi
if [[ "$#" -ne 0 ]]; then
  echo "usage: $0 [--base=<revision>]" >&2
  exit 2
fi
if [[ -z "${DIFF_COVER_BIN:-}" ]] && ! command -v uvx >/dev/null 2>&1; then
  echo "uvx is required; install uv from https://docs.astral.sh/uv/" >&2
  exit 2
fi

targets="$(${JULIA:-julia} --startup-file=no --project=. \
  test/ci/select_shards.jl --base="$base" --head=HEAD --format=lines)"
[[ -n "$targets" ]] || {
  echo "selector returned no targets" >&2
  exit 2
}

tmpdir="$(mktemp -d)"
trap 'rm -rf "$tmpdir"' EXIT
reports=()

while IFS= read -r target; do
  [[ -n "$target" ]] || continue
  if [[ "$target" == documentation-smoke ]]; then
    "${JULIA:-julia}" --startup-file=no --project=. \
      test/runtests.jl documentation-smoke
    continue
  fi
  report="$tmpdir/${target}.info"
  "${JULIA:-julia}" --startup-file=no --project=. \
    --code-coverage="$report" test/runtests.jl "shard:$target"
  reports+=("$report")
done <<<"$targets"

if [[ "${#reports[@]}" -eq 0 ]]; then
  exit 0
fi
test/ci/check_patch_coverage.sh --base="$base" "${reports[@]}"
