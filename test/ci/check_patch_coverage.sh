#!/usr/bin/env bash
set -euo pipefail

base="origin/main"
if [[ "${1:-}" == --base=* ]]; then
  base="${1#--base=}"
  shift
fi

if [[ "$#" -eq 0 ]]; then
  echo "no LCOV reports supplied" >&2
  exit 2
fi

git_root="$(git rev-parse --show-toplevel)"
git_root="${git_root%/}"
physical_git_root="$(cd "$git_root" && pwd -P)"
alternate_git_root=""
if [[ "$physical_git_root" == /private/var/* ]]; then
  alternate_git_root="/${physical_git_root#/private/}"
fi
tmpdir="$(mktemp -d)"
trap 'rm -rf "$tmpdir"' EXIT
reports=()
index=0
for report in "$@"; do
  index=$((index + 1))
  normalized_report="$tmpdir/report-$index.info"
  awk -v root1="$git_root" \
      -v root2="$physical_git_root" \
      -v root3="$alternate_git_root" '
    function strip_root(root) {
      if (root == "") {
        return 0
      }
      prefix = "SF:" root "/"
      if (index($0, prefix) == 1) {
        print "SF:" substr($0, length(prefix) + 1)
        return 1
      }
      return 0
    }
    strip_root(root1) || strip_root(root2) || strip_root(root3) { next }
    { print }
  ' "$report" > "$normalized_report"
  reports+=("$normalized_report")
done

if [[ -n "${DIFF_COVER_BIN:-}" ]]; then
  read -r -a diff_cover_command <<<"$DIFF_COVER_BIN"
elif command -v uvx >/dev/null 2>&1; then
  diff_cover_command=(uvx --from diff-cover==10.1.0 diff-cover)
else
  echo "uvx is required; install uv from https://docs.astral.sh/uv/" >&2
  exit 2
fi

"${diff_cover_command[@]}" "${reports[@]}" \
  --compare-branch="$base" \
  --include='src/**' \
  --show-uncovered \
  --fail-under=99
