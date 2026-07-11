#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
tmpdir="$(mktemp -d)"
trap 'rm -rf "$tmpdir"' EXIT

cd "$tmpdir"
git init -q
git config user.name "CI Coverage Test"
git config user.email "ci-coverage@example.invalid"
mkdir -p src
printf 'x = 1\n' > src/example.jl
git add src/example.jl
git commit -qm "base"
printf 'x = 2\n' > src/example.jl
git add src/example.jl
git commit -qm "change"

source_file="$tmpdir/src/example.jl"
printf 'TN:\nSF:%s\nDA:1,1\nend_of_record\n' "$source_file" > covered.info
printf 'TN:\nSF:%s\nDA:1,0\nend_of_record\n' "$source_file" > uncovered.info

DIFF_COVER_BIN="uvx --from diff-cover==10.1.0 diff-cover" \
  "$repo_root/test/ci/check_patch_coverage.sh" --base=HEAD~1 covered.info

if DIFF_COVER_BIN="uvx --from diff-cover==10.1.0 diff-cover" \
  "$repo_root/test/ci/check_patch_coverage.sh" --base=HEAD~1 uncovered.info; then
  echo "expected uncovered patch to fail" >&2
  exit 1
fi

fake_julia="$tmpdir/fake-julia"
cat > "$fake_julia" <<'FAKE_JULIA'
#!/usr/bin/env bash
case " $* " in
  *" test/ci/select_shards.jl "*)
    printf 'documentation-smoke\n'
    ;;
  *" test/runtests.jl documentation-smoke"*)
    ;;
  *)
    echo "unexpected fake julia invocation: $*" >&2
    exit 1
    ;;
esac
FAKE_JULIA
chmod +x "$fake_julia"

PATH="/usr/bin:/bin" DIFF_COVER_BIN=true JULIA="$fake_julia" \
  "$repo_root/test/ci/coverage_changed.sh" --base=HEAD
