#!/usr/bin/env bash
set -euo pipefail
# usage: ./uninstall.sh [/path/to/prefix]
REPO_ROOT="$(cd -P "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SAVED_PREFIX=$(cat "$REPO_ROOT/.install_prefix" 2>/dev/null || echo "")
PREFIX="${1:-${SAVED_PREFIX:-$HOME/.local/bin}}"

shift || true
while [[ $# -gt 0 ]]; do
  case "$1" in
    -h|--help)
      echo "Usage: ./uninstall.sh [/path/to/prefix]"
      exit 0
      ;;
    *)
      echo "[ERROR] Unknown option: $1"
      exit 1
      ;;
  esac
  shift
done

SYMLINKS=("v2t-step1" "v2t-step2" "v2t-sra")

echo "[INFO] Removing symlinks from: $PREFIX"
for name in "${SYMLINKS[@]}"; do
  target="$PREFIX/$name"
  if [[ -L "$target" ]]; then
    rm -f "$target"
    echo "[OK] Removed $target"
  elif [[ -e "$target" ]]; then
    echo "[WARN] $target exists but is not a symlink. Skipping."
  else
    echo "[INFO] $target not found."
  fi
done

removed_any=false
comp_paths=(
  "$HOME/.local/share/bash-completion/completions/virus2tree"
  "/usr/share/bash-completion/completions/virus2tree"
)
for comp in "${comp_paths[@]}"; do
  if [[ -e "$comp" ]]; then
    rm -f "$comp"
    echo "[OK] Removed $comp"
    removed_any=true
  fi
done
echo "[DONE] Uninstall complete."
