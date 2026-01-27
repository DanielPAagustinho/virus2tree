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

comp_symlinks=("v2t-step1" "v2t-step2" "v2t-sra")

for base_path in "$HOME/.local/share/bash-completion/completions" "/usr/share/bash-completion/completions"; do
  # Remover archivo principal
  if [[ -e "$base_path/virus2tree" ]]; then
    rm -f "$base_path/virus2tree"
    echo "[OK] Removed $base_path/virus2tree"
    removed_any=true
  fi
  
  for link in "${comp_symlinks[@]}"; do
    if [[ -L "$base_path/$link" ]]; then
      rm -f "$base_path/$link"
      echo "[OK] Removed completion symlink $base_path/$link"
      removed_any=true
    fi
  done
done

echo "[DONE] Uninstall complete."