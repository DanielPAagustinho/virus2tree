#!/usr/bin/env bash
set -euo pipefail
# usage: ./install.sh [/path/to/prefix]
PREFIX="${1:-/usr/local/bin}"  #prefix is the place where symlinks will be created 
BASENAMES=("virus2tree_step1.sh:v2t-step1"
           "virus2tree_step2.sh:v2t-step2"
           "download_sra_reads.sh:v2t-sra")

# Resolve install.sh absolute directory, following symlinks
SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do
  DIR="$(cd -P "$(dirname "$SOURCE")" && pwd)"
  SOURCE="$(readlink "$SOURCE")"
  [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE"
done
REPO_ROOT="$(cd -P "$(dirname "$SOURCE")" && pwd)"
SCRIPTS_DIR="$REPO_ROOT/scripts" #here will be the scripts to be installed

# Check the three targets exist and are executable
for pair in "${BASENAMES[@]}"; do
  src="${pair%%:*}"
  target="$SCRIPTS_DIR/$src"

  # Must exist
  [[ -f "$target" ]] || { echo "[ERROR] Missing file: $target"; exit 1; }

  # Ensure executable bit; try to fix if not
  if [[ ! -x "$target" ]]; then
    if [[ -w "$target" ]]; then
      chmod u+x "$target" || true
    fi
  fi

  # Final check
  [[ -x "$target" ]] || {
    echo "[ERROR] Not executable and could not fix perms: $target"
    echo "        Try: chmod +x \"$target\""
    exit 1
  }
done

# Create the prefix if it's a user-owned location
if ! mkdir -p "$PREFIX" 2>/dev/null; then
  echo "[ERROR] Can't create $PREFIX (need permissions?)."
  echo "        Try: sudo ./install.sh   or   ./install.sh \"\$HOME/.local/bin\""
  exit 1
fi
if [[ ! -w "$PREFIX" ]]; then
  echo "[ERROR] $PREFIX is not writable."
  echo "        Try: sudo ./install.sh   or   ./install.sh \"\$HOME/.local/bin\""
  exit 1
fi

# Create absolute symlinks
for pair in "${BASENAMES[@]}"; do
  src="${pair%%:*}"
  dst="${pair##*:}"
  # refuse to clobber a regular file
  if [[ -e "$PREFIX/$dst" && ! -L "$PREFIX/$dst" ]]; then
    echo "[ERROR] $PREFIX/$dst exists and is not a symlink. Remove it first."
    exit 1
  fi
  ln -sfn "$SCRIPTS_DIR/$src" "$PREFIX/$dst"
  echo "[OK] $PREFIX/$dst -> $SCRIPTS_DIR/$src"
done

# Install bash completions
COMP_SRC="$REPO_ROOT/share/bash-completion/completions/virus2tree"
COMP_INSTALLED=false
# install bash completions 
if [[ -f "$COMP_SRC" ]]; then
  if [[ "${EUID:-$(id -u)}" -eq 0 ]]; then
    COMP_SHARE="/usr/share/bash-completion/completions"
  else
    COMP_SHARE="$HOME/.local/share/bash-completion/completions"
  fi
  
  mkdir -p "$COMP_SHARE" 2>/dev/null || true
  
  if [[ -d "$COMP_SHARE" && -w "$COMP_SHARE" ]]; then
    # Copiar el archivo principal
    if cp -f "$COMP_SRC" "$COMP_SHARE/virus2tree" 2>/dev/null; then
      echo "[OK] Bash completions installed to $COMP_SHARE/virus2tree"
      COMP_INSTALLED=true
      
      # Crear symlinks para que bash-completion lo encuentre por cada comando
      for cmd in v2t-step1 v2t-step2 v2t-sra; do
        ln -sf virus2tree "$COMP_SHARE/$cmd" 2>/dev/null || true
      done
      echo "[OK] Created completion symlinks: v2t-step1, v2t-step2, v2t-sra"
    fi
  fi
fi

# Warn if the prefix is not in PATH
case ":$PATH:" in
  *":$PREFIX:"*) echo "[OK] $PREFIX is in PATH.";;
  *) echo "[WARN] $PREFIX is not in PATH. Add this to your shell rc:"
     echo "export PATH=\"$PREFIX:\$PATH\"";;
esac

echo "$PREFIX" > "$REPO_ROOT/.install_prefix"

echo "[DONE] Installed executables: v2t-step1, v2t-step2, v2t-sra"
if [[ "$COMP_INSTALLED" == true ]]; then
  echo "[INFO] To enable completions now, run:"
  echo "  source ~/.bashrc"
fi
echo "[INFO] To uninstall, run: ./uninstall.sh \"$PREFIX\""
