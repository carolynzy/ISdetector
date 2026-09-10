#!/usr/bin/env bash
# Install into a dedicated environment and publish a user-local launcher.
set -Eeuo pipefail
trap 'status=$?; echo "ERROR: installation failed (exit $status) at line $LINENO." >&2; exit "$status"' ERR

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
INSTALL_PREFIX="${HOME}/.local/share/isdetector/env"
LAUNCHER_DIR="${HOME}/.local/bin"
EDITABLE=false
usage() {
    cat <<'EOF'
Usage: ./install.sh [--editable] [--prefix PATH] [--launcher-dir PATH]

Always create a new dedicated Conda/Mamba environment with Python, BWA,
BLAST+, Samtools >=1.3, and ISdetector. No shell activation is needed to run
its launcher. Conda or Mamba must already be installed.

  --editable          Link Python source to this repository for development.
  --prefix PATH       New environment (default: ~/.local/share/isdetector/env).
  --launcher-dir PATH Launcher directory (default: ~/.local/bin).
  --isolated          Accepted for compatibility; isolation is always enabled.
  --help              Show this help.

Existing environment directories are never overwritten. To reinstall, choose
another --prefix; a launcher created by this installer is updated on success.
Machine-specific tool paths are stored inside the environment, outside Git.
EOF
}
while (($#)); do
    case "$1" in
        --editable) EDITABLE=true; shift ;;
        --isolated) shift ;;
        --prefix|--launcher-dir)
            if (($# < 2)) || [[ -z "$2" ]]; then
                echo "ERROR: $1 requires a path." >&2; exit 2
            fi
            if [[ "$1" == --prefix ]]; then INSTALL_PREFIX=$2; else LAUNCHER_DIR=$2; fi
            shift 2 ;;
        --help|-h) usage; exit 0 ;;
        *) echo "ERROR: unknown option: $1" >&2; usage; exit 2 ;;
    esac
done
[[ "$INSTALL_PREFIX" == /* ]] || INSTALL_PREFIX="$PWD/$INSTALL_PREFIX"
[[ "$LAUNCHER_DIR" == /* ]] || LAUNCHER_DIR="$PWD/$LAUNCHER_DIR"
MANAGER=$(command -v mamba || command -v conda || true)
if [[ -z "$MANAGER" ]]; then
    echo 'ERROR: Install Conda or Mamba first, then rerun this installer.' >&2; exit 1
fi
[[ "$MANAGER" == /* ]] || MANAGER="$PWD/$MANAGER"
if [[ -e "$INSTALL_PREFIX" || -L "$INSTALL_PREFIX" ]]; then
    echo "ERROR: $INSTALL_PREFIX already exists. Choose a new --prefix; nothing was changed." >&2; exit 1
fi
LAUNCHER="$LAUNCHER_DIR/isdetector"
if [[ -e "$LAUNCHER" || -L "$LAUNCHER" ]]; then
    if [[ -L "$LAUNCHER" ]] || ! grep -q '^# ISdetector managed launcher$' "$LAUNCHER"; then
        echo "ERROR: Refusing to replace unrelated launcher $LAUNCHER. Choose --launcher-dir." >&2; exit 1
    fi
fi
# Resolve existing parent paths and reject private generated files inside the repo.
# realpath -m also handles the not-yet-created environment directory.
INSTALL_PREFIX=$(realpath -m -- "$INSTALL_PREFIX")
LAUNCHER_DIR=$(realpath -m -- "$LAUNCHER_DIR")
case "$INSTALL_PREFIX/" in "$SCRIPT_DIR/"*) echo 'ERROR: Choose an environment outside the repository.' >&2; exit 1 ;; esac
case "$LAUNCHER_DIR/" in "$SCRIPT_DIR/"*) echo 'ERROR: Choose a launcher directory outside the repository.' >&2; exit 1 ;; esac
LAUNCHER="$LAUNCHER_DIR/isdetector"
echo "Creating dedicated environment: $INSTALL_PREFIX"
"$MANAGER" create --yes --prefix "$INSTALL_PREFIX" --override-channels \
    --strict-channel-priority --channel conda-forge --channel bioconda \
    'python>=3.10' pip bwa 'samtools>=1.3' blast
export PATH="$INSTALL_PREFIX/bin:$PATH"
PYTHON_BIN="$INSTALL_PREFIX/bin/python"
PIP_ARGS=(install)
if "$EDITABLE"; then PIP_ARGS+=(-e); fi
"$PYTHON_BIN" -m pip "${PIP_ARGS[@]}" "$SCRIPT_DIR"
# Only publish configuration and launcher after successful installation/validation.
"$PYTHON_BIN" "$SCRIPT_DIR/scripts/configure_installation.py" "$INSTALL_PREFIX" "$MANAGER" "$LAUNCHER"
echo 'Installation complete. Run: isdetector'
printf 'Launcher: %s\n' "$LAUNCHER"
if [[ "::${PATH}:" != *":$LAUNCHER_DIR:"* ]]; then
    echo 'Add the launcher directory to PATH in your shell startup file:'
    printf '  export PATH=%q:"$PATH"\n' "$LAUNCHER_DIR"
fi
echo 'If another isdetector comes first in PATH, use the launcher path above or put its directory first.'
