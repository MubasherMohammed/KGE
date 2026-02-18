#!/bin/bash
#
# KGE Installer — MAGeCK with interactive HTML reports
#
# Usage:
#   bash install.sh            # Create new 'mageckenv' conda environment
#   bash install.sh myenv      # Create environment with custom name
#
set -euo pipefail

ENV_NAME="${1:-mageckenv}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "============================================"
echo "  KGE Installer — MAGeCK + Interactive HTML"
echo "============================================"
echo ""

# ---- Step 1: Check / Install conda ----
if ! command -v conda &> /dev/null; then
    echo "[0/4] conda not found — installing Miniconda..."

    # Detect OS and architecture
    OS="$(uname -s)"
    ARCH="$(uname -m)"
    case "${OS}" in
        Linux)  PLATFORM="Linux" ;;
        Darwin) PLATFORM="MacOSX" ;;
        *)      echo "ERROR: Unsupported OS: ${OS}"; exit 1 ;;
    esac
    case "${ARCH}" in
        x86_64)  MARCH="x86_64" ;;
        aarch64|arm64) MARCH="aarch64" ;;
        *)       echo "ERROR: Unsupported architecture: ${ARCH}"; exit 1 ;;
    esac

    MINICONDA_URL="https://repo.anaconda.com/miniconda/Miniconda3-latest-${PLATFORM}-${MARCH}.sh"
    MINICONDA_INSTALLER="/tmp/miniconda_installer.sh"

    echo "  Downloading: ${MINICONDA_URL}"
    curl -fsSL -o "${MINICONDA_INSTALLER}" "${MINICONDA_URL}"
    chmod +x "${MINICONDA_INSTALLER}"

    MINICONDA_DIR="${HOME}/miniconda3"
    echo "  Installing to: ${MINICONDA_DIR}"
    bash "${MINICONDA_INSTALLER}" -b -p "${MINICONDA_DIR}"
    rm -f "${MINICONDA_INSTALLER}"

    # Initialize conda for both bash and zsh (macOS defaults to zsh)
    eval "$("${MINICONDA_DIR}/bin/conda" shell.bash hook)"
    conda init bash 2>/dev/null || true
    conda init zsh  2>/dev/null || true

    echo "  Miniconda installed successfully."
    echo "  Shell integration configured for bash and zsh."
    echo "  (Restart your terminal for conda to be available in new sessions.)"
    echo ""
fi

# ---- Step 1b: Verify conda platform compatibility with MAGeCK ----
CONDA_PYTHON="$(conda info --base)/bin/python"
CONDA_PLATFORM="$(conda info --json 2>/dev/null | "${CONDA_PYTHON}" -c "import sys,json; print(json.load(sys.stdin)['platform'])" 2>/dev/null || echo "unknown")"
echo "  Detected conda platform: ${CONDA_PLATFORM}"

# MAGeCK 0.5.x from bioconda is available for linux-64 and osx-64.
# On osx-arm64 (Apple Silicon), we need CONDA_SUBDIR=osx-64 to install
# the x86_64 build under Rosetta 2 emulation.
case "${CONDA_PLATFORM}" in
    linux-64)
        echo "  Platform linux-64: MAGeCK supported natively."
        ;;
    osx-64)
        echo "  Platform osx-64: MAGeCK supported natively."
        ;;
    osx-arm64)
        echo "  Platform osx-arm64 (Apple Silicon): MAGeCK not available natively."
        echo "  Switching to osx-64 (Rosetta 2 emulation) for environment creation."
        export CONDA_SUBDIR=osx-64
        ;;
    linux-aarch64)
        echo "  Platform linux-aarch64: MAGeCK not available natively."
        echo "  Switching to linux-64 emulation for environment creation."
        export CONDA_SUBDIR=linux-64
        ;;
    *)
        echo "WARNING: Unrecognised platform '${CONDA_PLATFORM}'."
        echo "  Proceeding — environment creation may fail if MAGeCK is unavailable."
        ;;
esac

# ---- Step 2: Create conda environment ----
echo "[1/4] Creating conda environment '${ENV_NAME}'..."
if conda env list | grep -q "^${ENV_NAME} "; then
    echo "  Environment '${ENV_NAME}' already exists."
    read -p "  Remove and recreate? [y/N] " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        conda env remove -n "${ENV_NAME}" -y
    else
        echo "  Keeping existing environment. Will overlay modifications."
    fi
fi

if ! conda env list | grep -q "^${ENV_NAME} "; then
    conda env create -f "${SCRIPT_DIR}/environment.yml" -n "${ENV_NAME}"
fi

# ---- Step 3: Find site-packages path ----
echo "[2/4] Locating site-packages..."
SITE_PACKAGES=$(conda run -n "${ENV_NAME}" python -c "import site; print(site.getsitepackages()[0])" 2>/dev/null)
MAGECK_DIR="${SITE_PACKAGES}/mageck"

if [ ! -d "${MAGECK_DIR}" ]; then
    echo "ERROR: mageck package not found at ${MAGECK_DIR}"
    exit 1
fi
echo "  Found: ${MAGECK_DIR}"

# ---- Step 4: Overlay modified Python source ----
echo "[3/4] Installing KGE modifications..."
cp "${SCRIPT_DIR}"/mageck/*.py "${MAGECK_DIR}/"
cp "${SCRIPT_DIR}"/mageck/*.gmt "${MAGECK_DIR}/"
cp "${SCRIPT_DIR}"/mageck/*.txt "${MAGECK_DIR}/"
cp "${SCRIPT_DIR}"/mageck/*.Rnw "${MAGECK_DIR}/" 2>/dev/null || true
cp "${SCRIPT_DIR}"/mageck/*.Rmd "${MAGECK_DIR}/" 2>/dev/null || true
cp "${SCRIPT_DIR}"/mageck/*.RTemplate "${MAGECK_DIR}/" 2>/dev/null || true

# Clear bytecode cache so Python picks up new source
rm -rf "${MAGECK_DIR}/__pycache__"

echo "  Copied modified files to ${MAGECK_DIR}"

# ---- Step 5: Verify ----
echo "[4/4] Verifying installation..."
conda run -n "${ENV_NAME}" python -c "
from mageck.version import __version__
from mageck.htmlReport import mageck_report_main
import plotly
import decoupler
print(f'  MAGeCK version:    {__version__}')
print(f'  Plotly version:    {plotly.__version__}')
print(f'  Decoupler version: {decoupler.__version__}')
print(f'  HTML report:       OK')
" 2>/dev/null

echo ""
echo "============================================"
echo "  Installation complete!"
echo "============================================"
echo ""
echo "  Activate:  conda activate ${ENV_NAME}"
echo ""
echo "  Usage:"
echo "    mageck count -l library.txt -n demo --fastq *.fastq --html-report"
echo "    mageck test  -k counts.txt -t treat -c ctrl -n demo --html-report"
echo "    mageck mle   -k counts.txt -d design.txt -n demo --html-report"
echo "    mageck report -n demo -k counts.txt"
echo ""
