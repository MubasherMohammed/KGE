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

# ---- Step 1: Check conda ----
if ! command -v conda &> /dev/null; then
    echo "ERROR: conda not found. Install Miniconda/Anaconda first."
    echo "  https://docs.conda.io/en/latest/miniconda.html"
    exit 1
fi

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
