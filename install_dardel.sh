#!/bin/bash
#
# KGE Installer for Dardel (PDC/KTH HPC)
#
# Loads required modules, creates a conda environment under Klemming,
# then overlays the modified MAGeCK source files.
#
# Usage:
#   bash install_dardel.sh              # Create 'mageckenv' environment
#   bash install_dardel.sh myenv        # Create environment with custom name
#
set -euo pipefail

ENV_NAME="${1:-mageckenv}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "============================================"
echo "  KGE Installer — Dardel (PDC/KTH)"
echo "============================================"
echo ""

# ---- Step 1: Load Dardel modules ----
echo "[1/5] Loading Dardel modules..."
if ! command -v ml &> /dev/null && ! command -v module &> /dev/null; then
    echo "ERROR: module system not found. Are you on Dardel?"
    echo "  This installer is for PDC Dardel. Use install.sh for local machines."
    exit 1
fi

# Load PDC and miniconda3 modules
# Adjust versions if needed for your Dardel allocation
ml PDC/23.12
ml miniconda3/25.3.1-1-cpeGNU-24.11

echo "  Loaded PDC and miniconda3 modules"

# Verify conda is available after module load
if ! command -v conda &> /dev/null; then
    echo "ERROR: conda not available after loading modules."
    echo "  Check available versions with: ml spider miniconda3"
    exit 1
fi

# ---- Step 2: Create conda environment ----
echo "[2/5] Creating conda environment '${ENV_NAME}'..."
if conda env list 2>/dev/null | grep -q "^${ENV_NAME} "; then
    echo "  Environment '${ENV_NAME}' already exists."
    read -p "  Remove and recreate? [y/N] " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        conda env remove -n "${ENV_NAME}" -y
    else
        echo "  Keeping existing environment. Will overlay modifications."
    fi
fi

if ! conda env list 2>/dev/null | grep -q "^${ENV_NAME} "; then
    conda env create -f "${SCRIPT_DIR}/environment.yml" -n "${ENV_NAME}"
fi

# ---- Step 3: Activate and find site-packages ----
echo "[3/5] Locating site-packages..."

# On Dardel, 'conda run' may not work reliably — use source activate
# Save current state so we can provide activation instructions
SITE_PACKAGES=$(conda run -n "${ENV_NAME}" python -c "import site; print(site.getsitepackages()[0])" 2>/dev/null || true)

# Fallback: try source activate if conda run failed
if [ -z "${SITE_PACKAGES}" ] || [ ! -d "${SITE_PACKAGES}" ]; then
    echo "  conda run failed, trying source activate..."
    source activate "${ENV_NAME}" 2>/dev/null || conda activate "${ENV_NAME}" 2>/dev/null || true
    SITE_PACKAGES=$(python -c "import site; print(site.getsitepackages()[0])" 2>/dev/null)
fi

MAGECK_DIR="${SITE_PACKAGES}/mageck"

if [ ! -d "${MAGECK_DIR}" ]; then
    echo "ERROR: mageck package not found at ${MAGECK_DIR}"
    echo "  Try activating manually: source activate ${ENV_NAME}"
    exit 1
fi
echo "  Found: ${MAGECK_DIR}"

# ---- Step 4: Overlay modified Python source ----
echo "[4/5] Installing KGE modifications..."
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
echo "[5/5] Verifying installation..."
conda run -n "${ENV_NAME}" python -c "
from mageck.version import __version__
from mageck.htmlReport import mageck_report_main
import plotly
import decoupler
print(f'  MAGeCK version:    {__version__}')
print(f'  Plotly version:    {plotly.__version__}')
print(f'  Decoupler version: {decoupler.__version__}')
print(f'  HTML report:       OK')
" 2>/dev/null || {
    # Fallback verification with source activate
    source activate "${ENV_NAME}" 2>/dev/null || true
    python -c "
from mageck.version import __version__
from mageck.htmlReport import mageck_report_main
import plotly
import decoupler
print(f'  MAGeCK version:    {__version__}')
print(f'  Plotly version:    {plotly.__version__}')
print(f'  Decoupler version: {decoupler.__version__}')
print(f'  HTML report:       OK')
"
}

echo ""
echo "============================================"
echo "  Installation complete!"
echo "============================================"
echo ""
echo "  To use on Dardel, first load modules then activate:"
echo ""
echo "    ml PDC/23.12"
echo "    ml miniconda3/25.3.1-1-cpeGNU-24.11"
echo "    source activate ${ENV_NAME}"
echo ""
echo "  Usage:"
echo "    mageck count -l library.txt -n demo --fastq *.fastq --html-report"
echo "    mageck test  -k counts.txt -t treat -c ctrl -n demo --html-report"
echo "    mageck mle   -k counts.txt -d design.txt -n demo --html-report"
echo "    mageck report -n demo -k counts.txt"
echo ""
