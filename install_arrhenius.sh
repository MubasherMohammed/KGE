#!/bin/bash
#
#SBATCH -n 1
#SBATCH -c 8
#SBATCH -t 01:00:00
#SBATCH -A naiss2026-4-1258-cpu
#SBATCH -p cpu
#SBATCH --job-name=kge-install
#SBATCH --output=kge-install-%j.out
#
# KGE Installer for Arrhenius HPC (NAISS)
#
# Runs either interactively on a login node or as a SLURM job:
#   sbatch install_arrhenius.sh [env_name]
#
# Unlike Dardel, Arrhenius does NOT provide a miniconda3 module.
# This installer:
#   1. Tries the Miniforge module (provided by Lmod on Arrhenius).
#   2. If that is unavailable, downloads a standalone Miniforge
#      (conda-compatible) build matched to the HPC platform
#      (Linux x86_64 / EPYC 9005) into $HOME/miniforge3.
#   3. Creates a conda environment (under ENV_PREFIX) that contains
#      MAGeCK AND the Illumina demultiplexer bcl2fastq2.
#   4. Overlays the modified KGE MAGeCK source files.
#
# Usage (cd into the KGE folder first, it contains install_arrhenius.sh,
# environment_arrhenius.yml and mageck/):
#   sbatch install_arrhenius.sh              # Create 'mageckenv' environment
#   sbatch install_arrhenius.sh myenv        # Create environment with custom name
#
set -euo pipefail

ENV_NAME="${1:-mageckenv}"

# Install the environment under the user's nobackup project space.
# Override with: ENV_PREFIX=/some/path sbatch install_arrhenius.sh
# (CONDA_ENVS_PATH / CONDA_PKGS_DIRS are exported AFTER the Miniforge module
#  is loaded, because the module sets its own HPC-specific conda paths.)
ENV_PREFIX="${ENV_PREFIX:-/nobackup/proj/disk/crispr-fun/personal/mubasher}"

# Resolve the directory containing the installer + environment yml + mageck/ overlay.
# Under sbatch, BASH_SOURCE[0] resolves into SLURM's private per-job directory
# (/var/lib/slurm/slurmd/job<id>/), which does NOT contain the companion files --
# so prefer SLURM_SUBMIT_DIR, i.e. cd into the KGE folder before submitting.
if [ -n "${SLURM_SUBMIT_DIR:-}" ]; then
    SCRIPT_DIR="${SLURM_SUBMIT_DIR}"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi

if [ ! -f "${SCRIPT_DIR}/environment_arrhenius.yml" ]; then
    echo "ERROR: environment_arrhenius.yml not found next to this script (${SCRIPT_DIR})."
    echo "  cd into the KGE folder that also contains environment_arrhenius.yml and mageck/,"
    echo "  then: sbatch install_arrhenius.sh"
    exit 1
fi

echo "============================================"
echo "  KGE Installer — Arrhenius HPC"
echo "============================================"
echo ""

# ---- Step 1: Get conda (Miniforge module or local install) ----
echo "[1/5] Setting up conda (Miniforge)..."

if ! command -v ml &> /dev/null && ! command -v module &> /dev/null; then
    echo "ERROR: module system not found. Are you on Arrhenius?"
    echo "  This installer is for Arrhenius HPC. Use install.sh for local machines."
    exit 1
fi

_has_conda() { command -v conda &> /dev/null; }

has_conda=false

# 1a. Try the Miniforge modules provided on Arrhenius
for mod in "Miniforge/26.3.2-2-eb" "Miniforge/recommendation" "Miniforge3/26.3.2-2"; do
    if ml show "${mod}" &> /dev/null || module show "${mod}" &> /dev/null; then
        echo "  Loading module: ${mod}"
        ml "${mod}" 2>/dev/null || module load "${mod}" 2>/dev/null || true
        if _has_conda; then
            has_conda=true
            break
        fi
    fi
done

# 1b. Fallback: download Miniforge compatible with the HPC platform
if ! $has_conda; then
    echo "  No conda-capable module found — downloading Miniforge..."

    OS="$(uname -s)"
    ARCH="$(uname -m)"
    case "${OS}" in
        Linux)  PLATFORM="Linux" ;;
        *)      echo "ERROR: Unsupported OS: ${OS}"; exit 1 ;;
    esac
    case "${ARCH}" in
        x86_64)  MARCH="x86_64" ;;
        aarch64|arm64) MARCH="aarch64" ;;
        *)       echo "ERROR: Unsupported architecture: ${ARCH}"; exit 1 ;;
    esac

    MINIFORGE_URL="https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-${PLATFORM}-${MARCH}.sh"
    MINIFORGE_INSTALLER="/tmp/miniforge_installer.sh"

    echo "  Downloading: ${MINIFORGE_URL}"
    curl -fsSL -o "${MINIFORGE_INSTALLER}" "${MINIFORGE_URL}"
    chmod +x "${MINIFORGE_INSTALLER}"

    MINIFORGE_DIR="${HOME}/miniforge3"
    echo "  Installing to: ${MINIFORGE_DIR}"
    bash "${MINIFORGE_INSTALLER}" -b -p "${MINIFORGE_DIR}"
    rm -f "${MINIFORGE_INSTALLER}"

    eval "$("${MINIFORGE_DIR}/bin/conda" shell.bash hook)"
    conda init bash 2>/dev/null || true
fi

if ! _has_conda; then
    echo "ERROR: conda still not available. Run the installer again in a fresh shell."
    exit 1
fi

echo "  conda base: $(conda info --base)"

# Point envs/pkgs into the nobackup project space AFTER module load.
CONDA_ENVS_PATH="${CONDA_ENVS_PATH:-${ENV_PREFIX}/conda_envs}"
CONDA_PKGS_DIRS="${CONDA_PKGS_DIRS:-${ENV_PREFIX}/conda_pkgs}"
mkdir -p "${CONDA_ENVS_PATH}" "${CONDA_PKGS_DIRS}"
export CONDA_ENVS_PATH CONDA_PKGS_DIRS

# ---- Step 1b: Ensure linux-64 subdir (Arrhenius CPU nodes are x86_64/EPYC) ----
CONDA_PLATFORM="$(conda info --json 2>/dev/null | "$(conda info --base)/bin/python" -c "import sys,json; print(json.load(sys.stdin)['platform'])" 2>/dev/null || echo "unknown")"
echo "  Detected conda platform: ${CONDA_PLATFORM}"
if [ "${CONDA_PLATFORM}" != "linux-64" ]; then
    export CONDA_SUBDIR=linux-64
    echo "  Setting CONDA_SUBDIR=linux-64 for MAGeCK/bcl2fastq2 compatibility."
fi

# ---- Step 2: Create conda environment (with bcl2fastq2) ----
echo "[2/5] Creating conda environment '${ENV_NAME}' (MAGeCK + bcl2fastq2)..."
if conda env list 2>/dev/null | grep -q "^${ENV_NAME} "; then
    echo "  Environment '${ENV_NAME}' already exists."
    if [ "${KGE_FORCE_RECREATE:-0}" = "1" ]; then
        echo "  KGE_FORCE_RECREATE=1 — removing and recreating."
        conda env remove -n "${ENV_NAME}" -y
    elif [ -t 0 ]; then
        read -p "  Remove and recreate? [y/N] " -n 1 -r
        echo
        if [[ $REPLY =~ ^[Yy]$ ]]; then
            conda env remove -n "${ENV_NAME}" -y
        else
            echo "  Keeping existing environment. Will overlay modifications."
        fi
    else
        echo "  Non-interactive (batch): keeping existing environment."
        echo "  Set KGE_FORCE_RECREATE=1 to recreate:"
        echo "    KGE_FORCE_RECREATE=1 sbatch install_arrhenius.sh"
    fi
fi

if ! conda env list 2>/dev/null | grep -q "^${ENV_NAME} "; then
    conda env create -f "${SCRIPT_DIR}/environment_arrhenius.yml" -n "${ENV_NAME}"
fi

# ---- Step 3: Activate and find site-packages ----
echo "[3/5] Locating site-packages..."

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
if ! compgen -G "${SCRIPT_DIR}/mageck/*.py" > /dev/null; then
    echo "ERROR: KGE modified source not found at ${SCRIPT_DIR}/mageck/"
    echo "  Copy the mageck/ folder next to this script:"
    echo "    scp -r mageck mubasher@login.hpc.arrhenius.naiss.se:${SCRIPT_DIR}/"
    exit 1
fi
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

# Verify bcl2fastq2 (Illumina demultiplexer)
if conda run -n "${ENV_NAME}" which bcl2fastq &>/dev/null; then
    BCL2FASTQ_VER="$(conda run -n "${ENV_NAME}" bcl2fastq --version 2>/dev/null | head -n1 || true)"
    echo "  bcl2fastq:           ${BCL2FASTQ_VER:-found}"
elif command -v bcl2fastq &>/dev/null; then
    echo "  bcl2fastq:           found (already on PATH)"
else
    echo "  bcl2fastq:           NOT found — check that bih-cubi channel resolved correctly"
fi

echo ""
echo "============================================"
echo "  Installation complete!"
echo "============================================"
echo ""
echo "  To use on Arrhenius:"
echo ""
echo "    module load Miniforge/26.3.2-2-eb   # or the downloaded Miniforge after re-login"
echo "    export CONDA_ENVS_PATH=${CONDA_ENVS_PATH}"
echo "    export CONDA_PKGS_DIRS=${CONDA_PKGS_DIRS}"
echo "    source activate ${ENV_NAME}"
echo ""
echo "  Demultiplexing with bcl2fastq2:"
echo "    bcl2fastq --runfolder-dir <Run> --output-dir <Out> --barcode-mismatches 1"
echo ""
echo "  Usage:"
echo "    mageck count -l library.txt -n demo --fastq *.fastq --html-report"
echo "    mageck test  -k counts.txt -t treat -c ctrl -n demo --html-report"
echo "    mageck mle   -k counts.txt -d design.txt -n demo --html-report"
echo "    mageck report -n demo -k counts.txt"
echo ""