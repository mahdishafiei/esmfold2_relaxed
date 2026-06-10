#!/bin/bash
# One-time setup: creates venv and downloads model weights (~26 GB)
# Usage: bash setup.sh
set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(dirname "$SCRIPT_DIR")"
VENV_DIR="$SCRIPT_DIR/env"
WEIGHTS_DIR="$SCRIPT_DIR/weights"

echo "=== ESMFold2 Antibody Complex — Setup ==="
echo "Repo:    $REPO_DIR"
echo "Venv:    $VENV_DIR"
echo "Weights: $WEIGHTS_DIR"
echo ""

# Check Python 3.12
PY=$(python3 --version 2>&1)
if [[ "$PY" != *"3.12"* ]]; then
  echo "ERROR: Python 3.12 required (found: $PY)"
  echo "Install with: conda install python=3.12 or use pyenv"
  exit 1
fi

# 1. Create virtualenv
echo "[1/4] Creating Python 3.12 virtualenv..."
python3 -m venv "$VENV_DIR"
source "$VENV_DIR/bin/activate"
pip install --upgrade pip -q

# 2. Install the ESM package from this repo
echo "[2/4] Installing ESM package and dependencies..."
pip install -e "$REPO_DIR" -q

# 3. Build xformers from source with CUDA ops (requires nvcc)
echo "[3/4] Building xformers from source with CUDA ops (~15 min)..."
CUDA_HOME=${CUDA_HOME:-$(dirname $(dirname $(which nvcc 2>/dev/null)) 2>/dev/null)} \
FORCE_CUDA=1 \
TORCH_CUDA_ARCH_LIST="8.0;8.6;8.9;9.0" \
MAX_JOBS=8 \
  pip install --no-build-isolation --force-reinstall \
    git+https://github.com/facebookresearch/xformers.git@main -q || \
  echo "  Warning: xformers build failed — predictions still work via PyTorch fallback."

# 4. Download model weights
echo "[4/4] Downloading model weights (~26 GB)..."
mkdir -p "$WEIGHTS_DIR"

WEIGHTS_DIR="$WEIGHTS_DIR" python3 - << 'PYEOF'
import os
from huggingface_hub import snapshot_download

weights_dir = os.environ["WEIGHTS_DIR"]

print("  Downloading ESMFold2-Fast (~900 MB)...")
snapshot_download("biohub/ESMFold2-Fast",
                  local_dir=os.path.join(weights_dir, "ESMFold2-Fast"))

print("  Downloading ESMC-6B backbone (~25 GB)...")
snapshot_download("biohub/ESMC-6B",
                  local_dir=os.path.join(weights_dir, "ESMC-6B"))
PYEOF

echo ""
echo "=== Setup complete ==="
echo ""
echo "Quick test:"
echo "  cd $SCRIPT_DIR"
echo "  ./fold.sh --test"
echo ""
echo "Run predictions:"
echo "  ./fold.sh --heavy VH_SEQ --light VL_SEQ --antigen AG_SEQ"
