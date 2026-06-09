#!/bin/bash
# One-time setup: creates venv and downloads model weights
# Run once before first use: bash setup.sh
set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(dirname "$SCRIPT_DIR")"
VENV_DIR="$SCRIPT_DIR/env"

echo "=== ESMFold2 Antibody Complex — Setup ==="
echo ""

# 1. Create virtualenv
echo "[1/4] Creating Python virtualenv..."
python3 -m venv "$VENV_DIR"
source "$VENV_DIR/bin/activate"
pip install --upgrade pip -q

# 2. Install the esm package from this repo
echo "[2/4] Installing ESM package..."
pip install -e "$REPO_DIR" -q

# 3. Install xformers (required for numerical accuracy)
echo "[3/4] Installing xformers..."
pip install xformers -q

# 4. Download weights
echo "[4/4] Downloading model weights (~26 GB)..."
python - << 'PYEOF'
from huggingface_hub import snapshot_download
import os

weights_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "weights")
os.makedirs(weights_dir, exist_ok=True)

print("  Downloading ESMFold2-Fast (~900 MB)...")
snapshot_download("biohub/ESMFold2-Fast", local_dir=f"{weights_dir}/ESMFold2-Fast")

print("  Downloading ESMC-6B backbone (~25 GB)...")
snapshot_download("biohub/ESMC-6B", local_dir=f"{weights_dir}/ESMC-6B")

print("  Done.")
PYEOF

echo ""
echo "=== Setup complete ==="
echo "Run predictions with: ./fold.sh --heavy SEQ --light SEQ --antigen SEQ"
