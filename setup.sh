#!/bin/bash
# setup.sh — one-time bootstrap, per machine.  bash setup.sh
#
# Builds:
#   1. the fold venv   ($HOME/.ef2_venv)        torch + this repo's esm  (the validated stack)
#   2. the DockQ venv  ($HOME/.ef2_dockq_venv)  DockQ — MUST be separate, it pins numpy<2
#   3. the weights     ($HOME/.ef2_hf_cache)    biohub/ESMFold2-Fast + biohub/ESMC-6B, ~26 GB
#
# Idempotent: anything already working is left alone, and the weights are never
# re-downloaded if a cache already has them.
#
# Overrides (export before running):
#   EF2_VENV_DIR    fold venv path        (default $HOME/.ef2_venv)
#   EF2_DOCKQ_VENV  DockQ venv path       (default $HOME/.ef2_dockq_venv)
#   EF2_HF_HOME     weights cache path    (default $HOME/.ef2_hf_cache, or <repo>/hf_cache if present)
#   EF2_PYTHON      python for uv         (default cpython-3.12.13; esm requires 3.12)
#   EF2_TORCH       torch pin             (default torch==2.13.0, validated with cu130)
#   EF2_FULL_MODEL=1  also fetch biohub/ESMFold2 (the 48-layer model; not used by the recipe)
set -e
REPO="$(cd "$(dirname "$0")" && pwd)"
VENV="${EF2_VENV_DIR:-$HOME/.ef2_venv}"
DOCKQ_VENV="${EF2_DOCKQ_VENV:-$HOME/.ef2_dockq_venv}"
PYSPEC="${EF2_PYTHON:-cpython-3.12.13}"
TORCH_PIN="${EF2_TORCH:-torch==2.13.0}"

# --- where the weights go: reuse any cache that already has them ---
HF_TARGET=""
for h in "$EF2_HF_HOME" "$HOME/.ef2_hf_cache" "$REPO/hf_cache"; do
  if [ -n "$h" ] && [ -d "$h/hub" ]; then HF_TARGET="$h"; break; fi
done
HF_TARGET="${HF_TARGET:-${EF2_HF_HOME:-$HOME/.ef2_hf_cache}}"

# --- uv (fast, and gives us a managed CPython 3.12); plain venv+pip is the fallback ---
USE_UV=1
if ! command -v uv >/dev/null 2>&1; then
  curl -LsSf https://astral.sh/uv/install.sh | sh >/dev/null 2>&1 || USE_UV=0
  export PATH="$HOME/.local/bin:$PATH"
  command -v uv >/dev/null 2>&1 || USE_UV=0
fi
export UV_LINK_MODE=copy   # repo and venv may be on different filesystems

# --- 1. fold venv -------------------------------------------------------------------
# esm is installed from a temp copy of the source: it keeps the checkout free of build
# artefacts and works even when the repo lives on read-only/network storage.
if [ -x "$VENV/bin/python" ] && "$VENV/bin/python" -c "import esm, torch" 2>/dev/null; then
  echo "[setup] fold venv ok: $VENV"
else
  echo "[setup] building fold venv at $VENV ..."
  ESMTMP="$(mktemp -d)"
  cp -a "$REPO/esm" "$REPO/pyproject.toml" "$REPO/README.md" "$REPO/LICENSE.md" "$ESMTMP/"
  rm -rf "$ESMTMP"/*.egg-info 2>/dev/null || true
  if [ "$USE_UV" = 1 ]; then
    uv venv --python "$PYSPEC" --python-preference only-managed "$VENV"
    uv pip install --python "$VENV/bin/python" "$TORCH_PIN" --torch-backend=auto
    uv pip install --python "$VENV/bin/python" "$ESMTMP"
  else
    PY_OK=$(python3 -c 'import sys; print(int(sys.version_info[:2]==(3,12)))' 2>/dev/null || echo 0)
    [ "$PY_OK" = 1 ] || { echo "ERROR: need python3.12 (esm requires >=3.12,<3.13), or install uv."; exit 1; }
    python3 -m venv "$VENV"
    "$VENV/bin/pip" install -q --upgrade pip
    "$VENV/bin/pip" install -q "$TORCH_PIN"
    "$VENV/bin/pip" install -q "$ESMTMP"
  fi
  rm -rf "$ESMTMP"
  "$VENV/bin/python" -c "import torch, esm, transformers, numpy; \
print('[setup] fold venv: torch', torch.__version__, '| transformers', transformers.__version__, \
'| numpy', numpy.__version__)"
fi

# --- 2. DockQ venv (optional but recommended; without it the DockQ columns are skipped) ---
if [ -x "$DOCKQ_VENV/bin/DockQ" ]; then
  echo "[setup] DockQ venv ok: $DOCKQ_VENV"
else
  echo "[setup] building DockQ venv at $DOCKQ_VENV ..."
  if [ "$USE_UV" = 1 ]; then
    uv venv --python "$PYSPEC" --python-preference only-managed "$DOCKQ_VENV"
    uv pip install --python "$DOCKQ_VENV/bin/python" DockQ
  else
    python3 -m venv "$DOCKQ_VENV"
    "$DOCKQ_VENV/bin/pip" install -q DockQ
  fi
  "$DOCKQ_VENV/bin/DockQ" --help >/dev/null 2>&1 \
    && echo "[setup] DockQ ok" || echo "[setup] WARN: DockQ did not install cleanly"
fi

# --- 3. weights ---------------------------------------------------------------------
export HF_HOME="$HF_TARGET"
mkdir -p "$HF_TARGET"
need_dl=0
for m in ESMFold2-Fast ESMC-6B; do
  [ -d "$HF_TARGET/hub/models--biohub--$m" ] || need_dl=1
done
if [ "$need_dl" = 0 ]; then
  echo "[setup] weights present: $HF_TARGET"
else
  echo "[setup] downloading weights -> $HF_TARGET (~26 GB, once) ..."
  "$VENV/bin/hf" download biohub/ESMFold2-Fast     # the folding model the recipe uses
  "$VENV/bin/hf" download biohub/ESMC-6B           # its language-model backbone (~25 GB)
fi
# esm loads a CCD chemistry dictionary at fold time — one 178 MB file that lives in the
# full ESMFold2 repo. Fetch it now, or the first prediction stalls on a download (and
# fails outright on a compute node without internet).
if ! ls "$HF_TARGET"/hub/models--biohub--ESMFold2/snapshots/*/ccd.pkl >/dev/null 2>&1; then
  echo "[setup] fetching ccd.pkl (178 MB) ..."
  "$VENV/bin/hf" download biohub/ESMFold2 ccd.pkl
fi
if [ "${EF2_FULL_MODEL:-0}" = 1 ]; then "$VENV/bin/hf" download biohub/ESMFold2; fi

cat <<EOF

[setup] done.
  source $REPO/env.sh
  bash $REPO/check.sh                      # confirm this machine is ready
  bash $REPO/predict.sh --target <dir with heavy.txt/light.txt/antigen.txt> --gpu 0
EOF
