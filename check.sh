#!/bin/bash
# check.sh — is this machine ready to run predictions?  bash check.sh
REPO="$(cd "$(dirname "$0")" && pwd)"
ok=1
echo "host : $(hostname)"
echo "repo : $REPO"
echo

# --- fold venv ---
V=""
for c in "$EF2_VENV_DIR" "$HOME/.ef2_venv" "$REPO/.venv"; do
  [ -n "$c" ] && [ -x "$c/bin/python" ] && { V="$c"; break; }
done
if [ -n "$V" ] && "$V/bin/python" -c "import esm, torch" 2>/dev/null; then
  echo "✅ fold venv : $V  ($("$V/bin/python" -c 'import torch;print("torch",torch.__version__)'))"
else
  echo "❌ fold venv : none on this machine"; ok=0
fi

# --- weights ---
H=""
for c in "$EF2_HF_HOME" "$HOME/.ef2_hf_cache" "$REPO/hf_cache"; do
  [ -n "$c" ] && [ -d "$c/hub/models--biohub--ESMFold2-Fast" ] && { H="$c"; break; }
done
if [ -n "$H" ]; then
  echo "✅ weights   : $H"
  [ -d "$H/hub/models--biohub--ESMC-6B" ] || { echo "❌ ESMC-6B backbone missing under $H"; ok=0; }
  ls "$H"/hub/models--biohub--ESMFold2/snapshots/*/ccd.pkl >/dev/null 2>&1 \
    || echo "⚠️  ccd.pkl   : missing — the first fold will download 178 MB (fails with no internet; 'bash setup.sh' fetches it)"
else
  echo "❌ weights   : biohub/ESMFold2-Fast not found in any cache"; ok=0
fi

# --- DockQ (optional) ---
D=""
for c in "$EF2_DOCKQ_VENV/bin/DockQ" "$HOME/.ef2_dockq_venv/bin/DockQ" "$(command -v DockQ)"; do
  [ -n "$c" ] && [ -x "$c" ] && { D="$c"; break; }
done
[ -n "$D" ] && echo "✅ DockQ     : $D" \
            || echo "⚠️  DockQ     : not installed (scores.csv will have no DockQ columns)"

# --- GPU ---
if command -v nvidia-smi >/dev/null 2>&1; then
  echo "✅ GPUs      : $(nvidia-smi --query-gpu=index,name,memory.total --format=csv,noheader | tr '\n' ' | ')"
else
  echo "❌ GPUs      : nvidia-smi not found"; ok=0
fi

echo
if [ "$ok" = 1 ] && [ -n "$D" ]; then
  echo "READY ✅  →  bash predict.sh --target <target dir> --gpu 0    (see targets/README.md)"
elif [ "$ok" = 1 ]; then
  echo "MOSTLY READY  →  usable now; run 'bash setup.sh' to add the DockQ venv."
else
  echo "NOT READY  →  run:  bash setup.sh    (builds the venvs; skips the weights if a cache already has them)"
fi
