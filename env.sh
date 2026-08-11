# env.sh — activate the ESMFold2 environment. `source` this file; do not execute it.
#
#   source /path/to/esmfold2_relaxed/env.sh
#
# It activates the fold venv (torch + esm), points HF_HOME at the model weights, and
# exports $EF2_PY (fold interpreter), $EF2_DOCKQ (DockQ binary), $EF2_VENV, $EF2_REPO.
# Build everything once per machine with:  bash setup.sh
#
# Override any of these by exporting them before sourcing:
#   EF2_VENV_DIR    fold venv          (default $HOME/.ef2_venv)
#   EF2_DOCKQ_VENV  DockQ venv         (default $HOME/.ef2_dockq_venv — separate on purpose:
#                                       DockQ pins numpy<2, which conflicts with esm)
#   EF2_HF_HOME     weights cache      (default $HOME/.ef2_hf_cache, else <repo>/hf_cache)
#
# On multi-node setups (e.g. JupyterHub named servers) keep the weights on the shared
# filesystem — symlink or point EF2_HF_HOME at them — and build one venv per node: a venv
# on network storage is pathologically slow to build and to import.

EF2_REPO="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"; export EF2_REPO

# --- fold venv: first one that exists wins ---
export EF2_VENV=""
for _v in "$EF2_VENV_DIR" "$HOME/.ef2_venv" "$EF2_REPO/.venv"; do
  if [ -n "$_v" ] && [ -x "$_v/bin/python" ]; then export EF2_VENV="$_v"; break; fi
done
if [ -z "$EF2_VENV" ]; then
  echo "ERROR: no fold venv on this machine. Build it once:  bash $EF2_REPO/setup.sh" >&2
  return 1 2>/dev/null || exit 1
fi

# --- weights (HF_HOME): prefer a node-local cache, fall back to one inside the repo ---
for _h in "$EF2_HF_HOME" "$HOME/.ef2_hf_cache" "$EF2_REPO/hf_cache"; do
  if [ -n "$_h" ] && [ -d "$_h/hub" ]; then export HF_HOME="$_h"; break; fi
done

# --- DockQ: separate venv, optional (without it the DockQ columns are simply skipped) ---
for _d in "$EF2_DOCKQ_VENV/bin/DockQ" "$HOME/.ef2_dockq_venv/bin/DockQ" "$(command -v DockQ)"; do
  if [ -n "$_d" ] && [ -x "$_d" ]; then export EF2_DOCKQ="$_d"; break; fi
done

export EF2_PY="$EF2_VENV/bin/python"
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
export TRITON_CACHE_DIR="${TRITON_CACHE_DIR:-/tmp/ef2_triton_$(id -u)}"
mkdir -p "$TRITON_CACHE_DIR" 2>/dev/null || true

# shellcheck disable=SC1091
source "$EF2_VENV/bin/activate"
echo "[env.sh] venv=$EF2_VENV  HF_HOME=${HF_HOME:-<unset — weights not found>}  DockQ=${EF2_DOCKQ:-<none>}"
