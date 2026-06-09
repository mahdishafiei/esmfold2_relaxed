#!/bin/bash
# Entry point — activates venv and runs prediction
# Usage: ./fold.sh --heavy SEQ --light SEQ --antigen SEQ
#        ./fold.sh --test   (smoke test to verify setup)
set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
VENV="$SCRIPT_DIR/env/bin/activate"

if [ ! -f "$VENV" ]; then
  echo "ERROR: virtualenv not found at $VENV"
  echo "Run setup first: bash $SCRIPT_DIR/setup.sh"
  exit 1
fi

source "$VENV"
cd "$SCRIPT_DIR"

# Smoke test mode
if [ "$1" = "--test" ]; then
  echo "Running smoke test (1 seed, short sequences)..."
  python predict_local.py \
    --heavy   "EVQLVESGGGLVKPGGSLRLSCAASGFTFS" \
    --antigen "DQICIGYHANNSTEQVDTIMEK" \
    --seeds 1 --loops 5 --diff-steps 10 \
    --out /tmp/esmfold2_test.cif
  echo "Smoke test passed. Output: /tmp/esmfold2_test.cif"
  exit 0
fi

python predict_local.py "$@"
