#!/bin/bash
# Entry point — activates venv and runs prediction
# Usage: ./fold.sh --heavy SEQ --light SEQ --antigen SEQ
set -e
cd "$(dirname "$0")"
source ../env/bin/activate 2>/dev/null || source env/bin/activate 2>/dev/null || true
python predict_local.py "$@"
