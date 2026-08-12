#!/bin/bash
# predict.sh — fold + score in one command, without remembering to `source env.sh`.
#
#   bash predict.sh --target targets/my_ab --gpu 0
#   bash predict.sh --heavy EVQ... --light DIV... --antigen DQI... --tag myAb --gpu 1
#   bash predict.sh --target targets/my_ab --gpu 0 --num_seeds 2 --num_loops 5   # quick smoke test
#
# Every flag is passed straight through to fold.py (`python fold.py --help`).
# Defaults are the validated recipe: ESMFold2-Fast, single-sequence, 20 loops,
# 100 sampling steps, lm_dropout 0.3, fp32, seeds 0-24, all seeds saved and scored.
set -e
REPO="$(cd "$(dirname "$0")" && pwd)"
# shellcheck disable=SC1091
source "$REPO/env.sh"
exec python "$REPO/fold.py" "$@"
