#!/usr/bin/env bash
set -euo pipefail

if [ $# -lt 2 ]; then
  echo "Usage: $0 <patient_id> <results_dir>"
  exit 1
fi

PATIENT="$1"
RESULTS_DIR="$2"

QMD="${PATIENT}.qmd"
HTML="${PATIENT}.html"

cp template.qmd "$QMD"

quarto render "$QMD" \
  -P patient:"$PATIENT" \
  -P results_dir:"$RESULTS_DIR" \
  --output "$HTML"

echo "Report written to $HTML"
