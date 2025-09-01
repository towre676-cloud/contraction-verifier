#!/usr/bin/env bash
set -euo pipefail

# Resolve repo root from script location (works no matter where you call it)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${REPO_ROOT}"

# Prefer project venv Python; else fallback to current python
if [ -x ".venv/Scripts/python.exe" ]; then
  PY=".venv/Scripts/python.exe"
elif [ -x "$HOME/contraction-mapping-verification/.venv/Scripts/python.exe" ]; then
  PY="$HOME/contraction-mapping-verification/.venv/Scripts/python.exe"
else
  PY="$(python - <<'PY'
import sys; print(sys.executable)
PY
)"
fi

echo "Using Python: $PY"
echo "Repo root: $REPO_ROOT"

echo "=== Core Banach ==="
"$PY" -m contraction_verifier.cli examples/core_banach/example_core.json --type core_banach

echo "=== Ledger ==="
"$PY" -m contraction_verifier.cli examples/ledger/example_ledger.json --type ledger

echo "=== Coupled Streams ==="
"$PY" -m contraction_verifier.cli examples/coupled_streams/example_coupled.json --type coupled_stream

echo "=== Gamma_V (SVD exact) ==="
"$PY" -m contraction_verifier.cli examples/gamma_v/example_gammaV.json --type core_banach --use-gammaV --svd-exact
