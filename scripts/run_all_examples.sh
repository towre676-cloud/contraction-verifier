#!/bin/bash
set -e

echo "=== Core Banach ==="
python src/verify_gaps.py examples/core_banach/example_core.json --type core_banach

echo "=== Ledger ==="
python src/verify_gaps.py examples/ledger/example_ledger.json --type ledger

echo "=== Coupled Streams ==="
python src/verify_gaps.py examples/coupled_streams/example_coupled.json --type coupled_stream

echo "=== Gamma_V (SVD exact) ==="
python src/verify_gaps.py examples/gamma_v/example_gammaV.json --type core_banach --use-gammaV --svd-exact
