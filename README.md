cat > README.md <<'MD'
# contraction-verifier

Production-grade verifier for contraction mappings:

- **Core Banach:** `γ(1 + t c_P) ≤ 1 − ε`
- **Ledger (discrete):** `γ_k(1 + t_k c_P) e^{-σ Δτ_k} ≤ 1 − ε` per step
- **Coupled streams (n-stream):** spectral radius `ρ(L) < 1`
- **γ_V bound:** use `‖R Γ‖₂` (projection-aware tightening)
- **Verified numerics:** outward rounding; strict margins
- **Interval arithmetic:** propagate measurement uncertainty

## Install (editable)
```bash
python -m venv .venv
source .venv/Scripts/activate   # Windows Git Bash
pip install -e .[dev]
