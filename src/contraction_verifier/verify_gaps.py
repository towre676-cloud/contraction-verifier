#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# verify_gaps.py - Advanced contraction verification script

import argparse
import json
import math
import sys
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np
from scipy.linalg import qr
from scipy.sparse.linalg import eigs


def load_json(path: str) -> Dict[str, Any]:
    try:
        with open(path, 'r', encoding='utf-8') as f:
            return json.load(f)
    except FileNotFoundError:
        print(f'Error: file not found: {path}')
        sys.exit(1)
    except json.JSONDecodeError as e:
        print(f'Error: invalid JSON in {path}: {e}')
        sys.exit(1)


def outward_round(val: float, direction: str) -> float:
    import numpy as _np
    if direction == 'up':
        return float(_np.nextafter(val, _np.inf))
    else:
        return float(_np.nextafter(val, -_np.inf))


@dataclass
class Interval:
    lo: float
    hi: float

    def __post_init__(self):
        if self.lo > self.hi:
            self.lo, self.hi = self.hi, self.lo

    def __add__(self, other: Union['Interval', float]) -> 'Interval':
        if isinstance(other, Interval):
            return Interval(self.lo + other.lo, self.hi + other.hi)
        return Interval(self.lo + other, self.hi + other)

    def __mul__(self, other: Union['Interval', float]) -> 'Interval':
        if isinstance(other, Interval):
            a, b, c, d = self.lo, self.hi, other.lo, other.hi
            prods = [a * c, a * d, b * c, b * d]
            return Interval(min(prods), max(prods))
        else:
            if other >= 0:
                return Interval(self.lo * other, self.hi * other)
            else:
                return Interval(self.hi * other, self.lo * other)

    def exp(self) -> 'Interval':
        return Interval(math.exp(self.lo), math.exp(self.hi))

    def __repr__(self) -> str:
        return f'[{self.lo:.6g}, {self.hi:.6g}]'


def collatz_wielandt_bounds(L: np.ndarray, v: Optional[np.ndarray] = None, iters: int = 100, tol: float = 1e-12) -> Tuple[float, float]:
    n = L.shape[0]
    if v is None:
        v = np.ones(n, dtype=float)
    v = v / np.linalg.norm(v)
    for _ in range(iters):
        w = L @ v
        nw = np.linalg.norm(w)
        if nw < tol:
            break
        v = w / nw
    w = L @ v
    ratios = []
    for i in range(n):
        if abs(v[i]) > tol:
            ratios.append(w[i] / v[i])
    if not ratios:
        bound = float(np.linalg.norm(L, ord=np.inf))
        return bound, bound
    r_min = float(np.min(ratios))
    r_max = float(np.max(ratios))
    return r_min, r_max


def compute_gamma_V(Gamma: np.ndarray, R: Optional[np.ndarray] = None, V_basis: Optional[np.ndarray] = None, use_svd: bool = False) -> float:
    if R is None and V_basis is None:
        raise ValueError('compute_gamma_V: provide R or V_basis')
    if R is None:
        Q, _ = qr(V_basis, mode='economic')
        R = Q @ Q.T
    RG = R @ Gamma
    if use_svd:
        return float(np.linalg.norm(RG, ord=2))
    try:
        vals, _ = eigs(RG.T @ RG, k=1, which='LM', maxiter=200, tol=1e-8)
        lam = max(0.0, float(np.real(vals[0])))
        return math.sqrt(lam)
    except Exception:
        return float(np.linalg.norm(RG, ord=2))


def epsilon_from_witness(alpha: float, beta: float, sigma: float, delta_tau_min: float) -> float:
    return alpha ** beta * math.exp(-sigma * delta_tau_min)


def apply_verified_margins(lhs: float, rhs: float, args) -> Tuple[float, float]:
    if getattr(args, 'verified_strict', False):
        rel_k = getattr(args, 'rel_k', 1e-15)
        abs_k = getattr(args, 'abs_k', 0.0)
        lhs = lhs * (1 + rel_k) + abs_k
        rhs = rhs * (1 - rel_k) - abs_k
    if getattr(args, 'verified', False) or getattr(args, 'verified_strict', False):
        lhs = outward_round(lhs, 'up')
        rhs = outward_round(rhs, 'down')
    return lhs, rhs


def verify_core_banach(data: Dict[str, Any], args) -> bool:
    alpha = float(data.get('alpha', 0.95))
    beta = float(data.get('beta', 1.0))
    sigma = float(data.get('sigma', 0.5))
    delta_tau_min = float(data.get('delta_tau_min', 0.1))
    c_P = float(data.get('c_P', 1.0))

    eps = float(args.eps) if args.eps is not None else epsilon_from_witness(alpha, beta, sigma, delta_tau_min)

    core = data.get('core_params', {})
    if args.interval:
        g_lo = float(core.get('gamma_lo', core.get('gamma', 0.85)))
        g_hi = float(core.get('gamma_hi', core.get('gamma', 0.85)))
        t_lo = float(core.get('t_lo', core.get('t', 0.02)))
        t_hi = float(core.get('t_hi', core.get('t', 0.02)))
        gamma_I = Interval(g_lo, g_hi)
        t_I = Interval(t_lo, t_hi)
        cP_I = Interval(c_P, c_P)
        one = Interval(1.0, 1.0)
        lhs_I = gamma_I * (one + t_I * cP_I)
        lhs = lhs_I.hi
    else:
        gamma = float(core.get('gamma', 0.85))
        t = float(core.get('t', 0.02))
        if args.use_gammaV and ('Gamma_matrix' in data):
            Gamma = np.array(data['Gamma_matrix'], dtype=float)
            if 'R_matrix' in data:
                R = np.array(data['R_matrix'], dtype=float)
                gamma = compute_gamma_V(Gamma, R=R, V_basis=None, use_svd=args.svd_exact)
            elif 'V_basis' in data:
                V_basis = np.array(data['V_basis'], dtype=float)
                gamma = compute_gamma_V(Gamma, R=None, V_basis=V_basis, use_svd=args.svd_exact)
            else:
                raise ValueError('use_gammaV requires R_matrix or V_basis in JSON')
            print(f'gamma_V = ||R*Gamma||_2 = {gamma:.9g}')
        lhs = gamma * (1.0 + t * c_P)

    rhs = 1.0 - eps
    lhs, rhs = apply_verified_margins(lhs, rhs, args)

    passed = lhs <= rhs
    margin = rhs - lhs
    print(f'eps={eps:.9g}  c_P={c_P:.9g}  type=core_banach')
    print(f'Core Banach: {"OK" if passed else "FAIL"} (margin = {margin:+.3e})')

    if getattr(args, 'sensitivity', False) and not args.interval:
        gamma_print = float(core.get('gamma', 0.85))
        t_print = float(core.get('t', 0.02))
        d_gamma = -(1.0 + t_print * c_P)
        d_t = -gamma_print * c_P
        d_sigma = eps * delta_tau_min
        print(f'Sensitivity: dG/dgamma={d_gamma:.6g}, dG/dt={d_t:.6g}, dG/dsigma={d_sigma:.6g}')
    return passed


def verify_ledger(data: Dict[str, Any], args) -> bool:
    alpha = float(data.get('alpha', 0.98))
    beta = float(data.get('beta', 1.0))
    sigma = float(data.get('sigma', 1.2))
    delta_tau_min = float(data.get('delta_tau_min', 0.05))
    c_P = float(data.get('c_P', 1.0))

    eps = float(args.eps) if args.eps is not None else epsilon_from_witness(alpha, beta, sigma, delta_tau_min)

    steps = data.get('steps', [])
    if not isinstance(steps, list) or not steps:
        print('Error: ledger: "steps" must be a non-empty list')
        return False

    all_ok = True
    worst_margin = float('+inf')
    worst_i = -1

    print(f'eps={eps:.9g}  c_P={c_P:.9g}  type=ledger')
    print(f'Verifying {len(steps)} steps...')

    for i, st in enumerate(steps):
        if args.interval:
            gI = Interval(float(st.get('gamma_lo', st.get('gamma', 0.9))),
                          float(st.get('gamma_hi', st.get('gamma', 0.9))))
            tI = Interval(float(st.get('t_lo', st.get('t', 0.01))),
                          float(st.get('t_hi', st.get('t', 0.01))))
            dI = Interval(float(st.get('delta_tau_lo', st.get('delta_tau', delta_tau_min))),
                          float(st.get('delta_tau_hi', st.get('delta_tau', delta_tau_min))))
            one = Interval(1.0, 1.0)
            cPI = Interval(c_P, c_P)
            facI = gI * (one + tI * cPI)
            expI = Interval(-sigma, -sigma) * dI
            lhs = (facI * expI.exp()).hi
        else:
            g = float(st.get('gamma', 0.9))
            t = float(st.get('t', 0.01))
            d = float(st.get('delta_tau', delta_tau_min))
            lhs = g * (1.0 + t * c_P) * math.exp(-sigma * d)

        rhs = 1.0 - eps
        lhs_adj, rhs_adj = apply_verified_margins(lhs, rhs, args)

        ok = lhs_adj <= rhs_adj
        margin = rhs_adj - lhs_adj
        all_ok &= ok

        if margin < worst_margin:
            worst_margin = margin
            worst_i = i

        if i < 10 or not ok:
            print(f'  step {i:03d}: {lhs_adj:.9g} <= {rhs_adj:.9g}  {"OK" if ok else "FAIL"}')

    print(f'Ledger: {"OK" if all_ok else "FAIL"}; worst margin {worst_margin:+.3e} at step {worst_i}')
    if worst_margin < 1e-6:
        print('WARNING: near-critical margin (<= 1e-6)')
    return all_ok


def verify_coupled_stream(data: Dict[str, Any], args) -> bool:
    if 'kappa' in data and 'eta' in data:
        kappa = np.array(data['kappa'], dtype=float).ravel()
        eta = np.array(data['eta'], dtype=float)
        n = kappa.size
        if eta.shape != (n, n):
            raise ValueError('eta must be n x n to match len(kappa)')
        L = np.array(eta, dtype=float)
        for i in range(n):
            L[i, i] = 1.0 - float(kappa[i])
    else:
        sp = data.get('stream_params', {})
        k1 = float(sp.get('kappa1', 0.4))
        k2 = float(sp.get('kappa2', 0.5))
        e1 = float(sp.get('eta1', 0.1))
        e2 = float(sp.get('eta2', 0.08))
        L = np.array([[1.0 - k1, e1], [e2, 1.0 - k2]], dtype=float)

    print('type=coupled_stream')
    print('Lipschitz matrix L:')
    for row in L:
        print('  ' + ' '.join(f'{x: .6g}' for x in row))

    if getattr(args, 'svd_exact', False) and L.shape[0] <= 5:
        vals = np.linalg.eigvals(L)
        rho = float(np.max(np.abs(vals)))
        r_min = r_max = rho
        print(f'rho(L) = {rho:.9g} (exact)')
    else:
        r_min, r_max = collatz_wielandt_bounds(L, iters=getattr(args, 'power_iters', 150), tol=getattr(args, 'tol', 1e-12))
        print(f'rho(L) in [{r_min:.9g}, {r_max:.9g}] (enclosure)')

    tol = float(getattr(args, 'tol', 1e-12))
    if r_max < 1.0 - tol:
        print('Coupled streams: OK')
        print(f'Effective epsilon = {1.0 - r_max:.9g}')
        return True
    elif r_min >= 1.0 - tol:
        print('Coupled streams: FAIL')
        return False
    else:
        print('Coupled streams: Inconclusive (near-critical)')
        return False


def create_example_manifest(example_type: str, output_file: str) -> None:
    if example_type == 'core':
        obj = {
            'verification_type': 'core_banach',
            'alpha': 0.95, 'beta': 1.0, 'sigma': 0.5, 'delta_tau_min': 0.1, 'c_P': 1.0,
            'core_params': {'gamma': 0.85, 't': 0.02}
        }
    elif example_type == 'ledger':
        obj = {
            'verification_type': 'ledger',
            'alpha': 0.98, 'beta': 1.0, 'sigma': 1.2, 'delta_tau_min': 0.05, 'c_P': 1.0,
            'steps': [
                {'gamma': 0.92, 't': 0.01, 'delta_tau': 0.08},
                {'gamma': 0.89, 't': 0.015, 'delta_tau': 0.06},
                {'gamma': 0.91, 't': 0.012, 'delta_tau': 0.07}
            ]
        }
    elif example_type == 'coupled':
        obj = {
            'verification_type': 'coupled_stream',
            'stream_params': {'kappa1': 0.4, 'kappa2': 0.5, 'eta1': 0.1, 'eta2': 0.08}
        }
    elif example_type == 'gammaV':
        obj = {
            'verification_type': 'core_banach',
            'alpha': 0.95, 'beta': 1.0, 'sigma': 0.5, 'delta_tau_min': 0.1, 'c_P': 1.0,
            'core_params': {'gamma': 0.95, 't': 0.02},
            'Gamma_matrix': [[0.9, 0.1], [0.0, 0.8]],
            'R_matrix': [[1.0, 0.0], [0.0, 0.6]]
        }
    else:
        raise ValueError(f'Unknown example type: {example_type}')
    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(obj, f, indent=2)
    print(f'Wrote example to {output_file}')


def main(argv: Optional[List[str]] = None) -> int:
    p = argparse.ArgumentParser(description='Verify contraction mapping conditions')
    p.add_argument('json_file', nargs='?', help='Manifest JSON (omit only with --create-example)')
    p.add_argument('--type', choices=['core_banach', 'ledger', 'coupled_stream'], help='Verification type')
    p.add_argument('--eps', type=float, help='Override epsilon')
    p.add_argument('--tol', type=float, default=1e-12, help='Numerical tolerance (default 1e-12)')
    p.add_argument('--verified', action='store_true', help='Outward rounding on final comparison')
    p.add_argument('--verified-strict', action='store_true', help='Inflate lhs/deflate rhs by margins before outward rounding')
    p.add_argument('--rel-k', type=float, default=1e-15, help='Relative margin for strict mode')
    p.add_argument('--abs-k', type=float, default=0.0, help='Absolute margin for strict mode')
    p.add_argument('--interval', action='store_true', help='Propagate measurement intervals')
    p.add_argument('--use-gammaV', action='store_true', help='Use gamma_V = ||R * Gamma||_2 if Gamma_matrix present')
    p.add_argument('--svd-exact', action='store_true', help='Use exact SVD/eig instead of power iteration for gamma_V / rho')
    p.add_argument('--power-iters', type=int, default=150, help='Power iteration steps for spectral enclosure')
    p.add_argument('--sensitivity', action='store_true', help='Print sensitivities (core_banach)')
    p.add_argument('--create-example', choices=['core', 'ledger', 'coupled', 'gammaV'], help='Write a sample manifest')
    p.add_argument('--out', default='example.json', help='Output for --create-example')

    args = p.parse_args(argv)

    if args.create_example:
        create_example_manifest(args.create_example, args.out)
        return 0

    if not args.json_file or not args.type:
        print('Error: provide json_file and --type (or use --create-example).')
        return 2

    data = load_json(args.json_file)
    if args.type == 'core_banach':
        ok = verify_core_banach(data, args)
    elif args.type == 'ledger':
        ok = verify_ledger(data, args)
    elif args.type == 'coupled_stream':
        ok = verify_coupled_stream(data, args)
    else:
        print(f'Unknown --type: {args.type}')
        return 2
    return 0 if ok else 1
