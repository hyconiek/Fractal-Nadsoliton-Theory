#!/usr/bin/env python3
"""Bounded repository intake; source artifacts are read, never rewritten.

Static bounds use exact-decimal mpmath.iv inputs and rational ANOVA.
Dynamics are independently assembled and checked by Poisson uniformization;
their floating arithmetic remains numerical, not interval-certified.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import itertools
import json
import math
import re
from fractions import Fraction
from pathlib import Path

import mpmath as mp
import numpy as np
from scipy.sparse import coo_matrix, eye
from scipy.special import logsumexp, softmax
from scipy.stats import poisson


ROOT = Path(__file__).resolve().parents[1]
PACKAGE = ROOT / "fin_generative_campaign_20261014"
Q = 12
G = 5.145228719489142


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def binary_fraction(endpoint):
    sign, mantissa, exponent, _ = endpoint
    n = -mantissa if sign else mantissa
    return Fraction(n << exponent, 1) if exponent >= 0 else Fraction(n, 1 << -exponent)


def public_decimal(value, digits=15, upper=True):
    scale = 10 ** digits
    n, d = value.numerator * scale, value.denominator
    rounded = -((-n) // d) if upper else n // d
    sign = "-" if rounded < 0 else ""
    a, b = divmod(abs(rounded), scale)
    return f"{sign}{a}.{b:0{digits}d}"


def fraction_iv(value):
    return mp.iv.mpf(value.numerator) / mp.iv.mpf(value.denominator)


def nominal_anova(F):
    f = np.array([Fraction.from_float(float(x)) for x in F.flat], dtype=object).reshape(F.shape)
    mean = sum(f.flat, Fraction()) / (Q ** 3)
    singles = []
    for axis in range(3):
        singles.append([sum((f[idx] for idx in np.ndindex(f.shape) if idx[axis] == a), Fraction()) / (Q ** 2) for a in range(Q)])
    pairs = []
    for a, b in ((0, 1), (0, 2), (1, 2)):
        table = np.empty((Q, Q), dtype=object)
        hidden = 3 - a - b
        for i, j in itertools.product(range(Q), repeat=2):
            indices = []
            for k in range(Q):
                idx = [0, 0, 0]
                idx[a], idx[b], idx[hidden] = i, j, k
                indices.append(tuple(idx))
            table[i, j] = sum((f[idx] for idx in indices), Fraction()) / Q
        pairs.append(table)
    residual = np.empty(f.shape, dtype=object)
    pair = np.empty(f.shape)
    for a, b, c in np.ndindex(f.shape):
        fp = pairs[0][a, b] + pairs[1][a, c] + pairs[2][b, c] - singles[0][a] - singles[1][b] - singles[2][c] + mean
        residual[a, b, c] = f[a, b, c] - fp
        pair[a, b, c] = float(fp)
    return max(residual.flat) - min(residual.flat), pair


def exact_input_static(F, A):
    mp.iv.dps = 55
    I = mp.iv
    omega, phi, eta, g = (I.mpf(s) for s in ("0.18575", "0.1625", "1.8", "5.145228719489142"))
    w = [I.mpf(0)]
    for d in range(1, Q):
        distance = min(d, Q - d)
        w.append(I.cos(omega * distance + phi) / (1 + I.power(I.mpf(distance), eta)))
    lambdas = {k: sum((w[d] * (1 - I.cos(2 * I.pi * k * d / Q)) for d in range(1, Q)), I.mpf(0)) for k in (3, 4, 5, 6)}
    row = [sum((lambdas[k] * I.cos(2 * I.pi * k * d / Q) / 6 for k in (3, 4, 5)), I.mpf(0)) + lambdas[6] * ((-1) ** d) / 12 for d in range(Q)]
    matrix_error = Fraction()
    for i, j in itertools.product(range(Q), repeat=2):
        for endpoint in row[(j - i) % Q]._mpi_:
            matrix_error = max(matrix_error, abs(binary_fraction(endpoint) - Fraction.from_float(float(A[i, j]))))
    epsilon = Fraction()
    selected = {}
    for state in itertools.product(range(Q), repeat=3):
        quadratic = sum((row[(b - a) % Q] for a in state for b in state), I.mpf(0)) * g / 8
        fields = [sum((row[(j - a) % Q] for a in state), I.mpf(0)) * g / 4 for j in range(Q)]
        val = quadratic + I.log(sum((I.exp(h) for h in fields), I.mpf(0)))
        nominal = Fraction.from_float(float(F[state]))
        lo, hi = map(binary_fraction, val._mpi_)
        epsilon = max(epsilon, abs(lo - nominal), abs(hi - nominal))
        if state[0] in (0, 1) and state[1] in (0, 2) and state[2] in (0, 4):
            selected[state] = val
    contrast = I.mpf(0)
    for bits in itertools.product((0, 1), repeat=3):
        state = (bits[0], 2 * bits[1], 4 * bits[2])
        contrast += ((-1) ** sum(bits)) * selected[state]
    contrast_lo, contrast_hi = map(binary_fraction, contrast._mpi_)
    assert contrast_lo > 0
    nominal_osc, pair = nominal_anova(F)
    osc_lo = max(Fraction(), nominal_osc - 16 * epsilon)
    osc_hi = nominal_osc + 16 * epsilon
    e = I.exp(fraction_iv(osc_hi) / 2)
    tv_hi = binary_fraction(((e - 1) / (e + 1))._mpi_[1])
    jump_rates = []
    for hidden in (0, 1):
        logits = [g * (row[j] + row[(j - hidden) % Q]) / 3 for j in range(Q)]
        weights = [I.exp(h) for h in logits]
        jump_rates.append(weights[1] / sum(weights, I.mpf(0)))
    rate_difference = jump_rates[1] - jump_rates[0]
    rate_lo, rate_hi = map(binary_fraction, rate_difference._mpi_)
    assert rate_lo > 0
    return {
        "arithmetic": "mpmath.iv exact-decimal inputs; exact rational ANOVA; public decimals rounded outward",
        "conditional_backend": "mpmath.iv correctness, not proof-assistant formalization",
        "all_labelled_states": Q ** 3,
        "matrix_discrepancy_upper": public_decimal(matrix_error, 25),
        "logweight_discrepancy_upper": public_decimal(epsilon, 25),
        "third_contrast_lower": public_decimal(contrast_lo, 15, False),
        "third_contrast_upper": public_decimal(contrast_hi, 15),
        "residual_oscillation_lower": public_decimal(osc_lo, 15, False),
        "residual_oscillation_upper": public_decimal(osc_hi, 15),
        "pair_TV_upper": public_decimal(tv_hi, 15),
        "pair_TV_public_ceiling": "0.051301995505",
        "nonlumpability_local_rate_difference_lower": public_decimal(rate_lo, 15, False),
        "nonlumpability_local_rate_difference_upper": public_decimal(rate_hi, 15),
        "proof_note": "||product_i(I-E_i)||_infinity <= 8, hence oscillation error <= 16 epsilon",
    }, pair


def build_labelled_parent(A):
    states = np.array(list(itertools.product(range(Q), repeat=3)), dtype=int)
    counts = np.eye(Q, dtype=int)[states].sum(axis=1)
    logpi = G * np.einsum("si,ij,sj->s", counts, A, counts) / 6
    pi = np.exp(logpi - logsumexp(logpi))
    rows, cols, data = [], [], []
    indices = np.arange(len(states))
    strides = (Q ** 2, Q, 1)
    for axis in range(3):
        other = counts - np.eye(Q, dtype=int)[states[:, axis]]
        rates = softmax((G / 3) * (other @ A), axis=1)
        for new in range(Q):
            mask = states[:, axis] != new
            ix = indices[mask]
            dest = ix + (new - states[mask, axis]) * strides[axis]
            rows.extend(ix); cols.extend(dest); data.extend(rates[mask, new])
    generator = coo_matrix((data, (rows, cols)), shape=(len(states), len(states))).tocsr()
    generator -= coo_matrix((np.asarray(generator.sum(axis=1)).ravel(), (indices, indices)), shape=generator.shape).tocsr()
    return states, pi, generator


def uniformized(generator, t, initial, rate):
    P = eye(generator.shape[0], format="csr") + generator / rate
    term = initial.copy()
    coefficient = math.exp(-rate * t)
    result = coefficient * term
    cutoff = 0
    while poisson.sf(cutoff, rate * t) > 1e-15:
        cutoff += 1
        term = P.T @ term
        coefficient *= rate * t / cutoff
        result += coefficient * term
    return result, float(poisson.sf(cutoff, rate * t)), cutoff


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    manifest = []
    for line in (PACKAGE / "MANIFEST.sha256").read_text().splitlines():
        expected, relative = line.split(None, 1)
        path = PACKAGE / relative
        assert path.is_file() and sha(path) == expected, relative
        manifest.append(relative)
    core = ROOT / "fin_physical_bridge_campaign_01/PHYS-002/safe_seed.py"
    spec = importlib.util.spec_from_file_location("intake_core", core)
    model = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(model)
    _, _, _, _, A = model.build_rank7()
    contract = json.loads((PACKAGE / "RG-C01/MODEL_CONTRACT.json").read_text())
    matrix_error = float(np.max(np.abs(A - np.asarray(contract["matrix_A7"]))))
    assert matrix_error < 1e-13
    assert contract["contract"]["gain"] == model.G_FROZEN
    states3 = np.array(list(itertools.product(range(Q), repeat=3)))
    counts3 = np.eye(Q, dtype=int)[states3].sum(axis=1)
    F = (G / 8) * np.einsum("si,ij,sj->s", counts3, A, counts3) + logsumexp((G / 4) * (counts3 @ A), axis=1)
    F = F.reshape((Q, Q, Q))
    certified, pair = exact_input_static(F, A)
    assert Fraction(certified["pair_TV_upper"]) <= Fraction(certified["pair_TV_public_ceiling"])
    p = np.exp(F - logsumexp(F))
    approx = np.exp(pair - logsumexp(pair))
    static_tv = float(np.abs(p - approx).sum() / 2)
    assert abs(static_tv - 0.0233588550889189) < 1e-12
    event_exact = float(sum(p[i, j, i] for i, j in itertools.product(range(Q), repeat=2)))
    event_approx = float(sum(approx[i, j, i] for i, j in itertools.product(range(Q), repeat=2)))
    # Invariant direct parent control retains M=4 throughout elimination.
    states4 = np.array(list(itertools.product(range(Q), repeat=4)))
    counts4 = np.eye(Q, dtype=int)[states4].sum(axis=1)
    log4 = (G / 8) * np.einsum("si,ij,sj->s", counts4, A, counts4)
    parent4 = np.exp(log4 - logsumexp(log4)).reshape((Q,) * 4)
    marginal3 = parent4.sum(axis=3)
    marginal_error = float(np.max(np.abs(marginal3 - p)))
    assert marginal_error < 1e-14
    states, pi, L = build_labelled_parent(A)
    y = states[:, 0] * Q + states[:, 1]
    py = np.bincount(y, weights=pi, minlength=Q ** 2)
    B = np.zeros((Q ** 2, Q ** 2))
    coo = L.tocoo()
    for i, j, value in zip(coo.row, coo.col, coo.data):
        if y[i] != y[j]:
            B[y[i], y[j]] += (pi[i] / py[y[i]]) * value
    np.fill_diagonal(B, -B.sum(axis=1))
    rowsum = float(np.max(np.abs(np.asarray(L.sum(axis=1)))))
    flux = L.multiply(pi[:, None])
    difference = (flux - flux.T).tocoo()
    db = float(np.max(np.abs(difference.data)))
    assert rowsum < 1e-12 and db < 1e-14
    initial = np.where(y == 0, pi / py[0], 0)
    coarse_initial = np.zeros(Q ** 2); coarse_initial[0] = 1
    records = []
    supplied = json.loads((PACKAGE / "RG-C03/PROCESS_MAP.json").read_text())["FIN_A7"]
    for t, source in zip((0.15, 0.5, 1.0, 1.5), supplied["finite_window"]):
        fine, tail, cutoff = uniformized(L, t, initial, 3.0)
        observed = np.bincount(y, weights=fine, minlength=Q ** 2)
        coarse, coarse_tail, _ = uniformized(coo_matrix(B).tocsr(), t, coarse_initial, 2.0)
        value = float(np.abs(observed - coarse).sum() / 2)
        assert abs(value - source["TV_exact_vs_equilibrium_generator"]) < 1e-12
        records.append({"t": t, "TV_to_averaged_generator": value, "Poisson_tail_numerical": tail, "coarse_Poisson_tail_numerical": coarse_tail, "series_cutoff": cutoff})
    f = (y == 0).astype(float)
    Lf = L @ f
    projected = np.bincount(y, weights=pi * Lf, minlength=Q ** 2) / py
    v = Lf - projected[y]
    K0 = float(pi @ (v * v))
    Kprime0 = float((pi * v) @ (L @ v))
    memory = supplied["memory_quadratic_channel"]
    assert abs(K0 - memory["M0_coupled"]) < 1e-13
    assert abs(Kprime0 - memory["M1_derivative"]) < 1e-13
    jump0 = float(L[0, Q ** 2])
    jump1 = float(L[1, Q ** 2 + 1])
    defect = abs(jump0 - jump1)
    assert abs(defect - supplied["strong_lumpability_witness"]["maximum_single_visible_refresh_rate_difference"]) < 1e-13
    a, _, _ = uniformized(L, 0.5, np.eye(1, len(pi), 0).ravel(), 3.0)
    b, _, _ = uniformized(L, 0.5, np.eye(1, len(pi), 1).ravel(), 3.0)
    prep_tv = float(np.abs(np.bincount(y, weights=a, minlength=Q ** 2) - np.bincount(y, weights=b, minlength=Q ** 2)).sum() / 2)
    assert abs(prep_tv - supplied["hidden_preparation_dependence"][1]["TV_hidden_0_vs_hidden_1"]) < 1e-12
    results = {
        "date": "2026-10-14",
        "manifest": {"entries": len(manifest), "mismatch": 0},
        "source_hashes": {str(core.relative_to(ROOT)): sha(core), "AGENTS.md": sha(ROOT / "AGENTS.md"), "FIN_GENERATIVE_NEXT_CAMPAIGN_20261014.md": sha(ROOT / "FIN_GENERATIVE_NEXT_CAMPAIGN_20261014.md"), "campaign/run_campaign.py": sha(PACKAGE / "run_campaign.py"), "campaign/interval_certify.py": sha(PACKAGE / "interval_certify.py"), "campaign/SOURCE_PLAN_20261014.md": sha(PACKAGE / "SOURCE_PLAN_20261014.md")},
        "campaign_matrix_vs_versioned_core": matrix_error,
        "plan_byte_identity": sha(ROOT / "FIN_GENERATIVE_NEXT_CAMPAIGN_20261014.md") == sha(PACKAGE / "SOURCE_PLAN_20261014.md"),
        "plan_format_normalized_identity": re.sub(r"\s+", "", (ROOT / "FIN_GENERATIVE_NEXT_CAMPAIGN_20261014.md").read_text().replace("\\", "").replace("*", "")) == re.sub(r"\s+", "", (PACKAGE / "SOURCE_PLAN_20261014.md").read_text().replace("\\", "").replace("*", "")),
        "static_certificate": certified,
        "static_numerical": {"actual_pair_TV": static_tv, "event_exact": event_exact, "event_pair": event_approx, "event_error": abs(event_exact - event_approx), "parent_to_formula_max_error": marginal_error},
        "dynamic_numerical": {"row_sum_max": rowsum, "detailed_balance_max": db, "lumpability_defect": defect, "K_at_zero": K0, "K_derivative_at_zero": Kprime0, "integrated_memory_moments_computed": False, "uniformization_comparison": records, "hidden_preparation_TV_at_t0p5": prep_tv},
        "scope": ["No dynamic interval propagation certificate", "K(0) and K'(0) are not integrated M0/M1", "RG-C05 is a static consistency check, not the requested dynamic finite-window validation", "No physical source/calibration result"]
    }
    if args.output:
        args.output.write_text(json.dumps(results, indent=2) + "\n")
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
