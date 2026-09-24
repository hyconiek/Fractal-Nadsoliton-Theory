#!/usr/bin/env python3
"""Independent core checks for FIN_MISSING_LAWS_POST_FLAT_HANDOFF_20260924.

This intentionally does not import or execute the handoff replay scripts.  It
reconstructs the algebra used by the strongest POST-01..08 claims from the
declared q=12 shell weights and standard finite-dimensional linear algebra.

The script does not attempt to certify POST-06 coexistence scouts, and it does
not claim to reproduce the exact POST-07 gap/KL numbers because the handoff
does not provide the equilibrium/proposal data needed for that replay.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np


Q = 12
STRICT_W = np.array(
    [
        0.4699856726450201,
        0.1920435516901028,
        0.09142861427792495,
        0.0470291687456504,
        0.02413122336363006,
        0.011070817321442113,
    ],
    dtype=float,
)
SPARSE_W = np.array([0.53385079, 0.19757896, 0.05169472, 0.00175987, 0.0, 0.0])


def shell_laplacian(weights: np.ndarray, q: int = Q) -> np.ndarray:
    a = np.zeros((q, q), dtype=float)
    for i in range(q):
        for j in range(q):
            if i == j:
                continue
            d = min((j - i) % q, (i - j) % q)
            a[i, j] = -weights[d - 1]
        a[i, i] = -a[i].sum()
    return a


def complete_graph_incidence_and_covariance(weights: np.ndarray, q: int = Q):
    edges = [(i, j) for i in range(q) for j in range(i + 1, q)]
    b = np.zeros((q, len(edges)), dtype=float)
    dvals = np.empty(len(edges), dtype=float)
    for e, (i, j) in enumerate(edges):
        b[i, e] = 1.0
        b[j, e] = -1.0
        shell = min((j - i) % q, (i - j) % q)
        dvals[e] = weights[shell - 1]
    return edges, b, np.diag(dvals)


def moment_rows(q: int = Q) -> np.ndarray:
    theta = 2.0 * np.pi * np.arange(q) / q
    return np.vstack(
        [np.cos(theta), np.sin(theta), np.cos(2.0 * theta), np.sin(2.0 * theta)]
    )


def conditioned_edge_covariance(b: np.ndarray, d: np.ndarray, c: np.ndarray) -> np.ndarray:
    f = c @ b
    return d - d @ f.T @ np.linalg.inv(f @ d @ f.T) @ f @ d


def top_rank_part(a: np.ndarray, rank: int) -> np.ndarray:
    vals, vecs = np.linalg.eigh(a)
    idx = np.argsort(vals)[-rank:]
    return (vecs[:, idx] * vals[idx]) @ vecs[:, idx].T


def permutation_matrix(kind: str, shift: int, q: int = Q) -> np.ndarray:
    p = np.zeros((q, q), dtype=float)
    for j in range(q):
        target = (j + shift) % q if kind == "rotation" else (shift - j) % q
        p[target, j] = 1.0
    return p


def visible_projector(q: int = Q) -> np.ndarray:
    j = np.arange(q)
    cols = []
    for k in (3, 4, 5):
        cols.extend(
            [np.cos(2.0 * np.pi * k * j / q), np.sin(2.0 * np.pi * k * j / q)]
        )
    cols.append((-1.0) ** j)
    u, _ = np.linalg.qr(np.column_stack(cols))
    return u @ u.T


def d12_hom_cycle_to_visible_dimension(q: int = Q) -> float:
    """Character inner product for Hom_D12(cycle(K_q), visible k=3..6)."""
    pvis = visible_projector(q)
    accum = 0.0
    for kind in ("rotation", "reflection"):
        for shift in range(q):
            p = permutation_matrix(kind, shift, q)
            chi_v = float(np.trace(p))
            chi_v_g2 = float(np.trace(p @ p))
            # Complete oriented edge space is Lambda^2 of the vertex permutation rep.
            chi_edge = 0.5 * (chi_v * chi_v - chi_v_g2)
            # Incidence image is the mean-zero vertex representation.
            chi_cycle = chi_edge - (chi_v - 1.0)
            chi_visible = float(np.trace(pvis @ p))
            accum += chi_cycle * chi_visible
    return accum / (2.0 * q)


def ctmc_generator(pi: np.ndarray, edges, rule: str) -> np.ndarray:
    n = len(pi)
    g = np.zeros((n, n), dtype=float)
    for i, j in edges:
        for a, b in ((i, j), (j, i)):
            ratio = pi[b] / pi[a]
            if rule == "metropolis":
                rate = min(1.0, ratio)
            elif rule == "barker":
                rate = ratio / (1.0 + ratio)
            else:
                raise ValueError(rule)
            g[a, b] = rate
    g[np.diag_indices(n)] = -g.sum(axis=1)
    return g


def spectral_gap(generator: np.ndarray) -> float:
    vals = np.linalg.eigvals(-generator)
    vals = np.sort(np.real_if_close(vals).real)
    return float(vals[1])


def stationary_path_kl_rate(pi: np.ndarray, q: np.ndarray, r: np.ndarray) -> float:
    total = 0.0
    for i in range(len(pi)):
        for j in range(len(pi)):
            if i == j or q[i, j] <= 0.0:
                continue
            total += pi[i] * (
                q[i, j] * math.log(q[i, j] / r[i, j]) - q[i, j] + r[i, j]
            )
    return float(total)


def main() -> None:
    edges, b, d = complete_graph_incidence_and_covariance(STRICT_W)
    c = moment_rows()
    a = shell_laplacian(STRICT_W)
    dc = conditioned_edge_covariance(b, d, c)
    a_cond = b @ dc @ b.T
    a7 = top_rank_part(a, 7)

    post01 = {
        "q": Q,
        "edges": len(edges),
        "rank_B": int(np.linalg.matrix_rank(b, tol=1e-10)),
        "rank_CB": int(np.linalg.matrix_rank(c @ b, tol=1e-10)),
        "cycle_dimension": len(edges) - int(np.linalg.matrix_rank(b, tol=1e-10)),
        "moment_allowed_edge_dimension": int(np.linalg.matrix_rank(dc, tol=1e-10)),
        "visible_quotient_dimension": int(np.linalg.matrix_rank(a_cond, tol=1e-10)),
        "A_minus_BDBt_fro": float(np.linalg.norm(a - b @ d @ b.T)),
        "conditioned_covariance_minus_A7_fro": float(np.linalg.norm(a_cond - a7)),
    }
    assert post01["edges"] == 66
    assert post01["cycle_dimension"] == 55
    assert post01["moment_allowed_edge_dimension"] == 62
    assert post01["visible_quotient_dimension"] == 7
    assert post01["conditioned_covariance_minus_A7_fro"] < 1e-10

    # POST-02: direct SPD infimal-composition check plus the K(l)=G/l family.
    rng = np.random.default_rng(20260924)
    x = rng.normal(size=(7, 7))
    g = x.T @ x + np.eye(7)
    lengths = rng.uniform(0.2, 2.0, size=8)
    resistances = [ell * np.linalg.inv(g) for ell in lengths]
    k_series = np.linalg.inv(sum(resistances))
    k_closed = g / lengths.sum()
    post02 = {
        "segments": len(lengths),
        "relative_series_error": float(np.linalg.norm(k_series - k_closed) / np.linalg.norm(k_closed)),
        "total_length": float(lengths.sum()),
    }
    assert post02["relative_series_error"] < 1e-12

    # POST-03: independent pseudoinverse/minimal-lift check on a deterministic h in range(A7).
    h = a_cond @ rng.normal(size=Q)
    a_plus = np.linalg.pinv(a_cond, rcond=1e-12)
    dc_plus = np.linalg.pinv(dc, rcond=1e-12)
    qstar = dc @ b.T @ a_plus @ h
    post03 = {
        "reconstruction_error": float(np.linalg.norm(b @ qstar - h)),
        "moment_constraint_error": float(np.linalg.norm(c @ b @ qstar)),
        "edge_cost": float(qstar @ dc_plus @ qstar),
        "visible_cost": float(h @ a_plus @ h),
        "cost_abs_error": float(abs(qstar @ dc_plus @ qstar - h @ a_plus @ h)),
        "cycle_dimension": post01["cycle_dimension"],
    }
    assert post03["reconstruction_error"] < 1e-10
    assert post03["moment_constraint_error"] < 1e-10
    assert post03["cost_abs_error"] < 1e-10

    post04_hom = d12_hom_cycle_to_visible_dimension()
    post04 = {"D12_Hom_cycle_to_visible_dimension": post04_hom}
    assert abs(post04_hom - 16.0) < 1e-10

    strict_lam = np.fft.fft(a[0]).real
    sparse_a = shell_laplacian(SPARSE_W)
    sparse_lam = np.fft.fft(sparse_a[0]).real
    lam_error = float(np.linalg.norm(sparse_lam[3:7] - strict_lam[3:7]))
    multiplicity = np.array([2, 2, 2, 2, 2, 1], dtype=float)
    p_strict = multiplicity * STRICT_W / float(multiplicity @ STRICT_W)
    p_sparse = multiplicity * SPARSE_W / float(multiplicity @ SPARSE_W)
    tv = float(0.5 * np.abs(p_strict - p_sparse).sum())
    soft = {
        str(tau): [float(lam / (1.0 + tau * lam)) for lam in strict_lam[1:3]]
        for tau in (1, 10, 100)
    }
    post05 = {
        "strict_lambda12": strict_lam[1:3].tolist(),
        "strict_lambda3456": strict_lam[3:7].tolist(),
        "sparse_lambda12": sparse_lam[1:3].tolist(),
        "lambda3456_l2_error": lam_error,
        "embedded_jump_shell_TV": tv,
        "strict_soft_values": soft,
    }
    assert abs(lam_error - 1.4705563783716023e-08) < 1e-12
    assert abs(tv - 0.13441699427176207) < 1e-12

    # POST-07: independent witness of the analytic nonuniqueness theorem.
    # These values are intentionally not the package's unreplayable q=12 numbers.
    pi = np.array([1.0, 2.0, 4.0, 8.0])
    pi /= pi.sum()
    ring_edges = [(0, 1), (1, 2), (2, 3), (3, 0)]
    metro = ctmc_generator(pi, ring_edges, "metropolis")
    barker = ctmc_generator(pi, ring_edges, "barker")
    detailed_balance_m = max(
        abs(pi[i] * metro[i, j] - pi[j] * metro[j, i]) for i, j in ring_edges
    )
    detailed_balance_b = max(
        abs(pi[i] * barker[i, j] - pi[j] * barker[j, i]) for i, j in ring_edges
    )
    conductance_ratios = [
        float((pi[i] * barker[i, j]) / (pi[i] * metro[i, j])) for i, j in ring_edges
    ]
    post07 = {
        "witness_pi": pi.tolist(),
        "detailed_balance_error_metropolis": float(detailed_balance_m),
        "detailed_balance_error_barker": float(detailed_balance_b),
        "conductance_ratio_range": [min(conductance_ratios), max(conductance_ratios)],
        "metropolis_gap": spectral_gap(metro),
        "barker_gap": spectral_gap(barker),
        "KL_M_to_B": stationary_path_kl_rate(pi, metro, barker),
        "KL_B_to_M": stationary_path_kl_rate(pi, barker, metro),
    }
    assert detailed_balance_m < 1e-14 and detailed_balance_b < 1e-14
    assert max(conductance_ratios) - min(conductance_ratios) > 1e-3
    assert post07["KL_M_to_B"] > 0.0 and post07["KL_B_to_M"] > 0.0

    # POST-08: min-plus elimination and Gaussian covariance addition share M/(t1+t2)
    # as endpoint quadratic precision, despite representing different history algebras.
    y = rng.normal(size=(5, 5))
    m = y.T @ y + np.eye(5)
    t1, t2 = 0.7, 1.9
    k1, k2 = m / t1, m / t2
    k_minplus = np.linalg.inv(np.linalg.inv(k1) + np.linalg.inv(k2))
    k_expected = m / (t1 + t2)
    cov1, cov2 = t1 * np.linalg.inv(m), t2 * np.linalg.inv(m)
    k_sumproduct = np.linalg.inv(cov1 + cov2)
    post08 = {
        "minplus_precision_error": float(np.linalg.norm(k_minplus - k_expected)),
        "gaussian_sumproduct_precision_error": float(np.linalg.norm(k_sumproduct - k_expected)),
    }
    assert post08["minplus_precision_error"] < 1e-12
    assert post08["gaussian_sumproduct_precision_error"] < 1e-12

    result = {
        "POST-01": post01,
        "POST-02": post02,
        "POST-03": post03,
        "POST-04": post04,
        "POST-05": post05,
        "POST-06": {
            "independent_numeric_replay": False,
            "reason": "handoff contains no replay code for coexistence/UV scouts",
        },
        "POST-07": post07,
        "POST-08": post08,
    }
    out = Path(__file__).with_name("verification_results.json")
    out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
