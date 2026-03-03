#!/usr/bin/env python3
"""
QW-2041: Canonical vs refrozen kernel reparameterization audit.

Question:
Can pre-QW1700 canonical kernel semantics be recovered from the QW-2039
refrozen kernel via a physically admissible monotone reparameterization
of the distance coordinate?

Transform family:
    d_can = a * d_ref**p + b
with a>0, p>0 and d_can>0 on the analysis domain.
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw2041_canonical_refrozen_reparameterization_audit.json"
OUT_MD = ROOT / "RAPORT_QW2041_CANONICAL_REFROZEN_REPARAMETERIZATION_AUDIT.md"


@dataclass(frozen=True)
class Kernel:
    omega: float
    phi: float
    beta: float
    eta: float


def kernel_eta(k: Kernel, d: np.ndarray) -> np.ndarray:
    return np.cos(k.omega * d + k.phi) / (1.0 + k.beta * (d ** k.eta))


def r2_score(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    ss_res = float(np.sum((y_true - y_pred) ** 2))
    ss_tot = float(np.sum((y_true - float(np.mean(y_true))) ** 2))
    if ss_tot <= 1e-15:
        return 0.0
    return float(1.0 - ss_res / ss_tot)


def safe_corr(a: np.ndarray, b: np.ndarray) -> float:
    sa = float(np.std(a))
    sb = float(np.std(b))
    if sa <= 1e-15 or sb <= 1e-15:
        return 0.0
    return float(np.corrcoef(a, b)[0, 1])


def affine_r2(y_true: np.ndarray, x: np.ndarray) -> Dict[str, float]:
    m = np.column_stack([x, np.ones(len(x), dtype=float)])
    coef, *_ = np.linalg.lstsq(m, y_true, rcond=None)
    y_hat = m @ coef
    return {
        "coef_a": float(coef[0]),
        "coef_b": float(coef[1]),
        "r2": r2_score(y_true, y_hat),
        "rmse": float(np.sqrt(np.mean((y_true - y_hat) ** 2))),
    }


def zeros_of_kernel(k: Kernel, n: int) -> np.ndarray:
    out = []
    for idx in range(n):
        out.append((math.pi / 2.0 - k.phi + idx * math.pi) / max(k.omega, 1e-12))
    return np.array(out, dtype=float)


def mapped_node_errors(
    canonical: Kernel,
    refrozen: Kernel,
    a: float,
    b: float,
    p: float,
    n_nodes: int,
) -> Dict[str, float]:
    z_ref = zeros_of_kernel(refrozen, n=n_nodes)
    z_can = zeros_of_kernel(canonical, n=n_nodes)

    x = (z_can - b) / max(a, 1e-15)
    pos = x > 0.0
    if int(np.sum(pos)) < 3:
        return {
            "n_compared": 0,
            "median_rel": float("inf"),
            "q95_rel": float("inf"),
            "max_rel": float("inf"),
        }

    d_pred = x[pos] ** (1.0 / max(p, 1e-15))
    n = int(min(len(d_pred), len(z_ref)))
    rel = np.abs(d_pred[:n] - z_ref[:n]) / np.maximum(np.abs(z_ref[:n]), 1e-15)

    return {
        "n_compared": n,
        "median_rel": float(np.median(rel)),
        "q95_rel": float(np.quantile(rel, 0.95)),
        "max_rel": float(np.max(rel)),
    }


def evaluate_candidate(
    canonical: Kernel,
    refrozen: Kernel,
    d_ref: np.ndarray,
    k_ref: np.ndarray,
    a: float,
    b: float,
    p: float,
    sign: float,
) -> Dict[str, object]:
    d_can = a * (d_ref ** p) + b
    if bool(np.any(d_can <= 0.0)):
        return {
            "a": a,
            "b": b,
            "p": p,
            "sign": sign,
            "invalid": True,
        }

    k_can = kernel_eta(canonical, d_can)
    k_model = sign * k_can

    rmse = float(np.sqrt(np.mean((k_ref - k_model) ** 2)))
    corr = safe_corr(k_ref, k_model)
    r2 = r2_score(k_ref, k_model)
    aff = affine_r2(k_ref, k_model)
    nodes = mapped_node_errors(canonical, refrozen, a=a, b=b, p=p, n_nodes=8)

    flags = {
        "strict_corr_ge_0p95": bool(corr >= 0.95),
        "strict_r2_ge_0p90": bool(r2 >= 0.90),
        "affine_r2_ge_0p95": bool(float(aff["r2"]) >= 0.95),
        "node_median_rel_le_0p10": bool(float(nodes["median_rel"]) <= 0.10),
        "node_q95_rel_le_0p25": bool(float(nodes["q95_rel"]) <= 0.25),
        "physical_p_band_0p5_2p0": bool(0.5 <= p <= 2.0),
        "physical_shift_absb_le_1": bool(abs(b) <= 1.0),
    }

    pass_count = int(sum(1 for v in flags.values() if v))
    total_flags = int(len(flags))

    return {
        "a": float(a),
        "b": float(b),
        "p": float(p),
        "sign": float(sign),
        "rmse": rmse,
        "corr": corr,
        "r2": r2,
        "affine": aff,
        "node_errors": nodes,
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "invalid": False,
    }


def scan_reparameterizations(canonical: Kernel, refrozen: Kernel) -> Dict[str, object]:
    d_ref = np.linspace(1.0, 24.0, 512, dtype=float)
    k_ref = kernel_eta(refrozen, d_ref)

    a_grid = np.linspace(0.03, 3.20, 90, dtype=float)
    b_grid = np.linspace(-1.20, 2.00, 90, dtype=float)
    p_grid = np.linspace(0.45, 2.60, 48, dtype=float)

    candidates: List[Dict[str, float]] = []

    for p in p_grid:
        dp = d_ref ** p
        for b in b_grid:
            x = a_grid[:, None] * dp[None, :] + b
            valid = np.all(x > 0.0, axis=1)
            if not bool(np.any(valid)):
                continue

            x_valid = x[valid]
            k_can = kernel_eta(canonical, x_valid)

            mse_pos = np.mean((k_can - k_ref[None, :]) ** 2, axis=1)
            mse_neg = np.mean((-k_can - k_ref[None, :]) ** 2, axis=1)
            use_neg = mse_neg < mse_pos
            mse = np.where(use_neg, mse_neg, mse_pos)

            k_take = min(3, int(len(mse)))
            if k_take <= 0:
                continue
            idx_local = np.argpartition(mse, k_take - 1)[:k_take]
            idx_global = np.where(valid)[0][idx_local]

            for local_i, global_i in zip(idx_local.tolist(), idx_global.tolist()):
                candidates.append(
                    {
                        "a": float(a_grid[global_i]),
                        "b": float(b),
                        "p": float(p),
                        "sign": -1.0 if bool(use_neg[local_i]) else 1.0,
                        "mse": float(mse[local_i]),
                    }
                )

        if len(candidates) > 3000:
            candidates = sorted(candidates, key=lambda r: r["mse"])[:700]

    candidates = sorted(candidates, key=lambda r: r["mse"])[:180]

    evaluated: List[Dict[str, object]] = []
    for c in candidates:
        ev = evaluate_candidate(
            canonical=canonical,
            refrozen=refrozen,
            d_ref=d_ref,
            k_ref=k_ref,
            a=float(c["a"]),
            b=float(c["b"]),
            p=float(c["p"]),
            sign=float(c["sign"]),
        )
        if not bool(ev.get("invalid", False)):
            evaluated.append(ev)

    if not evaluated:
        raise RuntimeError("No valid candidate reparameterizations found.")

    best_wave = min(evaluated, key=lambda r: float(r["rmse"]))
    best_gate = sorted(
        evaluated,
        key=lambda r: (
            -int(r["pass_count"]),
            float(r["rmse"]),
            -float(r["corr"]),
        ),
    )[0]

    return {
        "grid": {
            "a_min_max_n": [float(a_grid[0]), float(a_grid[-1]), int(len(a_grid))],
            "b_min_max_n": [float(b_grid[0]), float(b_grid[-1]), int(len(b_grid))],
            "p_min_max_n": [float(p_grid[0]), float(p_grid[-1]), int(len(p_grid))],
            "d_ref_range": [1.0, 24.0],
            "d_ref_points": int(len(d_ref)),
            "candidate_pool_evaluated": int(len(evaluated)),
        },
        "best_waveform": best_wave,
        "best_gate_candidate": best_gate,
        "top5_by_gate": sorted(
            evaluated,
            key=lambda r: (
                -int(r["pass_count"]),
                float(r["rmse"]),
            ),
        )[:5],
    }


def main() -> None:
    d2039 = json.loads((ROOT / "report_qw2039_derivation_compatible_refrozen_kernel_gate.json").read_text(encoding="utf-8"))

    canonical = Kernel(
        omega=float(math.pi / 4.0),
        phi=float(math.pi / 6.0),
        beta=0.01,
        eta=1.0,
    )
    refrozen = Kernel(
        omega=float(d2039["selected_kernel"]["omega"]),
        phi=float(d2039["selected_kernel"]["phi"]),
        beta=float(d2039["selected_kernel"]["beta"]),
        eta=float(d2039["selected_kernel"]["eta"]),
    )

    result = scan_reparameterizations(canonical=canonical, refrozen=refrozen)
    gate = result["best_gate_candidate"]
    pass_count = int(gate["pass_count"])
    total_flags = int(gate["total_flags"])

    if pass_count == total_flags:
        verdict = "CANONICAL_REFROZEN_REPARAMETERIZATION_PASS"
        readiness = "CANONICAL_SEMANTICS_RECOVERABLE_BY_MONOTONE_MAP"
    elif pass_count >= 5:
        verdict = "CANONICAL_REFROZEN_REPARAMETERIZATION_PARTIAL"
        readiness = "PARTIAL_RECOVERY_ONLY"
    else:
        verdict = "CANONICAL_REFROZEN_REPARAMETERIZATION_FAIL"
        readiness = "CANONICAL_SEMANTIC_DRIFT_CONFIRMED_INTERNAL"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "canonical_kernel": {
            "omega": canonical.omega,
            "phi": canonical.phi,
            "beta": canonical.beta,
            "eta": canonical.eta,
        },
        "refrozen_kernel_qw2039": {
            "omega": refrozen.omega,
            "phi": refrozen.phi,
            "beta": refrozen.beta,
            "eta": refrozen.eta,
        },
        "scan": result,
        "verdict": verdict,
        "readiness": readiness,
        "required_next_step": (
            "DEFINE_EXPLICIT_BRIDGE_OPERATOR_BETWEEN_CANONICAL_AND_EFFECTIVE_SEMANTICS"
            if verdict != "CANONICAL_REFROZEN_REPARAMETERIZATION_PASS"
            else "PROMOTE_BRIDGE_AS_CANONICAL_EQUIVALENCE_LAYER"
        ),
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    bw = result["best_waveform"]
    bg = result["best_gate_candidate"]
    lines = [
        "# RAPORT QW-2041: CANONICAL vs REFROZEN REPARAMETERIZATION AUDIT",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- Readiness: **{readiness}**",
        "",
        "## Kernels",
        (
            "- Canonical TeX omega/phi/beta/eta: "
            f"{canonical.omega:.6f} / {canonical.phi:.6f} / {canonical.beta:.6f} / {canonical.eta:.6f}"
        ),
        (
            "- Refrozen QW-2039 omega/phi/beta/eta: "
            f"{refrozen.omega:.6f} / {refrozen.phi:.6f} / {refrozen.beta:.6f} / {refrozen.eta:.6f}"
        ),
        "",
        "## Best Waveform Candidate (min RMSE)",
        f"- a/b/p/sign: {bw['a']:.6f} / {bw['b']:.6f} / {bw['p']:.6f} / {int(bw['sign'])}",
        f"- rmse: {bw['rmse']:.6f}",
        f"- corr: {bw['corr']:.6f}",
        f"- r2: {bw['r2']:.6f}",
        f"- affine_r2: {bw['affine']['r2']:.6f}",
        (
            "- node errors median/q95/max rel: "
            f"{bw['node_errors']['median_rel']:.6f} / {bw['node_errors']['q95_rel']:.6f} / {bw['node_errors']['max_rel']:.6f}"
        ),
        "",
        "## Best Gate Candidate (max hard flags)",
        f"- pass_count: {bg['pass_count']}/{bg['total_flags']}",
        f"- a/b/p/sign: {bg['a']:.6f} / {bg['b']:.6f} / {bg['p']:.6f} / {int(bg['sign'])}",
        f"- rmse: {bg['rmse']:.6f}",
        f"- corr: {bg['corr']:.6f}",
        f"- r2: {bg['r2']:.6f}",
        f"- affine_r2: {bg['affine']['r2']:.6f}",
        (
            "- node errors median/q95/max rel: "
            f"{bg['node_errors']['median_rel']:.6f} / {bg['node_errors']['q95_rel']:.6f} / {bg['node_errors']['max_rel']:.6f}"
        ),
        "",
        "## Gate Flags (best gate candidate)",
    ]
    for k, v in bg["flags"].items():
        lines.append(f"- {k}: {v}")

    lines.extend(
        [
            "",
            "## Required Next Step",
            f"- {out['required_next_step']}",
            "",
            "## Artifacts",
            f"- JSON: `{OUT_JSON.name}`",
        ]
    )

    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2041] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2041] Saved MD:   {OUT_MD.name}")
    print(f"[QW-2041] verdict={verdict} pass_count={pass_count}/{total_flags}")


if __name__ == "__main__":
    main()
