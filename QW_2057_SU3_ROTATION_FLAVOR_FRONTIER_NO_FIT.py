#!/usr/bin/env python3
"""
QW-2057: SU(3)-like rotation flavor frontier (no fit).

Deterministic test:
- build sector unitary operators from kernel-driven rotation angles,
- no numerical fitting,
- evaluate CKM/PMNS closure envelope.
"""

from __future__ import annotations

import itertools
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw2057_su3_rotation_flavor_frontier_no_fit.json"
OUT_MD = ROOT / "RAPORT_QW2057_SU3_ROTATION_FLAVOR_FRONTIER_NO_FIT.md"


CKM_REF = np.array(
    [
        [0.97401, 0.22650, 0.00361],
        [0.22636, 0.97320, 0.04053],
        [0.00854, 0.03978, 0.99917],
    ],
    dtype=float,
)
PMNS_REF = np.array(
    [
        [0.821, 0.550, 0.150],
        [0.432, 0.582, 0.693],
        [0.378, 0.598, 0.707],
    ],
    dtype=float,
)


def kernel_fn(d: np.ndarray, omega: float, phi: float, beta: float, eta: float) -> np.ndarray:
    return np.cos(omega * d + phi) / (1.0 + beta * (d**eta))


def cyclic_distance_matrix(q: np.ndarray, modulus: int = 24) -> np.ndarray:
    dq = np.abs(q[:, None] - q[None, :]) % float(modulus)
    return np.minimum(dq, float(modulus) - dq)


def rel_mean_pct(pred: np.ndarray, ref: np.ndarray) -> float:
    return float(100.0 * np.mean(np.abs(pred - ref) / np.clip(ref, 1e-12, None)))


def r12(theta: float) -> np.ndarray:
    c, s = float(np.cos(theta)), float(np.sin(theta))
    return np.array([[c, s, 0.0], [-s, c, 0.0], [0.0, 0.0, 1.0]], dtype=float)


def r13(theta: float) -> np.ndarray:
    c, s = float(np.cos(theta)), float(np.sin(theta))
    return np.array([[c, 0.0, s], [0.0, 1.0, 0.0], [-s, 0.0, c]], dtype=float)


def r23(theta: float) -> np.ndarray:
    c, s = float(np.cos(theta)), float(np.sin(theta))
    return np.array([[1.0, 0.0, 0.0], [0.0, c, s], [0.0, -s, c]], dtype=float)


def sector_unitary_from_mode(qvec: np.ndarray, kernel: Dict[str, float], mode: str) -> np.ndarray:
    d = 1.0 + cyclic_distance_matrix(qvec, modulus=24)
    kabs = np.abs(kernel_fn(d, kernel["omega"], kernel["phi"], kernel["beta"], kernel["eta"]))
    ksgn = np.sign(kernel_fn(d, kernel["omega"], kernel["phi"], kernel["beta"], kernel["eta"]) + 1e-12)
    alpha = 1.0 + abs(kernel["phi"])

    if mode == "base":
        t12 = np.arctan2(kabs[0, 1], d[0, 1] ** alpha)
        t13 = np.arctan2(kabs[0, 2], d[0, 2] ** alpha)
        t23 = np.arctan2(kabs[1, 2], d[1, 2] ** alpha)
    elif mode == "enhanced":
        t12 = np.arctan2(kabs[0, 1] * (1.0 + kernel["omega"]), d[0, 1] ** (0.7 + abs(kernel["phi"])))
        t13 = np.arctan2(kabs[0, 2] * (1.0 + 0.5 * kernel["omega"]), d[0, 2] ** (0.7 + abs(kernel["phi"])))
        t23 = np.arctan2(kabs[1, 2] * (1.0 + 1.5 * kernel["omega"]), d[1, 2] ** (0.7 + abs(kernel["phi"])))
    elif mode == "signed":
        t12 = np.arctan2(ksgn[0, 1] * kabs[0, 1], d[0, 1] ** alpha)
        t13 = np.arctan2(ksgn[0, 2] * kabs[0, 2], d[0, 2] ** alpha)
        t23 = np.arctan2(ksgn[1, 2] * kabs[1, 2], d[1, 2] ** alpha)
    elif mode == "boost_local":
        t12 = np.arctan2(2.5 * kabs[0, 1], d[0, 1] ** (0.5 + 0.5 * abs(kernel["phi"])))
        t13 = np.arctan2(2.5 * kabs[0, 2], d[0, 2] ** (0.5 + 0.5 * abs(kernel["phi"])))
        t23 = np.arctan2(3.0 * kabs[1, 2], d[1, 2] ** (0.5 + 0.5 * abs(kernel["phi"])))
    else:
        raise ValueError(f"Unknown mode: {mode}")

    return r23(float(t23)) @ r13(float(t13)) @ r12(float(t12))


def qeff_from_mass(m_mev: float, gamma: float) -> float:
    return float(-(4.0 / gamma) * (np.log(m_mev / 173_000.0) / np.log(4.0)))


def main() -> None:
    d2049 = json.loads((ROOT / "report_qw2049_spectral_micro_stagec_intersection_gate.json").read_text(encoding="utf-8"))
    kernel = {k: float(v) for k, v in d2049["stagec_pool"]["selected_kernel"].items()}

    r1961 = json.loads((ROOT / "report_qw1961_noncircular_gamma_q_derivation_matrix.json").read_text(encoding="utf-8"))
    best = r1961["summary"]["best_noncircular"]
    q_assign_name = str(best["q_assignment"])
    q_assign = r1961["inputs"]["q_assignments"][q_assign_name]
    gamma = float(best["gamma_value"])
    delta_info = float(r1961["info_split_source"]["delta_info"])

    q_proxy_up = np.array([q_assign["Top"], q_assign["Charm"], q_assign["Muon"]], dtype=float)
    q_proxy_down = np.array([q_assign["Bottom"], q_assign["Tau"], q_assign["Muon"]], dtype=float)

    mass_quark_mev = {"u": 2.16, "d": 4.67, "s": 93.4, "c": 1270.0, "b": 4180.0, "t": 173000.0}
    q_quark_up = np.array(
        [
            qeff_from_mass(mass_quark_mev["u"], gamma) - delta_info,
            qeff_from_mass(mass_quark_mev["c"], gamma) - delta_info,
            qeff_from_mass(mass_quark_mev["t"], gamma) - delta_info,
        ],
        dtype=float,
    )
    q_quark_down = np.array(
        [
            qeff_from_mass(mass_quark_mev["d"], gamma) - delta_info,
            qeff_from_mass(mass_quark_mev["s"], gamma) - delta_info,
            qeff_from_mass(mass_quark_mev["b"], gamma) - delta_info,
        ],
        dtype=float,
    )

    q_nu = np.array([0.0, 1.0, 2.0], dtype=float)
    q_lep = np.array([q_assign["Electron"], q_assign["Muon"], q_assign["Tau"]], dtype=float)

    q_schemes = {
        "proxy_old": {"q_up": q_proxy_up, "q_down": q_proxy_down},
        "quark_mass_inversion": {"q_up": q_quark_up, "q_down": q_quark_down},
    }
    modes = ["base", "enhanced", "signed", "boost_local"]
    thresholds = {"ckm_mean_rel_pct_max": 15.0, "pmns_mean_rel_pct_max": 15.0}

    rows: List[Dict[str, object]] = []
    for q_name, qs in q_schemes.items():
        for mode_up, mode_down, mode_nu, mode_lep in itertools.product(modes, modes, modes, modes):
            u_up = sector_unitary_from_mode(qs["q_up"], kernel, mode_up)
            u_down = sector_unitary_from_mode(qs["q_down"], kernel, mode_down)
            u_nu = sector_unitary_from_mode(q_nu, kernel, mode_nu)
            u_lep = sector_unitary_from_mode(q_lep, kernel, mode_lep)

            ckm_pred = np.abs(u_up.T @ u_down)
            pmns_pred = np.abs(u_nu.T @ u_lep)
            ckm_err = rel_mean_pct(ckm_pred, CKM_REF)
            pmns_err = rel_mean_pct(pmns_pred, PMNS_REF)

            flags = {
                "ckm_mean_rel_pct_le_max": bool(ckm_err <= thresholds["ckm_mean_rel_pct_max"]),
                "pmns_mean_rel_pct_le_max": bool(pmns_err <= thresholds["pmns_mean_rel_pct_max"]),
            }
            rows.append(
                {
                    "q_scheme": q_name,
                    "modes": {"up": mode_up, "down": mode_down, "nu": mode_nu, "lep": mode_lep},
                    "metrics": {"ckm_mean_rel_pct": ckm_err, "pmns_mean_rel_pct": pmns_err},
                    "flags": flags,
                    "pass_count": int(sum(1 for v in flags.values() if v)),
                    "total_flags": int(len(flags)),
                    "all_pass": bool(all(flags.values())),
                }
            )

    rows_sorted = sorted(
        rows,
        key=lambda r: (
            float(max(0.0, r["metrics"]["ckm_mean_rel_pct"] - thresholds["ckm_mean_rel_pct_max"])),
            float(max(0.0, r["metrics"]["pmns_mean_rel_pct"] - thresholds["pmns_mean_rel_pct_max"])),
            float(r["metrics"]["ckm_mean_rel_pct"] + r["metrics"]["pmns_mean_rel_pct"]),
        ),
    )
    any_all_pass = any(bool(r["all_pass"]) for r in rows_sorted)
    best_row = rows_sorted[0]
    best_ckm = min(rows, key=lambda r: float(r["metrics"]["ckm_mean_rel_pct"]))
    best_pmns = min(rows, key=lambda r: float(r["metrics"]["pmns_mean_rel_pct"]))

    verdict = (
        "SU3_ROTATION_FLAVOR_FRONTIER_HAS_CLOSURE_CANDIDATE"
        if any_all_pass
        else "SU3_ROTATION_FLAVOR_FRONTIER_FAILS_TO_CLOSE_BOTH_CKM_PMNS"
    )
    required_next = (
        "FREEZE_SU3_ROTATION_CANDIDATE_AND_RUN_STRICT_TRIAD_GATE"
        if any_all_pass
        else "MOVE_TO_DERIVATIONAL_NONABELIAN_FLAVOR_GENERATOR_WITH_EXPLICIT_CKM_PMNS_STRUCTURE"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "kernel_source": "report_qw2049_spectral_micro_stagec_intersection_gate.json:stagec_pool.selected_kernel",
        "kernel": kernel,
        "mass_chain_source": {
            "file": "report_qw1961_noncircular_gamma_q_derivation_matrix.json",
            "q_assignment_name": q_assign_name,
            "gamma_value": gamma,
            "delta_info": delta_info,
        },
        "thresholds": thresholds,
        "rows_sorted": rows_sorted,
        "best_row": best_row,
        "best_ckm": best_ckm,
        "best_pmns": best_pmns,
        "any_all_pass": any_all_pass,
        "verdict": verdict,
        "required_next_step": required_next,
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2057: SU3 ROTATION FLAVOR FRONTIER (NO FIT)",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- any_all_pass: {any_all_pass}",
        "",
        "## Best Closure Row",
        f"- q_scheme: {best_row['q_scheme']}",
        f"- modes up/down/nu/lep: {best_row['modes']['up']} / {best_row['modes']['down']} / {best_row['modes']['nu']} / {best_row['modes']['lep']}",
        f"- CKM mean rel%: {best_row['metrics']['ckm_mean_rel_pct']:.3f}",
        f"- PMNS mean rel%: {best_row['metrics']['pmns_mean_rel_pct']:.3f}",
        "",
        "## Best CKM Row",
        f"- q_scheme: {best_ckm['q_scheme']}",
        f"- modes up/down/nu/lep: {best_ckm['modes']['up']} / {best_ckm['modes']['down']} / {best_ckm['modes']['nu']} / {best_ckm['modes']['lep']}",
        f"- CKM mean rel%: {best_ckm['metrics']['ckm_mean_rel_pct']:.3f}",
        f"- PMNS mean rel%: {best_ckm['metrics']['pmns_mean_rel_pct']:.3f}",
        "",
        "## Best PMNS Row",
        f"- q_scheme: {best_pmns['q_scheme']}",
        f"- modes up/down/nu/lep: {best_pmns['modes']['up']} / {best_pmns['modes']['down']} / {best_pmns['modes']['nu']} / {best_pmns['modes']['lep']}",
        f"- CKM mean rel%: {best_pmns['metrics']['ckm_mean_rel_pct']:.3f}",
        f"- PMNS mean rel%: {best_pmns['metrics']['pmns_mean_rel_pct']:.3f}",
        "",
        "## Top 10 Rows",
    ]
    for i, row in enumerate(rows_sorted[:10], start=1):
        lines.append(
            (
                f"- {i}. {row['q_scheme']} | {row['modes']['up']}/{row['modes']['down']}/{row['modes']['nu']}/{row['modes']['lep']} | "
                f"CKM={row['metrics']['ckm_mean_rel_pct']:.3f} | PMNS={row['metrics']['pmns_mean_rel_pct']:.3f}"
            )
        )

    lines.extend(
        [
            "",
            "## Required Next Step",
            f"- {required_next}",
            "",
            "## Artifacts",
            f"- JSON: `{OUT_JSON.name}`",
        ]
    )
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2057] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2057] Saved MD:   {OUT_MD.name}")
    print(f"[QW-2057] verdict={verdict} any_all_pass={any_all_pass}")


if __name__ == "__main__":
    main()
