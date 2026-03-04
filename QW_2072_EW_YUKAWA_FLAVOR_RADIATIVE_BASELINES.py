#!/usr/bin/env python3
"""
QW-2072: EW + Yukawa + flavor RGE radiative baselines (non-closing).

Purpose:
- implement transparent baseline blocks for missing radiative channels,
- keep strict scientific labeling: baseline implemented, closure not yet achieved.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw2072_ew_yukawa_flavor_radiative_baselines.json"
OUT_MD = ROOT / "RAPORT_QW2072_EW_YUKAWA_FLAVOR_RADIATIVE_BASELINES.md"


def load_json(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def get_ref(groups: Dict, param_id: str):
    for _, items in groups.items():
        for item in items:
            if item["id"] == param_id:
                return item["value"]
    raise KeyError(f"Missing parameter in registry: {param_id}")


def solve_mw_from_alpha_gf_mz(alpha: float, g_f: float, m_z: float, delta_r: float) -> float:
    # M_W^2 (1 - M_W^2 / M_Z^2) = (pi alpha)/(sqrt(2) G_F) * 1/(1-delta_r)
    if not (0.0 < alpha and 0.0 < g_f and 0.0 < m_z):
        raise ValueError("alpha, G_F, M_Z must be positive.")
    if delta_r >= 1.0:
        raise ValueError("delta_r must satisfy delta_r < 1.")

    a = (math.pi * alpha) / (math.sqrt(2.0) * g_f) * (1.0 / (1.0 - delta_r))
    mz2 = m_z * m_z
    disc = 1.0 - 4.0 * a / mz2
    if disc < 0.0:
        raise ValueError("No real MW solution for given inputs.")
    x = 0.5 * mz2 * (1.0 + math.sqrt(disc))
    return float(math.sqrt(max(x, 0.0)))


def delta_r_required(alpha: float, g_f: float, m_w: float, m_z: float) -> float:
    a = (math.pi * alpha) / (math.sqrt(2.0) * g_f)
    denom = (m_w * m_w) * (1.0 - (m_w * m_w) / (m_z * m_z))
    if abs(denom) < 1e-18:
        raise ValueError("Invalid denominator for delta_r computation.")
    return float(1.0 - a / denom)


def normalized_rowwise(a: np.ndarray) -> np.ndarray:
    b = np.array(a, dtype=float)
    for i in range(b.shape[0]):
        s = float(np.linalg.norm(b[i, :]))
        if s > 0.0:
            b[i, :] /= s
    return b


def mean_rel_pct(pred: np.ndarray, ref: np.ndarray) -> float:
    rel = np.abs(pred - ref) / np.clip(np.abs(ref), 1e-12, None)
    return float(100.0 * np.mean(rel))


def main() -> None:
    reg = load_json("report_qw2068_sm_gr_parameter_registry.json")
    groups = reg["groups"]

    alpha_inv = float(get_ref(groups, "alpha_em_inv_mz"))
    alpha = 1.0 / alpha_inv
    sin2 = float(get_ref(groups, "sin2_theta_w_mz"))
    g_f = float(get_ref(groups, "g_f"))
    m_w_ref = float(get_ref(groups, "m_w"))
    m_z_ref = float(get_ref(groups, "m_z"))
    v_ref = float(get_ref(groups, "v_higgs"))

    sinw = math.sqrt(max(sin2, 0.0))
    cosw = math.sqrt(max(1.0 - sin2, 0.0))
    e_coupling = math.sqrt(4.0 * math.pi * alpha)
    g_coupling = e_coupling / max(sinw, 1e-15)
    gp_coupling = e_coupling / max(cosw, 1e-15)

    v_from_gf = (math.sqrt(2.0) * g_f) ** (-0.5)
    m_w_tree = 0.5 * g_coupling * v_from_gf
    m_z_tree = 0.5 * math.sqrt(g_coupling * g_coupling + gp_coupling * gp_coupling) * v_from_gf
    sin2_tree_from_masses = 1.0 - (m_w_tree * m_w_tree) / (m_z_tree * m_z_tree)
    m_w_delta_r0 = solve_mw_from_alpha_gf_mz(alpha=alpha, g_f=g_f, m_z=m_z_ref, delta_r=0.0)
    delta_r_req = delta_r_required(alpha=alpha, g_f=g_f, m_w=m_w_ref, m_z=m_z_ref)

    # Yukawa baselines from masses at reference v.
    mass_ids = [
        "m_top",
        "m_bottom",
        "m_charm",
        "m_tau",
        "m_muon",
        "m_electron",
        "m_up",
        "m_down",
        "m_strange",
    ]
    yukawa = {}
    for pid in mass_ids:
        m = float(get_ref(groups, pid))
        yukawa[pid.replace("m_", "y_")] = float(math.sqrt(2.0) * m / v_ref)

    # Flavor RGE toy baseline (non-closing): first-order drift proxy.
    ckm_ref = np.array(get_ref(groups, "ckm_matrix_abs"), dtype=float)
    pmns_ref = np.array(get_ref(groups, "pmns_matrix_abs"), dtype=float)
    mu_low = 1.0
    mu_high = 1000.0
    t = math.log(mu_high / mu_low) / (16.0 * math.pi * math.pi)

    yu = np.diag([yukawa["y_up"] ** 2, yukawa["y_charm"] ** 2, yukawa["y_top"] ** 2])
    yd = np.diag([yukawa["y_down"] ** 2, yukawa["y_strange"] ** 2, yukawa["y_bottom"] ** 2])
    ckm_high = normalized_rowwise(ckm_ref + t * (yu @ ckm_ref - ckm_ref @ yd))

    yl = np.diag([yukawa["y_electron"] ** 2, yukawa["y_muon"] ** 2, yukawa["y_tau"] ** 2])
    yn = np.diag([1.0e-12, 2.0e-12, 3.0e-12])  # placeholder neutrino-yukawa proxy
    pmns_high = normalized_rowwise(pmns_ref + t * (yl @ pmns_ref - pmns_ref @ yn))

    ckm_drift = mean_rel_pct(ckm_high, ckm_ref)
    pmns_drift = mean_rel_pct(pmns_high, pmns_ref)

    channel_updates: List[Dict] = [
        {
            "id": "electroweak_oblique_and_delta_r",
            "status": "implemented_baseline_nonclosing",
            "strict_level": "model_radiative_baseline",
            "closure_ready": False,
            "notes": "Tree+Delta_r baseline implemented; full EW loop closure still missing.",
        },
        {
            "id": "yukawa_running_and_threshold_matching",
            "status": "implemented_baseline_nonclosing",
            "strict_level": "model_radiative_baseline",
            "closure_ready": False,
            "notes": "Reference-scale Yukawa map implemented; full threshold/RG matching still missing.",
        },
        {
            "id": "ckm_pmns_rge_transport",
            "status": "implemented_baseline_nonclosing",
            "strict_level": "model_radiative_baseline",
            "closure_ready": False,
            "notes": "First-order transport proxy implemented; not a full precision RGE closure.",
        },
    ]

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "registry": "report_qw2068_sm_gr_parameter_registry.json",
        },
        "electroweak_baseline": {
            "alpha_em_mz": alpha,
            "sin2_theta_w_input": sin2,
            "e": e_coupling,
            "g": g_coupling,
            "g_prime": gp_coupling,
            "v_from_gf": v_from_gf,
            "m_w_tree": m_w_tree,
            "m_z_tree": m_z_tree,
            "sin2_from_tree_masses": sin2_tree_from_masses,
            "m_w_delta_r0_from_alpha_gf_mz": m_w_delta_r0,
            "delta_r_required_for_mw_ref": delta_r_req,
            "mw_ref": m_w_ref,
            "mz_ref": m_z_ref,
        },
        "yukawa_baseline": yukawa,
        "flavor_rge_baseline": {
            "mu_low_gev": mu_low,
            "mu_high_gev": mu_high,
            "t_log_over_16pi2": t,
            "ckm_mean_drift_rel_pct": ckm_drift,
            "pmns_mean_drift_rel_pct": pmns_drift,
            "ckm_high": ckm_high.tolist(),
            "pmns_high": pmns_high.tolist(),
        },
        "channel_updates": channel_updates,
        "verdict": "EW_YUKAWA_FLAVOR_RADIATIVE_BASELINES_IMPLEMENTED_NONCLOSING",
        "required_next_step": "UPGRADE_TO_FULL_PRECISION_EW_YUKAWA_CKM_PMNS_RGE_CLOSURE",
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2072: EW + YUKAWA + FLAVOR RADIATIVE BASELINES",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{out['verdict']}**",
        "",
        "## Electroweak Baseline",
        f"- v_from_gf: {v_from_gf:.6f} GeV",
        f"- m_w_tree: {m_w_tree:.6f} GeV",
        f"- m_w(delta_r=0): {m_w_delta_r0:.6f} GeV",
        f"- delta_r_required_for_mw_ref: {delta_r_req:.6f}",
        "",
        "## Flavor Transport Baseline",
        f"- CKM mean drift rel% (1->1000 GeV): {ckm_drift:.6f}",
        f"- PMNS mean drift rel% (1->1000 GeV): {pmns_drift:.6f}",
        "",
        "## Channel Updates",
    ]
    for c in channel_updates:
        lines.append(
            f"- {c['id']}: {c['status']} (closure_ready={c['closure_ready']})"
        )

    lines.extend(["", "## Artifact", f"- JSON: `{OUT_JSON.name}`"])
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2072] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2072] Saved MD:   {OUT_MD.name}")
    print(
        f"[QW-2072] verdict={out['verdict']} "
        f"delta_r_req={delta_r_req:.6f} ckm_drift={ckm_drift:.6f}"
    )


if __name__ == "__main__":
    main()
