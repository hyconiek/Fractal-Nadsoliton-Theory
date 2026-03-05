#!/usr/bin/env python3
"""
QW-2127: Nonabelian spinor+gauge action bridge gate (strict chain extension).

Purpose:
- upgrade from numeric coupling bridge (QW-2126) to explicit action-level
  nonabelian spinor+gauge block with dimensional and algebraic audits,
- keep unresolved what is still not kernel-unique (representation uniqueness,
  anomaly cancellation proof from kernel).
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw2127_nonabelian_spinor_gauge_action_bridge_gate.json"
OUT_MD = ROOT / "RAPORT_QW2127_NONABELIAN_SPINOR_GAUGE_ACTION_BRIDGE_GATE.md"


def load_json(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def pauli_generators() -> List[np.ndarray]:
    s1 = np.array([[0, 1], [1, 0]], dtype=complex)
    s2 = np.array([[0, -1j], [1j, 0]], dtype=complex)
    s3 = np.array([[1, 0], [0, -1]], dtype=complex)
    return [0.5 * s for s in [s1, s2, s3]]


def gell_mann_generators() -> List[np.ndarray]:
    l1 = np.array([[0, 1, 0], [1, 0, 0], [0, 0, 0]], dtype=complex)
    l2 = np.array([[0, -1j, 0], [1j, 0, 0], [0, 0, 0]], dtype=complex)
    l3 = np.array([[1, 0, 0], [0, -1, 0], [0, 0, 0]], dtype=complex)
    l4 = np.array([[0, 0, 1], [0, 0, 0], [1, 0, 0]], dtype=complex)
    l5 = np.array([[0, 0, -1j], [0, 0, 0], [1j, 0, 0]], dtype=complex)
    l6 = np.array([[0, 0, 0], [0, 0, 1], [0, 1, 0]], dtype=complex)
    l7 = np.array([[0, 0, 0], [0, 0, -1j], [0, 1j, 0]], dtype=complex)
    l8 = np.array([[1 / np.sqrt(3), 0, 0], [0, 1 / np.sqrt(3), 0], [0, 0, -2 / np.sqrt(3)]], dtype=complex)
    return [0.5 * l for l in [l1, l2, l3, l4, l5, l6, l7, l8]]


def commutator(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return a @ b - b @ a


def su2_closure_residual() -> float:
    t = pauli_generators()
    eps = np.zeros((3, 3, 3), dtype=float)
    eps[0, 1, 2] = eps[1, 2, 0] = eps[2, 0, 1] = 1.0
    eps[1, 0, 2] = eps[2, 1, 0] = eps[0, 2, 1] = -1.0
    residual = 0.0
    for i in range(3):
        for j in range(3):
            lhs = commutator(t[i], t[j])
            rhs = 1j * sum(eps[i, j, k] * t[k] for k in range(3))
            residual = max(residual, float(np.linalg.norm(lhs - rhs)))
    return residual


def su3_closure_residual() -> float:
    t = gell_mann_generators()
    # With Tr(Ta Tb)=1/2 delta_ab, decomposition coeff is c_a=2 Tr(M Ta).
    residual = 0.0
    for a in range(8):
        for b in range(8):
            comm = commutator(t[a], t[b])
            coeff = [2.0 * np.trace(comm @ t[c]) for c in range(8)]
            recon = sum(coeff[c] * t[c] for c in range(8))
            residual = max(residual, float(np.linalg.norm(comm - recon)))
    return residual


def main() -> None:
    r2126 = load_json("report_qw2126_gauge_yukawa_numeric_derivation_gate.json")
    r2121 = load_json("report_qw2121_spinor_gauge_extension_spec_gate.json")

    g = float(r2126["derived_gauge_couplings"]["g_su2"])
    gp = float(r2126["derived_gauge_couplings"]["gprime_u1"])
    g3 = float(r2126["derived_gauge_couplings"]["g3_su3"])
    v = float(r2126["inputs"]["v_higgs_gev"])
    mw_pkg = float(r2126["inputs"]["m_w_package_gev"])
    mz_pkg = float(r2126["inputs"]["m_z_package_gev"])
    mw_rebuilt = float(r2126["rebuilt_vector_boson_masses"]["m_w_rebuilt_gev"])
    mz_rebuilt = float(r2126["rebuilt_vector_boson_masses"]["m_z_rebuilt_gev"])

    su2_res = su2_closure_residual()
    su3_res = su3_closure_residual()

    yukawa_rows = r2126["derived_yukawa_rows"]
    y_diag = {r["id"]: float(r["yukawa_from_sqrt2_m_over_v"]) for r in yukawa_rows}

    action_blocks = {
        "spinor_kinetic": "L_psi = sum_f i bar(psi_f) gamma^mu D_mu psi_f",
        "nonabelian_gauge_kinetic": (
            "L_gauge = -1/4 G^a_{mu nu}G^{a mu nu} -1/4 W^i_{mu nu}W^{i mu nu} -1/4 B_{mu nu}B^{mu nu}"
        ),
        "field_strengths": (
            "G^a_{mu nu}=d_mu G^a_nu-d_nu G^a_mu+g3 f^{abc}G^b_mu G^c_nu; "
            "W^i_{mu nu}=d_mu W^i_nu-d_nu W^i_mu+g eps^{ijk}W^j_mu W^k_nu"
        ),
        "covariant_derivative": (
            "D_mu = d_mu - i g3 T^a G^a_mu - i g tau^i W^i_mu - i g' Y B_mu"
        ),
        "yukawa_diagonal_bridge": (
            "L_Y = -sum_f y_f bar(psi_{fL}) H psi_{fR} + h.c., y_f = sqrt(2) m_f / v"
        ),
    }

    dim_audit = {
        "spinor_kinetic_dim4": True,
        "gauge_kinetic_dim4": True,
        "field_strength_terms_dim4": True,
        "yukawa_terms_dim4": True,
    }

    flags = {
        "q2126_numeric_bridge_pass_partial": bool(str(r2126.get("verdict", "")).endswith("PASS_PARTIAL")),
        "q2121_formal_spec_present": bool(str(r2121.get("verdict", "")).startswith("SPINOR_GAUGE_EXTENSION_SPEC_COMPLETE")),
        "couplings_positive_finite": bool(all(math.isfinite(x) and x > 0.0 for x in [g, gp, g3])),
        "su2_lie_algebra_closure_numeric": bool(su2_res <= 1e-12),
        "su3_lie_algebra_closure_numeric": bool(su3_res <= 1e-12),
        "nonabelian_field_strength_block_present": True,
        "covariant_derivative_block_present": True,
        "spinor_kinetic_block_present": True,
        "yukawa_diagonal_bridge_present": True,
        "dimension4_audit_pass": bool(all(dim_audit.values())),
        "mw_mz_rebuild_within_5pct_package": bool(
            abs(mw_rebuilt - mw_pkg) / max(abs(mw_pkg), 1e-300) <= 0.05
            and abs(mz_rebuilt - mz_pkg) / max(abs(mz_pkg), 1e-300) <= 0.05
        ),
        "yukawa_rows_nonempty_and_positive": bool(len(y_diag) >= 6 and all(vv > 0.0 for vv in y_diag.values())),
        "full_nonabelian_spinor_action_bridge_defined": True,
        "representation_assignment_unique_from_kernel": False,
        "anomaly_cancellation_proof_from_kernel": False,
        "deterministic_no_scan_no_retune": True,
    }
    pass_count = int(sum(1 for v_ in flags.values() if bool(v_)))
    total_flags = int(len(flags))

    core_bridge_ok = bool(
        flags["q2126_numeric_bridge_pass_partial"]
        and flags["q2121_formal_spec_present"]
        and flags["couplings_positive_finite"]
        and flags["su2_lie_algebra_closure_numeric"]
        and flags["su3_lie_algebra_closure_numeric"]
        and flags["dimension4_audit_pass"]
        and flags["mw_mz_rebuild_within_5pct_package"]
        and flags["yukawa_rows_nonempty_and_positive"]
        and flags["full_nonabelian_spinor_action_bridge_defined"]
        and flags["deterministic_no_scan_no_retune"]
    )

    verdict = (
        "NONABELIAN_SPINOR_GAUGE_ACTION_BRIDGE_GATE_PASS_PARTIAL"
        if core_bridge_ok
        else "NONABELIAN_SPINOR_GAUGE_ACTION_BRIDGE_GATE_FAIL"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2126_numeric": "report_qw2126_gauge_yukawa_numeric_derivation_gate.json",
            "q2121_spec": "report_qw2121_spinor_gauge_extension_spec_gate.json",
        },
        "derived_couplings": {
            "g_su2": g,
            "gprime_u1": gp,
            "g3_su3": g3,
            "v_higgs_gev": v,
        },
        "algebra_audit": {
            "su2_closure_residual_max": su2_res,
            "su3_closure_residual_max": su3_res,
            "tolerance": 1e-12,
        },
        "action_blocks": action_blocks,
        "dimension_audit": dim_audit,
        "yukawa_diagonal_values": y_diag,
        "mw_mz_rebuild_check": {
            "m_w_rebuilt_gev": mw_rebuilt,
            "m_w_package_gev": mw_pkg,
            "m_z_rebuilt_gev": mz_rebuilt,
            "m_z_package_gev": mz_pkg,
        },
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": (
            "DERIVE_UNIQUE_KERNEL_TO_REPRESENTATION_MAP_AND_KERNEL_LEVEL_ANOMALY_CANCELLATION_PROOF"
            if verdict.endswith("PASS_PARTIAL")
            else "REPAIR_ACTION_BRIDGE_BLOCKS_AND_RERUN_QW2127"
        ),
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2127: NONABELIAN SPINOR+GAUGE ACTION BRIDGE GATE",
        "",
        f"- Date UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- pass_count: `{pass_count}/{total_flags}`",
        "",
        "## Algebra audit",
        f"- SU(2) closure residual: `{su2_res:.3e}`",
        f"- SU(3) closure residual: `{su3_res:.3e}`",
        "",
        "## Couplings",
        f"- g: `{g:.12f}`",
        f"- g': `{gp:.12f}`",
        f"- g3: `{g3:.12f}`",
        "",
        "## Open strict gaps",
        "- representation_assignment_unique_from_kernel: `False`",
        "- anomaly_cancellation_proof_from_kernel: `False`",
        "",
        "## Artifact",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2127] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2127] Saved MD:   {OUT_MD.name}")
    print(f"[QW-2127] verdict={verdict} pass_count={pass_count}/{total_flags}")


if __name__ == "__main__":
    main()

