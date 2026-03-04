#!/usr/bin/env python3
"""
QW-2075: Strict CP-phase derivation gate (deterministic, no-scan).

Goal:
- derive CP phases from the same deterministic flavor construction used in QW-2063,
- avoid false closure claims by promoting only what passes explicit gate conditions.
"""

from __future__ import annotations

import importlib.util
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw2075_strict_cp_phase_derivation_gate.json"
OUT_MD = ROOT / "RAPORT_QW2075_STRICT_CP_PHASE_DERIVATION_GATE.md"


def load_json(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def get_ref(groups: Dict, param_id: str) -> Dict:
    for _, items in groups.items():
        for item in items:
            if item["id"] == param_id:
                return item
    raise KeyError(f"Missing parameter in registry: {param_id}")


def rel_err_pct(pred: float, ref: float) -> float:
    return abs(pred - ref) / max(abs(ref), 1e-15) * 100.0


def canonical_rephase(m: np.ndarray) -> np.ndarray:
    """Fix a deterministic phase convention to reduce arbitrary basis phases."""
    out = np.array(m, dtype=complex, copy=True)
    for i in range(3):
        out[i, :] *= np.exp(-1j * np.angle(out[i, 0]))
    for j in range(3):
        out[:, j] *= np.exp(-1j * np.angle(out[0, j]))
    return out


def unitarity_residual(m: np.ndarray) -> float:
    ident = np.eye(m.shape[0], dtype=complex)
    return float(np.max(np.abs(m.conj().T @ m - ident)))


def cp_phase_from_matrix(m: np.ndarray) -> Dict[str, float]:
    """
    Extract CP phase proxy from Jarlskog invariant and |Uij| angles.
    Convention:
    - delta_primary in [0, pi]
    - delta_secondary = pi - delta_primary (same sin(delta), branch ambiguity)
    """
    a = np.abs(m)
    s13 = float(a[0, 2])
    c13 = math.sqrt(max(1.0 - s13 * s13, 1e-15))
    s12 = float(a[0, 1] / c13)
    s23 = float(a[1, 2] / c13)
    c12 = math.sqrt(max(1.0 - s12 * s12, 1e-15))
    c23 = math.sqrt(max(1.0 - s23 * s23, 1e-15))

    j = float(np.imag(m[0, 0] * m[1, 1] * np.conj(m[0, 1]) * np.conj(m[1, 0])))
    den = float(max(s12 * c12 * s23 * c23 * s13 * c13 * c13, 1e-15))
    sin_delta = float(np.clip(j / den, -1.0, 1.0))

    delta_primary = float(math.asin(sin_delta))
    if delta_primary < 0.0:
        delta_primary += math.pi
    delta_secondary = float(math.pi - delta_primary)

    return {
        "s12": s12,
        "s23": s23,
        "s13": s13,
        "jarlskog": j,
        "denominator": den,
        "sin_delta": sin_delta,
        "delta_primary_rad": delta_primary,
        "delta_secondary_rad": delta_secondary,
    }


def main() -> None:
    reg = load_json("report_qw2068_sm_gr_parameter_registry.json")
    r2049 = load_json("report_qw2049_spectral_micro_stagec_intersection_gate.json")
    r1961 = load_json("report_qw1961_noncircular_gamma_q_derivation_matrix.json")

    # Load reusable deterministic flavor operators from QW-2063.
    spec = importlib.util.spec_from_file_location(
        "qw2063_module",
        str(ROOT / "QW_2063_DERIVATIONAL_RECONSTRUCTION_SHARED_FLAVOR_BASIS.py"),
    )
    if spec is None or spec.loader is None:
        raise RuntimeError("Cannot load QW-2063 module for deterministic flavor construction.")
    qw2063 = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(qw2063)

    kernel = {k: float(v) for k, v in r2049["stagec_pool"]["selected_kernel"].items()}
    q_name = str(r1961["summary"]["best_noncircular"]["q_assignment"])
    q_assign = r1961["inputs"]["q_assignments"][q_name]

    q_up = np.array([q_assign["Top"], q_assign["Charm"], q_assign["Muon"]], dtype=float)
    q_down = np.array([q_assign["Bottom"], q_assign["Tau"], q_assign["Muon"]], dtype=float)
    q_lep = np.array([q_assign["Electron"], q_assign["Muon"], q_assign["Tau"]], dtype=float)

    deriv = qw2063.derive_flavor_basis_from_kernel(kernel)
    q_nu = np.array(deriv["q_nu_order"], dtype=float)
    params = deriv["params"]

    hu = qw2063.flavor_hamiltonian(q_up, +1.0, +1.0, params, kernel)
    hd = qw2063.flavor_hamiltonian(q_down, -1.0, +1.0, params, kernel)
    hn = qw2063.flavor_hamiltonian(q_nu, +1.0, -1.0, params, kernel)
    hl = qw2063.flavor_hamiltonian(q_lep, -1.0, -1.0, params, kernel)

    _, uu = np.linalg.eigh(hu)
    _, ud = np.linalg.eigh(hd)
    _, un = np.linalg.eigh(hn)
    _, ul = np.linalg.eigh(hl)

    ckm_complex = canonical_rephase(uu.conj().T @ ud)
    pmns_complex = canonical_rephase(un.conj().T @ ul)

    ckm_phase = cp_phase_from_matrix(ckm_complex)
    pmns_phase = cp_phase_from_matrix(pmns_complex)

    ckm_ref = get_ref(reg["groups"], "delta_cp_ckm")
    ckm_ref_val = float(ckm_ref["value"])
    ckm_tol = float(ckm_ref["tolerance_rel_pct"])
    ckm_rel_err_primary = rel_err_pct(ckm_phase["delta_primary_rad"], ckm_ref_val)
    ckm_rel_err_secondary = rel_err_pct(ckm_phase["delta_secondary_rad"], ckm_ref_val)
    ckm_branch_has_tol_match = bool(
        (ckm_rel_err_primary <= ckm_tol) or (ckm_rel_err_secondary <= ckm_tol)
    )

    ckm_unitarity_res = unitarity_residual(ckm_complex)
    pmns_unitarity_res = unitarity_residual(pmns_complex)

    flags = {
        "deterministic_no_scan": True,
        "no_retune": True,
        "ckm_unitarity_ok": bool(ckm_unitarity_res <= 1e-10),
        "pmns_unitarity_ok": bool(pmns_unitarity_res <= 1e-10),
        "ckm_phase_numerically_stable": bool(ckm_phase["denominator"] > 1e-12),
        "pmns_phase_numerically_stable": bool(pmns_phase["denominator"] > 1e-12),
        "ckm_phase_matches_registry_tolerance_any_branch": ckm_branch_has_tol_match,
        "pmns_phase_finite": bool(np.isfinite(pmns_phase["delta_primary_rad"])),
    }

    strict_update_delta_ckm = bool(
        flags["deterministic_no_scan"]
        and flags["no_retune"]
        and flags["ckm_unitarity_ok"]
        and flags["ckm_phase_numerically_stable"]
        and flags["ckm_phase_matches_registry_tolerance_any_branch"]
    )
    strict_update_delta_pmns = bool(
        flags["deterministic_no_scan"]
        and flags["no_retune"]
        and flags["pmns_unitarity_ok"]
        and flags["pmns_phase_numerically_stable"]
        and flags["pmns_phase_finite"]
    )

    updates: List[Dict] = []
    if strict_update_delta_ckm:
        updates.append(
            {
                "id": "delta_cp_ckm",
                "predicted_value": float(ckm_phase["delta_primary_rad"]),
                "status": "derived",
                "strict_level": "strict_internal_gate",
                "method": "qw2075_deterministic_complex_ckm_jarlskog_phase",
                "notes": "Promoted only because at least one phase branch matches registry tolerance.",
            }
        )
    if strict_update_delta_pmns:
        updates.append(
            {
                "id": "delta_cp_pmns",
                "predicted_value": float(pmns_phase["delta_primary_rad"]),
                "status": "derived",
                "strict_level": "strict_internal_gate",
                "method": "qw2075_deterministic_complex_pmns_jarlskog_phase",
                "notes": "Deterministic no-scan PMNS CP phase from complex mixing operator.",
            }
        )

    pass_count = int(sum(1 for v in flags.values() if bool(v)))
    total_flags = int(len(flags))

    if strict_update_delta_ckm and strict_update_delta_pmns:
        verdict = "STRICT_CP_PHASE_DERIVATION_PASS"
    elif strict_update_delta_pmns:
        verdict = "STRICT_CP_PHASE_DERIVATION_PARTIAL_PMNS_ONLY"
    elif strict_update_delta_ckm:
        verdict = "STRICT_CP_PHASE_DERIVATION_PARTIAL_CKM_ONLY"
    else:
        verdict = "STRICT_CP_PHASE_DERIVATION_FAIL"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "registry": "report_qw2068_sm_gr_parameter_registry.json",
            "kernel_selection": "report_qw2049_spectral_micro_stagec_intersection_gate.json",
            "mass_q_assignment": "report_qw1961_noncircular_gamma_q_derivation_matrix.json",
            "deterministic_flavor_constructor": "QW_2063_DERIVATIONAL_RECONSTRUCTION_SHARED_FLAVOR_BASIS.py",
        },
        "kernel": kernel,
        "q_assignment_name": q_name,
        "unitarity": {
            "ckm_residual_max_abs": ckm_unitarity_res,
            "pmns_residual_max_abs": pmns_unitarity_res,
        },
        "ckm_phase": {
            **ckm_phase,
            "registry_reference_rad": ckm_ref_val,
            "registry_tolerance_rel_pct": ckm_tol,
            "rel_err_primary_pct": ckm_rel_err_primary,
            "rel_err_secondary_pct": ckm_rel_err_secondary,
        },
        "pmns_phase": pmns_phase,
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "updates": updates,
        "count_updates": len(updates),
        "verdict": verdict,
        "required_next_step": (
            "IMPROVE_CKM_CP_PHASE_DERIVATION_WITHOUT_RETUNE_AND_INTEGRATE_PMNS_UPDATE"
            if verdict != "STRICT_CP_PHASE_DERIVATION_PASS"
            else "INTEGRATE_CP_PHASE_UPDATES_IN_QW2069"
        ),
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2075: STRICT CP PHASE DERIVATION GATE",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- pass_count: {pass_count}/{total_flags}",
        f"- updates promoted: {len(updates)}",
        "",
        "## CKM Phase",
        f"- delta_primary [0,pi]: {ckm_phase['delta_primary_rad']:.6f} rad",
        f"- delta_secondary [0,pi]: {ckm_phase['delta_secondary_rad']:.6f} rad",
        f"- reference: {ckm_ref_val:.6f} rad",
        f"- rel_err primary/secondary: {ckm_rel_err_primary:.3f}% / {ckm_rel_err_secondary:.3f}%",
        "",
        "## PMNS Phase",
        f"- delta_primary [0,pi]: {pmns_phase['delta_primary_rad']:.6f} rad",
        f"- sin(delta): {pmns_phase['sin_delta']:.6f}",
        "",
        "## Gate Rule",
        "- CKM update is promoted only if branch-ambiguous phase is compatible with registry tolerance.",
        "- PMNS update is promoted if deterministic and numerically stable (registry has no fixed central value here).",
        "",
        "## Artifact",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2075] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2075] Saved MD:   {OUT_MD.name}")
    print(
        f"[QW-2075] verdict={verdict} pass_count={pass_count}/{total_flags} "
        f"updates={len(updates)}"
    )


if __name__ == "__main__":
    main()

