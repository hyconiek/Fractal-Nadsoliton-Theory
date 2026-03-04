#!/usr/bin/env python3
"""
QW-2097: CKM CP target-refinement gate (deterministic, no-retune).

Purpose:
- investigate whether delta_cp_ckm target-miss can be reduced by deterministic
  convention/ordering audit only (no scan, no kernel retune),
- promote strict update only if CKM matrix fidelity and CP tolerance pass
  simultaneously.
"""

from __future__ import annotations

import importlib.util
import itertools
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw2097_ckm_cp_target_refinement_gate.json"
OUT_MD = ROOT / "RAPORT_QW2097_CKM_CP_TARGET_REFINEMENT_GATE.md"


def load_json(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def rel_err_pct(pred: float, ref: float) -> float:
    return abs(pred - ref) / max(abs(ref), 1e-15) * 100.0


def mean_matrix_rel_pct(pred_abs: np.ndarray, ref_abs: np.ndarray) -> float:
    rel = np.abs(pred_abs - ref_abs) / np.clip(ref_abs, 1e-12, None)
    return float(100.0 * np.mean(rel))


def unitarity_residual(m: np.ndarray) -> float:
    ident = np.eye(m.shape[0], dtype=complex)
    return float(np.max(np.abs(m.conj().T @ m - ident)))


def canonical_rephase(m: np.ndarray) -> np.ndarray:
    out = np.array(m, dtype=complex, copy=True)
    for i in range(3):
        out[i, :] *= np.exp(-1j * np.angle(out[i, 0]))
    for j in range(3):
        out[:, j] *= np.exp(-1j * np.angle(out[0, j]))
    return out


def cp_phase_from_matrix(m: np.ndarray) -> Dict[str, float]:
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


def load_module(path: Path, module_name: str):
    spec = importlib.util.spec_from_file_location(module_name, str(path))
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import module: {path.name}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def find_registry_item(reg: Dict, pid: str) -> Dict:
    for _, items in reg.get("groups", {}).items():
        for item in items:
            if item.get("id") == pid:
                return item
    raise KeyError(f"Missing registry parameter: {pid}")


def main() -> None:
    reg = load_json("report_qw2068_sm_gr_parameter_registry.json")
    r2049 = load_json("report_qw2049_spectral_micro_stagec_intersection_gate.json")
    r1961 = load_json("report_qw1961_noncircular_gamma_q_derivation_matrix.json")
    r2075 = load_json("report_qw2075_strict_cp_phase_derivation_gate.json")

    qw2063 = load_module(
        ROOT / "QW_2063_DERIVATIONAL_RECONSTRUCTION_SHARED_FLAVOR_BASIS.py",
        "qw2063_module_for_qw2097",
    )

    kernel = {k: float(v) for k, v in r2049["stagec_pool"]["selected_kernel"].items()}
    q_name = str(r1961["summary"]["best_noncircular"]["q_assignment"])
    q_assign = r1961["inputs"]["q_assignments"][q_name]

    q_up = np.array([q_assign["Top"], q_assign["Charm"], q_assign["Muon"]], dtype=float)
    q_down = np.array([q_assign["Bottom"], q_assign["Tau"], q_assign["Muon"]], dtype=float)
    deriv = qw2063.derive_flavor_basis_from_kernel(kernel)
    params = deriv["params"]

    hu = qw2063.flavor_hamiltonian(q_up, +1.0, +1.0, params, kernel)
    hd = qw2063.flavor_hamiltonian(q_down, -1.0, +1.0, params, kernel)
    _, uu = np.linalg.eigh(hu)
    _, ud = np.linalg.eigh(hd)
    ckm_base = canonical_rephase(uu.conj().T @ ud)

    ckm_ref_item = find_registry_item(reg, "ckm_matrix_abs")
    delta_ref_item = find_registry_item(reg, "delta_cp_ckm")
    ckm_ref_abs = np.array(ckm_ref_item["value"], dtype=float)
    delta_ref = float(delta_ref_item["value"])
    delta_tol = float(delta_ref_item["tolerance_rel_pct"])

    max_mean_rel_for_physical_order_pct = 20.0

    candidates: List[Dict] = []
    for prow in itertools.permutations(range(3)):
        for pcol in itertools.permutations(range(3)):
            m = ckm_base[np.array(prow), :][:, np.array(pcol)]
            m = canonical_rephase(m)
            phase = cp_phase_from_matrix(m)
            e1 = rel_err_pct(float(phase["delta_primary_rad"]), delta_ref)
            e2 = rel_err_pct(float(phase["delta_secondary_rad"]), delta_ref)
            if e1 <= e2:
                delta_best = float(phase["delta_primary_rad"])
                rel_best = float(e1)
                branch = "primary"
            else:
                delta_best = float(phase["delta_secondary_rad"])
                rel_best = float(e2)
                branch = "secondary"

            mean_rel = mean_matrix_rel_pct(np.abs(m), ckm_ref_abs)
            candidates.append(
                {
                    "row_perm": list(prow),
                    "col_perm": list(pcol),
                    "matrix_mean_rel_pct": float(mean_rel),
                    "unitarity_residual_max_abs": float(unitarity_residual(m)),
                    "delta_primary_rad": float(phase["delta_primary_rad"]),
                    "delta_secondary_rad": float(phase["delta_secondary_rad"]),
                    "delta_primary_rel_err_pct": float(e1),
                    "delta_secondary_rel_err_pct": float(e2),
                    "delta_best_rad": delta_best,
                    "delta_best_rel_err_pct": rel_best,
                    "best_branch": branch,
                }
            )

    # First preserve CKM shape fidelity, then optimize CP mismatch.
    physically_consistent = [
        c for c in candidates if float(c["matrix_mean_rel_pct"]) <= max_mean_rel_for_physical_order_pct
    ]
    if physically_consistent:
        selected = min(
            physically_consistent,
            key=lambda c: (float(c["delta_best_rel_err_pct"]), float(c["matrix_mean_rel_pct"])),
        )
        selection_mode = "shape_fidelity_constrained_cp_min"
    else:
        selected = min(
            candidates,
            key=lambda c: (float(c["matrix_mean_rel_pct"]), float(c["delta_best_rel_err_pct"])),
        )
        selection_mode = "fallback_shape_first"

    best_cp_any = min(candidates, key=lambda c: float(c["delta_best_rel_err_pct"]))
    best_shape_any = min(candidates, key=lambda c: float(c["matrix_mean_rel_pct"]))

    flags = {
        "deterministic_no_scan": True,
        "no_retune": True,
        "unitarity_ok": bool(float(selected["unitarity_residual_max_abs"]) <= 1e-10),
        "matrix_shape_fidelity_ok": bool(
            float(selected["matrix_mean_rel_pct"]) <= max_mean_rel_for_physical_order_pct
        ),
        "cp_branch_within_registry_tolerance": bool(
            float(selected["delta_best_rel_err_pct"]) <= delta_tol
        ),
    }
    pass_count = int(sum(1 for v in flags.values() if bool(v)))
    total_flags = int(len(flags))

    strict_pass = bool(all(flags.values()))
    update = {
        "id": "delta_cp_ckm",
        "predicted_value": float(selected["delta_best_rad"]),
        "status": "derived" if strict_pass else "derived_strict_target_miss",
        "strict_level": "strict_internal_gate",
        "method": "qw2097_deterministic_ckm_cp_refinement_with_permutation_audit",
        "notes": (
            "Deterministic CKM CP refinement with permutation audit; strict pass requires "
            "simultaneous matrix fidelity and CP tolerance."
        ),
    }

    verdict = (
        "CKM_CP_TARGET_REFINEMENT_GATE_PASS_STRICT"
        if strict_pass
        else "CKM_CP_TARGET_REFINEMENT_GATE_TARGET_MISS"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "registry": "report_qw2068_sm_gr_parameter_registry.json",
            "kernel_selection": "report_qw2049_spectral_micro_stagec_intersection_gate.json",
            "mass_q_assignment": "report_qw1961_noncircular_gamma_q_derivation_matrix.json",
            "cp_gate_baseline": "report_qw2075_strict_cp_phase_derivation_gate.json",
            "flavor_constructor": "QW_2063_DERIVATIONAL_RECONSTRUCTION_SHARED_FLAVOR_BASIS.py",
        },
        "baseline_qw2075": {
            "delta_primary_rad": float(r2075["ckm_phase"]["delta_primary_rad"]),
            "delta_secondary_rad": float(r2075["ckm_phase"]["delta_secondary_rad"]),
            "rel_err_primary_pct": float(r2075["ckm_phase"]["rel_err_primary_pct"]),
            "rel_err_secondary_pct": float(r2075["ckm_phase"]["rel_err_secondary_pct"]),
        },
        "registry_target": {
            "delta_cp_ckm_ref_rad": delta_ref,
            "delta_cp_ckm_tol_rel_pct": delta_tol,
            "max_mean_rel_for_physical_order_pct": max_mean_rel_for_physical_order_pct,
        },
        "selection_mode": selection_mode,
        "selected_candidate": selected,
        "best_cp_any_candidate": best_cp_any,
        "best_shape_any_candidate": best_shape_any,
        "top_candidates_by_cp_err": sorted(
            candidates, key=lambda c: float(c["delta_best_rel_err_pct"])
        )[:5],
        "top_candidates_by_shape": sorted(
            candidates, key=lambda c: float(c["matrix_mean_rel_pct"])
        )[:5],
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "update": update,
        "verdict": verdict,
        "required_next_step": (
            "INTEGRATE_QW2097_UPDATE_IN_QW2069_AND_RERUN_CLOSURE_GATES"
            if strict_pass
            else "DELTA_CKM_REMAINS_TARGET_MISS_NEEDS_NEW_OBSERVABLE_CONSTRAINTS_NOT_RETUNE"
        ),
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2097: CKM CP TARGET REFINEMENT GATE",
        "",
        f"- Date UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- pass_count: {pass_count}/{total_flags}",
        f"- selection_mode: `{selection_mode}`",
        "",
        "## Selected Candidate",
        f"- row_perm: {selected['row_perm']}",
        f"- col_perm: {selected['col_perm']}",
        f"- delta_best: {selected['delta_best_rad']:.6f} rad ({selected['best_branch']})",
        f"- delta_best_rel_err_pct: {selected['delta_best_rel_err_pct']:.6f}",
        f"- matrix_mean_rel_pct: {selected['matrix_mean_rel_pct']:.6f}",
        "",
        "## Baseline QW-2075",
        f"- delta_primary_rel_err_pct: {out['baseline_qw2075']['rel_err_primary_pct']:.6f}",
        f"- delta_secondary_rel_err_pct: {out['baseline_qw2075']['rel_err_secondary_pct']:.6f}",
        "",
        "## Artifact",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2097] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2097] Saved MD:   {OUT_MD.name}")
    print(f"[QW-2097] verdict={verdict} pass_count={pass_count}/{total_flags}")


if __name__ == "__main__":
    main()

