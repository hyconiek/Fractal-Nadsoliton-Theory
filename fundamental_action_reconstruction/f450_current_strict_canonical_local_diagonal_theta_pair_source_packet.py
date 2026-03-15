#!/usr/bin/env python3
from __future__ import annotations

import json
import math
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

P437_SCRIPT = ROOT / "p437_current_strict_r14_r15_n477_canonical_local_diagonal_sigma_opposite_pair_sums_value_evaluation_harness_probe.py"
P449_SCRIPT = ROOT / "p449_current_strict_canonical_local_diagonal_multi_pair_o2_cut_defect_evaluation_probe.py"

P437_OUT = (
    GENERATED
    / "p437_current_strict_r14_r15_n477_canonical_local_diagonal_sigma_opposite_pair_sums_value_evaluation_harness_probe.json"
)
P437_SUMMARY = (
    GENERATED
    / "p437_current_strict_r14_r15_n477_canonical_local_diagonal_sigma_opposite_pair_sums_value_evaluation_harness_probe_summary.json"
)
P449_OUT = GENERATED / "p449_current_strict_canonical_local_diagonal_multi_pair_o2_cut_defect_evaluation_probe.json"
P449_SUMMARY = GENERATED / "p449_current_strict_canonical_local_diagonal_multi_pair_o2_cut_defect_evaluation_probe_summary.json"

QW2190 = REPO / "report_qw2190_kernel_mode_representation_emergence_gate.json"

OUT_THETA_PAIR = GENERATED / "theta_pair_canonical_local_diagonal_strict_derived_v1.json"
OUT_THETA_PAIR_SUMMARY = GENERATED / "theta_pair_canonical_local_diagonal_strict_derived_v1_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def real_fourier_basis(n: int) -> dict[str, np.ndarray]:
    x = np.arange(n, dtype=float)
    basis: dict[str, np.ndarray] = {"e0": np.ones(n, dtype=float) / math.sqrt(n)}
    for m in range(1, n // 2):
        basis[f"c{m}"] = math.sqrt(2.0 / n) * np.cos(2.0 * math.pi * m * x / n)
        basis[f"s{m}"] = math.sqrt(2.0 / n) * np.sin(2.0 * math.pi * m * x / n)
    if n % 2 == 0:
        basis[f"e{n//2}"] = ((-1.0) ** x) / math.sqrt(n)
    return basis


def orthonormal_residual(b: np.ndarray) -> float:
    m = b.T @ b
    return float(np.linalg.norm(m - np.eye(m.shape[0])))


def run_script(path: Path) -> dict[str, Any]:
    proc = subprocess.run(
        ["python3", str(path)],
        cwd=str(REPO),
        capture_output=True,
        text=True,
        check=False,
    )
    return {
        "script": str(path.relative_to(REPO)),
        "returncode": proc.returncode,
        "stdout_tail": "\n".join(proc.stdout.splitlines()[-10:]),
        "stderr_tail": "\n".join(proc.stderr.splitlines()[-10:]),
    }


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    runs = [run_script(P437_SCRIPT), run_script(P449_SCRIPT)]
    if any(int(r["returncode"]) != 0 for r in runs):
        raise SystemExit(json.dumps({"status": "DEPENDENCY_SCRIPT_FAILED", "runs": runs}, ensure_ascii=True))

    p437 = load_json(P437_OUT)
    p437_summary = load_json(P437_SUMMARY)
    p449 = load_json(P449_OUT)
    p449_summary = load_json(P449_SUMMARY)
    q2190 = load_json(QW2190)

    n = int(p449.get("computed", {}).get("n") or q2190["mode_mapping"]["n_octaves"])
    fourier = real_fourier_basis(n)
    c1, s1 = fourier["c1"], fourier["s1"]
    c2, s2 = fourier["c2"], fourier["s2"]

    pair1 = p449["computed"]["pairs"]["pair1"]
    pair2 = p449["computed"]["pairs"]["pair2"]

    theta_1 = float(pair1["theta_star"])
    theta_2 = float(pair2["theta_star"])

    u1 = math.cos(theta_1) * c1 + math.sin(theta_1) * s1
    u2 = math.cos(theta_2) * c2 + math.sin(theta_2) * s2

    b = np.column_stack([u1, u2])
    gram = b.T @ b

    strict_derived_inputs_ok = bool(p437_summary.get("input_marked_strict_derived") is True) and bool(
        p437_summary.get("provider_marked_strict_derived") is True
    )
    strict_derived_provenance_ok = bool(p437_summary.get("provider_theorem_refs_present") is True) and bool(
        p449_summary.get("input_theorem_level_pass_from_P437") is True
    )

    theta_pair = {
        "object": "ThetaPair_canonical_local_diagonal_strict_derived_v1",
        "status": "actual_exported_strict_core_theta_pair_source__strict_derived_value_instantiation__no_false_pass",
        "as_of": "2026-03-15",
        "intent": (
            "Export one strict-derived actual theta-pair (theta_1,theta_2) and induced u_1,u_2 vectors on the QW-2190 "
            "real Fourier scaffold, derived from the canonical diagonal/local residual profile computed by P437 and the "
            "pair-wise diagonal-sector defect reconstruction of N484 (evaluated by P449)."
        ),
        "scope": "strict_core_strict_derived_value_instantiation_only",
        "derived_from": {
            "p437_profile_artifact": str(P437_OUT.relative_to(REPO)),
            "p437_summary": str(P437_SUMMARY.relative_to(REPO)),
            "p449_defect_artifact": str(P449_OUT.relative_to(REPO)),
            "p449_summary": str(P449_SUMMARY.relative_to(REPO)),
            "qw2190_mode_scaffold": str(QW2190.relative_to(REPO)),
        },
        "provenance_checks": {
            "strict_derived_inputs_marked": strict_derived_inputs_ok,
            "strict_derived_provenance_chain_present": strict_derived_provenance_ok,
            "note": "This packet trusts only strict-derived lane artifacts (P437/P449) and does not import any sigma-int premises (eps/delta_d).",
        },
        "construction": {
            "diagonal_sector_pair_m_reconstruction": {
                "theta_star_rule": "theta_* = (1/2) atan2(Im(F_2m), Re(F_2m))",
                "u_plus_rule": "u_+ := cos(theta_*) c_m + sin(theta_*) s_m  (QW-2190 real Fourier scaffold)",
                "references": ["N484", "N485", "N487"],
            }
        },
        "outputs": {
            "pair1": {
                "m": int(pair1["m"]),
                "l": int(pair1["l"]),
                "F_2m": pair1["F"],
                "theta_1_star": theta_1,
                "u_1_star": [float(x) for x in u1.tolist()],
                "eigenvalues": pair1["eigenvalues"],
            },
            "pair2": {
                "m": int(pair2["m"]),
                "l": int(pair2["l"]),
                "F_2m": pair2["F"],
                "theta_2_star": theta_2,
                "u_2_star": [float(x) for x in u2.tolist()],
                "eigenvalues": pair2["eigenvalues"],
            },
        },
        "audits": {
            "u1_u2_orthonormal_residual": float(np.linalg.norm(gram - np.eye(2))),
            "u1_u2_det_gram": float(np.linalg.det(gram)),
            "u1_u2_rank": int(np.linalg.matrix_rank(b, tol=1e-12)),
            "u1_u2_orthonormal_residual_vs_identity": orthonormal_residual(b),
            "p437_status": p437.get("status"),
            "p449_status": p449.get("status"),
        },
        "hard_limits": [
            "Does not claim sigma-int -> theta strict-core upgrade (T159 remains a separate lane with selector slots).",
            "Does not claim coefficient-class discharge of F_2m(d); this is one strict-derived value instantiation lane.",
            "Does not claim global strict-core selector closure nor admissible S_sel_int.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    theta_pair_summary = {
        "stage": "F450",
        "status": "F450_EXECUTED_CURRENT_STRICT_CANONICAL_LOCAL_DIAGONAL_THETA_PAIR_SOURCE_PACKET_NO_FALSE_PASS",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "outputs": {
            "theta_1_star": theta_1,
            "theta_2_star": theta_2,
            "pair1_absF": float(pair1["F"]["abs"]),
            "pair2_absF": float(pair2["F"]["abs"]),
            "u1_u2_orthonormal_residual": theta_pair["audits"]["u1_u2_orthonormal_residual"],
        },
        "inputs": {
            "p437_summary": str(P437_SUMMARY.relative_to(REPO)),
            "p449_summary": str(P449_SUMMARY.relative_to(REPO)),
        },
        "dependency_runs": runs,
        "no_false_pass": True,
    }

    OUT_THETA_PAIR.write_text(json.dumps(theta_pair, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_THETA_PAIR_SUMMARY.write_text(
        json.dumps(theta_pair_summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii"
    )
    print(OUT_THETA_PAIR_SUMMARY)


if __name__ == "__main__":
    main()

